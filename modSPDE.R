# this script fits SPDE model the data and generates predictions

# generate default priors for SPDE model
# from Lindgren Rue (2015) "Bayesian Spatial Modelling with R-INLA"
# sigma0: field standard deviation
# NOTE: by default, this constructs a spde prior with unit median marginal variance 
#       and median effective range equal to a fifth of the spatial range 
# or use inla.spde2.pcmatern (possibly allow (1/4,4) variance here rather than (1/2,2))
# U and alpha are the threshold and probability of crossing the threshold for the variance prior
getSPDEPrior = function(mesh, U=1, alpha=0.01, medianRange=NULL) {
  size <- min(c(diff(range(mesh$loc[, 1])), diff(range(mesh$loc[, 2])))) # 1444.772
  # sigma0=1
  # range0 <- size/5
  # kappa0 <- sqrt(8)/range0
  # tau0 <- 1/(sqrt(4 * pi) * kappa0 * sigma0)
  # spde <- inla.spde2.matern(mesh, B.tau = cbind(log(tau0), -1, +1),
  #                           B.kappa = cbind(log(kappa0), 0, -1), theta.prior.mean = c(0, 0),
  #                           theta.prior.prec = c(0.1, 1))
  
  if(is.null(medianRange))
    range0 <- size/5
  else
    range0 = medianRange
  spde = inla.spde2.pcmatern(mesh, prior.range=c(range0, 0.5), prior.sigma = c(U, alpha))
  spde
}

# get a reasonable default mesh triangulation for the SPDE model for [-1,1] x [-1,1] spatial domain
getSPDEMesh = function(locs=cbind(c(-1, -1, 1, 1), c(-1, 1, -1, 1)), n=3500, max.n=5000, doPlot=TRUE, max.edge=c(.01, .1), 
                            offset=-.08, cutoff=.005) {
  
  
  # generate mesh on R2
  mesh = inla.mesh.2d(loc.domain=locs, n=n, max.n=max.n, offset=offset, cutoff=cutoff, max.edge=max.edge)
  
  # plot the mesh if user wants
  if(doPlot) {
    plot(mesh)
  }
  
  mesh
}

# get a reasonable default mesh triangulation for the SPDE model for the Kenya data
getSPDEMeshKenya = function(locs=NULL, n=5000, max.n=5000, doPlot=FALSE, max.edge=c(7, 200), 
                            offset=-.08, cutoff=4) {
  
  if(is.null(locs)) {
    out = load("../U5MR/kenyaDataEd.RData")
    locs=cbind(ed$east, ed$north)
  }
  
  # generate mesh on R2
  mesh = inla.mesh.2d(loc=locs, n=n, max.n=max.n, offset=offset, cutoff=cutoff, max.edge=max.edge)
  
  # plot the mesh if user wants
  if(doPlot) {
    plot(mesh)
    plotMapDat(project=TRUE, border="blue")
  }
  
  mesh
}

# use the fitSPDE function to fit SPDE model to binomial data within Kenya
fitSPDEKenyaDat = function(dat=NULL, dataType=c("mort", "ed"), 
                           mesh=getSPDEMeshKenya(), prior=getSPDEPrior(mesh), 
                           significanceCI=.8, int.strategy="ccd", strategy="gaussian", 
                           nPostSamples=1000, verbose=TRUE, link=1, seed=123, 
                           urbanEffect=TRUE, clusterEffect=TRUE, 
                           leaveOutRegionName=NULL, kmres=5, doValidation=FALSE, previousFit=NULL, 
                           family=c("binomial", "betabinomial"), leaveOutI=NULL) {
  # load observations
  family = match.arg(family)
  dataType = match.arg(dataType)
  if(is.null(dat)) {
    if(dataType == "mort") {
      out = load("../U5MR/kenyaData.RData")
      dat = mort
    }
    else {
      out = load("../U5MR/kenyaDataEd.RData")
      dat = ed
    }
  }
  
  # construct prediction locations
  if(!is.null(leaveOutRegionName) || !is.null(leaveOutI)) {
    if(!is.null(leaveOutRegionName) && !is.null(leaveOutI))
      stop("Neither leaveOutRegionName nor leaveOutI are NULL. At least one of them must be")
    
    # determine what observations will be left out for validation sake
    regionNames = countyToRegion(dat$admin1)
    if(!is.null(leaveOutRegionName))
      leaveOutI = regionNames == leaveOutRegionName
    
    # modify prediction locations and covariates based on left out region
    obsCoords = cbind(dat$east, dat$north)
    predPts = obsCoords[leaveOutI,]
    predsUrban = dat$urban[leaveOutI]
    predClusterI = rep(TRUE, nrow(predPts))
    
    # modify observation locations and covariates based on left out region
    dat = dat[!leaveOutI,]
  } else {
    # make prediction coordinates on a fine grid for plotting, and add coarser grid of testing points
    out = load(paste0("dataPointsKenya.RData"))
    xRangeDat = dataPointsKenya$xRange
    yRangeDat = dataPointsKenya$yRange
    
    # mx = 100
    # my = 100
    # predPts = make.surface.grid(list(x=seq(xRangeDat[1], xRangeDat[2], l=mx), y=seq(yRangeDat[1], yRangeDat[2], l=my)))
    
    # # remove grid points outside of Kenya national boundaries
    # load("../U5MR/adminMapData.RData")
    # polys = adm0@polygons
    # kenyaPoly = polys[[1]]@Polygons[[77]]@coords
    # kenyaPolyProj = projKenya(kenyaPoly)
    # inKenya = in.poly(predPts, kenyaPolyProj)
    # predPts = predPts[inKenya,]
    
    # get prediction locations from population grid (population density adjustment 
    # for target population doesn't matter here, since we only need the prediction 
    # grid at this point)
    if(kmres == 5) {
      load("../U5MR/popGrid.RData")
    }
    else
      popGrid = makeInterpPopGrid(kmres, FALSE)
    
    predPts = cbind(popGrid$east, popGrid$north)
    predsUrban = popGrid$urban
    
    # # add other testing locations to matrix of prediction locations and remember which 
    # plotGridI = 1:sum(inKenya)
    # gridTestI = (max(plotGridI) + 1):(max(plotGridI) + length(simulationData$xGrid))
    # predPts = rbind(predPts, cbind(simulationData$xGrid, simulationData$yGrid))
    
    predClusterI = rep(FALSE, nrow(predPts))
  }
  xPred = matrix(rep(1, nrow(predPts)), ncol=1)
  
  # set observations
  obsValues = dat$y
  obsCoords = cbind(dat$east, dat$north)
  obsNs = dat$n
  xObs = matrix(rep(1, length(obsValues)), ncol=1)
  obsUrban = dat$urban
  
  # add urban effect to the design matrix if necessary
  if(urbanEffect) {
    # predsUrban = getUrbanRural(predPts)
    if(any(is.na(predsUrban))) {
      nans = is.na(predsUrban)
      goodPoints = which(!nans)
      predsUrban = predsUrban[goodPoints]
      predPts = predPts[goodPoints,]
      xPred = matrix(xPred[goodPoints,], ncol=ncol(xPred))
      predClusterI = predClusterI[goodPoints]
    }
    if(any(is.na(obsUrban))) {
      stop("Some observations not in any counties")
    }
    
    xObs = cbind(xObs, obsUrban)
    xPred = cbind(xPred, predsUrban)
  }
  
  c(fitSPDE(obsCoords, obsValues, xObs, 
            predPts, xPred, 
            mesh, prior, 
            significanceCI, int.strategy, strategy, 
            nPostSamples, verbose, link, seed, 
            family, obsNs, clusterEffect, predClusterI, doValidation, previousFit), 
    list(obsCoords=obsCoords, obsValues=obsValues, xObs=xObs, xPred=xPred, obsNs=obsNs, obsUrban=obsUrban, 
         predPts=predPts, predClusterI=predClusterI, kmres=kmres)
    )
}

# use the fitSPDEKenyaDat function to validate SPDE model to binomial data within Kenya with leave one 
# region out validation, and prediction at cluster level
validateSPDEKenyaDat = function(dat=NULL, dataType=c("mort", "ed"), 
                                mesh=getSPDEMeshKenya(), prior=getSPDEPrior(mesh), 
                                significanceCI=.8, int.strategy="ccd", strategy="gaussian", 
                                nPostSamples=1000, verbose=FALSE, link=1, seed=123, 
                                urbanEffect=TRUE, clusterEffect=FALSE, kmres=5, 
                                loadPreviousFit=TRUE, saveResults=TRUE, family=c("binomial", "betabinomial"), 
                                sampleTable=NULL, stratifiedValidation=TRUE, loadPreviousResults=FALSE) {
  if(!is.null(seed))
    set.seed(seed)
  family = match.arg(family)
  
  # load observations
  dataType = match.arg(dataType)
  if(is.null(dat)) {
    if(dataType == "mort") {
      out = load("../U5MR/kenyaData.RData")
      dat = mort
    }
    else {
      out = load("../U5MR/kenyaDataEd.RData")
      dat = ed
    }
  }
  if(dataType == "mort")
    fileNameRoot = "Mort"
  else
    fileNameRoot = "Ed"
  familyText = "Bin"
  if(family == "betabinomial")
    familyText = "BBin"
  
  # first fit the full model (we will use this to initialize the model during the validation fits for each left out county)
  
  fileName = paste0("savedOutput/validation/resultsSPDE", fileNameRoot, "ValidationFull", "_", familyText, "_clust", clusterEffect, 
                    "_urb", urbanEffect, ".RData")
  if(!loadPreviousFit || !file.exists(fileName)) {
    print("Fitting full model")
    time = system.time(fit <- fitSPDEKenyaDat(dat, dataType, mesh, prior, significanceCI, int.strategy, strategy, nPostSamples, 
                                              verbose, link, NULL, urbanEffect, clusterEffect,
                                              kmres=kmres, doValidation=TRUE, family=family))
    
    # get observations and prediction summary statistics
    truth = (dat$y / dat$n)
    obsUrban = dat$urban
    est = fit$obsPreds
    vars = fit$obsSDs^2
    # lower = fit$obsLower
    # upper = fit$obsUpper
    lower = NULL
    upper = NULL
    estMat = fit$obsMat
    estMatBinomial = addBinomialVar(estMat, dat$n)
    
    cpo = fit$mod$cpo$cpo
    cpoFailure = fit$mod$cpo$failure
    dic = fit$mod$dic$dic
    waic = fit$mod$waic$waic
    modelFit = fit$mod
    
    # calculate validation scoring rules
    print("Pooled scores:")
    fullPooledScoresBinomial = data.frame(c(getScores(truth, est, vars, lower, upper, estMatBinomial, doRandomReject=TRUE), WAIC=waic, DIC=dic, CPO=mean(cpo, na.rm=TRUE), Time=time[3]))
    print(fullPooledScoresBinomial)
    print("Rural scores:")
    fullRuralScoresBinomial = data.frame(c(getScores(truth[!obsUrban], est[!obsUrban], vars[!obsUrban], lower[!obsUrban], upper[!obsUrban], estMatBinomial[!obsUrban,], doRandomReject=TRUE), WAIC=NA, DIC=NA, CPO=mean(cpo[!obsUrban], na.rm=TRUE), Time=time[3]))
    print(fullRuralScoresBinomial)
    print("Urban scores:")
    fullUrbanScoresBinomial = data.frame(c(getScores(truth[obsUrban], est[obsUrban], vars[obsUrban], lower[obsUrban], upper[obsUrban], estMatBinomial[obsUrban,], doRandomReject=TRUE), WAIC=NA, DIC=NA, CPO=mean(cpo[obsUrban], na.rm=TRUE), Time=time[3]))
    print(fullUrbanScoresBinomial)
    
    fullPooledScores = data.frame(c(getScores(truth, est, vars, lower, upper, estMat), WAIC=waic, DIC=dic, CPO=mean(cpo, na.rm=TRUE), Time=time[3]))
    fullRuralScores = data.frame(c(getScores(truth[!obsUrban], est[!obsUrban], vars[!obsUrban], lower[!obsUrban], upper[!obsUrban], estMat[!obsUrban,]), WAIC=NA, DIC=NA, CPO=mean(cpo[!obsUrban], na.rm=TRUE), Time=time[3]))
    fullUrbanScores = data.frame(c(getScores(truth[obsUrban], est[obsUrban], vars[obsUrban], lower[obsUrban], upper[obsUrban], estMat[obsUrban,]), WAIC=NA, DIC=NA, CPO=mean(cpo[obsUrban], na.rm=TRUE), Time=time[3]))
    
    if(saveResults)
      save(time=time, fit=fit, fullPooledScoresBinomial=fullPooledScoresBinomial, fullRuralScoresBinomial=fullRuralScoresBinomial, fullUrbanScoresBinomial=fullUrbanScoresBinomial, 
           fullPooledScores=fullPooledScores, fullRuralScores=fullRuralScores, fullUrbanScores=fullUrbanScores, file=fileName)
  }
  else {
    print("Loading previous full model fit")
    load(fileName)
  }
  previousFit = fit
  
  # set up sample table of indices if using stratified validation
  if(stratifiedValidation && is.null(sampleTable))
    sampleTable = getValidationI(dat=dat, dataType=dataType)
  
  # get region names
  allRegions = countyToRegion(dat$admin1)
  regions = sort(unique(allRegions))
  if(!stratifiedValidation)
    nFold = length(regions)
  else
    nFold = 8
  
  # calculate bins for nearest neighbor distances
  distanceMax = 0
  for(i in 1:nFold) {
    if(!stratifiedValidation) {
      thisRegion = regions[i]
      thisSampleI = allRegions == thisRegion
      leaveOutI = NULL
    } else {
      thisRegion = NULL
      leaveOutI = sampleTable[,i]
      thisSampleI = leaveOutI
    }
    
    ##### Break scores down by distance
    predPts = cbind(dat$east, dat$north)[thisSampleI,]
    obsCoords = cbind(dat$east, dat$north)[!thisSampleI,]
    predUrban = dat$urban[thisSampleI]
    obsUrban = dat$urban[!thisSampleI]
    
    # first calculate all distances, broken down by urban, rural, and all aggregated observations
    distMatuu = rdist(obsCoords[!obsUrban,], predPts[!predUrban,])
    distMatuU = rdist(obsCoords[!obsUrban,], predPts[predUrban,])
    distMatUu = rdist(obsCoords[obsUrban,], predPts[!predUrban,])
    distMatUU = rdist(obsCoords[obsUrban,], predPts[predUrban,])
    distMatAu = rdist(obsCoords, predPts[!predUrban,])
    distMatAU = rdist(obsCoords, predPts[predUrban,])
    distMatuA = rdist(obsCoords[!obsUrban,], predPts)
    distMatUA = rdist(obsCoords[obsUrban,], predPts)
    distMatAA = rdist(obsCoords, predPts)
    
    # now calculate nearest distances
    nndistsuu = apply(distMatuu, 2, function(x) {min(x[x != 0])})
    nndistsuU = apply(distMatuU, 2, function(x) {min(x[x != 0])})
    nndistsUu = apply(distMatUu, 2, function(x) {min(x[x != 0])})
    nndistsUU = apply(distMatUU, 2, function(x) {min(x[x != 0])})
    nndistsAu = apply(distMatAu, 2, function(x) {min(x[x != 0])})
    nndistsAU = apply(distMatAU, 2, function(x) {min(x[x != 0])})
    nndistsuA = apply(distMatuA, 2, function(x) {min(x[x != 0])})
    nndistsUA = apply(distMatUA, 2, function(x) {min(x[x != 0])})
    nndistsAA = apply(distMatAA, 2, function(x) {min(x[x != 0])})
    tempMax = c(nndistsuu, nndistsuU, nndistsUu, nndistsUU, nndistsAu, nndistsAU, nndistsuA, nndistsuA, nndistsUA, nndistsUA, nndistsAA, nndistsAA)
    distanceMax = max(distanceMax, tempMax)
  }
  distanceBreaks = seq(0, distanceMax+1, l=20)
  
  # set up sample table of indices if using stratified validation
  if(stratifiedValidation && is.null(sampleTable))
    sampleTable = getValidationI(dat=dat, dataType=dataType)
  
  completeScoreTableBinomial = c()
  pooledScoreTableBinomial = c()
  urbanScoreTableBinomial = c()
  ruralScoreTableBinomial = c()
  completeScoreTable = c()
  pooledScoreTable = c()
  urbanScoreTable = c()
  ruralScoreTable = c()
  binnedScoringRulesuuAll = list()
  binnedScoringRulesuUAll = list()
  binnedScoringRulesUuAll = list()
  binnedScoringRulesUUAll = list()
  binnedScoringRulesAuAll = list()
  binnedScoringRulesAUAll = list()
  binnedScoringRulesuAAll = list()
  binnedScoringRulesUAAll = list()
  binnedScoringRulesAAAll = list()
  binnedScoringRulesuuBinomialAll = list()
  binnedScoringRulesuUBinomialAll = list()
  binnedScoringRulesUuBinomialAll = list()
  binnedScoringRulesUUBinomialAll = list()
  binnedScoringRulesAuBinomialAll = list()
  binnedScoringRulesAUBinomialAll = list()
  binnedScoringRulesuABinomialAll = list()
  binnedScoringRulesUABinomialAll = list()
  binnedScoringRulesAABinomialAll = list()
  singleScores = c()
  singleScoresBinomial = c()
  startFrom = 1
  
  # load previous results if necessary
  fileName = paste0("savedOutput/validation/resultsSPDE", fileNameRoot, "ValidationAllTemp", "_clust", clusterEffect, 
                    "_urb", urbanEffect, ".RData")
  if(loadPreviousResults && file.exists(fileName)) {
    load(fileName)
    startFrom = i+1
  }
  
  for(i in startFrom:nFold) {
    if(!stratifiedValidation) {
      thisRegion = regions[i]
      thisSampleI = allRegions == thisRegion
      print(paste0("Validating SPDE model with urban=", urbanEffect, ", cluster=", clusterEffect, 
                   ", region=", thisRegion, " (", i, "/", length(regions), " regions)"))
      leaveOutI = NULL
    } else {
      thisRegion = NULL
      print(paste0("Validating SPDE model with urban=", urbanEffect, ", cluster=", clusterEffect, 
                   ", (", i, "/", nFold, " folds)"))
      leaveOutI = sampleTable[,i]
      thisSampleI = leaveOutI
    }
    
    # fit the model
    time = system.time(fit <- fitSPDEKenyaDat(dat, dataType, mesh, prior, significanceCI, int.strategy, strategy, nPostSamples, 
                                              verbose, link, NULL, urbanEffect, clusterEffect, thisRegion,
                                              kmres, TRUE, previousFit, family=family, leaveOutI=leaveOutI))
    
    # get observations and prediction summary statistics
    truth = (dat$y / dat$n)[thisSampleI]
    obsUrban = dat$urban[thisSampleI]
    est = fit$preds
    vars = fit$sigmas^2
    # lower = fit$lower
    # upper = fit$upper
    lower = NULL
    upper = NULL
    estMat = fit$predMat
    estMatBinomial = addBinomialVar(estMat, dat$n[thisSampleI])
    
    # calculate validation scoring rules
    print("Pooled scores:")
    if(!stratifiedValidation)
      thisPooledScoresBinomial = data.frame(c(list(Region=thisRegion), getScores(truth, est, vars, lower, upper, estMatBinomial, doRandomReject=TRUE), Time=time[3]))
    else
      thisPooledScoresBinomial = data.frame(c(list(Fold=i), getScores(truth, est, vars, lower, upper, estMatBinomial, doRandomReject=TRUE), Time=time[3]))
    print(thisPooledScoresBinomial)
    
    if(stratifiedValidation || thisRegion != "Nairobi") {
      print("Rural scores:")
      if(!stratifiedValidation)
        thisRuralScoresBinomial = data.frame(c(list(Region=thisRegion), getScores(truth[!obsUrban], est[!obsUrban], vars[!obsUrban], lower[!obsUrban], upper[!obsUrban], estMatBinomial[!obsUrban,], doRandomReject=TRUE), Time=time[3]))
      else
        thisRuralScoresBinomial = data.frame(c(list(Fold=i), getScores(truth[!obsUrban], est[!obsUrban], vars[!obsUrban], lower[!obsUrban], upper[!obsUrban], estMatBinomial[!obsUrban,], doRandomReject=TRUE), Time=time[3]))
      print(thisRuralScoresBinomial)
      
      if(!stratifiedValidation)
        thisRuralScores = data.frame(c(list(Region=thisRegion), getScores(truth[!obsUrban], est[!obsUrban], vars[!obsUrban], lower[!obsUrban], upper[!obsUrban], estMat[!obsUrban,]), Time=time[3]))
      else
        thisRuralScores = data.frame(c(list(Fold=i), getScores(truth[!obsUrban], est[!obsUrban], vars[!obsUrban], lower[!obsUrban], upper[!obsUrban], estMat[!obsUrban,]), Time=time[3]))
    } else {
      thisRuralScoresBinomial = thisPooledScoresBinomial
      thisRuralScoresBinomial[,2:(ncol(thisRuralScoresBinomial)-1)] = NA
      thisRuralScores = thisRuralScoresBinomial
    }
    
    print("Urban scores:")
    if(!stratifiedValidation)
      thisUrbanScoresBinomial = data.frame(c(list(Region=thisRegion), getScores(truth[obsUrban], est[obsUrban], vars[obsUrban], lower[obsUrban], upper[obsUrban], estMatBinomial[obsUrban,], doRandomReject=TRUE), Time=time[3]))
    else
      thisUrbanScoresBinomial = data.frame(c(list(Fold=i), getScores(truth[obsUrban], est[obsUrban], vars[obsUrban], lower[obsUrban], upper[obsUrban], estMatBinomial[obsUrban,], doRandomReject=TRUE), Time=time[3]))
    print(thisUrbanScoresBinomial)
    
    if(!stratifiedValidation)
      thisPooledScores = data.frame(c(list(Region=thisRegion), getScores(truth, est, vars, lower, upper, estMat), Time=time[3]))
    else
      thisPooledScores = data.frame(c(list(Fold=i), getScores(truth, est, vars, lower, upper, estMat), Time=time[3]))
    if(!stratifiedValidation)
      thisUrbanScores = data.frame(c(list(Region=thisRegion), getScores(truth[obsUrban], est[obsUrban], vars[obsUrban], lower[obsUrban], upper[obsUrban], estMat[obsUrban,]), Time=time[3]))
    else
      thisUrbanScores = data.frame(c(list(Fold=i), getScores(truth[obsUrban], est[obsUrban], vars[obsUrban], lower[obsUrban], upper[obsUrban], estMat[obsUrban,]), Time=time[3]))
    
    # append scoring rule tables
    completeScoreTable = rbind(completeScoreTable, thisPooledScores)
    completeScoreTable = rbind(completeScoreTable, thisRuralScores)
    completeScoreTable = rbind(completeScoreTable, thisUrbanScores)
    
    pooledScoreTable = rbind(pooledScoreTable, thisPooledScores)
    ruralScoreTable = rbind(ruralScoreTable, thisRuralScores)
    urbanScoreTable = rbind(urbanScoreTable, thisUrbanScores)
    
    completeScoreTableBinomial = rbind(completeScoreTableBinomial, thisPooledScoresBinomial)
    completeScoreTableBinomial = rbind(completeScoreTableBinomial, thisRuralScoresBinomial)
    completeScoreTableBinomial = rbind(completeScoreTableBinomial, thisUrbanScoresBinomial)
    
    pooledScoreTableBinomial = rbind(pooledScoreTableBinomial, thisPooledScoresBinomial)
    ruralScoreTableBinomial = rbind(ruralScoreTableBinomial, thisRuralScoresBinomial)
    urbanScoreTableBinomial = rbind(urbanScoreTableBinomial, thisUrbanScoresBinomial)
    
    ##### Break scores down by distance
    predPts = fit$predPts
    obsCoords = fit$obsCoords
    predUrban = dat$urban[thisSampleI]
    obsUrban = dat$urban[!thisSampleI]
    
    # first calculate all distances, broken down by urban, rural, and all aggregated observations
    distMatuU = rdist(obsCoords[!obsUrban,], predPts[predUrban,])
    distMatUU = rdist(obsCoords[obsUrban,], predPts[predUrban,])
    distMatAU = rdist(obsCoords, predPts[predUrban,])
    distMatuA = rdist(obsCoords[!obsUrban,], predPts)
    distMatUA = rdist(obsCoords[obsUrban,], predPts)
    distMatAA = rdist(obsCoords, predPts)
    
    # now calculate nearest distances
    nndistsuU = apply(distMatuU, 2, function(x) {min(x[x != 0])})
    nndistsUU = apply(distMatUU, 2, function(x) {min(x[x != 0])})
    nndistsAU = apply(distMatAU, 2, function(x) {min(x[x != 0])})
    nndistsuA = apply(distMatuA, 2, function(x) {min(x[x != 0])})
    nndistsUA = apply(distMatUA, 2, function(x) {min(x[x != 0])})
    nndistsAA = apply(distMatAA, 2, function(x) {min(x[x != 0])})
    if(stratifiedValidation || thisRegion != "Nairobi") {
      distMatuu = rdist(obsCoords[!obsUrban,], predPts[!predUrban,])
      distMatUu = rdist(obsCoords[obsUrban,], predPts[!predUrban,])
      distMatAu = rdist(obsCoords, predPts[!predUrban,])
      
      nndistsuu = apply(distMatuu, 2, function(x) {min(x[x != 0])})
      nndistsUu = apply(distMatUu, 2, function(x) {min(x[x != 0])})
      nndistsAu = apply(distMatAu, 2, function(x) {min(x[x != 0])})
      
      binnedScoringRulesuu = getScores(truth[!predUrban], est[!predUrban], vars[!predUrban], lower[!predUrban], upper[!predUrban], distances=nndistsuu, breaks=distanceBreaks)$binnedResults
      binnedScoringRulesUu = getScores(truth[!predUrban], est[!predUrban], vars[!predUrban], lower[!predUrban], upper[!predUrban], distances=nndistsUu, breaks=distanceBreaks)$binnedResults
      binnedScoringRulesAu = getScores(truth[!predUrban], est[!predUrban], vars[!predUrban], lower[!predUrban], upper[!predUrban], distances=nndistsAu, breaks=distanceBreaks)$binnedResults
      
      binnedScoringRulesuuBinomial = getScores(truth[!predUrban], est[!predUrban], vars[!predUrban], estMat=estMatBinomial[!predUrban,], doRandomReject=TRUE, distances=nndistsuu, breaks=distanceBreaks)$binnedResults
      binnedScoringRulesUuBinomial = getScores(truth[!predUrban], est[!predUrban], vars[!predUrban], estMat=estMatBinomial[!predUrban,], doRandomReject=TRUE, distances=nndistsUu, breaks=distanceBreaks)$binnedResults
      binnedScoringRulesAuBinomial = getScores(truth[!predUrban], est[!predUrban], vars[!predUrban], estMat=estMatBinomial[!predUrban,], doRandomReject=TRUE, distances=nndistsAu, breaks=distanceBreaks)$binnedResults
    } else {
      binnedScoringRulesuu = NULL
      binnedScoringRulesUu = NULL
      binnedScoringRulesAu = NULL
      
      binnedScoringRulesuuBinomial = NULL
      binnedScoringRulesUuBinomial = NULL
      binnedScoringRulesAuBinomial = NULL
    }
    
    # calculate scores without accounting for binomial variation
    binnedScoringRulesuU = getScores(truth[predUrban], est[predUrban], vars[predUrban], lower[predUrban], upper[predUrban], distances=nndistsuU, breaks=distanceBreaks)$binnedResults
    binnedScoringRulesUU = getScores(truth[predUrban], est[predUrban], vars[predUrban], lower[predUrban], upper[predUrban], distances=nndistsUU, breaks=distanceBreaks)$binnedResults
    binnedScoringRulesAU = getScores(truth[predUrban], est[predUrban], vars[predUrban], lower[predUrban], upper[predUrban], distances=nndistsAU, breaks=distanceBreaks)$binnedResults
    binnedScoringRulesuA = getScores(truth, est, vars, lower, upper, distances=nndistsuA, breaks=distanceBreaks)$binnedResults
    binnedScoringRulesUA = getScores(truth, est, vars, lower, upper, distances=nndistsUA, breaks=distanceBreaks)$binnedResults
    binnedScoringRulesAA = getScores(truth, est, vars, lower, upper, distances=nndistsAA, breaks=distanceBreaks)$binnedResults
    
    # calculate scores accounting for binomial variation
    binnedScoringRulesuUBinomial = getScores(truth[predUrban], est[predUrban], vars[predUrban], estMat=estMatBinomial[predUrban,], doRandomReject=TRUE, distances=nndistsuU, breaks=distanceBreaks)$binnedResults
    binnedScoringRulesUUBinomial = getScores(truth[predUrban], est[predUrban], vars[predUrban], estMat=estMatBinomial[predUrban,], doRandomReject=TRUE, distances=nndistsUU, breaks=distanceBreaks)$binnedResults
    binnedScoringRulesAUBinomial = getScores(truth[predUrban], est[predUrban], vars[predUrban], estMat=estMatBinomial[predUrban,], doRandomReject=TRUE, distances=nndistsAU, breaks=distanceBreaks)$binnedResults
    binnedScoringRulesuABinomial = getScores(truth, est, vars, estMat=estMatBinomial, doRandomReject=TRUE, distances=nndistsuA, breaks=distanceBreaks)$binnedResults
    binnedScoringRulesUABinomial = getScores(truth, est, vars, estMat=estMatBinomial, doRandomReject=TRUE, distances=nndistsUA, breaks=distanceBreaks)$binnedResults
    binnedScoringRulesAABinomial = getScores(truth, est, vars, estMat=estMatBinomial, doRandomReject=TRUE, distances=nndistsAA, breaks=distanceBreaks)$binnedResults
    
    # concatenate binned scoring rule results
    binnedScoringRulesuuAll = c(binnedScoringRulesuuAll, list(binnedScoringRulesuu))
    binnedScoringRulesuUAll = c(binnedScoringRulesuUAll, list(binnedScoringRulesuU))
    binnedScoringRulesUuAll = c(binnedScoringRulesUuAll, list(binnedScoringRulesUu))
    binnedScoringRulesUUAll = c(binnedScoringRulesUUAll, list(binnedScoringRulesUU))
    binnedScoringRulesAuAll = c(binnedScoringRulesAuAll, list(binnedScoringRulesAu))
    binnedScoringRulesAUAll = c(binnedScoringRulesAUAll, list(binnedScoringRulesAU))
    binnedScoringRulesuAAll = c(binnedScoringRulesuAAll, list(binnedScoringRulesuA))
    binnedScoringRulesUAAll = c(binnedScoringRulesUAAll, list(binnedScoringRulesUA))
    binnedScoringRulesAAAll = c(binnedScoringRulesAAAll, list(binnedScoringRulesAA))
    binnedScoringRulesuuBinomialAll = c(binnedScoringRulesuuBinomialAll, list(binnedScoringRulesuuBinomial))
    binnedScoringRulesuUBinomialAll = c(binnedScoringRulesuUBinomialAll, list(binnedScoringRulesuUBinomial))
    binnedScoringRulesUuBinomialAll = c(binnedScoringRulesUuBinomialAll, list(binnedScoringRulesUuBinomial))
    binnedScoringRulesUUBinomialAll = c(binnedScoringRulesUUBinomialAll, list(binnedScoringRulesUUBinomial))
    binnedScoringRulesAuBinomialAll = c(binnedScoringRulesAuBinomialAll, list(binnedScoringRulesAuBinomial))
    binnedScoringRulesAUBinomialAll = c(binnedScoringRulesAUBinomialAll, list(binnedScoringRulesAUBinomial))
    binnedScoringRulesuABinomialAll = c(binnedScoringRulesuABinomialAll, list(binnedScoringRulesuABinomial))
    binnedScoringRulesUABinomialAll = c(binnedScoringRulesUABinomialAll, list(binnedScoringRulesUABinomial))
    binnedScoringRulesAABinomialAll = c(binnedScoringRulesAABinomialAll, list(binnedScoringRulesAABinomial))
    
    ##### Calculate individual scoring rules
    # calculate the scoring rules, and add nearest neighbor distances for each stratum
    if(!stratifiedValidation) {
      thisSingleScores = data.frame(c(list(Region=thisRegion, dataI=which(thisSampleI), NNDist=nndistsAA, NNDistU=nndistsUA, NNDistu=nndistsuA), getScores(truth, est, vars, lower, upper, estMat, getAverage=FALSE), Time=time[3]))
      thisSingleScoresBinomial = data.frame(c(list(Region=thisRegion, dataI=which(thisSampleI), NNDist=nndistsAA, NNDistU=nndistsUA, NNDistu=nndistsuA), getScores(truth, est, vars, lower, upper, estMatBinomial, getAverage=FALSE), Time=time[3]))
    }
    else {
      thisSingleScores = data.frame(c(list(Fold=i, dataI=which(thisSampleI), NNDist=nndistsAA, NNDistU=nndistsUA, NNDistu=nndistsuA), getScores(truth, est, vars, lower, upper, estMat, getAverage=FALSE), Time=time[3]))
      thisSingleScoresBinomial = data.frame(c(list(Fold=i, dataI=which(thisSampleI), NNDist=nndistsAA, NNDistU=nndistsUA, NNDistu=nndistsuA), getScores(truth, est, vars, lower, upper, estMatBinomial, getAverage=FALSE), Time=time[3]))
    }
    
    # concatenate the results
    singleScoresBinomial = rbind(singleScoresBinomial, thisSingleScoresBinomial)
    singleScores = rbind(singleScores, thisSingleScores)
    
    # save results so far
    save(completeScoreTable, pooledScoreTable, ruralScoreTable, urbanScoreTable, 
         completeScoreTableBinomial, pooledScoreTableBinomial, ruralScoreTableBinomial, urbanScoreTableBinomial, 
         binnedScoringRulesuuAll, binnedScoringRulesuUAll, binnedScoringRulesUuAll, binnedScoringRulesUUAll, 
         binnedScoringRulesAuAll, binnedScoringRulesAUAll, binnedScoringRulesuAAll, binnedScoringRulesUAAll, 
         binnedScoringRulesAAAll, 
         binnedScoringRulesuuBinomialAll, binnedScoringRulesuUBinomialAll, binnedScoringRulesUuBinomialAll, binnedScoringRulesUUBinomialAll, 
         binnedScoringRulesAuBinomialAll, binnedScoringRulesAUBinomialAll, binnedScoringRulesuABinomialAll, binnedScoringRulesUABinomialAll, 
         binnedScoringRulesAABinomialAll, 
         singleScores, singleScoresBinomial, 
         i, file=fileName)
  }
  
  list(completeScoreTable=completeScoreTable, 
       pooledScoreTable=pooledScoreTable, 
       ruralScoreTable=ruralScoreTable, 
       urbanScoreTable=urbanScoreTable, 
       inSamplePooledScores=fullPooledScores, 
       inSampleUrbanScores=fullUrbanScores, 
       inSampleRuralScores=fullRuralScores, 
       
       completeScoreTableBinomial=completeScoreTableBinomial, 
       pooledScoreTableBinomial=pooledScoreTableBinomial, 
       ruralScoreTableBinomial=ruralScoreTableBinomial, 
       urbanScoreTableBinomial=urbanScoreTableBinomial, 
       inSamplePooledScoresBinomial=fullPooledScoresBinomial, 
       inSampleUrbanScoresBinomial=fullUrbanScoresBinomial, 
       inSampleRuralScoresBinomial=fullRuralScoresBinomial, 
       
       binnedScoringRulesuuAll=averageBinnedScores(binnedScoringRulesuuAll), binnedScoringRulesuUAll=averageBinnedScores(binnedScoringRulesuUAll), 
       binnedScoringRulesUuAll=averageBinnedScores(binnedScoringRulesUuAll), binnedScoringRulesUUAll=averageBinnedScores(binnedScoringRulesUUAll), 
       binnedScoringRulesAuAll=averageBinnedScores(binnedScoringRulesAuAll), binnedScoringRulesAUAll=averageBinnedScores(binnedScoringRulesAUAll), 
       binnedScoringRulesuAAll=averageBinnedScores(binnedScoringRulesuAAll), binnedScoringRulesUAAll=averageBinnedScores(binnedScoringRulesUAAll), 
       binnedScoringRulesAAAll=averageBinnedScores(binnedScoringRulesAAAll), 
       binnedScoringRulesuuBinomialAll=averageBinnedScores(binnedScoringRulesuuBinomialAll), binnedScoringRulesuUBinomialAll=averageBinnedScores(binnedScoringRulesuUBinomialAll), 
       binnedScoringRulesUuBinomialAll=averageBinnedScores(binnedScoringRulesUuBinomialAll), binnedScoringRulesUUBinomialAll=averageBinnedScores(binnedScoringRulesUUBinomialAll), 
       binnedScoringRulesAuBinomialAll=averageBinnedScores(binnedScoringRulesAuBinomialAll), binnedScoringRulesAUBinomialAll=averageBinnedScores(binnedScoringRulesAUBinomialAll), 
       binnedScoringRulesuABinomialAll=averageBinnedScores(binnedScoringRulesuABinomialAll), binnedScoringRulesUABinomialAll=averageBinnedScores(binnedScoringRulesUABinomialAll), 
       binnedScoringRulesAABinomialAll=averageBinnedScores(binnedScoringRulesAABinomialAll), 
       
       singleScores=singleScores, singleScoresBinomial=singleScoresBinomial, 
       
       fullModelFit=previousFit)
}

# function for fitting the SPDE model to data
# predClusterI: whether or not cluster effects are included in the posterior draws at the prediction locations
fitSPDE = function(obsCoords, obsValues, xObs=matrix(rep(1, length(obsValues)), ncol=1), 
                   predCoords, xPred = matrix(rep(1, nrow(predCoords)), ncol=1), 
                   mesh=getSPDEMesh(obsCoords), prior=getSPDEPrior(mesh), 
                   significanceCI=.8, int.strategy="grid", strategy="laplace", 
                   nPostSamples=1000, verbose=TRUE, link=1, seed=NULL, 
                   family=c("normal", "binomial", "betabinomial"), obsNs=rep(1, length(obsValues)), clusterEffect=TRUE, 
                   predClusterI=rep(TRUE, nrow(predCoords)), doValidation=FALSE, previousFit=NULL) {
  family = match.arg(family)
  startTime = proc.time()[3]
  if(!is.null(seed))
    set.seed(seed)
  
  startTimeDefineModel = proc.time()[3]
  
  # set beta binomial prior if necessary
  if(family == "betabinomial") {
    if(clusterEffect)
      stop("cluster effect must not be set to TRUE for betaBinomial model")
    
    # The following code uses a PC prior for the beta overdose and perimeter that is 
    # a bit sketchy mathematically. We've therefore opted for a different prior
    # lambda = getLambdapcBeta(U=1, logitU=TRUE, alpha=0.01, p=.5, normalize=TRUE)
    # bbpcPriorTable = getpcBetaLogitTableForINLA(lambda, p=0.5, tailProb=1e-4, n=500)
    # control.family = list(hyper = list(rho = list(prior = bbpcPriorTable)))
    
    # set median at .04 and upper 97.5th pctile at 0.2
    mu = logit(0.04)
    prec = 1/((logit(.2)-logit(.04))/qnorm(.975))^2
    control.family = list(hyper = list(rho = list(prior="logtnormal", param=c(mu, prec))))
  } else {
    control.family = list()
  }
  
  # construct A matrix for observations
  m = nrow(obsCoords)
  AEst = inla.spde.make.A(mesh, loc = obsCoords)
  
  # construct A matrix for predictions
  APred = inla.spde.make.A(mesh, loc = predCoords)
  
  # make inla stack
  ys = obsValues
  n = ncol(AEst) # number of basis elements
  nObs = length(ys) # number of observations
  nPreds = nrow(predCoords)
  latticeInds = 1:n
  cluster = 1:length(ys)
  
  # construct the observation stack 
  if(family == "normal") {
    if(!is.null(xObs)) {
      stack.est = inla.stack(A =list(AEst, 1),
                             effects =list(field=latticeInds, X=xObs),
                             data =list(y=ys, link=1),
                             tag ="est",
                             remove.unused=FALSE)
    } else {
      stack.est = inla.stack(A =list(AEst),
                             effects =list(field=latticeInds),
                             data =list(y=ys, link=1),
                             tag ="est",
                             remove.unused=FALSE)
    }
  } else {
    if(!clusterEffect) {
      if(!is.null(xObs)) {
        stack.est = inla.stack(A =list(AEst, 1),
                               effects =list(field=latticeInds, X=xObs),
                               data =list(y=ys, Ntrials = obsNs, link=1),
                               tag ="est",
                               remove.unused=FALSE)
      } else {
        stack.est = inla.stack(A =list(AEst),
                               effects =list(field=latticeInds),
                               data =list(y=ys, Ntrials = obsNs, link=1),
                               tag ="est",
                               remove.unused=FALSE)
      }
    } else {
      if(!is.null(xObs)) {
        stack.est = inla.stack(A =list(AEst, 1, 1),
                               effects =list(field=latticeInds, X=xObs, cluster=cluster),
                               data =list(y=ys, Ntrials = obsNs, link=1),
                               tag ="est",
                               remove.unused=FALSE)
      } else {
        stack.est = inla.stack(A =list(AEst, 1),
                               effects =list(field=latticeInds, cluster=cluster),
                               data =list(y=ys, Ntrials = obsNs, link=1),
                               tag ="est",
                               remove.unused=FALSE)
      }
    }
  }
  endTimeDefineModel = proc.time()[3]
  totalTimeDefineModel = endTimeDefineModel - startTimeDefineModel
  
  # make mesh index
  mesh.index <- inla.spde.make.index(name = "field", n.spde = prior$n.spde)
  
  # fit model
  control.inla = list(strategy=strategy, int.strategy=int.strategy)
  modeControl = inla.set.control.mode.default()
  if(!is.null(previousFit)) {
    # initialize the fitting process based on a previous optimum
    # modeControl$result = previousFit
    modeControl$theta = previousFit$mode$theta
    modeControl$x = previousFit$mode$x
    modeControl$restart = TRUE
  }
  
  stack.full = stack.est
  stackDat = inla.stack.data(stack.full, spde=prior)
  allQuantiles = c(0.5, (1-significanceCI) / 2, 1 - (1-significanceCI) / 2)
  
  # see: inla.doc("loggamma")
  # shape=.1, scale=10 for unit mean, variance 100 prior
  startModelFitTime = proc.time()[3]
  if(family == "normal") {
    if(!is.null(xObs)) {
      mod = inla(y ~ - 1 + X + f(field, model=prior), 
                 data = stackDat, 
                 control.predictor=list(A=inla.stack.A(stack.full), compute=TRUE, link=stackDat$link, quantiles=allQuantiles), 
                 family=family, verbose=verbose, control.inla=control.inla, 
                 control.compute=list(config=TRUE, cpo=doValidation, dic=doValidation, waic=doValidation), 
                 control.mode=modeControl, 
                 control.fixed=list(quantiles=allQuantiles), 
                 control.family=list(hyper = list(prec = list(prior="loggamma", param=c(0.1,0.1)))))
    } else {
      mod = inla(y ~ - 1 + f(field, model=prior), 
                 data = stackDat, 
                 control.predictor=list(A=inla.stack.A(stack.full), compute=TRUE, link=stackDat$link, quantiles=allQuantiles), 
                 family=family, verbose=verbose, control.inla=control.inla, 
                 control.compute=list(config=TRUE, cpo=doValidation, dic=doValidation, waic=doValidation), 
                 control.mode=modeControl, 
                 control.fixed=list(quantiles=allQuantiles), 
                 control.family=list(hyper = list(prec = list(prior="loggamma", param=c(0.1,0.1)))))
    }
  } else if(family == "binomial" || family == "betabinomial") {
    clusterList = list(param=c(1, 0.01), prior="pc.prec")
    if(!is.null(xObs) && clusterEffect) {
      mod = inla(y ~ - 1 + X + f(field, model=prior) + f(cluster, model="iid", hyper = list(prec = clusterList)), 
                 data = stackDat, 
                 control.predictor=list(A=inla.stack.A(stack.full), compute=TRUE, link=stackDat$link, quantiles=allQuantiles), 
                 family=family, verbose=verbose, control.inla=control.inla, 
                 control.compute=list(config=TRUE, cpo=doValidation, dic=doValidation, waic=doValidation), 
                 control.mode=modeControl, Ntrials=stackDat$Ntrials, 
                 control.fixed=list(quantiles=allQuantiles), control.family = control.family)
    } else if(is.null(xObs) && clusterEffect) {
      mod = inla(y ~ - 1 + f(field, model=prior) + f(cluster, model="iid", hyper = list(prec = clusterList)), 
                 data = stackDat, 
                 control.predictor=list(A=inla.stack.A(stack.full), compute=TRUE, link=stackDat$link, quantiles=allQuantiles), 
                 family=family, verbose=verbose, control.inla=control.inla, 
                 control.compute=list(config=TRUE, cpo=doValidation, dic=doValidation, waic=doValidation), 
                 control.mode=modeControl, Ntrials=stackDat$Ntrials, 
                 control.fixed=list(quantiles=allQuantiles), control.family = control.family)
    } else if(!is.null(xObs) && !clusterEffect) {
      mod = inla(y ~ - 1 + X + f(field, model=prior), 
                 data = stackDat, 
                 control.predictor=list(A=inla.stack.A(stack.full), compute=TRUE, link=stackDat$link, quantiles=allQuantiles), 
                 family=family, verbose=verbose, control.inla=control.inla, 
                 control.compute=list(config=TRUE, cpo=doValidation, dic=doValidation, waic=doValidation), 
                 control.mode=modeControl, Ntrials=stackDat$Ntrials, 
                 control.fixed=list(quantiles=allQuantiles), control.family = control.family)
    } else if(is.null(xObs) && !clusterEffect) {
      mod = inla(y ~ - 1 + f(field, model=prior), 
                 data = stackDat, 
                 control.predictor=list(A=inla.stack.A(stack.full), compute=TRUE, link=stackDat$link, quantiles=allQuantiles), 
                 family=family, verbose=verbose, control.inla=control.inla, 
                 control.compute=list(config=TRUE, cpo=doValidation, dic=doValidation, waic=doValidation), 
                 control.mode=modeControl, Ntrials=stackDat$Ntrials, 
                 control.fixed=list(quantiles=allQuantiles), control.family = control.family)
    }
  } else {
    stop(paste0("Unsupported family: ", family))
  }
  endModelFitTime = proc.time()[3]
  totalModelFitTime = endModelFitTime - startModelFitTime
  
  # get predictive surface, SD, and data
  n = nrow(obsCoords)
  # obsInds = 1:n
  obsInds = inla.stack.index(stack.full, "est")$data
  predInds = inla.stack.index(stack.full, "pred")$data
  index = inla.stack.index(stack.full, "pred")$data
  linpreds = mod[["summary.fitted.values"]]$mean
  linpred.sd = mod[["summary.fitted.values"]]$sd
  
  preds = linpreds
  predSDs = linpred.sd
  
  # generate samples from posterior
  startTimePosteriorSampling = proc.time()[3]
  postSamples = inla.posterior.sample(nPostSamples, mod)
  endTimePosteriorSampling = proc.time()[3]
  totalTimePosteriorSampling = endTimePosteriorSampling - startTimePosteriorSampling
  latentMat = sapply(postSamples, function(x) {x$latent})
  # if(clusterEffect)
  #   clusterVars = sapply(postSamples, function(x) {1 / x$hyperpar[3]})
  if(family == "normal") {
    nuggetVars = sapply(postSamples, function(x) {1 / x$hyperpar[1]})
  } else if(family == "binomial" && clusterEffect) {
    nuggetVars = sapply(postSamples, function(x) {1 / x$hyperpar[3]})
  } else if(family == "betabinomial") {
    overdispersions = sapply(postSamples, function(x) {x$hyperpar[1]})
  }
  
  latentVarNames = rownames(postSamples[[1]]$latent)
  fieldIndices = which(grepl("field", latentVarNames))
  fixedIndices = which(grepl("X", latentVarNames))
  # if(clusterEffect)
  #   clustIndices = grepl("clust", latentVarNames)
  
  # generate (logit in binomial case) predictions: first without cluster effect then add the cluster effect in
  if(length(xPred) != 0)
    fixedPart = xPred  %*% latentMat[fixedIndices,]
  else
    fixedPart = 0
  predMat = fixedPart + APred %*% latentMat[fieldIndices,]
  
  # do the same for the observations
  if(length(xObs) != 0)
    fixedPart = xObs  %*% latentMat[fixedIndices,]
  else
    fixedPart = 0
  obsMat = fixedPart + AEst %*% latentMat[fieldIndices,]
  
  # add in cluster effect if necessary
  if((family == "binomial" && clusterEffect) || family == "normal") {
    # get cluster effect variance
    clusterVars = nuggetVars
    rhos = NULL
    predMatClustEffect = predMat + sweep(matrix(rnorm(length(predMat), sd=rep(sqrt(clusterVars), each=nrow(predMat))), nrow=nrow(predMat)), 1, predClusterI, "*")
    obsMatClustEffect = obsMat + matrix(rnorm(length(obsMat), sd=rep(sqrt(clusterVars), each=nrow(obsMat))), nrow=nrow(obsMat))
  } else if(family == "betabinomial") {
    # get cluster induced overdispersion
    rhos = overdispersions
    predMat = expit(predMat)
    as = sweep(predMat, 2, 1/rhos-1, "*")
    bs = sweep(1-predMat, 2, 1/rhos-1, "*")
    predMatClustEffect = matrix(rbeta(length(predMat), c(as.matrix(as)), c(as.matrix(bs))), nrow=nrow(predMat))
    obsMat = expit(obsMat)
    as = sweep(obsMat, 2, 1/rhos-1, "*")
    bs = sweep(1-obsMat, 2, 1/rhos-1, "*")
    obsMatClustEffect = matrix(rbeta(length(obsMat), c(as.matrix(as)), c(as.matrix(bs))), nrow=nrow(obsMat))
    clusterVars = NULL
  } else {
    clusterVars = NULL
    rhos = NULL
    predMatClustEffect = predMat
    obsMatClustEffect = obsMat
  }
  
  # transform predictions from logit to probability scale
  if(family == "binomial") {
    predMat = expit(predMat)
    predMatClustEffect = expit(predMatClustEffect)
    obsMat = expit(obsMat)
    obsMatClustEffect = expit(obsMatClustEffect)
  }
  
  # get summary statistics
  preds = rowMeans(predMatClustEffect)
  predSDs = apply(predMatClustEffect, 1, sd)
  lower = apply(predMatClustEffect, 1, quantile, probs=(1-significanceCI)/2)
  medians = apply(predMatClustEffect, 1, median)
  upper = apply(predMatClustEffect, 1, quantile, probs=1-(1-significanceCI)/2)
  obsPreds = rowMeans(obsMat)
  obsSDs = apply(obsMatClustEffect, 1, sd)
  obsLower = apply(obsMatClustEffect, 1, quantile, probs=(1-significanceCI)/2)
  obsMedian = apply(obsMat, 1, median)
  obsUpper = apply(obsMatClustEffect, 1, quantile, probs=1-(1-significanceCI)/2)
  
  if(length(xPred) != 0) {
    interceptSummary=mod$summary.fixed[1,1:5]
    fixedEffectSummary=mod$summary.fixed[,1:5]
  } 
  else {
    interceptSummary = matrix(rep(0, 5), nrow=1)
    fixedEffectSummary = mod$summary.fixed
  }
  rangeSummary=mod$summary.hyperpar[2,1:5]
  spatialSDSummary = mod$summary.hyperpar[3,1:5]
  
  # get posterior hyperparameter samples and transform them as necessary
  hyperMat = sapply(postSamples, function(x) {x$hyperpar})
  if(family == "normal") {
    clusterVarI = 1
    spatialRangeI = 2
    spatialSDI = 3
    mat = apply(hyperMat, 2, function(x) {c(totalVar=x[spatialSDI]^2+1/x[clusterVarI], spatialVar=x[spatialSDI]^2, errorVar=1/x[clusterVarI], 
                                            totalSD=sqrt(x[spatialSDI]^2+1/x[clusterVarI]), spatialSD=x[spatialSDI], errorSD=sqrt(1/x[clusterVarI]), 
                                            spatialRange=x[spatialRangeI])})
  } else if(family == "binomial") {
    spatialRangeI = 1
    spatialSDI = 2
    if(!clusterEffect) {
      clusterVarI = NULL
      mat = apply(hyperMat, 2, function(x) {c(totalVar=x[spatialSDI]^2, spatialVar=x[spatialSDI]^2, 
                                              totalSD=x[spatialSDI], spatialSD=x[spatialSDI], 
                                              spatialRange=x[spatialRangeI])})
    } else {
      clusterVarI = 3
      mat = apply(hyperMat, 2, function(x) {c(totalVar=x[spatialSDI]^2+1/x[clusterVarI], spatialVar=x[spatialSDI]^2, errorVar=1/x[clusterVarI], 
                                              totalSD=sqrt(x[spatialSDI]^2+1/x[clusterVarI]), spatialSD=x[spatialSDI], errorSD=sqrt(1/x[clusterVarI]), 
                                              spatialRange=x[spatialRangeI])})
    }
  } else if(family == "betabinomial") {
    overdispersionI = 1
    spatialRangeI = 2
    spatialSDI = 3
    clusterVarI = NULL
    mat = apply(hyperMat, 2, function(x) {c(totalVar=x[spatialSDI]^2, spatialVar=x[spatialSDI]^2, 
                                            totalSD=x[spatialSDI], spatialSD=x[spatialSDI], 
                                            spatialRange=x[spatialRangeI], overdispersion=x[overdispersionI])})
  }
  
  if(clusterEffect || family == "normal")
    hyperNames = c("totalVar", "spatialVar", "errorVar", "totalSD", "spatialSD", "errorSD", "spatialRange")
  else if(family == "binomial")
    hyperNames = c("totalVar", "spatialVar", "totalSD", "spatialSD", "spatialRange")
  else if(family == "betabinomial")
    hyperNames = c("totalVar", "spatialVar", "totalSD", "spatialSD", "spatialRange", "overdispersion")
  rownames(mat) = hyperNames
  
  getSummaryStatistics = function(draws) {
    c(Est=mean(draws), SD=sd(draws), 
      Qlower=quantile(probs=(1 - significanceCI) / 2, draws), 
      Q50=quantile(probs=0.5, draws), 
      Qupper=quantile(probs=1 - (1 - significanceCI) / 2, draws))
  }
  summaryNames = c("Est", "SD", "Qlower", "Q50", "Qupper")
  parameterSummaryTable = t(apply(mat, 1, getSummaryStatistics))
  colnames(parameterSummaryTable) = summaryNames
  
  # separate out default parameter summaries
  if(family == "normal" || clusterEffect) {
    sdSummary=parameterSummaryTable[6,]
    varSummary=parameterSummaryTable[3,]
    rangeSummary=parameterSummaryTable[7,]
  } else {
    sdSummary = matrix(rep(0, 5), nrow=1) # these are specifically error/cluster sd and var
    varSummary = matrix(rep(0, 5), nrow=1)
    rangeSummary=parameterSummaryTable[5,]
  }
  overdispersionSummary = matrix(rep(0, 5), nrow=1)
  if(family == "betabinomial")
    overdispersionSummary=parameterSummaryTable[6,]
  
  endTime = proc.time()[3]
  totalTime = endTime - startTime
  timings = data.frame(totalTime=totalTime, 
                       modelDefineTime=totalTimeDefineModel, 
                       modelFitTime=totalModelFitTime, 
                       posteriorSamplingTime=totalTimePosteriorSampling, 
                       otherTime=totalTime-(totalTimeDefineModel + totalModelFitTime + totalTimePosteriorSampling))
  timings$modelDefinePct = timings$modelDefineTime / timings$totalTime
  timings$modelFitTimePct = timings$modelFitTime / timings$totalTime
  timings$posteriorSamplingTimePct = timings$posteriorSamplingTime / timings$totalTime
  timings$otherTimePct = timings$otherTime / timings$totalTime
  
  list(mod=mod, preds=preds, sigmas=predSDs, lower=lower, median=medians, upper=upper, 
       obsPreds=obsPreds, obsSDs=obsSDs, obsLower=obsLower, obsMedian=obsMedian, obsUpper=obsUpper, 
       mesh=mesh, prior=prior, stack=stack.full, 
       interceptSummary=interceptSummary, fixedEffectSummary=fixedEffectSummary, rangeSummary=rangeSummary, 
       sdSummary=sdSummary, varSummary=varSummary, overdispersionSummary=overdispersionSummary, 
       parameterSummaryTable=parameterSummaryTable, 
       predMat=predMatClustEffect, obsMat=obsMatClustEffect, hyperMat=hyperMat, timings=timings, clusterVars=clusterVars, rhos=rhos)
}

# this function generates results for the simulation study for the SPDE model
# input arguments:
#   argument specifying the dataset type
resultsSPDE = function(randomSeeds=NULL, covType=c("exponential", "matern", "mixture"), rangeText=c("01", "05", "1", "mix"), 
                       maxDataSets=NULL) {
  # determine the type of covariance for the data set
  covType = match.arg(covType)
  
  # determine the spatial range for the data set. No range text means it's a mixture
  rangeText = match.arg(rangeText)
  
  # construct the file name for the desired data set and load it
  if(rangeText == "mix")
    dataText = paste0(covType, "DataSet.RData")
  else
    dataText = paste0(covType, rangeText, "DataSet.RData")
  out = load(dataText)
  dataSets = simulationData
  
  # construct the SPDE mesh using all of the locations from all data sets
  mesh = getSPDEMesh(cbind(c(dataSets$xTrain, dataSets$xTest), c(dataSets$yTrain, dataSets$yTest)))
  
  # generate results for all data sets and return results (TODO: otherVariableNames)
  resultsSPDE = fitModelToDataSets(fitSPDE, dataSets, randomSeeds=randomSeeds, otherArgs=list(mesh=mesh), 
                                   maxDataSets=maxDataSets)
  
  # save results
  fileName = paste0("resultsSPDE_cov", covType, "Range", rangeText, "maxDataSets", maxDataSets, ".RData")
  save(resultsSPDE, file=fileName)
  
  resultsSPDE
}







