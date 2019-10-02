# this script fits SPDE model the data and generates predictions

# generate default priors for SPDE model
# from Lindgren Rue (2015) "Bayesian Spatial Modelling with R-INLA"
# sigma0: field standard deviation
# NOTE: by default, this constructs a spde prior with unit median marginal variance 
#       and median effective range equal to a fifth of the spatial range 
# or use inla.spde2.pcmatern (possibly allow (1/4,4) variance here rather than (1/2,2))
# U and alpha are the threshold and probability of crossing the threshold for the variance prior
getSPDEPrior = function(mesh, sigma0=1, strictPrior=FALSE, U=3, alpha=0.01) {
  size <- min(c(diff(range(mesh$loc[, 1])), diff(range(mesh$loc[, 2])))) # 1277.237
  # range0 <- size/5
  # kappa0 <- sqrt(8)/range0
  # tau0 <- 1/(sqrt(4 * pi) * kappa0 * sigma0)
  # spde <- inla.spde2.matern(mesh, B.tau = cbind(log(tau0), -1, +1),
  #                           B.kappa = cbind(log(kappa0), 0, -1), theta.prior.mean = c(0, 0),
  #                           theta.prior.prec = c(0.1, 1))
  
  range0 <- size/5
  spde = inla.spde2.pcmatern(mesh, prior.range=c(range0, 0.5), prior.sigma = c(U, alpha))
  spde
}

# get a reasonable default mesh triangulation for the SPDE model
getSPDEMesh = function(locs, n=3500, max.n=5000, doPlot=TRUE, max.edge=c(.01, .1), 
                            offset=-.08, cutoff=.005) {
  
  
  # generate mesh on R2
  mesh = inla.mesh.2d(loc.domain=locs, n=n, max.n=max.n, offset=offset, cutoff=cutoff, max.edge=max.edge)
  
  # plot the mesh if user wants
  if(doPlot) {
    plot(mesh)
  }
  
  mesh
}

# function for fitting the SPDE model to gaussian data
fitSPDE = function(obsCoords, obsValues, xObs=matrix(rep(1, length(obsValues)), ncol=1), 
                   predCoords, xPred = matrix(rep(1, nrow(predCoords)), ncol=1), 
                   mesh=getSPDEMesh(obsCoords), prior=getSPDEPrior(mesh), 
                   significanceCI=.8, int.strategy="grid", strategy="laplace", 
                   nPostSamples=1000, verbose=TRUE, link=1, seed=NULL) {
  
  if(!is.null(seed))
    set.seed(seed)
  
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
  
  # construct the observation stack 
  stack.est = inla.stack(A =list(AEst, 1),
                         effects =list(field=latticeInds, X=xObs),
                         data =list(y=ys, link=1),
                         tag ="est",
                         remove.unused=FALSE)
  
  # make mesh index
  mesh.index <- inla.spde.make.index(name = "field", n.spde = prior$n.spde)
  
  # fit model
  control.inla = list(strategy=strategy, int.strategy=int.strategy)
  modeControl = inla.set.control.mode.default()
  
  stack.full = stack.est
  stackDat = inla.stack.data(stack.full, spde=prior)
  allQuantiles = c(0.5, (1-significanceCI) / 2, 1 - (1-significanceCI) / 2)
  
  mod = inla(y ~ - 1 + X + f(field, model=prior), 
             data = stackDat, 
             control.predictor=list(A=inla.stack.A(stack.full), compute=TRUE, link=stackDat$link, quantiles=allQuantiles), 
             family="gaussian", verbose=verbose, control.inla=control.inla, 
             control.compute=list(config=TRUE), 
             control.mode=modeControl)
  
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
  postSamples = inla.posterior.sample(nPostSamples, mod)
  latentMat = sapply(postSamples, function(x) {x$latent})
  # if(clusterEffect)
  #   clusterVars = sapply(postSamples, function(x) {1 / x$hyperpar[3]})
  latentVarNames = rownames(postSamples[[1]]$latent)
  fieldIndices = which(grepl("field", latentVarNames))
  fixedIndices = which(grepl("X", latentVarNames))
  # if(clusterEffect)
  #   clustIndices = grepl("clust", latentVarNames)
  
  # generate logit predictions (first without cluster effect then add the cluster effect in)
  if(length(xPred) != 0)
    fixedPart = xPred  %*% latentMat[fixedIndices,]
  else
    fixedPart = 0
  predMatNoClust = fixedPart + APred %*% latentMat[fieldIndices,]
  predMat <- predMatNoClust
  
  # if(clusterEffect) {
  #   # add in estimated cluster effects at the sampled enumeration areas
  #   predMat[clusterIndices,] = predMat[clusterIndices,] + latentMat[clustIndices,]
  #   
  #   # addend cluster effects we have no information about at the on sampled enumeration areas
  #   predMat[eaIndices[-clusterIndices],] = predMat[eaIndices[-clusterIndices],] + 
  #     rnorm(length(eaIndices[-clusterIndices]) * nPostSamples, sd = rep(sqrt(clusterVars), each=length(eaIndices[-clusterIndices])))
  # }
  
  lower = mod[['summary.pred']]
  
  list(mod=mod, preds=preds, SDs=predSDs, obsInds=obsInds, predInds=index, pixelInds=pixelIndices, 
       eaInds=eaIndices, mesh=mesh, prior=prior, stack=stack.full, 
       countyPreds=countyPreds, regionPreds=regionPreds, pixelPreds=pixelPreds, eaPreds=eaPreds)
}

# this function generates results for the simulation study for the SPDE model
# input arguments:
#   argument specifying the dataset type
resultsSPDE = function(randomSeeds=NULL, covType=c("exponential", "matern", "mixture"), rangeText=c("01", "05", "1", ""), 
                       maxDataSets=NULL) {
  # determine the type of covariance for the data set
  covType = match.arg(covType)
  
  # determine the spatial range for the data set. No range text means it's a mixture
  rangeText = match.arg(rangeText)
  
  # construct the file name for the desired data set and load it
  dataText = paste0(covType, rangeText, "DataSet.RData")
  out = load(dataText)
  dataSets = simulationData
  
  # construct the SPDE mesh using all of the locations from all data sets
  mesh = getSPDEMesh(cbind(c(dataSets$xTrain, dataSets$xTest), c(dataSets$yTrain, dataSets$yTest)))
  
  # generate results for all data sets and return results (TODO: otherVariableNames)
  fitModelToDataSets(fitSPDE, dataSets, randomSeeds=randomSeeds, otherArgs=list(mesh=mesh), 
                     maxDataSets=maxDataSets)
}







