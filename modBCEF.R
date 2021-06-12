# model for fitting BCEF dataset with flexible covariate model

modBCEF = function(dat, predCoords, predPTC, latInfo=NULL, 
                   nu=1.5, seed=1, nLayer=ifelse(is.null(latInfo), 2, length(latInfo)), 
                   NC=c(25, 100), nBuffer=5, priorPar=NULL, 
                   rwPrior=NULL, rwKnots=NULL, 
                   rwModel=c("rw1", "rw2"), nNonlinearBasis=20, 
                   normalize=TRUE, fastNormalize=TRUE, 
                   intStrategy="ccd", strategy="gaussian", 
                   significanceCI=.8, printVerboseTimings=FALSE, nPostSamples=1000, 
                   clusterEffect=TRUE, predictionType=c("mean", "median"), 
                   initialEffectiveRange=NULL, initialAlphas=rep(1/nLayer, nLayer-1), 
                   effRangeRange=NULL, separateRanges=TRUE, doValidation=FALSE, 
                   precomputedNormalizationFun=NULL, loadPrecomputationResults=FALSE, 
                   precomputationFileNameRoot="BCEFprecomputations", 
                   savePrecomputationResults=FALSE, family=c("normal"), 
                   leaveOutI=NULL, previousFit=NULL, verbose=TRUE, diagonal=0.0, 
                   rwConstr=TRUE) {
  rwModel = match.arg(rwModel)
  
  # construct lattice info if necessary
  if(is.null(latInfo)) {
    rotationAngle = 49.5 * (pi/180)
    center = colMeans(cbind(BCEF$x, BCEF$y))
    rotationMat = rbind(c(cos(rotationAngle), -sin(rotationAngle)), 
                        c(sin(rotationAngle), cos(rotationAngle)))
    
    # construct basis functions on the rotated centered coordinate system
    xRange = range(rotatedCentered[,1])
    yRange = range(rotatedCentered[,2])
    latInfo = makeLatGrids(xRange, yRange, NC=NC, nLayer=nLayer)
    
    # transform back to the original coordinate system, and set range of data in correct units
    for(i in 1:length(latInfo)) {
      latInfo[[i]]$latCoords = t(t(rotationMat) %*% t(latInfo[[i]]$latCoords))
      latInfo[[i]]$latCoords = sweep(latInfo[[i]]$latCoords, 2, center, "+")
      
      latInfo[[i]]$xRangeDat = range(dat$x)
      latInfo[[i]]$yRangeDat = range(dat$y)
    }
  }
  
  # setup
  obsCoords = cbind(dat$x, dat$y)
  obsValues = log(dat$FCH)
  xObs = cbind(1, dat$PTC)
  xPred = cbind(1, predPTC)
  
  if(is.null(priorPar)) {
    priorPar = getPCPrior(diff(range(BCEF$x))/5, .01, 5, nLayer=nLayer, separateRanges=separateRanges, latticeInfo=latInfo, useUrbanPrior=FALSE)
  }
  
  if(is.null(rwPrior)) {
    rwPrior = getRandomWalkPriors(obsValues, models=rwModel, priorType="pc.prec", paramList=NULL, family="gaussian", n=NULL)[1]
  }
  
  # set knots if necessary. set spline knots in the same way as we set random walk knots, except we must
  # subtract an amount to get 25 basis functions since we're making cubic/order 3 b splines by default
  if(is.null(rwKnots)) {
    rwKnots = make_knots(c(dat$PTC, predPTC), n=nNonlinearBasis)[[1]]
  }
  
  # must set each value of the covariates to the closest knot
  rwEffects = match_with_knots(dat$PTC, rwKnots)
  rwEffectsNew = match_with_knots(predPTC, rwKnots)
  rwEffectsInds = match(rwEffects, rwKnots)
  rwEffectsIndsNew = match(rwEffectsNew, rwKnots)
  
  startTime = proc.time()[3]
  set.seed(seed)
  
  # get the type of prediction the user wants
  predictionType = match.arg(predictionType)
  family = match.arg(family)
  
  # check to make sure inputs make sense
  if(!clusterEffect && family == "normal")
    stop("cluster effect must be included for the normal family")
  
  # generate lattice basis matrix
  AObs = makeA(obsCoords, latInfo)
  
  ## run precomputations
  print("running precomputations...")
  startTimePrecomputations = proc.time()[3]
  
  # either run the precomputations or load them in
  if(!loadPrecomputationResults) {
    precomputedMatrices = precomputationsQ2(latInfo)
    if(is.null(precomputedNormalizationFun)) {
      precomputedNormalizationFun = precomputeNormalization(saveResults=FALSE, latticeInfo=latInfo, effRangeRange=effRangeRange, 
                                                            plotNormalizationSplines=FALSE)
    }
  } else {
    load(paste0("savedOutput/precomputations/", precomputationFileNameRoot, ".RData"))
  }
  
  endTimePrecomputations = proc.time()[3]
  precomputationTime = endTimePrecomputations - startTimePrecomputations
  print(paste0("finished precomputations. Took ", round(precomputationTime / 60, 2), " minutes"))
  
  # save the precomputations if necessary
  if(savePrecomputationResults) {
    save(precomputedMatrices, precomputedNormalizationFun, precomputationTime, 
         file=paste0("savedOutput/precomputations/", precomputationFileNameRoot, ".RData"))
  }
  
  # define the LK model
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
  
  n = length(obsValues)
  obsNs = rep(1, n) # this won't do anything unless we are in binomial or betabinomial family
  if(separateRanges) {
    rgen = inla.rgeneric.define(model=inla.rgeneric.lk.model.full, latInfo=latInfo, ys=obsValues, 
                                prior=priorPar, normalize=normalize, precomputedMatrices=precomputedMatrices, 
                                X=xObs, nu=nu, datCoords=obsCoords, fastNormalize=fastNormalize, 
                                printVerboseTimings=printVerboseTimings, 
                                initialEffectiveRange=initialEffectiveRange, initialAlphas=initialAlphas, 
                                precomputedNormalizationFun=precomputedNormalizationFun, ns=obsNs)
  } else {
    rgen = inla.rgeneric.define(model=inla.rgeneric.lk.model.standard, latInfo=latInfo, ys=obsValues, 
                                prior=priorPar, normalize=normalize, precomputedMatrices=precomputedMatrices, 
                                X=xObs, nu=nu, datCoords=obsCoords, fastNormalize=fastNormalize, 
                                printVerboseTimings=printVerboseTimings, 
                                initialEffectiveRange=initialEffectiveRange, initialAlphas=initialAlphas, 
                                precomputedNormalizationFun=precomputedNormalizationFun, ns=obsNs)
  }
  
  # use these global variables for testing calls to inla.rgeneric.lk.model.simple
  # latInfo<<-latInfo; ys<<-obsValues; ns<<-obsNs;
  # prior<<-priorPar; normalize<<-normalize; precomputedMatrices<<-precomputedMatrices;
  # X<<-xObs; nu<<-nu; datCoords<<-obsCoords; fastNormalize<<-fastNormalize;
  # printVerboseTimings<<-printVerboseTimings; initialEffectiveRange<<-initialEffectiveRange;
  # initialAlphas<<-initialAlphas; precomputedNormalizationFun<<-precomputedNormalizationFun
  ## generate inla stack:
  # Stacked A matrix (A_s from notation of LR2015 Bayesian Spatial Modelling with R_INLA):
  # (AEst   0  )
  # ( 0   APred)
  # eta_s = (c^T c^T)^T
  # where c is the vector of lattice coefficients
  AEst = makeA(obsCoords, latInfo)
  APred = makeA(predCoords, latInfo)
  latticeInds = 1:ncol(AEst)
  rwInds = 1:nNonlinearBasis
  clust = 1:length(obsValues)
  
  # update the formula to include this random walk effect with the desired knots
  if(family == "normal") {
    if(!is.null(xObs)) {
      stack.est = inla.stack(A =list(AEst, 1, 1), 
                             effects =list(field=latticeInds, X=xObs, rw=rwEffects), 
                             data =list(y=obsValues, link=1), 
                             tag ="est", 
                             remove.unused=FALSE)
      stack.pred = inla.stack(A =list(matrix(APred[1,], nrow=1), 1, 1), 
                              effects =list(field=latticeInds, X=matrix(xPred[1,], nrow=1), rw=rwEffectsNew[1]), 
                              data =list(y=NA, link=1), 
                              tag ="pred", 
                              remove.unused=FALSE)
    } else {
      stop("must include covariates in this application...")
    }
  }
  else if(family == "binomial" || family == "betabinomial") {
    if(clusterEffect) {
      if(!is.null(xObs)) {
        stack.est = inla.stack(A =list(AEst, 1, 1, 1), 
                               effects =list(field=latticeInds, clust=clust, X=xObs, rw=rwEffects), 
                               data =list(y=obsValues, link=1, Ntrials = obsNs), 
                               tag ="est", 
                               remove.unused=FALSE)
      } else {
        stop("must include covariates in this application...")
      }
    } else {
      if(!is.null(xObs)) {
        stack.est = inla.stack(A =list(AEst, 1, 1), 
                               effects =list(field=latticeInds, X=xObs, rw=rwEffects), 
                               data =list(y=obsValues, link=1, Ntrials = obsNs), 
                               tag ="est", 
                               remove.unused=FALSE)
      } else {
        stop("must include covariates in this application...")
      }
    }
    if(!is.null(xObs)) {
      stack.pred = inla.stack(A =list(matrix(APred[1,], nrow=1), 1, 1), 
                              effects =list(field=latticeInds, X=matrix(xPred[1,], nrow=1), rw=rwEffectsNew[1]), 
                              data =list(y=NA, link=1, Ntrials=1, 
                                         tag ="pred"), 
                              remove.unused=FALSE)
    } else {
      stop("must include covariates in this application...")
    }
  }
  stack.full = inla.stack(stack.est, stack.pred, 
                          remove.unused=FALSE)
  dat = inla.stack.data(stack.full, rgen=rgen, remove.unused=FALSE)
  
  # initialize the fitting process based on a previous optimum if necessary
  modeControl = inla.set.control.mode.default()
  if(!is.null(previousFit)) {
    # modeControl$result = previousFit
    modeControl$theta = previousFit$mode$theta
    # modeControl$x = previousFit$mode$x
    modeControl$restart = TRUE
  }
  endTimeDefineModel = proc.time()[3]
  totalTimeDefineModel = endTimeDefineModel - startTimeDefineModel
  
  # fit the model
  # control.inla = list(cmin = 0, int.strategy=int.strategy) 
  # see: inla.doc("loggamma")
  # shape=.1, scale=10 for unit mean, variance 100 prior
  controls = list(strategy=strategy, int.strategy=intStrategy, diagonal=diagonal) 
  allQuantiles = c(0.5, (1-significanceCI) / 2, 1 - (1-significanceCI) / 2)
  startTimeFitModel = proc.time()[3]
  if(family == "normal") {
    if(!is.null(xObs)) {
      mod = inla(y ~ - 1 + X + f(field, model=rgen) + 
                   f(rw, model=rwModel, values=rwKnots, hyper=rwPrior, 
                     scale.model=TRUE, constr=rwConstr), 
                 data=dat, quantiles=allQuantiles, family=family, verbose=verbose, 
                 control.inla=controls, 
                 control.mode=modeControl, 
                 control.compute=list(config=TRUE, cpo=doValidation, dic=doValidation, waic=doValidation), 
                 control.predictor=list(A=inla.stack.A(stack.full), compute=FALSE), 
                 control.fixed=list(quantiles=allQuantiles), 
                 control.family=list(hyper = list(prec = list(prior="loggamma", param=c(0.1,0.1)))))
      # control.family=list(hyper = list(prec = list(param=c(1, 0.05), prior="pc.prec")))
    } else {
      stop("must include covariates in this application...")
    }
  } else if(family == "binomial" || family == "betabinomial") {
    if(clusterEffect) {
      # clusterList = list(param=c(.15, 0.01), prior="pc.prec")
      if(!is.null(xObs)) {
        mod = inla(y ~ - 1 + X + f(field, model=rgen) + 
                     f(clust, model="iid", hyper=list(prec=list(param=c(1, 0.01), prior="pc.prec"))) + 
                     f(rw, model=rwModel, values=rwKnots, hyper=rwPrior, 
                       scale.model=TRUE, constr=rwConstr), 
                   data=dat, quantiles=allQuantiles, family=family, verbose=verbose, 
                   control.inla=controls, Ntrials=dat$Ntrials, 
                   control.mode=modeControl, 
                   control.compute=list(config=TRUE, cpo=doValidation, dic=doValidation, waic=doValidation), 
                   control.predictor=list(A=inla.stack.A(stack.full), compute=FALSE), 
                   control.fixed=list(quantiles=allQuantiles), control.family = control.family)
      } else {
        stop("must include covariates in this application...")
      }
    } else {
      if(!is.null(xObs)) {
        mod = inla(y ~ - 1 + X + f(field, model=rgen) + 
                     f(rw, model=rwModel, values=rwKnots, hyper=rwPrior, 
                       scale.model=TRUE, constr=rwConstr), 
                   data=dat, quantiles=allQuantiles, family=family, verbose=verbose, 
                   control.inla=controls, Ntrials=dat$Ntrials, 
                   control.mode=modeControl, 
                   control.compute=list(config=TRUE, cpo=doValidation, dic=doValidation, waic=doValidation), 
                   control.predictor=list(A=inla.stack.A(stack.full), compute=FALSE), 
                   control.fixed=list(quantiles=allQuantiles), control.family = control.family)
      } else {
        stop("must include covariates in this application...")
      }
    }
  }
  endTimeFitModel = proc.time()[3]
  totalTimeFitModel = endTimeFitModel - startTimeFitModel
  
  print(paste0("finished fitting model. Took ", round(totalTimeFitModel / 60, 2), " minutes"))
  
  # # improve the approximation of the posterior for the hyperparameters if requested by the user
  # if(improveHyperpar) {
  #   browser()
  #   mod = inla.hyperpar(mod)
  # }
  
  # get predictive surface, SD, and data
  index = inla.stack.index(stack.full, "pred")$data
  
  obsInds = 1:n
  predInds = (n+1):(n+nrow(predCoords))
  if(predictionType == "mean") {
    linpreds = mod[["summary.linear.predictor"]]$mean
    # preds = linpreds[predInds]
    # obsPreds = linpreds[obsInds]
  } else {
    linpreds = mod[["summary.linear.predictor"]]$`0.5quant`
    # preds = linpreds[predInds]
    # obsPreds = linpreds[obsInds]
  }
  linpred.sd = mod[["summary.linear.predictor"]]$sd
  # predSDs = linpred.sd[predInds]
  # obsSDs = linpred.sd[obsInds]
  
  # generate samples from posterior
  print("Sampling from posterior")
  startTimePosteriorSampling = proc.time()[3]
  postSamples = inla.posterior.sample(nPostSamples, mod)
  endTimePosteriorSampling = proc.time()[3]
  totalTimePosteriorSampling = endTimePosteriorSampling - startTimePosteriorSampling
  print(paste0("finished sampling from the posterior. Took ", round(totalTimePosteriorSampling / 60, 2), " minutes"))
  
  print("Processing posterior samples...")
  startTimeSampleProcessing = proc.time()[3]
  # get posterior hyperparameter samples and transform them as necessary
  # interpretation of hyperparameters: 
  
  # 1: error precision
  # 2: log effective range
  # 3: log spatial variance
  # 4-(3 + nLayer - 1): multivariateLogit alpha
  # 3 + nLayer: rw Precision
  hyperMat = sapply(postSamples, function(x) {x$hyperpar})
  if(family == "normal") {
    if(!separateRanges) {
      mat = apply(hyperMat, 2, function(x) {c(totalVar=exp(x[3])+1/x[1]+1/x[3 + nLayer], spatialVar=exp(x[3]), errorVar=1/x[1], 
                                              totalSD=sqrt(exp(x[3])+1/x[1]), spatialSD=sqrt(exp(x[3])), errorSD=sqrt(1/x[1]), 
                                              spatialRange=exp(x[2]), alpha=multivariateExpit(x[4:(3 + nLayer - 1)]), 
                                              rwVar=1/x[3 + nLayer], rwSD=1/sqrt(x[3 + nLayer]))})
      mat = rbind(mat, alpha=1-colSums(mat[8:(7+nLayer-1),]))
      mat = mat[c(1:(7+nLayer-1), nrow(mat), (nrow(mat)-2):(nrow(mat)-1)),]
      hyperNames = c("totalVar", "spatialVar", "clusterVar", "totalSD", "spatialSD", "clusterSD", "spatialRange", 
                     paste0("alpha", 1:nLayer), "rwVar", "rwSD")
    } else {
      mat = apply(hyperMat, 2, function(x) {c(totalVar=exp(x[2+nLayer])+1/x[1]+1/x[3+2*nLayer-1], spatialVar=exp(x[2+nLayer]), errorVar=1/x[1], 
                                              totalSD=sqrt(exp(x[2+nLayer])+1/x[1]), spatialSD=sqrt(exp(x[2+nLayer])), errorSD=sqrt(1/x[1]), 
                                              spatialRange=exp(x[2:(1+nLayer)]), alpha=multivariateExpit(x[(3+nLayer):(3+2*nLayer-2)]), 
                                              rwVar=1/x[3+2*nLayer-1], rwSD=1/sqrt(x[3+2*nLayer-1]))})
      mat = rbind(mat, alpha=1-colSums(matrix(mat[(6+nLayer+1):(6+nLayer+1 + nLayer-2),], nrow=nLayer-1)))
      mat = mat[c(1:(6+nLayer+1 + nLayer-2), nrow(mat), (nrow(mat)-2):(nrow(mat)-1)),]
      hyperNames = c("totalVar", "spatialVar", "clusterVar", "totalSD", "spatialSD", "clusterSD", paste0("spatialRange", 1:nLayer), 
                     paste0("alpha", 1:nLayer), "rwVar", "rwSD")
    }
  } else if(family == "binomial") {
    stop("family must be normal for this example")
    if(clusterEffect) {
      if(!separateRanges) {
        logSpatialRangeI = 1
        logSpatialVarI = 2
        logitAlphaI = 3:(2 + nLayer - 1)
        precisionI = 3 + nLayer - 1
        mat = apply(hyperMat, 2, function(x) {c(totalVar=exp(x[logSpatialVarI])+1/x[precisionI], spatialVar=exp(x[logSpatialVarI]), clusterVar=1/x[precisionI], 
                                                totalSD=sqrt(exp(x[logSpatialVarI])+1/x[precisionI]), spatialSD=sqrt(exp(x[logSpatialVarI])), clusterSD=sqrt(1/x[precisionI]), 
                                                spatialRange=exp(x[logSpatialRangeI]), alpha=multivariateExpit(x[logitAlphaI]))})
        mat = rbind(mat, alpha=1-colSums(mat[8:(7+nLayer-1),]))
        hyperNames = c("totalVar", "spatialVar", "clusterVar", "totalSD", "spatialSD", "clusterSD", "spatialRange", 
                       paste0("alpha", 1:nLayer))
      } else {
        logSpatialRangeI = 1:nLayer
        logSpatialVarI = 1 + nLayer
        logitAlphaI = (logSpatialVarI+1):(logSpatialVarI+1 + nLayer-2)
        precisionI = max(logitAlphaI)+1
        mat = apply(hyperMat, 2, function(x) {c(totalVar=exp(x[logSpatialVarI])+1/x[precisionI], spatialVar=exp(x[logSpatialVarI]), clusterVar=1/x[precisionI], 
                                                totalSD=sqrt(exp(x[logSpatialVarI])+1/x[precisionI]), spatialSD=sqrt(exp(x[logSpatialVarI])), clusterSD=sqrt(1/x[precisionI]), 
                                                spatialRange=exp(x[logSpatialRangeI]), alpha=multivariateExpit(x[logitAlphaI]))})
        mat = rbind(mat, alpha=1-colSums(matrix(mat[(6+nLayer+1):(6+nLayer+1 + nLayer-2),], nrow=nLayer-1)))
        hyperNames = c("totalVar", "spatialVar", "clusterVar", "totalSD", "spatialSD", "clusterSD", paste0("spatialRange", 1:nLayer), 
                       paste0("alpha", 1:nLayer))
      }
    } else {
      if(!separateRanges) {
        logSpatialRangeI = 1
        logSpatialVarI = 2
        logitAlphaI = 3:(2 + nLayer - 1)
        mat = apply(hyperMat, 2, function(x) {c(spatialVar=exp(x[logSpatialVarI]), spatialSD=sqrt(exp(x[logSpatialVarI])), 
                                                spatialRange=exp(x[logSpatialRangeI]), alpha=multivariateExpit(x[logitAlphaI]))})
        mat = rbind(mat, alpha=1-colSums(matrix(mat[4:(3+nLayer-1),], nrow=nLayer-1)))
        hyperNames = c("spatialVar", "spatialSD", "spatialRange", paste0("alpha", 1:nLayer))
      } else {
        logSpatialRangeI = 1:nLayer
        logSpatialVarI = 1 + nLayer
        logitAlphaI = (logSpatialVarI+1):(logSpatialVarI+1 + nLayer-2)
        mat = apply(hyperMat, 2, function(x) {c(spatialVar=exp(x[logSpatialVarI]), spatialSD=sqrt(exp(x[logSpatialVarI])), 
                                                spatialRange=exp(x[logSpatialRangeI]), alpha=multivariateExpit(x[logitAlphaI]))})
        # mat = rbind(mat, alpha=1-colSums(mat[4:(3+nLayer-1),]))
        mat = rbind(mat, alpha=1-colSums(matrix(mat[(3+nLayer):(3+nLayer + nLayer-2),], nrow=nLayer-1)))
        hyperNames = c("spatialVar", "spatialSD", paste0("spatialRange", 1:nLayer), paste0("alpha", 1:nLayer))
      }
    }
  } else if(family == "betabinomial") {
    stop("family must be normal for this example")
    overdispersionI = 1
    if(!separateRanges) {
      logSpatialRangeI = 2
      logSpatialVarI = 3
      logitAlphaI = 4:(3 + nLayer - 1)
      mat = apply(hyperMat, 2, function(x) {c(spatialVar=exp(x[logSpatialVarI]), spatialSD=sqrt(exp(x[logSpatialVarI])), 
                                              spatialRange=exp(x[logSpatialRangeI]), alpha=multivariateExpit(x[logitAlphaI]), 
                                              overdispersion=x[overdispersionI])})
      mat = rbind(mat[-nrow(mat),], alpha=1-colSums(matrix(mat[4:(3+nLayer-1),], nrow=nLayer-1)), mat[nrow(mat),])
      hyperNames = c("spatialVar", "spatialSD", "spatialRange", paste0("alpha", 1:nLayer), "overdispersion")
    } else {
      logSpatialRangeI = 2:(1+nLayer)
      logSpatialVarI = 2 + nLayer
      logitAlphaI = (logSpatialVarI+1):(logSpatialVarI+1 + nLayer-2)
      mat = apply(hyperMat, 2, function(x) {c(spatialVar=exp(x[logSpatialVarI]), spatialSD=sqrt(exp(x[logSpatialVarI])), 
                                              spatialRange=exp(x[logSpatialRangeI]), alpha=multivariateExpit(x[logitAlphaI]), 
                                              overdispersion=x[overdispersionI])})
      # mat = rbind(mat, alpha=1-colSums(mat[4:(3+nLayer-1),]))
      mat = rbind(mat[-nrow(mat),], alpha=1-colSums(matrix(mat[(3+nLayer):(3+nLayer + nLayer-2),], nrow=nLayer-1)), mat[nrow(mat),])
      hyperNames = c("spatialVar", "spatialSD", paste0("spatialRange", 1:nLayer), paste0("alpha", 1:nLayer), "overdispersion")
    }
  }
  
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
  varI = which(grepl("Var", hyperNames))
  rangeI = which(grepl("Range", hyperNames))
  alphaI = which(grepl("alpha", hyperNames))
  varSummary=parameterSummaryTable[varI,]
  rangeSummary=parameterSummaryTable[rangeI,]
  alphaSummary=parameterSummaryTable[alphaI,]
  if(family == "normal" || clusterEffect) {
    sdI = which(grepl("SD", hyperNames))
    sdSummary=parameterSummaryTable[sdI,]
  } else {
    sdSummary = matrix(rep(0, 5), nrow=1)
  }
  overdispersionSummary = matrix(rep(0, 5), nrow=1)
  if(family == "betabinomial")
    overdispersionSummary=parameterSummaryTable[nrow(parameterSummaryTable),]
  
  # get samples of the latent field
  latentMat = sapply(postSamples, function(x) {x$latent})
  latentVarNames = rownames(postSamples[[1]]$latent)
  fieldIndices = which(grepl("field", latentVarNames))
  rwIndices = which(grepl("rw", latentVarNames))
  fixedIndices = which(grepl("X", latentVarNames))
  # if(clusterEffect)
  #   clustIndices = grepl("clust", latentVarNames)
  
  ## generate logit predictions (first without cluster effect then add the cluster effect in)
  # for prediction locations
  if(length(xPred) != 0) {
    fixedMat = latentMat[fixedIndices,]
    fixedPart = xPred %*% fixedMat
  }
  else {
    fixedMat = NULL
    fixedPart = 0
  }
  
  predMat = fixedPart + APred %*% latentMat[fieldIndices,] + latentMat[rwIndices,][rwEffectsIndsNew,]
  
  # for observation locations
  if(length(xObs) != 0)
    fixedPart = xObs  %*% fixedMat
  else
    fixedPart = 0
  obsMat = fixedPart + AObs %*% latentMat[fieldIndices,] + latentMat[rwIndices,][rwEffectsInds,]
  
  # get draws from basis function coefficients
  basisCoefMat = latentMat[fieldIndices,]
  rwCoefMat = latentMat[rwIndices,]
  
  # add in cluster effect if necessary
  predClusterI = rep(TRUE, nrow(predCoords))
  if((family == "binomial" && clusterEffect) || family == "normal") {
    # get betabinomial overdispersion parameter
    clusterVarI = which(grepl("clusterVar", hyperNames))
    clusterVars = mat[clusterVarI,]
    predMatClustEffect = predMat + sweep(matrix(rnorm(length(predMat), sd=rep(sqrt(clusterVars), each=nrow(predMat))), nrow=nrow(predMat)), 1, predClusterI, "*")
    obsMatClustEffect = obsMat + matrix(rnorm(length(obsMat), sd=rep(sqrt(clusterVars), each=nrow(obsMat))), nrow=nrow(obsMat))
    rhos = NULL
  } else if(family == "betabinomial") {
    # get cluster induced overdispersion
    overdispersionI = which(grepl("overdispersion", hyperNames))
    rhos = mat[overdispersionI,]
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
    rhos=NULL
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
  
  # compute predictive credible intervals
  if(is.null(predMatClustEffect)) {
    preds = rowMeans(predMat)
    predSDs = apply(predMat, 1, sd)
    lowerPreds = apply(predMat, 1, quantile, probs=(1-significanceCI)/2)
    medianPreds = apply(predMat, 1, median)
    upperPreds = apply(predMat, 1, quantile, probs=1-(1-significanceCI)/2)
    obsPreds = rowMeans(obsMat)
    obsSDs = apply(obsMat, 1, sd)
    lowerObs = apply(obsMat, 1, quantile, probs=(1-significanceCI)/2)
    medianObs = apply(obsMat, 1, median)
    upperObs = apply(obsMat, 1, quantile, probs=1-(1-significanceCI)/2)
  } else {
    if(family == "normal") {
      preds = rowMeans(predMat)
      medianPreds = apply(predMat, 1, median)
      obsPreds = rowMeans(obsMat)
      medianObs = apply(obsMat, 1, median)
    } else {
      preds = rowMeans(predMatClustEffect)
      medianPreds = apply(predMatClustEffect, 1, median)
      obsPreds = rowMeans(obsMatClustEffect)
      medianObs = apply(obsMatClustEffect, 1, median)
    }
    predSDs = apply(predMatClustEffect, 1, sd)
    lowerPreds = apply(predMatClustEffect, 1, quantile, probs=(1-significanceCI)/2)
    upperPreds = apply(predMatClustEffect, 1, quantile, probs=1-(1-significanceCI)/2)
    obsSDs = apply(obsMatClustEffect, 1, sd)
    lowerObs = apply(obsMatClustEffect, 1, quantile, probs=(1-significanceCI)/2)
    upperObs = apply(obsMatClustEffect, 1, quantile, probs=1-(1-significanceCI)/2)
  }
  
  if(!is.null(xObs) && all(xObs[,1]==1))
    interceptSummary=mod$summary.fixed[1,c(1, 2, 4, 3, 5)]
  else
    interceptSummary = matrix(rep(0, 5), nrow=1)
  
  if(!is.null(xObs))
    fixedEffectSummary = mod$summary.fixed[,c(1, 2, 4, 3, 5)]
  else
    fixedEffectSummary = mod$summary.fixed
  
  rwSummary = mod$summary.random$rw[,c(1, 2, 4, 3, 5, 6)]
  
  # compute basis function coefficient predictions and standard deviations
  nx = latInfo[[1]]$nx
  ny = latInfo[[1]]$ny
  # startI = n+nrow(predCoords)+1
  # endI = n+nrow(predCoords)+nx*ny
  startI = 1
  endI = nx*ny
  thesePreds = rowMeans(basisCoefMat)
  theseSDs = apply(basisCoefMat, 1, sd)
  coefPreds = list(layer1 = thesePreds[startI:endI])
  coefSDs = list(layer1 = theseSDs[startI:endI])
  if(nLayer >= 2) {
    for(i in 2:nLayer) {
      startI = endI + 1
      endI = startI + nrow(latInfo[[i]]$latCoords) - 1
      # coefPreds = c(coefPreds, list(linpreds[startI:endI]))
      # coefSDs = c(coefSDs, list(linpred.sd[startI:endI]))
      coefPreds = c(coefPreds, list(thesePreds[startI:endI]))
      coefSDs = c(coefSDs, list(theseSDs[startI:endI]))
    }
  }
  
  endTime = proc.time()[3]
  browser()
  totalTimeSampleProcessing = endTime - startTimeSampleProcessing
  print(paste0("finished processing samples. Took ", round(totalTimeSampleProcessing / 60, 2), " minutes"))
  totalTime = endTime - startTime
  timings = data.frame(totalTime=totalTime, 
                       precomputationTime=precomputationTime, 
                       modelDefineTime=totalTimeDefineModel, 
                       modelFitTime=totalTimeFitModel, 
                       posteriorSamplingTime=totalTimePosteriorSampling, 
                       sampleProcessingTime=totalTimeSampleProcessing, 
                       otherTime=totalTime-(precomputationTime + totalTimeDefineModel + totalTimeFitModel + totalTimePosteriorSampling + totalTimeSampleProcessing))
  timings$precomputationTimePct = timings$precomputationTime / timings$totalTime
  timings$modelDefinePct = timings$modelDefineTime / timings$totalTime
  timings$modelFitTimePct = timings$modelFitTime / timings$totalTime
  timings$posteriorSamplingTimePct = timings$posteriorSamplingTime / timings$totalTime
  timings$sampleProcessingTimePct = timings$sampleProcessingTime / timings$totalTime
  timings$otherTimePct = timings$otherTime / timings$totalTime
  
  list(preds=preds, sigmas=predSDs, lower=lowerPreds, median=medianPreds, upper=upperPreds, 
       obsPreds=obsPreds, obsSDs=obsSDs, obsLower=lowerObs, obsMedian=medianObs, obsUpper=upperObs, 
       mod=mod, latInfo=latInfo, coefPreds=coefPreds, coefSDs=coefSDs, 
       interceptSummary=interceptSummary, fixedEffectSummary=fixedEffectSummary, fixedMat=fixedMat, rangeSummary=rangeSummary, 
       sdSummary=sdSummary, varSummary=varSummary, overdispersionSummary=overdispersionSummary, parameterSummaryTable=parameterSummaryTable, 
       alphaSummary=alphaSummary, timings=timings, priorPar=priorPar, precomputedNormalizationFun=precomputedNormalizationFun, 
       # the rest of the outputs are saved to be used for spatial aggregations later on
       predMat=predMatClustEffect, obsMat=obsMatClustEffect, hyperMat=hyperMat, clusterVars=clusterVars, rhos=rhos, 
       rwSummary=rwSummary, rwKnots=rwKnots, rwMat=rwCoefMat, rwPrior=rwPrior
  )
}


# construct priors for random walk effects
## Inputs:
# y: the observations, only used lim constructing the priors for the gaussian family, as 
#    recommended by INLA
# models: either a single string, or a list of them with the same length as the number of 
#         nonlinear effects in the model. The strings describe the order of random walk of the 
#         nonlinear effect.
# priorType: either a single string, or a list of them with the same length as the number 
#             of nonlinear effects in the model. The strings describe the type of model for 
#             the prior of each nonlinear effect
# paramList: a list of param vectors, which are inputs to the random walk models in inla. 
#            each param vector contains 2 values, which have different interpretations 
#            depending on the type of prior used. For the log gamma prior, the first element 
#            param is the shape parameter, and the second is the inverse scale. For the 
#            pc.prec prior, the first element is the standard deviation threshold, and 
#            the second is the probability the standard deviation is over that threshold. 
#            See inla.doc("rw1"), inla.doc("rw2"), inla.doc("loggamma"), and 
#            inla.doc("pc.prec") for more details on priors for random walks. Must be same 
#            length as the desired value of n.
# n: the number of priors in the output list to generate
# family: the distribution family of the observations, as in a GLM
getRandomWalkPriors = function(y, models=c("rw1", "rw2"), priorType=c("pc.prec", "loggamma"), paramList=NULL, family="gaussian", n=NULL) {
  ##### get user inputs
  # determine the number of priors and make sure the number of priors specified by each input is consistent
  if(!is.null(paramList))
    n0 = length(paramList)
  else
    n0 = 1
  if(!identical(models, c("rw1", "rw2")))
    n1 = length(models)
  else {
    n1 = 1
    models = "rw1"
  }
  if(!identical(priorType, c("pc.prec", "loggamma")))
    n2 = length(priorType)
  else {
    n2 = 1
    priorType = "pc.prec"
  }
  nTemp = max(n0, n1, n2)
  if(nTemp != 1 && !is.null(n) && nTemp != n)
    stop("length of either paramList, models, or priorType is not 1 and does not match the specified number of random walk priors")
  if(is.null(n))
    n = nTemp
  
  # now expand input arguments to match the number of random walk priors
  if(length(priorType) == 1)
    priorType = rep(priorType, n)
  if(length(models) == 1)
    models = rep(models, n)
  if(length(paramList) == 1)
    paramList = rep(paramList, n)
  
  # if family is a string, get the actual family
  if(is.character(family))
    family = do.call(family, list())
  
  # set values of param in paramList
  if(is.null(paramList)) {
    paramList = list()
    
    if(any(priorType == "loggamma"))
      warning("INLA's default log gamma priors for random walk log precision are not recommended (shape=1, inverse scale=.01). Instead, the user should use a penalized complexity prior, or set their own the log gamma parameters")
    
    for(i in 1:length(priorType)) {
      prior = priorType[i]
      
      if(prior == "loggamma") {
        # these are the parameters listed in the example in inla.doc("rw1")
        u = 1
      }
      else if(prior == "pc.prec") {
        # these are the parameters recommended in inla.doc("rw1")
        if(family$family == "gaussian")
          u = sd(y)
        else if(family$family == "poisson" && family$link == "log")
          u = 1
        else if(family$family == "binomial" && family$link == "logit")
          u = .5
        else if(family$family == "binomial" && family$link == "probit")
          u = .33
        else
          stop(paste0("No default prior values for family '", family$family, "' and link '", family$link, "'. In this case, user must specify paramList"))
      }
      paramList = c(paramList, list(c(u, 0.01)))
    }
  }
  
  # construct list of priors
  rwPriors = list()
  for(i in 1:n) {
    rwPriors = c(rwPriors, list(prec=list(prior=priorType[i], param=paramList[[i]])))
  }
  
  # return results
  rwPriors
}

# Makes knots for random walk and b splines.  Either generates equally spaced knots or 
# non-equally spaced but instead based on quantiles from equally increasing probabilities.
# x: covariate values over which to construct default knot points. It can be a matrix, 
#    where each column corresponds to a covariate
# n: number of knots. Can be a vector with length equal to the number of columns of x
# quantiles: either used
# NOTE: for cubic b splines there are 3 more basis functions than knots
# NOTE2: fda package requires knots to cover range of x
make_knots = function(x, n=25, quantiles=FALSE, return.vector=FALSE) {
  # ensure x has the correct format
  if(prod(dim(x)) == 0) {
    return(NULL)
  }
  if(class(x) != "data.frame") {
    x = data.frame(x)
  }
  
  # ensure n is a vector with length equal to the number of columns of x
  if(length(n) == 1) {
    n = rep(n, ncol(x))
  }
  else if(length(n) != ncol(x)) {
    stop("n must be either scalar or have the same number of columns as x")
  }
  
  # generate the knots for each covariate
  make_knots_single = function(i) {
    rangeVals = range(x[, i])
    if(!quantiles)
      seq(rangeVals[1], rangeVals[2], l=n[i])
    else {
      ps = seq(0, 1, l=n[i])
      c(rangeVals[1], quantile(x[, i], probs=ps[2:(n[i]-1)]), rangeVals[2])
    }
  }
  
  knots = lapply(1:length(n), make_knots_single)
  if(return.vector) {
    unlist(knots)
  }
  else {
    knots
  }
}

# for each element of vals, rounds it to the nearest value in knots. Assumes knots is sorted and increasing
match_with_knots = function(vals, knots, returnIndices=FALSE) {
  # get basic info about the knots
  minVal = knots[1]
  maxVal = knots[length(knots)]
  
  # threshold values
  if(any(vals < maxVal)) {
    # warning("some coordinates larger than maximum of knots")
    vals[vals > maxVal] = maxVal
  }
  if(any(vals < minVal)) {
    # warning("some coordinates smaller than minimum of knots")
    vals[vals < minVal] = minVal
  }
  
  # round the values to the knots
  findIndex = function(i) {
    # first find the first knot the value is smaller than
    index = match(TRUE, vals[i] <= knots)
    
    # if NA, closest to the biggest knot.  If 1, closest to smallest knot
    if(is.na(index))
      return(length(knots))
    else if(index == 1)
      return(index)
    else {
      # otherwise find which of the 2 knots near the returned index the value is closer to
      dist1 = abs(vals[i] - knots[index])
      dist2 = abs(vals[i] - knots[index - 1])
      if(dist1 <= dist2)
        index
      else
        index - 1
    }
  }
  
  # get the indices and knots the values are closest to
  indices = sapply(1:length(vals), findIndex)
  result = knots[indices]
  
  # return results
  if(returnIndices)
    return(list(values=result, indices=indices))
  else
    return(result)
}








