# fit standard LK model given data and relevant parameters
# obsCoords: spatial coordinates of the data
# obsValues: observations at the coordinates
# predCoords: places at which to make predictions
# nu: Matern smoothness parameter from "simple" LKrig model
# xObs: observation design matrix (intercept and covariates at obs locations)
# xPred: prediction design matrix (intercept and covariates at pred locations)
# NC: number of coarse lattice points over the longest dimension of the data
# int.strategy: "auto" or "eb" for empirical Bayes
# predClusterI: whether or not to include cluster effects for which predictions
# priorPar default was getPrior(.1, .1, 10, nLayer=nLayer)
fitLKINLAStandard2 = function(obsCoords, obsValues, predCoords=obsCoords, nu=1.5, seed=1, nLayer=3, NC=5,
                              nBuffer=5, priorPar=NULL, 
                              xObs=cbind(1, obsCoords), xPred=cbind(1, predCoords), normalize=TRUE, 
                              intStrategy="grid", strategy="gaussian", fastNormalize=TRUE, 
                              predictionType=c("mean", "median"), significanceCI=0.8, 
                              printVerboseTimings=FALSE, nPostSamples=1000, family=c("normal", "binomial", "betabinomial"),
                              obsNs=rep(1, length(obsValues)), clusterEffect=TRUE, latInfo=NULL, 
                              initialEffectiveRange=NULL, initialAlphas=rep(1/nLayer, nLayer-1), 
                              effRangeRange=NULL, predClusterI=rep(TRUE, nrow(predCoords)), 
                              plotNormalizationSplines=FALSE, verbose=TRUE, separateRanges=FALSE, 
                              doValidation=FALSE, previousFit=NULL, precomputedNormalizationFun=NULL, 
                              useUrbanPrior=FALSE, savePrecomputationResults=FALSE, loadPrecomputationResults=FALSE, 
                              precomputationFileNameRoot="precomputationResults") {
  
  startTime = proc.time()[3]
  set.seed(seed)
  
  # get the type of prediction the user wants
  predictionType = match.arg(predictionType)
  family = match.arg(family)
  
  # check to make sure inputs make sense
  if(!clusterEffect && family == "normal")
    stop("cluster effect must be included for the normal family")
  
  # get lattice points, prediction points
  
  if(is.null(latInfo)) {
    xRangeDat = range(obsCoords[,1])
    yRangeDat = range(obsCoords[,2])
    latInfo = makeLatGrids(xRangeDat, yRangeDat, NC, nBuffer, nLayer)
  } else {
    xRangeDat = latInfo[[1]]$xRangeDat
    yRangeDat = latInfo[[1]]$yRangeDat
  }
  nx = latInfo[[1]]$nx
  ny = latInfo[[1]]$ny
  
  if(is.null(priorPar))
    priorPar = getPCPrior(1444.772/5, .01, 1, nLayer=nLayer, separateRanges=separateRanges, latticeInfo=latInfo, useUrbanPrior=useUrbanPrior) # 1444.772/5
  
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
                                                            plotNormalizationSplines=plotNormalizationSplines)
    }
  } else {
    load(paste0("savedOutput/precomputations/", precomputationFileNameRoot, ".RData"))
  }
  
  # save the precomputations if necessary
  if(savePrecomputationResults) {
    save(precomputedMatrices, precomputedNormalizationFun, 
         file=paste0("savedOutput/precomputations/", precomputationFileNameRoot, ".RData"))
  }
  
  endTimePrecomputations = proc.time()[3]
  precomputationTime = endTimePrecomputations - startTimePrecomputations
  print(paste0("finished precomputations. Took ", round(precomputationTime / 60, 2), " minutes"))
  
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
  clust = 1:length(obsValues)
  if(family == "normal") {
    if(!is.null(xObs)) {
      stack.est = inla.stack(A =list(AEst, 1), 
                             effects =list(field=latticeInds, X=xObs), 
                             data =list(y=obsValues, link=1), 
                             tag ="est", 
                             remove.unused=FALSE)
      stack.pred = inla.stack(A =list(APred, 1), 
                              effects =list(field=latticeInds, X=xPred), 
                              data =list(y=NA, link=1), 
                              tag ="pred", 
                              remove.unused=FALSE)
    } else {
      stack.est = inla.stack(A =list(AEst), 
                             effects =list(field=latticeInds), 
                             data =list(y=obsValues, link=1), 
                             tag ="est", 
                             remove.unused=FALSE)
      stack.pred = inla.stack(A =list(APred), 
                              effects =list(field=latticeInds), 
                              data =list(y=NA, link=1), 
                              tag ="pred", 
                              remove.unused=FALSE)
    }
  }
  else if(family == "binomial" || family == "betabinomial") {
    if(clusterEffect) {
      if(!is.null(xObs)) {
        stack.est = inla.stack(A =list(AEst, 1, 1), 
                               effects =list(field=latticeInds, clust=clust, X=xObs), 
                               data =list(y=obsValues, link=1, Ntrials = obsNs), 
                               tag ="est", 
                               remove.unused=FALSE)
      } else {
        stack.est = inla.stack(A =list(AEst, 1), 
                               effects =list(field=latticeInds, clust=clust), 
                               data =list(y=obsValues, link=1, Ntrials = obsNs), 
                               tag ="est", 
                               remove.unused=FALSE)
      }
    } else {
      if(!is.null(xObs)) {
        stack.est = inla.stack(A =list(AEst, 1), 
                               effects =list(field=latticeInds, X=xObs), 
                               data =list(y=obsValues, link=1, Ntrials = obsNs), 
                               tag ="est", 
                               remove.unused=FALSE)
      } else {
        stack.est = inla.stack(A =list(AEst), 
                               effects =list(field=latticeInds), 
                               data =list(y=obsValues, link=1, Ntrials = obsNs), 
                               tag ="est", 
                               remove.unused=FALSE)
      }
    }
    if(!is.null(xObs)) {
      stack.pred = inla.stack(A =list(APred, 1), 
                              effects =list(field=latticeInds, X=xPred), 
                              data =list(y=NA, link=1, Ntrials=rep(1, nrow(APred)), 
                                         tag ="pred"), 
                              remove.unused=FALSE)
    } else {
      stack.pred = inla.stack(A =list(APred), 
                              effects =list(field=latticeInds), 
                              data =list(y=NA, link=1, Ntrials=rep(1, nrow(APred)), 
                                         tag ="pred"), 
                              remove.unused=FALSE)
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
    modeControl$x = previousFit$mode$x
    modeControl$restart = TRUE
  }
  endTimeDefineModel = proc.time()[3]
  totalTimeDefineModel = endTimeDefineModel - startTimeDefineModel
  
  # fit the model
  # control.inla = list(cmin = 0, int.strategy=int.strategy) 
  # see: inla.doc("loggamma")
  # shape=.1, scale=10 for unit mean, variance 100 prior
  controls = list(strategy=strategy, int.strategy=intStrategy) 
  allQuantiles = c(0.5, (1-significanceCI) / 2, 1 - (1-significanceCI) / 2)
  startTimeFitModel = proc.time()[3]
  if(family == "normal") {
    if(!is.null(xObs)) {
      mod = inla(y ~ - 1 + X + f(field, model=rgen), 
                 data=dat, quantiles=allQuantiles, family=family, verbose=verbose, 
                 control.inla=controls, 
                 control.mode=modeControl, 
                 control.compute=list(config=TRUE, cpo=doValidation, dic=doValidation, waic=doValidation), 
                 control.predictor=list(A=inla.stack.A(stack.full), compute=FALSE), 
                 control.fixed=list(quantiles=allQuantiles), 
                 control.family=list(hyper = list(prec = list(param=c(1, 0.05), prior="pc.prec"))))
    } else {
      mod = inla(y ~ - 1 + f(field, model=rgen), 
                 data=dat, quantiles=allQuantiles, family=family, verbose=verbose, 
                 control.inla=controls, 
                 control.mode=modeControl, 
                 control.compute=list(config=TRUE, cpo=doValidation, dic=doValidation, waic=doValidation), 
                 control.predictor=list(A=inla.stack.A(stack.full), compute=FALSE), 
                 control.fixed=list(quantiles=allQuantiles), 
                 control.family=list(hyper = list(prec = list(param=c(1, 0.05), prior="pc.prec"))))
    }
  } else if(family == "binomial" || family == "betabinomial") {
    if(clusterEffect) {
      # clusterList = list(param=c(.15, 0.01), prior="pc.prec")
      if(!is.null(xObs)) {
        mod = inla(y ~ - 1 + X + f(field, model=rgen) + 
                     f(clust, model="iid", hyper = list(prec = list(param=c(1, 0.01), prior="pc.prec"))), 
                   data=dat, quantiles=allQuantiles, family=family, verbose=verbose, 
                   control.inla=controls, Ntrials=dat$Ntrials, 
                   control.mode=modeControl, 
                   control.compute=list(config=TRUE, cpo=doValidation, dic=doValidation, waic=doValidation), 
                   control.predictor=list(A=inla.stack.A(stack.full), compute=FALSE), 
                   control.fixed=list(quantiles=allQuantiles), control.family = control.family)
      } else {
        mod = inla(y ~ - 1 + f(field, model=rgen) + 
                     f(clust, model="iid", hyper = list(prec = list(param=c(1, 0.01), prior="pc.prec"))), 
                   data=dat, quantiles=allQuantiles, family=family, verbose=verbose, 
                   control.inla=controls, Ntrials=dat$Ntrials, 
                   control.mode=modeControl, 
                   control.compute=list(config=TRUE, cpo=doValidation, dic=doValidation, waic=doValidation), 
                   control.predictor=list(A=inla.stack.A(stack.full), compute=FALSE), 
                   control.fixed=list(quantiles=allQuantiles), control.family = control.family)
      }
    } else {
      if(!is.null(xObs)) {
        mod = inla(y ~ - 1 + X + f(field, model=rgen), 
                   data=dat, quantiles=allQuantiles, family=family, verbose=verbose, 
                   control.inla=controls, Ntrials=dat$Ntrials, 
                   control.mode=modeControl, 
                   control.compute=list(config=TRUE, cpo=doValidation, dic=doValidation, waic=doValidation), 
                   control.predictor=list(A=inla.stack.A(stack.full), compute=FALSE), 
                   control.fixed=list(quantiles=allQuantiles), control.family = control.family)
      } else {
        mod = inla(y ~ - 1 + f(field, model=rgen), 
                   data=dat, quantiles=allQuantiles, family=family, verbose=verbose, 
                   control.inla=controls, Ntrials=dat$Ntrials, 
                   control.mode=modeControl, 
                   control.compute=list(config=TRUE, cpo=doValidation, dic=doValidation, waic=doValidation), 
                   control.predictor=list(A=inla.stack.A(stack.full), compute=FALSE, quantiles=allQuantiles), 
                   control.fixed=list(quantiles=allQuantiles), control.family = control.family)
      }
    }
  }
  endTimeFitModel = proc.time()[3]
  totalTimeFitModel = endTimeFitModel - startTimeFitModel
  
  print(paste0("finished fitting model. Took ", round(totalTimeFitModel / 60, 2), " minutes"))
  
  # # improve the approximation of the posterior if requested by the user
  # if(improveHyperpar) {
  #   browser()
  #   mod = inla.hyperpar(mod)
  # }
  
  # browser()
  
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
  hyperMat = sapply(postSamples, function(x) {x$hyperpar})
  if(family == "normal") {
    if(!separateRanges) {
      mat = apply(hyperMat, 2, function(x) {c(totalVar=exp(x[3])+1/x[1], spatialVar=exp(x[3]), errorVar=1/x[1], 
                                              totalSD=sqrt(exp(x[3])+1/x[1]), spatialSD=sqrt(exp(x[3])), errorSD=sqrt(1/x[1]), 
                                              spatialRange=exp(x[2]), alpha=multivariateExpit(x[4:(3 + nLayer - 1)]))})
      mat = rbind(mat, alpha=1-colSums(mat[8:(7+nLayer-1),]))
      hyperNames = c("totalVar", "spatialVar", "clusterVar", "totalSD", "spatialSD", "clusterSD", "spatialRange", 
                     paste0("alpha", 1:nLayer))
    } else {
      mat = apply(hyperMat, 2, function(x) {c(totalVar=exp(x[2+nLayer])+1/x[1], spatialVar=exp(x[2+nLayer]), errorVar=1/x[1], 
                                              totalSD=sqrt(exp(x[2+nLayer])+1/x[1]), spatialSD=sqrt(exp(x[2+nLayer])), errorSD=sqrt(1/x[1]), 
                                              spatialRange=exp(x[2:(1+nLayer)]), alpha=multivariateExpit(x[(3+nLayer):(3+2*nLayer-2)]))})
      mat = rbind(mat, alpha=1-colSums(matrix(mat[(6+nLayer+1):(6+nLayer+1 + nLayer-2),], nrow=nLayer-1)))
      hyperNames = c("totalVar", "spatialVar", "clusterVar", "totalSD", "spatialSD", "clusterSD", paste0("spatialRange", 1:nLayer), 
                     paste0("alpha", 1:nLayer))
    }
  } else if(family == "binomial") {
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
  varI = which(grepl("clusterVar", hyperNames))
  rangeI = which(grepl("Range", hyperNames))
  alphaI = which(grepl("alpha", hyperNames))
  varSummary=parameterSummaryTable[varI,]
  rangeSummary=parameterSummaryTable[rangeI,]
  alphaSummary=parameterSummaryTable[alphaI,]
  if(family == "normal" || clusterEffect) {
    sdI = which(grepl("clusterSD", hyperNames))
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
  fixedIndices = which(grepl("X", latentVarNames))
  # if(clusterEffect)
  #   clustIndices = grepl("clust", latentVarNames)
  
  ## generate logit predictions (first without cluster effect then add the cluster effect in)
  # for prediction locations
  if(length(xPred) != 0)
    fixedPart = xPred  %*% latentMat[fixedIndices,]
  else
    fixedPart = 0
  predMat = fixedPart + APred %*% latentMat[fieldIndices,]
  
  # for observation locations
  if(length(xObs) != 0)
    fixedPart = xObs  %*% latentMat[fixedIndices,]
  else
    fixedPart = 0
  obsMat = fixedPart + AObs %*% latentMat[fieldIndices,]
  
  # get draws from basis function coefficients
  basisCoefMat = latentMat[fieldIndices,]
  
  # add in cluster effect if necessary
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
  
  # compute basis function coefficient predictions and standard deviations
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
       interceptSummary=interceptSummary, fixedEffectSummary=fixedEffectSummary, rangeSummary=rangeSummary, 
       sdSummary=sdSummary, varSummary=varSummary, overdispersionSummary=overdispersionSummary, parameterSummaryTable=parameterSummaryTable, 
       alphaSummary=alphaSummary, timings=timings, priorPar=priorPar, precomputedNormalizationFun=precomputedNormalizationFun, 
       # the rest of the outputs are saved to be used for spatial aggregations later on
       predMat=predMatClustEffect, obsMat=obsMatClustEffect, hyperMat=hyperMat, clusterVars=clusterVars, rhos=rhos
       )
}

# use the fitSPDE function to fit SPDE model to binomial data within Kenya
fitLKINLAKenyaDat = function(dat=NULL, dataType=c("mort", "ed"), 
                             nu=1.5, seed=1, nLayer=3, NC=14,
                             nBuffer=5, priorPar=NULL, 
                             normalize=TRUE, fastNormalize=TRUE, latInfo=NULL, 
                             intStrategy="ccd", strategy="gaussian", 
                             significanceCI=.8, printVerboseTimings=FALSE, nPostSamples=1000, 
                             urbanEffect=TRUE, clusterEffect=TRUE, predictionType=c("mean", "median"), 
                             initialEffectiveRange=NULL, initialAlphas=rep(1/nLayer, nLayer-1), 
                             effRangeRange=NULL, leaveOutRegionName=NULL, kmres=5, 
                             separateRanges=FALSE, doValidation=FALSE, previousFit=NULL, 
                             precomputedNormalizationFun=NULL, family=c("binomial", "betabinomial"), 
                             leaveOutI=NULL, verbose=TRUE, useUrbanPrior=TRUE) {
  # load observations
  dataType = match.arg(dataType)
  predictionType = match.arg(predictionType)
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
  
  # modify data and prediction information based on left out region
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
    # # make prediction coordinates on a fine grid
    # out = load(paste0("dataPointsKenya.RData"))
    # xRangeDat = dataPointsKenya$xRange
    # yRangeDat = dataPointsKenya$yRange
    # mx = 100
    # my = 100
    # predPts = make.surface.grid(list(x=seq(xRangeDat[1], xRangeDat[2], l=mx), y=seq(yRangeDat[1], yRangeDat[2], l=my)))
    # 
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
  
  if(length(NC) == 1) {
    # if(separateRanges)
    #   NC = c(30, 107) # by default, use two layers with the finest layer having resolution equal to 10km
    if(separateRanges)
      NC = c(30, 214) # by default, use two layers with the finest layer having resolution equal to 5km
  }
  if(separateRanges)
    nLayer = length(NC)
  
  c(fitLKINLAStandard2(obsCoords, obsValues, predPts, nu, seed, nLayer, NC, nBuffer, priorPar, 
                       xObs, xPred, normalize, intStrategy, strategy, fastNormalize, 
                       predictionType, significanceCI, printVerboseTimings, 
                       nPostSamples, family, obsNs, clusterEffect, latInfo, 
                       initialEffectiveRange, initialAlphas, effRangeRange, predClusterI, 
                       separateRanges=separateRanges, doValidation=doValidation, 
                       previousFit=previousFit, precomputedNormalizationFun=precomputedNormalizationFun, 
                       verbose=verbose, useUrbanPrior=useUrbanPrior), 
    list(obsCoords=obsCoords, obsValues=obsValues, xObs=xObs, xPred=xPred, obsNs=obsNs, obsUrban=obsUrban, 
         predPts=predPts, predClusterI=predClusterI, predsUrban=predsUrban, kmres=kmres)
  )
}

# use the fitLKINLAKenyaDat function to validate LKINLA model to binomial data within Kenya with leave one 
# region out validation, and prediction at cluster level
validateLKINLAKenyaDat = function(dat=NULL, dataType=c("mort", "ed"), 
                                  nu=1.5, seed=1, nLayer=3, NC=14,
                                  nBuffer=15, priorPar=NULL, 
                                  normalize=TRUE, fastNormalize=TRUE, latInfo=NULL, 
                                  intStrategy="ccd", strategy="gaussian", 
                                  significanceCI=.8, printVerboseTimings=FALSE, nPostSamples=1000, 
                                  urbanEffect=TRUE, clusterEffect=FALSE, predictionType=c("mean", "median"), 
                                  initialEffectiveRange=NULL, initialAlphas=rep(1/nLayer, nLayer-1), 
                                  effRangeRange=NULL, leaveOutRegionName=NULL, kmres=5, 
                                  separateRanges=FALSE, family=c("binomial", "betabinomial"), 
                                  loadPreviousFit=FALSE, saveResults=TRUE, 
                                  sampleTable=NULL, stratifiedValidation=TRUE, loadPreviousResults=FALSE) {
  if(!is.null(seed))
    set.seed(seed)
  
  if(separateRanges && is.null(latInfo)) {
    nLayer = 2
    NC = c(30, 214)
  }
  
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
  
  # set the type of problem we're working on
  if(dataType == "mort") {
    fileNameRoot = "Mort"
  }
  else {
    fileNameRoot = "Ed"
  }
  
  # first fit the full model (we will use this to initialize the model during the validation fits for each left out county)
  
  fileName = paste0("savedOutput/validation/resultsLKINLA", fileNameRoot, "ValidationFull", "_clust", clusterEffect, 
                    "_urb", urbanEffect, "_sep", separateRanges, ".RData")
  if(!loadPreviousFit || !file.exists(fileName)) {
    print("Fitting full model")
    time = system.time(fit <- fitLKINLAKenyaDat(dat, dataType, nu, seed, nLayer, NC,
                                                nBuffer, priorPar, 
                                                normalize, fastNormalize, latInfo, 
                                                intStrategy, strategy, 
                                                significanceCI, printVerboseTimings, nPostSamples, 
                                                urbanEffect, clusterEffect, predictionType, 
                                                initialEffectiveRange, initialAlphas, 
                                                effRangeRange, leaveOutRegionName, kmres, 
                                                separateRanges, doValidation=TRUE, family=family))
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
      save(time, fit, fullPooledScoresBinomial, fullRuralScoresBinomial, fullUrbanScoresBinomial, 
           fullPooledScores, fullRuralScores, fullUrbanScores, file=fileName)
  }
  else {
    print("Loading previous full model fit")
    load(fileName)
  }
  previousFit = fit
  latInfo = previousFit$latInfo
  precomputedNormalizationFun = previousFit$precomputedNormalizationFun
  
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
  fileName = paste0("savedOutput/validation/resultsLKINLA", fileNameRoot, "ValidationAllTemp", "_clust", clusterEffect, 
                    "_urb", urbanEffect, "_sep", separateRanges, ".RData")
  if(loadPreviousResults && file.exists(fileName)) {
    load(fileName)
    startFrom = i+1
  }
  for(i in startFrom:nFold) {
    if(i > nFold)
      break
    
    if(!stratifiedValidation) {
      thisRegion = regions[i]
      thisSampleI = allRegions == thisRegion
      print(paste0("Validating LKINLA model with urban=", urbanEffect, ", cluster=", clusterEffect, 
                   ", region=", thisRegion, " (", i, "/", length(regions), " regions)"))
      leaveOutI = NULL
    } else {
      thisRegion = NULL
      print(paste0("Validating LKINLA model with urban=", urbanEffect, ", cluster=", clusterEffect, 
                   ", (", i, "/", nFold, " folds)"))
      leaveOutI = sampleTable[,i]
      thisSampleI = leaveOutI
    }
    
    
    # fit the model
    time = system.time(fit <- fitLKINLAKenyaDat(dat, dataType, nu, seed, nLayer, NC,
                                                nBuffer, priorPar, 
                                                normalize, fastNormalize, latInfo, 
                                                intStrategy, strategy, 
                                                significanceCI, printVerboseTimings, nPostSamples, 
                                                urbanEffect, clusterEffect, predictionType, 
                                                initialEffectiveRange, initialAlphas, 
                                                effRangeRange, thisRegion, kmres, 
                                                separateRanges, previousFit=previousFit, doValidation=TRUE, 
                                                precomputedNormalizationFun=precomputedNormalizationFun, 
                                                leaveOutI=leaveOutI, family=family))
    
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

# fit standard LK model given data and relevant parameters
# obsCoords: spatial coordinates of the data
# obsValues: observations at the coordinates
# predCoords: places at which to make predictions
# nu: Matern smoothness parameter from "simple" LKrig model
# xObs: observation design matrix (intercept and covariates at obs locations)
# xPred: prediction design matrix (intercept and covariates at pred locations)
# NC: number of coarse lattice points over the longest dimension of the data
# int.strategy: "auto" or "eb" for empirical Bayes
fitLKINLAStandard = function(obsCoords, obsValues, predCoords=obsCoords, nu=1.5, seed=1, nLayer=3, NC=5,
                              nBuffer=5, priorPar=getPrior(.1, .1, 10), 
                              xObs=cbind(1, obsCoords), xPred=cbind(1, predCoords), normalize=TRUE, 
                              intStrategy="grid", strategy="laplace", fastNormalize=TRUE, 
                              predictionType=c("mean", "median"), significanceCI=0.8, 
                              printVerboseTimings=FALSE, nPostSamples=1000) {
  set.seed(seed)
  
  # get the type of prediction the user wants
  predictionType = match.arg(predictionType)
  
  # get lattice points, prediction points
  xRangeDat = range(obsCoords[,1])
  yRangeDat = range(obsCoords[,2])
  latInfo = makeLatGrids(xRangeDat, yRangeDat, NC, nBuffer, nLayer)
  nx = latInfo[[1]]$nx
  ny = latInfo[[1]]$ny
  
  # generate lattice basis matrix
  AObs = makeA(obsCoords, latInfo)
  
  # run matrix precomputations
  precomputedMatrices = precomputationsQ2(latInfo)
  
  # define the LK model
  n = length(obsValues)
  rgen = inla.rgeneric.define(model=inla.rgeneric.lk.model.standard, latInfo=latInfo, ys=obsValues, 
                              prior=priorPar, normalize=normalize, precomputedMatrices=precomputedMatrices, 
                              X=xObs, nu=nu, datCoords=obsCoords, fastNormalize=fastNormalize, 
                              printVerboseTimings=printVerboseTimings)
  # use these global variables for testing calls to inla.rgeneric.lk.model.simple
  # latInfo<<-latInfo; ys<<-obsValues;
  # prior<<-priorPar; normalize<<-normalize; precomputedMatrices<<-precomputedMatrices;
  # X<<-xObs; nu<<-nu; datCoords<<-obsCoords; fastNormalize<<-fastNormalize;
  # printVerboseTimings<<-printVerboseTimings
  ## generate inla stack:
  # Stacked A matrix (A_s from notation of LR2015 Bayesian Spatial Modelling with R_INLA):
  # (AEst   0  )
  # ( 0   APred)
  # eta_s = (c^T c^T)^T
  # where c is the vector of lattice coefficients
  AEst = makeA(obsCoords, latInfo)
  APred = makeA(predCoords, latInfo)
  latticeInds = 1:ncol(AEst)
  stack.est = inla.stack(A =list(AEst, 1), 
                         effects =list(field=latticeInds, X=xObs), 
                         data =list(y=obsValues), 
                         tag ="est", 
                         remove.unused=FALSE)
  stack.pred = inla.stack(A =list(APred[1:2,], 1), 
                          effects =list(field=latticeInds, X=xPred[1:2,]), 
                          data =list(y=NA), 
                          tag ="pred", 
                          remove.unused=FALSE)
  stack.full = inla.stack(stack.est, stack.pred, 
                          remove.unused=FALSE)
  dat = inla.stack.data(stack.full, rgen=rgen, remove.unused=FALSE)
  
  # fit the model
  # control.inla = list(cmin = 0, int.strategy=int.strategy) 
  # see: inla.doc("loggamma")
  # shape=.1, scale=10 for unit mean, variance 100 prior
  controls = list(strategy=strategy, int.strategy=intStrategy) 
  allQuantiles = c(0.5, (1-significanceCI) / 2, 1 - (1-significanceCI) / 2)
  mod = inla(y ~ - 1 + X + f(field, model=rgen), 
             data=dat, quantiles=allQuantiles, family="normal", verbose=TRUE, 
             control.inla=controls, 
             control.compute=list(config=TRUE), 
             control.predictor=list(A=inla.stack.A(stack.full), compute=TRUE, quantiles=allQuantiles), 
             control.fixed=list(quantiles=allQuantiles), 
             control.family=list(hyper = list(prec = list(prior="loggamma", param=c(0.1,0.1))))
  )
  
  # get predictive surface, SD, and data
  # index = inla.stack.index(stack.full, "pred")$data
  obsInds = 1:n
  # predInds = (n+1):(n+nrow(predCoords))
  if(predictionType == "mean") {
    linpreds = mod[["summary.linear.predictor"]]$mean
    # preds = linpreds[predInds]
    obsPreds = linpreds[obsInds]
  } else {
    linpreds = mod[["summary.linear.predictor"]]$`0.5quant`
    # preds = linpreds[predInds]
    obsPreds = linpreds[obsInds]
  }
  linpred.sd = mod[["summary.linear.predictor"]]$sd
  # predSDs = linpred.sd[predInds]
  obsSDs = linpred.sd[obsInds]
  
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
  predMat = fixedPart + APred %*% latentMat[fieldIndices,]
  
  # compute predictive credible intervals
  preds = rowMeans(predMat)
  predSDs = apply(predMat, 1, sd)
  lower = apply(predMat, 1, quantile, probs=(1-significanceCI)/2)
  medians = apply(predMat, 1, median)
  upper = apply(predMat, 1, quantile, probs=1-(1-significanceCI)/2)
  
  if(!is.null(xObs) && all(xObs[,1]==1))
    interceptSummary=mod$summary.fixed[,c(1, 2, 4, 3, 5)]
  else
    interceptSummary = matrix(rep(0, 5), nrow=1)
  
  # get posterior hyperparameter samples and transform them as necessary
  # interpretation of hyperparameters: 
  # 1: error precision
  # 2: log effective range
  # 3: log spatial variance
  # 4-(3 + nLayer - 1): multivariateLogit alpha
  hyperMat = sapply(postSamples, function(x) {x$hyperpar})
  mat = apply(hyperMat, 2, function(x) {c(totalVar=exp(x[3])+1/x[1], spatialVar=exp(x[3]), errorVar=1/x[1], 
                                          totalSD=sqrt(exp(x[3])+1/x[1]), spatialSD=sqrt(exp(x[3])), errorSD=sqrt(1/x[1]), 
                                          spatialRange=exp(x[2]), alpha=multivariateExpit(x[4:(3 + nLayer - 1)]))})
  mat = rbind(mat, alpha=1-colSums(mat[8:(7+nLayer-1),]))
  hyperNames = c("totalVar", "spatialVar", "errorVar", "totalSD", "spatialSD", "errorSD", "spatialRange", 
                 paste0("alpha", 1:nLayer))
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
  sdSummary=parameterSummaryTable[6,]
  varSummary=parameterSummaryTable[3,]
  rangeSummary=parameterSummaryTable[7,]
  alphaSummary=parameterSummaryTable[8:(8+nLayer-1),]
  
  # compute basis function coefficient predictions and standard deviations
  startI = n+2+1
  endI = n+2+nx*ny
  coefPreds = list(layer1 = linpreds[startI:endI])
  coefSDs = list(layer1 = linpred.sd[startI:endI])
  if(nLayer >= 2) {
    for(i in 2:nLayer) {
      startI = endI + 1
      endI = startI + nrow(latInfo[[i]]$latCoords) - 1
      coefPreds = c(coefPreds, list(linpreds[startI:endI]))
      coefSDs = c(coefSDs, list(linpred.sd[startI:endI]))
    }
  }
  
  list(mod=mod, preds=preds, sigmas=predSDs, latInfo=latInfo, obsPreds=obsPreds, 
       obsSDs=obsSDs, coefPreds=coefPreds, coefSDs=coefSDs, 
       interceptSummary=interceptSummary, rangeSummary=rangeSummary, 
       sdSummary=sdSummary, varSummary=varSummary, parameterSummaryTable=parameterSummaryTable, 
       alphaSummary=alphaSummary)
}

# this function generates results for the simulation study for the LKINLA (standard) model
# input arguments:
#   argument specifying the dataset type
resultsLKINLA = function(randomSeeds=NULL, covType=c("exponential", "matern", "mixture"), rangeText=c("01", "05", "1", "mix"), 
                         maxDataSets=NULL, NC=5, nLayer=3, normalize=TRUE, nBuffer=5) {
  
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
  
  # generate results for all data sets and return results (TODO: otherVariableNames)
  resultsLKINLA = fitModelToDataSets(fitLKINLAStandard, dataSets, randomSeeds=randomSeeds, 
                                     otherArgs=list(NC=NC, nLayer=nLayer, normalize=normalize, nBuffer=nBuffer), 
                                     maxDataSets=maxDataSets)
  
  # save results
  fileName = paste0("resultsLKINLA_cov", covType, "Range", rangeText, "NC", NC, "nLayer", nLayer, 
                    "maxDataSets", maxDataSets, ".RData")
  save(resultsLKINLA, file=fileName)
  
  resultsLKINLA
}









