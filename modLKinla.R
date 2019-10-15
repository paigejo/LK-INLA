# fit standard LK model given data and relevant parameters
# obsCoords: spatial coordinates of the data
# obsValues: observations at the coordinates
# predCoords: places at which to make predictions
# nu: Matern smoothness parameter from "simple" LKrig model
# xObs: observation design matrix (intercept and covariates at obs locations)
# xPred: prediction design matrix (intercept and covariates at pred locations)
# NC: number of coarse lattice points over the longest dimension of the data
# int.strategy: "auto" or "eb" for empirical Bayes
fitLKINLAStandard2 = function(obsCoords, obsValues, predCoords=obsCoords, nu=1.5, seed=1, nLayer=3, NC=5,
                             nBuffer=5, priorPar=getPrior(.1, .1, 10), 
                             xObs=cbind(1, obsCoords), xPred=cbind(1, predCoords), normalize=TRUE, 
                             intStrategy="grid", strategy="laplace", fastNormalize=TRUE, 
                             predictionType=c("mean", "median"), significanceCI=0.8, 
                             printVerboseTimings=FALSE, nPostSamples=10000, family=c("normal", "binomial"),
                             obsNs=rep(1, length(obsValues)), clusterEffect=TRUE, improveHyperpar=TRUE) {
  set.seed(seed)
  
  # get the type of prediction the user wants
  predictionType = match.arg(predictionType)
  family = match.arg(family)
  
  # check to make sure inputs make sense
  if(!clusterEffect && family == "normal")
    stop("cluster effect must be included for the normal family")
  
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
  clust = 1:length(obsValues)
  if(family == "normal") {
    stack.est = inla.stack(A =list(AEst, 1), 
                           effects =list(field=latticeInds, X=xObs), 
                           data =list(y=obsValues, link=1, Ntrials = obsNs), 
                           tag ="est", 
                           remove.unused=FALSE)
    stack.pred = inla.stack(A =list(APred, 1), 
                            effects =list(field=latticeInds, X=xPred), 
                            data =list(y=NA, link=1), 
                            tag ="pred", 
                            remove.unused=FALSE)
  }
  else if(family == "binomial") {
    if(clusterEffect) {
      stack.est = inla.stack(A =list(AEst, 1, 1), 
                             effects =list(field=latticeInds, clust=clust, X=xObs), 
                             data =list(y=obsValues, link=1, Ntrials = obsNs), 
                             tag ="est", 
                             remove.unused=FALSE)
    } else {
      stack.est = inla.stack(A =list(AEst, 1), 
                             effects =list(field=latticeInds, X=xObs), 
                             data =list(y=obsValues, link=1, Ntrials = obsNs), 
                             tag ="est", 
                             remove.unused=FALSE)
    }
    stack.pred = inla.stack(A =list(APred, 1), 
                            effects =list(field=latticeInds, X=xPred), 
                            data =list(y=NA, link=1, Ntrials=rep(1, nrow(APred)), 
                            tag ="pred"), 
                            remove.unused=FALSE)
  }
  stack.full = inla.stack(stack.est, stack.pred, 
                          remove.unused=FALSE)
  dat = inla.stack.data(stack.full, rgen=rgen, remove.unused=FALSE)
  
  # fit the model
  # control.inla = list(cmin = 0, int.strategy=int.strategy) 
  # see: inla.doc("loggamma")
  # shape=.1, scale=10 for unit mean, variance 100 prior
  controls = list(strategy=strategy, int.strategy=intStrategy) 
  allQuantiles = c(0.5, (1-significanceCI) / 2, 1 - (1-significanceCI) / 2)
  if(family == "normal") {
    mod = inla(y ~ - 1 + X + f(field, model=rgen), 
               data=dat, quantiles=allQuantiles, family=family, verbose=TRUE, 
               control.inla=controls, 
               control.compute=list(config=TRUE), 
               control.predictor=list(A=inla.stack.A(stack.full), compute=TRUE, quantiles=allQuantiles), 
               control.fixed=list(quantiles=allQuantiles), 
               control.family=list(hyper = list(prec = list(prior="loggamma", param=c(0.1,0.1)))))
  }
  else if(family == "binomial") {
    if(clusterEffect) {
      clusterList = list(param=c(.15, 0.01), prior="pc.prec")
      mod = inla(y ~ - 1 + X + f(field, model=rgen) + 
                   f(clust, model="iid", hyper = list(prec = clusterList)), 
                 data=dat, quantiles=allQuantiles, family=family, verbose=TRUE, 
                 control.inla=controls, Ntrials=dat$Ntrials, 
                 control.compute=list(config=TRUE), 
                 control.predictor=list(A=inla.stack.A(stack.full), compute=TRUE, quantiles=allQuantiles), 
                 control.fixed=list(quantiles=allQuantiles))
    } else {
      mod = inla(y ~ - 1 + X + f(field, model=rgen), 
                 data=dat, quantiles=allQuantiles, family=family, verbose=TRUE, 
                 control.inla=controls, Ntrials=dat$Ntrials, 
                 control.compute=list(config=TRUE), 
                 control.predictor=list(A=inla.stack.A(stack.full), compute=TRUE, quantiles=allQuantiles), 
                 control.fixed=list(quantiles=allQuantiles))
    }
  }
  
  # improve the approximation of the posterior if requested by the user
  if(improveHyperpar) {
    test = inla.hyperpar(mod)
    browser()
    # mod = inla.hyperpar(mod)
  }
  
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
  postSamples = inla.posterior.sample(nPostSamples, mod)
  
  # get posterior hyperparameter samples and transform them as necessary
  # interpretation of hyperparameters: 
  # 1: error precision
  # 2: log effective range
  # 3: log spatial variance
  # 4-(3 + nLayer - 1): multivariateLogit alpha
  hyperMat = sapply(postSamples, function(x) {x$hyperpar})
  if(family == "normal") {
    mat = apply(hyperMat, 2, function(x) {c(totalVar=exp(x[3])+1/x[1], spatialVar=exp(x[3]), errorVar=1/x[1], 
                                            totalSD=sqrt(exp(x[3])+1/x[1]), spatialSD=sqrt(exp(x[3])), errorSD=sqrt(1/x[1]), 
                                            spatialRange=exp(x[2]), alpha=multivariateExpit(x[4:(3 + nLayer - 1)]))})
    mat = rbind(mat, alpha=1-colSums(mat[8:(7+nLayer-1),]))
    hyperNames = c("totalVar", "spatialVar", "errorVar", "totalSD", "spatialSD", "errorSD", "spatialRange", 
                   paste0("alpha", 1:nLayer))
  } else {
    if(clusterEffect) {
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
      logSpatialRangeI = 1
      logSpatialVarI = 2
      logitAlphaI = 3:(2 + nLayer - 1)
      mat = apply(hyperMat, 2, function(x) {c(spatialVar=exp(x[logSpatialVarI]), spatialSD=sqrt(exp(x[logSpatialVarI])), 
                                              spatialRange=exp(x[logSpatialRangeI]), alpha=multivariateExpit(x[logitAlphaI]))})
      mat = rbind(mat, alpha=1-colSums(mat[4:(3+nLayer-1),]))
      hyperNames = c("spatialVar", "spatialSD", "spatialRange", paste0("alpha", 1:nLayer))
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
  sdI = which(grepl("clusterSD", hyperNames))
  rangeI = which(grepl("Range", hyperNames))
  alphaI = which(grepl("alpha", hyperNames))
  sdSummary=parameterSummaryTable[sdI,]
  varSummary=parameterSummaryTable[varI,]
  rangeSummary=parameterSummaryTable[rangeI,]
  alphaSummary=parameterSummaryTable[alphaI,]
  
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
  
  # add in cluster effect if necessary
  if((family == "binomial" && clusterEffect) || family == "normal") {
    # get cluster effect variance
    clusterVarI = which(grepl("clusterVar", hyperNames))
    clusterVars = mat[clusterVarI,]
    predMatClustEffect = predMat + matrix(rnorm(length(predMat), sd=rep(clusterVars, each=nrow(predMat))), nrow=nrow(predMat))
    obsMatClustEffect = obsMat + matrix(rnorm(length(obsMat), sd=rep(clusterVars, each=nrow(obsMat))), nrow=nrow(obsMat))
  } else {
    clusterVars = NULL
    predMatClustEffect = NULL
    obsMatClustEffect = NULL
  }
  
  # compute predictive credible intervals
  if(is.null(predMatClustEffect)) {
    preds = rowMeans(expit(predMat))
    predSDs = apply(expit(predMat), 1, sd)
    lowerPreds = apply(expit(predMat), 1, quantile, probs=(1-significanceCI)/2)
    medianPreds = apply(expit(predMat), 1, median)
    upperPreds = apply(expit(predMat), 1, quantile, probs=1-(1-significanceCI)/2)
    obsPreds = rowMeans(expit(obsMat))
    obsSDs = apply(expit(obsMat), 1, sd)
    lowerObs = apply(expit(obsMat), 1, quantile, probs=(1-significanceCI)/2)
    medianObs = apply(expit(obsMat), 1, median)
    upperObs = apply(expit(obsMat), 1, quantile, probs=1-(1-significanceCI)/2)
  } else {
    preds = rowMeans(expit(predMatClustEffect))
    predSDs = apply(expit(predMatClustEffect), 1, sd)
    lowerPreds = apply(expit(predMatClustEffect), 1, quantile, probs=(1-significanceCI)/2)
    medianPreds = apply(expit(predMatClustEffect), 1, median)
    upperPreds = apply(expit(predMatClustEffect), 1, quantile, probs=1-(1-significanceCI)/2)
    obsPreds = rowMeans(expit(obsMatClustEffect))
    obsSDs = apply(expit(obsMatClustEffect), 1, sd)
    lowerObs = apply(expit(obsMatClustEffect), 1, quantile, probs=(1-significanceCI)/2)
    medianObs = apply(expit(obsMatClustEffect), 1, median)
    upperObs = apply(expit(obsMatClustEffect), 1, quantile, probs=1-(1-significanceCI)/2)
  }
  
  interceptSummary=mod$summary.fixed[1,c(1, 2, 4, 3, 5)]
  
  # compute basis function coefficient predictions and standard deviations
  startI = n+nrow(predCoords)+1
  endI = n+nrow(predCoords)+nx*ny
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
  
  list(preds=preds, sigmas=predSDs, lower=lowerPreds, median=medianPreds, upper=upperPreds, 
       obsPreds=obsPreds, obsSDs=obsSDs, obsLower=lowerObs, medianObs=medianObs, upperObs=medianObs, 
       mod=mod, latInfo=latInfo, coefPreds=coefPreds, coefSDs=coefSDs, 
       interceptSummary=interceptSummary, rangeSummary=rangeSummary, 
       sdSummary=sdSummary, varSummary=varSummary, parameterSummaryTable=parameterSummaryTable, 
       alphaSummary=alphaSummary, 
       # the rest of the outputs are saved to be used for spatial aggregations later on
       noClustPredMat=predMat, hyperMat=hyperMat, clusterVars=clusterVars
       )
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
                              printVerboseTimings=FALSE, nPostSamples=10000) {
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
  
  interceptSummary=mod$summary.fixed[,c(1, 2, 4, 3, 5)]
  
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









