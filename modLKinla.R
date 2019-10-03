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
                             intStrategy="auto", strategy="gaussian", fastNormalize=FALSE, 
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
  rgen = inla.rgeneric.define(model=inla.rgeneric.lk.model.standard, latInfo=latInfo, obsValues=obsValues, 
                              prior=priorPar, normalize=normalize, precomputedMatrices=precomputedMatrices, 
                              X=xObs, nu=nu, datCoords=obsCoords, fastNormalize=fastNormalize, 
                              printVerboseTimings=printVerboseTimings)
  # use these global variables for testing calls to inla.rgeneric.lk.model.simple
  # latInfo<<-latInfo; obsValues<<-obsValues;
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
  stack.pred = inla.stack(A =list(APred, 1), 
                          effects =list(field=latticeInds, X=xPred), 
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
  mod = inla(y ~ - 1 + X + f(field, model=rgen), data=dat, quantiles=allQuantiles, 
             control.predictor=list(A=inla.stack.A(stack.full), compute=TRUE, quantiles=allQuantiles), 
             family="normal", verbose=TRUE, control.inla=controls, control.fixed=list(quantiles=allQuantiles), 
             control.family=list(hyper = list(prec = list(prior="loggamma", param=c(0.1,0.1))))
  )
  
  # get predictive surface, SD, and data
  index = inla.stack.index(stack.full, "pred")$data
  obsInds = 1:n
  predInds = (n+1):(n+nrow(predCoords))
  if(predictionType == "mean") {
    linpreds = mod[["summary.linear.predictor"]]$mean
    preds = linpreds[predInds]
    obsPreds = linpreds[obsInds]
  } else {
    linpreds = mod[["summary.linear.predictor"]]$`0.5quant`
    preds = linpreds[predInds]
    obsPreds = linpreds[obsInds]
  }
  linpred.sd = mod[["summary.linear.predictor"]]$sd
  predSDs = linpred.sd[predInds]
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
  lower = apply(predMat, 1, quantile, probs=(1-significanceCI)/2)
  medians = apply(predMat, 1, median)
  upper = apply(predMat, 1, quantile, probs=1-(1-significanceCI)/2)
  
  interceptSummary=mod$summary.fixed[,1:5]
  rangeSummary=mod$summary.random[2,1:5]
  spatialSDSummary = mod$summary.random[3,1:5]
  
  # get posterior hyperparameter samples and transform them as necessary
  hyperMat = sapply(postSamples, function(x) {x$hyperpar})
  mat = apply(hyperMat, 2, function(x) {c(totalVar=x[3]^2+1/x[1], spatialVar=x[3]^2, errorVar=1/x[1], 
                                          totalSD=sqrt(x[3]^2+1/x[1]), spatialSD=x[3], errorSD=sqrt(1/x[1]), 
                                          spatialRange=x[2])})
  hyperNames = c("totalVar", "spatialVar", "errorVar", "totalSD", "spatialSD", "errorSD", "spatialRange")
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
  
  list(mod=mod, preds=preds, SDs=predSDs, latInfo=latInfo, latWidth=latticeWidth, obsPreds=obsPreds, 
       obsSDs=obsSDs, coefPreds=coefPreds, coefSDs=coefSDs, 
       interceptSummary=interceptSummary, rangeSummary=rangeSummary, 
       sdSummary=sdSummary, varSummary=varSummary, parameterSummaryTable=parameterSummaryTable)
}

# this function generates results for the simulation study for the LKINLA (standard) model
# input arguments:
#   argument specifying the dataset type
resultsLKINLA = function(randomSeeds=NULL, covType=c("exponential", "matern", "mixture"), rangeText=c("01", "05", "1", ""), 
                         maxDataSets=NULL, NC=5, nLayer=3, normalize=TRUE, nBuffer=5) {
  
  # determine the type of covariance for the data set
  covType = match.arg(covType)
  
  # determine the spatial range for the data set. No range text means it's a mixture
  rangeText = match.arg(rangeText)
  
  # construct the file name for the desired data set and load it
  dataText = paste0(covType, rangeText, "DataSet.RData")
  out = load(dataText)
  dataSets = simulationData
  
  # generate results for all data sets and return results (TODO: otherVariableNames)
  fitModelToDataSets(fitLKINLAStandard, dataSets, randomSeeds=randomSeeds, 
                     otherArgs=list(NC=NC, nLayer=nLayer, normalize=normalize, nBuffer=nBuffer), 
                     maxDataSets=maxDataSets)
}









