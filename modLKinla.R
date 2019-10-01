# fit standard LK model given data and relevant parameters
# coords: spatial coordinates of the data
# ys: observations at the coordinates
# predPts: places at which to make predictions
# nu: Matern smoothness parameter from "simple" LKrig model
# X: observation design matrix (intercept and covariates at obs locations)
# XPred: prediction design matrix (intercept and covariates at pred locations)
# NC: number of coarse lattice points over the longest dimension of the data
# int.strategy: "auto" or "eb" for empirical Bayes
fitLKINLAStandard = function(coords, ys, predPts=coords, nu=1.5, seed=1, nLayer=3, NC=5,
                             nBuffer=5, priorPar=getPrior(.1, .1, 10), 
                             X=cbind(1, coords), XPred=cbind(1, predPts), normalize=TRUE, 
                             intStrategy="auto", strategy="gaussian", fastNormalize=FALSE, 
                             predictionType=c("median", "mean"), printVerboseTimings=FALSE) {
  set.seed(seed)
  
  # get the type of prediction the user wants
  predictionType = match.arg(predictionType)
  
  # get lattice points, prediction points
  xRangeDat = range(coords[,1])
  yRangeDat = range(coords[,2])
  latInfo = makeLatGrids(xRangeDat, yRangeDat, NC, nBuffer, nLayer)
  nx = latInfo[[1]]$nx
  ny = latInfo[[1]]$ny
  
  # generate lattice basis matrix
  AObs = makeA(coords, latInfo)
  
  # run matrix precomputations
  precomputedMatrices = precomputationsQ2(latInfo)
  
  # define the LK model
  n = length(ys)
  rgen = inla.rgeneric.define(model=inla.rgeneric.lk.model.standard, latInfo=latInfo, ys=ys, 
                              prior=priorPar, normalize=normalize, precomputedMatrices=precomputedMatrices, 
                              X=X, nu=nu, datCoords=coords, fastNormalize=fastNormalize, 
                              printVerboseTimings=printVerboseTimings)
  # use these global variables for testing calls to inla.rgeneric.lk.model.simple
  # latInfo<<-latInfo; ys<<-ys;
  # prior<<-priorPar; normalize<<-normalize; precomputedMatrices<<-precomputedMatrices;
  # X<<-X; nu<<-nu; datCoords<<-coords; fastNormalize<<-fastNormalize;
  # printVerboseTimings<<-printVerboseTimings
  ## generate inla stack:
  # Stacked A matrix (A_s from notation of LR2015 Bayesian Spatial Modelling with R_INLA):
  # (AEst   0  )
  # ( 0   APred)
  # eta_s = (c^T c^T)^T
  # where c is the vector of lattice coefficients
  AEst = makeA(coords, latInfo)
  APred = makeA(predPts, latInfo)
  latticeInds = 1:ncol(AEst)
  stack.est = inla.stack(A =list(AEst, 1), 
                         effects =list(field=latticeInds, X=X), 
                         data =list(y=ys), 
                         tag ="est", 
                         remove.unused=FALSE)
  stack.pred = inla.stack(A =list(APred, 1), 
                          effects =list(field=latticeInds, X=XPred), 
                          data =list(y=NA), 
                          tag ="pred", 
                          remove.unused=FALSE)
  stack.full = inla.stack(stack.est, stack.pred, 
                          remove.unused=FALSE)
  dat = inla.stack.data(stack.full, rgen=rgen, remove.unused=FALSE)
  
  # fit the model
  # control.inla = list(cmin = 0, int.strategy=int.strategy) 
  # see: inla.doc("loggamma")
  controls = list(strategy=strategy, int.strategy=intStrategy) 
  mod = inla(y ~ - 1 + X + f(field, model=rgen), data=dat, 
             control.predictor=list(A=inla.stack.A(stack.full), compute=TRUE), 
             family="normal", verbose=TRUE, control.inla=controls, 
             control.family=list(hyper = list(prec = list(prior="loggamma", param=c(0.1,0.1))))
  )
  
  # get predictive surface, SD, and data
  index = inla.stack.index(stack.full, "pred")$data
  obsInds = 1:n
  predInds = (n+1):(n+nrow(predPts))
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
  
  startI = n+nrow(predPts)+1
  endI = n+nrow(predPts)+nx*ny
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
  
  # get true effective range and marginal variance:
  latticeWidth = latInfo[[1]]$latWidth
  
  list(mod=mod, preds=preds, SDs=predSDs, latInfo=latInfo, latWidth=latticeWidth, obsPreds=obsPreds, 
       obsSDs=obsSDs, coefPreds=coefPreds, coefSDs=coefSDs)
}