# given a dataset, fit LatticeKrig with package
# coords and obs: coordinates (nx2 matrix) and observed values (n-vector)
# nDat: resolution of coarsest lattice in largest direction over the domain 
# of the data (including the buffer there are more lattice points though)
# nLayer: number of layers to include
# NOTE: includes a first order linear polynomial in spatial coordinates, so do not include that in the covariates
fitLK = function(coords, obs, predCoords=coords, XObs=NULL, XPred=NULL, NC=5, nLayer=3, simpleMod=TRUE, normalize=FALSE, 
                 nBuffer=5, nu=1.5, verbose=TRUE, lambdaStart=.1, a.wghtStart=5, maxit=15, doSEs=TRUE, significanceCI=.8) {
  # if(!simpleMod)
  #   stop("non simple models not yet supported for LatticeKrig prediction function")
  
  # set up (the lattice, the arguments to LatticeKrig)
  LKinfo = LKrigSetup(coords, nlevel=nLayer, nu=nu, NC=NC, normalize=normalize, 
                      lambda=lambdaStart, a.wght=a.wghtStart)
  
  # if are missing the predictive covariates, predict them using the observation covariates
  if(!is.null(XObs) && is.null(XPred)) {
    if(!is.matrix(XObs))
      XObs = as.matrix(XObs)
    XPred = matrix(nrow=nrow(predCoords), ncol=ncol(XObs))
    
    # make LatticeKrig predictive model for each covariate
    print(paste0("predicting covariates..."))
    for(i in 1:ncol(XObs)) {
      print(paste0("predicting covariate ", i, "/", ncol(XObs)))
      thisCov = XObs[,i]
      out = fitLK(coords, thisCov, NC=NC, nLayer=nLayer, simpleMod=simpleMod, normalize=normalize, 
                  nBuffer=nBuffer, nu=nu, verbose=verbose, lambdaStart=lambdaStart, 
                  a.wghtStart=a.wghtStart, maxit=maxit, doSEs=FALSE)
      XPred[,i] = out$preds
    }
  }
  
  # Maximum Likelihood
  if(verbose)
    print("Beginning LatticeKrig MLE fit...")
  # LKMLE = LKrig.MLE(coords, obs, LKinfo=LKinfo)
  if(is.null(XObs) || is.null(XPred))
    LKMLE = LKrigFindLambdaAwght(coords, obs, LKinfo=LKinfo, verbose=verbose, maxit=maxit)
  else
    LKMLE = LKrigFindLambdaAwght(coords, obs, LKinfo=LKinfo, verbose=verbose, maxit=maxit, Z=XObs)
  if(verbose)
    print(LKMLE$summary)
  
  # final fit
  if(is.null(XObs) || is.null(XPred)) {
    out = LKrig(coords, obs, LKinfo=LKMLE$LKinfo)
    preds = predict.LKrig(out, predCoords)
    if(doSEs)
      predSEs = predictSE.LKrig(out, predCoords)
    else
      predSEs = NULL
  }
  else {
    out = LKrig(coords, obs, LKinfo=LKMLE$LKinfo, Z=XObs)
    preds = predict.LKrig(out, predCoords, Znew=XPred)
    if(doSEs)
      predSEs = predictSE.LKrig(out, predCoords, Znew=XPred)
    else
      predSEs=NULL
  }
  
  # calculate confidence intervals
  lower = preds + qnorm((1 - significanceCI) / 2, sd=predSEs)
  upper = preds + qnorm(1 - (1 - significanceCI) / 2, sd=predSEs)
  
  ## preds
  ## sigmas
  ## lower
  ## upper
  ## interceptSummary
  ## rangeSummary
  ## sdSummary
  ## varSummary
  return(list(mod=out, preds=preds, sigmas=predSEs, lower=lower, upper=upper, LKinfo=LKinfo))
}

# given a dataset, fit LatticeKrig optimizing over alpha and other parameters 
# coords and obs: coordinates (nx2 matrix) and observed values (n-vector)
# nDat: resolution of coarsest lattice in largest direction over the domain 
# of the data (including the buffer there are more lattice points though)
# nLayer: number of layers to include
# NOTE: includes a first order linear polynomial in spatial coordinates, so do not include that in the covariates
# NOTE: lambda is sigma^2 / rho
fitLKStandard = function(coords, obs, predCoords=coords, XObs=NULL, XPred=NULL, NC=5, nLayer=3, normalize=TRUE, 
                 nBuffer=5, nu=1.5, verbose=TRUE, lambdaStart=.1, a.wghtStart=5, doSEs=TRUE, significanceCI=.8, 
                 lowerBoundLogLambda =-16,
                 upperBoundLogLambda = 4,
                 lowerBoundLogitAlpha = rep(-10, nLayer-1),
                 upperBoundLogitAlpha= rep(10, nLayer-1),
                 lowerBoundOmega = -3,
                 upperBoundOmega = .75,
                 factr=1e7,
                 pgtol=1e-1,
                 maxit=15, 
                 nsimConditional=100, 
                 fixedFunctionArgs = list(m = 1)) {
  # if(!simpleMod)
  #   stop("non simple models not yet supported for LatticeKrig prediction function")
  
  # if are missing the predictive covariates, predict them using the observation covariates
  if(!is.null(XObs) && is.null(XPred)) {
    if(!is.matrix(XObs))
      XObs = as.matrix(XObs)
    XPred = matrix(nrow=nrow(predCoords), ncol=ncol(XObs))
    
    # make LatticeKrig predictive model for each covariate
    print(paste0("predicting covariates..."))
    for(i in 1:ncol(XObs)) {
      print(paste0("predicting covariate ", i, "/", ncol(XObs)))
      thisCov = XObs[,i]
      out = fitLKStandard(coords, thisCov, NC=NC, nLayer=nLayer, simpleMod=simpleMod, normalize=normalize, 
                  nBuffer=nBuffer, nu=nu, verbose=verbose, lambdaStart=lambdaStart, 
                  a.wghtStart=a.wghtStart, maxit=maxit, doSEs=FALSE)
      XPred[,i] = out$preds
    }
  }
  
  # do initial latticeKrig fit
  # set up the lattice, the arguments to LatticeKrig
  LKinfoStart = LKrigSetup(coords, nlevel=nLayer, nu=1, NC=NC, normalize=normalize, 
                           lambda=lambdaStart, a.wght=a.wghtStart, fixedFunctionArgs=fixedFunctionArgs)
  
  # make a function to convert from a vector of parameters to a corresponding set of different, named parameters
  getParameters = function(parameters) {
    # omega =  log( a.wght -4)/2
    # transform from optimized parameters to probabilities summing to 1 to get alphas
    if(nLayer != 1) {
      alphas = multivariateExpit(parameters[1:(nLayer-1)])
      alphas = c(alphas, 1 - sum(alphas))
      log.lambda = parameters[nLayer-1 + 1]
      omega = parameters[nLayer-1 + 2]
      list(alphas=alphas, log.lambda=log.lambda, lambda=exp(log.lambda), omega=omega, a.wght=omega2Awght(omega, LKinfoStart))
    }
    else {
      log.lambda = parameters[1]
      omega = parameters[2]
      list(alphas=1, log.lambda=log.lambda, lambda=exp(log.lambda), omega=omega, a.wght=omega2Awght(omega, LKinfoStart))
    }
  }
  
  # make a wrapper function around LKrig in order to optimize overall parameters including alpha
  outerFun = function(parameters, thisVerbose=verbose) {
    parameterList = getParameters(parameters)
    alphas = parameterList$alphas
    log.lambda = parameterList$log.lambda
    lambda = parameterList$lambda
    omega = parameterList$omega
    a.wght = parameterList$a.wght
    
    # set up the lattice, the arguments to LatticeKrig
    LKinfo = LKrigSetup(coords, nlevel=nLayer, nu=NULL, NC=NC, normalize=normalize, 
                        lambda=lambda, a.wght=a.wght, alpha=alphas, fixedFunctionArgs=fixedFunctionArgs)
    
    # if(is.null(XObs) || is.null(XPred))
    #   LKMLE = LKrig(coords, obs, LKinfo=LKinfo, verbose=verbose, maxit=maxit)
    # else
    #   LKMLE = LKrig(coords, obs, LKinfo=LKinfo, verbose=verbose, maxit=maxit, Z=XObs)
    out = LKrig(coords, obs, LKinfo=LKinfo, verbose=thisVerbose, Z=XObs)
    out$lnProfileLike.FULL
  }
  
  # get initial parameters and optimization bounds
  if(nLayer == 1) {
    init = c(log.lambda=log(lambdaStart), omega=Awght2Omega(a.wghtStart, LKinfoStart))
    lower = c(lowerBoundLogLambda, lowerBoundOmega)
    upper = c(upperBoundLogLambda, upperBoundOmega)
  }
  else {
    init = c(logit.alphas=multivariateLogit(rep(1 / nLayer, nLayer-1)), log.lambda=log(lambdaStart), omega=Awght2Omega(a.wghtStart, LKinfoStart))
    lower = c(lowerBoundLogitAlpha, lowerBoundLogLambda, lowerBoundOmega)
    upper = c(upperBoundLogitAlpha, upperBoundLogLambda, upperBoundOmega)
  }
  
  # Maximum Likelihood
  if(verbose)
    print("Beginning LatticeKrig MLE fit...")
  
  result <- try(optim(init,
                      outerFun, 
                      lower=lower, 
                      upper=upper, 
                      method="L-BFGS-B",
                      #                      method="BFGS",
                      control=list(fnscale = -1,factr=factr,
                                   pgtol=pgtol, maxit=maxit,
                                   ndeps = rep(.05,length(init)))
  ))
  if(verbose){
    cat("Results from optimize:", fill=TRUE)
    print( result )
  }
  
  # final fit
  parameterList = getParameters(result$par)
  alphasMLE = parameterList$alphas
  log.lambdaMLE = parameterList$log.lambda
  lambdaMLE = parameterList$lambda
  omegaMLE = parameterList$omega
  a.wghtMLE = parameterList$a.wght
  
  # # set up the lattice, the arguments to LatticeKrig
  # LKinfo = LKrigSetup(coords, nlevel=nLayer, nu=NULL, NC=NC, normalize=normalize, 
  #                     lambda=lambdaMLE, a.wght=a.wghtMLE, alpha=alphasMLE, fixedFunctionArgs=fixedFunctionArgs)
  # if(is.null(XObs) || is.null(XPred)) {
  #   mod = LKrig(coords, obs, LKinfo=LKinfo)
  #   preds = predict.LKrig(mod, predCoords)
  #   if(doSEs)
  #     predSimulations = LKrig.sim.conditional(mod, x.grid=predCoords, M=nsimConditional)
  #   else
  #     predSimulations = NULL
  # }
  # else {
  #   mod = LKrig(coords, obs, LKinfo=LKMLE$LKinfo, Z=XObs)
  #   preds = predict.LKrig(mod, predCoords, Znew=XPred)
  #   if(doSEs)
  #     predSimulations = LKrig.sim.conditional(mod, x.grid=predCoords, Z.grid=XPred, M=nsimConditional)
  #   else
  #     predSimulations=NULL
  # }
  # 
  # # calculate predictive standard errors
  # predSEs = predSimulations$SE
  
  # # calculate confidence intervals
  # lower = preds + qnorm((1 - significanceCI) / 2, sd=predSEs)
  # upper = preds + qnorm(1 - (1 - significanceCI) / 2, sd=predSEs)
  # 
  # ## now we calculate uncertainty intervals for all parameters
  # # intercept
  # interceptSummary = c(Est=mod$d.coef[1], SD=sd(predSimulations$d.coef.draw[1,]), 
  #                      Qlower=quantile(probs=(1 - significanceCI) / 2, predSimulations$d.coef.draw[1,]), 
  #                      Q50=quantile(probs=.5, predSimulations$d.coef.draw[1,]), 
  #                      Qupper=quantile(probs=1 - (1 - significanceCI) / 2, predSimulations$d.coef.draw[1,]), )
  
  # to calculate summaries for the parameters, must calculate inverse of negative hessian
  print("Calculating hessian...")
  hess = hessian(outerFun, result$par, thisVerbose=FALSE)
  parSigma = solve(-hess)
  
  # simulate possible parameter values
  L = t(chol(parSigma))
  zSim = matrix(rnorm(nsimConditional * nrow(L)), nrow=nrow(L))
  parSim = L %*% zSim
  parSim = sweep(parSim, 1, result$par, "+")
  
  # for each simulated set of parameters, call LKrig and calculate the reparameterizations
  getMLEs = function(parameters) {
    i = parameters[1]
    print(paste0("Generating draw ", i, "/", nsimConditional, " from conditional distribution"))
    parameters = parameters[-1]
    parameterList = getParameters(parameters)
    alphas = parameterList$alphas
    log.lambda = parameterList$log.lambda
    lambda = parameterList$lambda
    omega = parameterList$omega
    a.wght = parameterList$a.wght
    
    # set up the lattice, the arguments to LatticeKrig
    LKinfo = LKrigSetup(coords, nlevel=nLayer, nu=NULL, NC=NC, normalize=normalize, 
                        lambda=lambda, a.wght=a.wght, alpha=alphas, fixedFunctionArgs=fixedFunctionArgs)
    
    # fit the model to the parameters
    out = LKrig(coords, obs, LKinfo=LKinfo, verbose=FALSE, Z=XObs)
    
    # calculate the effective range (the distance at which there is only 10% correlation)
    test = LKrig.cov.plot(LKinfo)
    sortI = sort(test$d, index.return=TRUE)$ix
    firstI = match(TRUE, (test$cov[sortI] * (1 / max(test$cov))) <= .1)
    effectiveRange = test$d[sortI][firstI]
    if(is.na(effectiveRange)) {
      effectiveRange = max(test$d)
      warning("Effective range at least as large as the grid width. Setting it to the grid width")
    }
    
    # get relevant MLEs
    c(intercept=out$d.coef[1], range=effectiveRange, sqrtrho=sqrt(out$rho.MLE), rho=out$rho.MLE, sigma=out$shat.MLE, sigmasq=out$shat.MLE^2, alphas=alphas, a.wght=a.wght)
  }
  
  parameterDrawTable = apply(rbind(1:nsimConditional, parSim), 2, getMLEs)
  
  getSummaryStatistics = function(draws) {
    c(Est=mean(draws), SD=sd(draws), 
      Qlower=quantile(probs=(1 - significanceCI) / 2, draws), 
      Q50=quantile(probs=0.5, draws), 
      Qupper=quantile(probs=1 - (1 - significanceCI) / 2, draws))
  }
  
  parameterSummaryTable = apply(parameterDrawTable, 1, getSummaryStatistics)
  
  ## preds
  ## sigmas
  ## lower
  ## upper
  ## interceptSummary
  ## rangeSummary
  ## sdSummary
  ## varSummary
  return(list(mod=mod, preds=preds, sigmas=predSEs, lower=lower, upper=upper, parameterSummaryTable=parameterSummaryTable, LKinfo=LKinfo))
}