## create the rgeneric model for inla
# library(MCMCpack)

# fixed constants that must be included in the inla.rgeneric.define call:
# n: number of observations
# xRange: range of x coordinates
# yRange: range of y coordinates
# nx: number of basis function lattice points in x directions
# NOTE: ny is determined automatically to match scale of x lattice points
# corRange: log(corScale) ~ log(Exp^(-1))(rate=corRange) (high rate -> longer range correlation)
# alpha1, alpha2: log(variance) ~ log(InvGamma(alpha1, rate=alpha2)) (high rate -> high variance)
# ys: observations
# first.time: is first time evaluating function.  User should always set to FALSE
inla.rgeneric.lk.model.test = function(
  cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
          "log.prior", "quit"),
  theta = NULL)
{
  if(first.time) {
    # load relevant external functions
    library(spam)
    library(Matrix)
    library(fields)
    source("~/Google Drive/UW/Wakefield/WakefieldShared/LK-INLA/code/LKinla.R")
    
    first.time <<- FALSE
  }
  
  # theta is of the form:
  # c(exp(kappa), exp(rho), )
  interpret.theta = function()
  {
    ## internal helper-function to map the parameters from the internal-scale to the
    ## user-scale
    
    # precomputations: get lattice grid cell width, convert parameters from effective correlation 
    # and marginal variance to kappa and rho
    latticeWidth = (xRange[2]-xRange[1])/(nx-1)
    effectiveCor = exp(theta[1L])
    sigmaSq = exp(theta[2L])
    kappa = 2.3/effectiveCor * latticeWidth
    rho = sigmaSq * 4*pi * kappa^2
    
    list(effectiveCor = effectiveCor, 
         sigmaSq = sigmaSq, 
         kappa = kappa, 
         rho = rho)
  }
  
  if(cmd != "initial") {
    pars = interpret.theta()
    effectiveCor = pars$effectiveCor
    sigmaSq = pars$sigmaSq
    kappa = pars$kappa
    rho = pars$rho
  }
  
  # returns matrix just like Q except ony zero and 1 for nonzero elements
  graph = function(){
    # first compute ny
    xLength = xRange[2] - xRange[1]
    yLength = yRange[2] - yRange[1]
    latticeWidth = xLength/(nx-1)
    ny = floor(yLength/latticeWidth) + 1
    
    # now compute the graph/adjacency matrix
    makeGraph(nx, ny)
  }
  
  # compute the precision matrix
  Q = function() {
    # first compute ny and determine rounded yRange to have lattice pts at edges
    xLength = xRange[2] - xRange[1]
    yLength = yRange[2] - yRange[1]
    latticeWidth = xLength/(nx-1)
    ny = floor(yLength/latticeWidth) + 1
    yRangeMod = (ny-1) * latticeWidth
    yRangeMod = c(yRange[1], yRange[1] + yRangeMod)
    
    # compute precision matrix
    return(makeQ(kappa, rho, xRange, yRangeMod, nx, ny, rawKnotRange=TRUE))
  }
  
  # get mean of each observation
  mu = function() {
    # first compute ny
    xLength = xRange[2] - xRange[1]
    yLength = yRange[2] - yRange[1]
    latticeWidth = xLength/(nx-1)
    ny = floor(yLength/latticeWidth) + 1
    
    # return(rep(0, n))
    return(rep(0, nx*ny))
  }
  
  log.norm.const = function()
  {
    ## let INLA compute it as -n/2 log(2pi) + 1/2 * log(|Q|)
    return (numeric(0))
  }
  
  log.prior = function() {
    require(invgamma)
    
    # get prior (note the Jacobian factors)
    dinvexp(effectiveCor, rate=corRange, log=TRUE) + log(effectiveCor) + 
      dinvgamma(sigmaSq, shape=alpha1, rate=alpha2, log=TRUE) + log(sigmaSq)
  }
  
  # NOTE: not yet tested for being reasonable
  initial = function() {
    
    # initialize kappa based on LR2011 Eq at top of page 6, Besag (1981)
    # theorCov = (1/(2*pi)) * besselK(corLInit*kappa, nu=0)
    # > besselK(.0256, nu=0)/(2*pi)
    # [1] 0.6019045
    # > besselK(1.24172, nu=0)*2
    # [1] 0.601908
    # > besselK(2.687, nu=0)*2
    # [1] 0.100024
    # > Matern(1, range=1, smoothness=1)
    # [1] 0.6019072
    # > Matern(3.21, range=1, smoothness=1)
    # [1] 0.1003779
    
    effectiveRangeInit = ((xRange[2] - xRange[1]) + (yRange[2] - yRange[1]))/2
    xLength = xRange[2] - xRange[1]
    squareWidth = xLength/(nx-1)
    sigmaSq = sum(ys - 0)^2/n
    
    return(c(log(effectiveRangeInit), log(sigmaSq)))
    
    ## don't do this stuff, since we're reparameterizing to correlation scale, marginal variance:
    # # kappaInit = .0256/corLInit
    # # kappaInit = 1.24172/corLInit
    # kappaInit = 2.687/effectiveRangeInit*squareWidth
    # 
    # # initialize rho based on inline formula of LR2011 middle of p. 6
    # # theorMargVar = rho/(4*pi*(4+kappa^2-4))
    # 
    # rhoInit = sigmaSq * 4 * pi * kappaInit^2
    # 
    # return(c(kappaInit, rhoInit))
  }
  
  quit = function() {
    return(invisible())
  }
  
  val = do.call(match.arg(cmd), args = list())
  return (val)
}

# test regeneric model
# buffer: buffer distance between domain edge of basis lattice and domain edge of data.
# n: number of observations
# xRange: range of x coordinates
# yRange: range of y coordinates
# nx: number of basis function lattice points in x directions
# NOTE: ny is determined automatically to match scale of x lattice points
# Xmat: design matrix
# ys: observations
# first.time: is first time evaluating function.  User should always set to FALSE
testLKModelBase = function(buffer=1) {
  
  # generate lattice and simulate observations
  nx=30
  ny=30
  n=50
  coords = matrix(runif(2*n), ncol=2)
  ys = sin(coords[,2]*2*pi)*2 + rnorm(n)
  
  # plot the observations
  par(mfrow=c(1,1))
  quilt.plot(coords, ys)
  
  # get lattice points, prediction points
  xRangeBasis = c(0-buffer, 1+buffer)
  yRangeBasis = c(0-buffer, 1+buffer)
  gridPts = make.surface.grid(list(xs=seq(xRangeBasis[1], xRangeBasis[2], l=nx), 
                                   ys=seq(yRangeBasis[1], yRangeBasis[2], l=ny)))
  
  # make prediction coordinates on a grid
  xRange=c(0,1)
  yRange=c(0,1)
  mx = 20
  my = 20
  predPts = make.surface.grid(list(x=seq(xRange[1], xRange[2], l=mx), y=seq(yRange[1], yRange[2], l=my)))
  
  # generate hyperparameters based on median and quantiles of inverse exponential and inverse gamma
  priorPar = getPrior(.1, .1, 10)
  
  # define the LK model
  first.time <<- TRUE
  rgen = inla.rgeneric.define(model=inla.rgeneric.lk.model.test, n=n, xRange=xRangeBasis, yRange=yRangeBasis, nx=nx, 
                              ys=ys, corRange=priorPar$corScalePar, alpha1=priorPar$varPar1, 
                              alpha2=priorPar$varPar2, first.time=TRUE)
  
  ## generate inla stack:
  # Stacked A matrix (A_s from notation of LR2015 Bayesian Spatial Modelling with R_INLA):
  # (AEst   0  )
  # ( 0   APred)
  # eta_s = (c^T c^T)^T
  # where c is the vector of lattice coefficients
  AEst = makeA(coords, xRangeBasis, nx, yRangeBasis, ny)
  APred = makeA(predPts, xRangeBasis, nx, yRangeBasis, ny)
  latticeInds = 1:ncol(AEst)
  stack.est = inla.stack(A =list(AEst), 
                         effects =list(field=latticeInds), 
                         data =list(y=ys), 
                         tag ="est", 
                         remove.unused=FALSE)
  stack.pred = inla.stack(A =list(APred), 
                         effects =list(field=latticeInds), 
                         data =list(y=NA), 
                         tag ="pred", 
                         remove.unused=FALSE)
  stack.full = inla.stack(stack.est, stack.pred, 
                          remove.unused=FALSE)
  dat = inla.stack.data(stack.full, rgen=rgen, remove.unused=FALSE)
  
  # show priors on effective correlation and marginal variance:
  xs1 = seq(.01, 1.5, l=100)
  plot(xs1, dinvexp(xs1, rate=priorPar$corScalePar), type="l", col="blue", 
       xlab="Effective Correlation", main="Effective Correlation Prior", 
       ylab="Prior Density")
  abline(v=qinvexp(.5, rate=priorPar$corScalePar), col="red")
  
  xs2 = seq(.01, 15, l=100)
  plot(xs2, dinvgamma(xs2, shape=priorPar$varPar1, rate=priorPar$varPar2), type="l", col="blue", 
       xlab="Marginal Variance", main="Marginal Variance Prior", 
       ylab="Prior Density")
  abline(v=qinvgamma(.1, shape=priorPar$varPar1, rate=priorPar$varPar2), col="red")
  abline(v=qinvgamma(.9, shape=priorPar$varPar1, rate=priorPar$varPar2), col="red")
  
  # fit the model
  mod = inla(y ~ - 1 + f(field, model=rgen), data=dat, 
             control.predictor=list(A=inla.stack.A(stack.full), compute=TRUE), 
             family="normal", verbose=TRUE)
  
  # show a model summary
  print(summary(mod))
  
  # show predictive surface, SD, and data
  index = inla.stack.index(stack.full, "pred")$data
  linpred.mean = mod[["summary.linear.predictor"]]$mean
  linpred.sd = mod[["summary.linear.predictor"]]$sd
  
  par(mfrow=c(2,3))
  colRange = range(ys)
  colRangeSD = range(linpred.sd)
  obsInds = 1:n
  predInds = (n+1):(n+mx*my)
  coefInds = (n+mx*my+1):(n+mx*my+nx*ny)
  quilt.plot(coords, ys, xlim=xRangeBasis, ylim=yRangeBasis, main="Data", zlim=colRange)
  quilt.plot(gridPts[,1], gridPts[,2], linpred.mean[coefInds], main="Basis Coefficient Mean", zlim=colRange)
  quilt.plot(gridPts[,1], gridPts[,2], linpred.sd[coefInds], main="Basis Coefficient SD", zlim=colRangeSD)
  
  quilt.plot(coords, linpred.mean[obsInds], xlim=xRangeBasis, ylim=yRangeBasis, main="Observation Mean", zlim=colRange)
  quilt.plot(predPts[,1], predPts[,2], linpred.mean[predInds], main="Prediction Mean", zlim=colRange, 
             xlim=xRangeBasis, ylim=yRangeBasis)
  quilt.plot(predPts[,1], predPts[,2], linpred.sd[predInds], main="Prediction SD", zlim=colRangeSD, 
             xlim=xRangeBasis, ylim=yRangeBasis)
  
  # plot marginals on interpretable scale (effective range, marginal variance)
  effRangeMarg = inla.tmarginal(exp, mod$marginals.hyperpar$`Theta1 for field`)
  varMarg = inla.tmarginal(exp, mod$marginals.hyperpar$`Theta2 for field`)
  
  par(mfrow=c(1,1))
  plot(effRangeMarg, type="l", main="Marginal for effective range")
  # plot(mod$marginals.hyperpar$`Theta1 for field`, type="l", main="Marginal for log range")
  plot(varMarg, type="l", main="Marginal for process variance")
  # plot(mod$marginals.hyperpar$`Theta2 for field`, type="l", main="Marginal for log variance")
  
  # do the same for kappa, rho
  # in order to get distribution for rho, must sample from joint hyperparameters
  latticeWidth = (xRangeBasis[2] - xRangeBasis[1])/(nx-1)
  kappaMarg = inla.tmarginal(function(x) {2.3/exp(x) * latticeWidth}, mod$marginals.hyperpar$`Theta1 for field`)
  thetasToRho = function(xs) {
    logCor = xs[2]
    logVar = xs[3]
    kappa = 2.3/exp(logCor) * latticeWidth
    sigma2 = exp(logVar)
    sigma2 * 4*pi * kappa^2
  }
  samples = inla.hyperpar.sample(50000, mod, TRUE)
  rhos = apply(samples, 1, thetasToRho)
  plot(kappaMarg, type="l", xlab="kappa", main="Marginal for kappa")
  hist(rhos, xlab="rho", main="Marginal for Rho", breaks=1000, freq=F, xlim=c(0, quantile(probs=.95, rhos)))
}

# test regeneric model using data simulated from LK model
# buffer: buffer distance between domain edge of basis lattice and domain edge of data.
# n: number of observations
# xRange: range of x coordinates
# yRange: range of y coordinates
# nx: number of basis function lattice points in x directions
# NOTE: ny is determined automatically to match scale of x lattice points
# Xmat: design matrix
# ys: observations
# first.time: is first time evaluating function.  User should always set to FALSE
testLKModelBase2 = function(buffer=1, kappa=1, rho=1, seed=1) {
  set.seed(seed)
  
  # get lattice points, prediction points
  nx=30
  ny=30
  xRangeBasis = c(0-buffer, 1+buffer)
  yRangeBasis = c(0-buffer, 1+buffer)
  gridPts = make.surface.grid(list(xs=seq(xRangeBasis[1], xRangeBasis[2], l=nx), 
                                   ys=seq(yRangeBasis[1], yRangeBasis[2], l=ny)))
  
  # generate lattice and simulate observations
  n=50
  coords = matrix(runif(2*n), ncol=2)
  AObs = makeA(coords, xRangeBasis, nx, yRangeBasis, ny)
  Q = makeQ(kappa=kappa, rho=rho, xRange=xRangeBasis, yRange=yRangeBasis, nx=nx, ny=ny)
  L = t(chol(solve(Q)))
  zsims = matrix(rnorm(nx*ny), ncol=1)
  fieldSims = L %*% zsims
  ys = as.numeric(AObs %*% fieldSims)
  
  # plot the observations
  par(mfrow=c(1,1))
  quilt.plot(coords, ys)
  
  # make prediction coordinates on a grid
  xRange=c(0,1)
  yRange=c(0,1)
  mx = 20
  my = 20
  predPts = make.surface.grid(list(x=seq(xRange[1], xRange[2], l=mx), y=seq(yRange[1], yRange[2], l=my)))
  
  # generate hyperparameters based on median and quantiles of inverse exponential and inverse gamma
  priorPar = getPrior(.1, .1, 10)
  
  # define the LK model
  first.time <<- TRUE
  rgen = inla.rgeneric.define(model=inla.rgeneric.lk.model.test, n=n, xRange=xRangeBasis, yRange=yRangeBasis, nx=nx, 
                              ys=ys, corRange=priorPar$corScalePar, alpha1=priorPar$varPar1, 
                              alpha2=priorPar$varPar2, first.time=TRUE)
  
  ## generate inla stack:
  # Stacked A matrix (A_s from notation of LR2015 Bayesian Spatial Modelling with R_INLA):
  # (AEst   0  )
  # ( 0   APred)
  # eta_s = (c^T c^T)^T
  # where c is the vector of lattice coefficients
  AEst = makeA(coords, xRangeBasis, nx, yRangeBasis, ny)
  APred = makeA(predPts, xRangeBasis, nx, yRangeBasis, ny)
  latticeInds = 1:ncol(AEst)
  stack.est = inla.stack(A =list(AEst), 
                         effects =list(field=latticeInds), 
                         data =list(y=ys), 
                         tag ="est", 
                         remove.unused=FALSE)
  stack.pred = inla.stack(A =list(APred), 
                          effects =list(field=latticeInds), 
                          data =list(y=NA), 
                          tag ="pred", 
                          remove.unused=FALSE)
  stack.full = inla.stack(stack.est, stack.pred, 
                          remove.unused=FALSE)
  dat = inla.stack.data(stack.full, rgen=rgen, remove.unused=FALSE)
  
  # show priors on effective correlation and marginal variance:
  xs1 = seq(.01, 1.5, l=100)
  plot(xs1, dinvexp(xs1, rate=priorPar$corScalePar), type="l", col="blue", 
       xlab="Effective Correlation", main="Effective Correlation Prior", 
       ylab="Prior Density")
  abline(v=qinvexp(.5, rate=priorPar$corScalePar), col="red")
  
  xs2 = seq(.01, 15, l=100)
  plot(xs2, dinvgamma(xs2, shape=priorPar$varPar1, rate=priorPar$varPar2), type="l", col="blue", 
       xlab="Marginal Variance", main="Marginal Variance Prior", 
       ylab="Prior Density")
  abline(v=qinvgamma(.1, shape=priorPar$varPar1, rate=priorPar$varPar2), col="red")
  abline(v=qinvgamma(.9, shape=priorPar$varPar1, rate=priorPar$varPar2), col="red")
  
  # fit the model
  mod = inla(y ~ - 1 + f(field, model=rgen), data=dat, 
             control.predictor=list(A=inla.stack.A(stack.full), compute=TRUE), 
             family="normal", verbose=TRUE)
  
  # show a model summary
  print(summary(mod))
  
  # show predictive surface, SD, and data
  index = inla.stack.index(stack.full, "pred")$data
  linpred.mean = mod[["summary.linear.predictor"]]$mean
  linpred.sd = mod[["summary.linear.predictor"]]$sd
  
  par(mfrow=c(2,3))
  colRange = range(ys)
  colRangeSD = range(linpred.sd)
  obsInds = 1:n
  predInds = (n+1):(n+mx*my)
  coefInds = (n+mx*my+1):(n+mx*my+nx*ny)
  quilt.plot(coords, ys, xlim=xRangeBasis, ylim=yRangeBasis, main="Data", zlim=colRange)
  quilt.plot(gridPts[,1], gridPts[,2], linpred.mean[coefInds], main="Basis Coefficient Mean", zlim=colRange)
  quilt.plot(gridPts[,1], gridPts[,2], linpred.sd[coefInds], main="Basis Coefficient SD", zlim=colRangeSD)
  
  quilt.plot(coords, linpred.mean[obsInds], xlim=xRangeBasis, ylim=yRangeBasis, main="Observation Mean", zlim=colRange)
  quilt.plot(predPts[,1], predPts[,2], linpred.mean[predInds], main="Prediction Mean", zlim=colRange, 
             xlim=xRangeBasis, ylim=yRangeBasis)
  quilt.plot(predPts[,1], predPts[,2], linpred.sd[predInds], main="Prediction SD", zlim=colRangeSD, 
             xlim=xRangeBasis, ylim=yRangeBasis)
  
  # calculate true effective range and marginal variance:
  latticeWidth = (xRangeBasis[2] - xRangeBasis[1])/(nx-1)
  effRange = 2.3/kappa * latticeWidth
  marginalVar = rho/(4*pi * kappa^2)
  
  # plot marginals on interpretable scale (effective range, marginal variance)
  effRangeMarg = inla.tmarginal(exp, mod$marginals.hyperpar$`Theta1 for field`)
  varMarg = inla.tmarginal(exp, mod$marginals.hyperpar$`Theta2 for field`)
  
  par(mfrow=c(1,1))
  plot(effRangeMarg, type="l", main="Marginal for effective range")
  abline(v=effRange, col="green")
  # plot(mod$marginals.hyperpar$`Theta1 for field`, type="l", main="Marginal for log range")
  plot(varMarg, type="l", main="Marginal for process variance")
  abline(v=marginalVar, col="green")
  # plot(mod$marginals.hyperpar$`Theta2 for field`, type="l", main="Marginal for log variance")
  
  # do the same for kappa, rho
  # in order to get distribution for rho, must sample from joint hyperparameters
  kappaMarg = inla.tmarginal(function(x) {2.3/exp(x) * latticeWidth}, mod$marginals.hyperpar$`Theta1 for field`)
  thetasToRho = function(xs) {
    logCor = xs[2]
    logVar = xs[3]
    kappa = 2.3/exp(logCor) * latticeWidth
    sigma2 = exp(logVar)
    sigma2 * 4*pi * kappa^2
  }
  samples = inla.hyperpar.sample(50000, mod, TRUE)
  rhos = apply(samples, 1, thetasToRho)
  plot(kappaMarg, type="l", xlab="kappa", main="Marginal for kappa")
  abline(v=kappa, col="green")
  hist(rhos, xlab="rho", main="Marginal for Rho", breaks=100, freq=F, xlim=c(0, quantile(probs=.95, rhos)))
  abline(v=rho, col="green")
}

## function for making prior on corelation scale and marginal variance
# corScaleMed: median correlation scale
# var10quant: 10th percentile of the marginal variance
# var90quant: 90th percentile of the marginal variance
getPrior = function(corScaleMed, var10quant, var90quant) {
  ## set the correlation scale hyperparameter
  corScalePar = corScaleMed*log(2)
  
  ## set the marginal variance hyperparameters
  # define objective function for quantile matching
  quantMatch = function(pars) {
    require(invgamma)
    alpha1 = pars[1]
    alpha2 = pars[2]
    q10 = qinvgamma(.1, shape=alpha1, rate=alpha2)
    q90 = qinvgamma(.9, shape=alpha1, rate=alpha2)
    (var10quant - q10)^2 + (var90quant - q90)^2
  }
  out = optim(c(1,1), quantMatch)
  alphas = out$par
  varPar1 = alphas[1]
  varPar2 = alphas[2]
  
  list(corScalePar=corScalePar, varPar1=varPar1, varPar2=varPar2)
}








##### now make "simple" model for multi-layer LK-INLA:

# fixed constants that must be included in the inla.rgeneric.define call:
# n: number of observations
# nLayer: number of lattice layers
# nx: number of basis function lattice points in x directions for coarsest lattice
# NOTE: ny is determined automatically to match scale of x lattice points
# gamma: log(corScale) ~ log(Exp^(-1))(rate=gamma) (high rate -> longer range correlation)
# alpha1, alpha2: log(variance) ~ log(InvGamma(alpha1, rate=alpha2)) (high rate -> high variance)
# ys: observations
# datCoords: observation locations
# first.time: is first time evaluating function.  User should always set to FALSE
# X: design matrices for observations locations (only needed for initialization)
# nu: fixed smoothness parameter
# nBuffer: number of basis functions to use as buffer outside knot range
# NC: number of basis elements along longest dimension of data range
# normalize: normalize precision matrix so marginal variance is 1 at center and all coeffcients have same variance. slower
# fastNormalize: simple normalization, marginal variance is 1 at center, but not all coefficients have same variance.  only slightly slower
inla.rgeneric.lk.model.simple = function(
  cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
          "log.prior", "quit"),
  theta = NULL)
{
  if(first.time) {
    # load relevant external functions
    library(spam)
    library(Matrix)
    library(fields)
    source("~/Google Drive/UW/Wakefield/WakefieldShared/LK-INLA/code/LKinla.R")
    source("~/Google Drive/UW/Wakefield/WakefieldShared/LK-INLA/code/LKinla_rgeneric.R")
    browser()
    ## build marginal variance to kappa transformation (make it on a log scale):
    # first we need to know over what scale to build the original estimator.  Base this on 
    # an initial guess for marginal variance
    
    # # get range of data
    # xRangeDat = range(datCoords[,1])
    # yRangeDat = range(datCoords[,2])
    # 
    # # get lattice knots
    # out = makeLatGrid(xRangeDat, yRangeDat, NC, nBuffer)
    # xRange = out$xRangeKnots
    # yRange = out$yRangeKnots
    # nx = out$nx
    # ny = out$ny
    # xLength = diff(xRange)
    # yLength = diff(yRange)
    # latWidth = xLength/(nx-1)
    # knotPts = make.surface.grid(list(x=seq(xRange[1], xRange[2], l=nx),
    #                                  y=seq(yRange[1], yRange[2], l=ny)))
    # 
    # # initialize covariance parameters
    # effectiveRangeInit = ((xRange[2] - xRange[1]) + (yRange[2] - yRange[1]))/2
    
    # use global assignment for testing:
    # get range of data
    xRangeDat <<- range(datCoords[,1])
    yRangeDat <<- range(datCoords[,2])

    # get lattice knots
    out = makeLatGrid(xRangeDat, yRangeDat, NC, nBuffer)
    xRange <<- out$xRangeKnots
    yRange <<- out$yRangeKnots
    nx <<- out$nx
    ny <<- out$ny
    xLength <<- diff(xRange)
    yLength <<- diff(yRange)
    latWidth <<- xLength/(nx-1)
    knotPts <<- make.surface.grid(list(x=seq(xRange[1], xRange[2], l=nx),
                                       y=seq(yRange[1], yRange[2], l=ny)))

    # initialize covariance parameters
    effectiveRangeInit <<- ((xRange[2] - xRange[1]) + (yRange[2] - yRange[1]))/2
    
    # # convert to theoretical kappa estimate
    # kappaEst = 2.3/effectiveRangeInit * latWidth
    # 
    # # now to the important part: get an initial estimate of marginal variance and a range (go with /100 and *100 of initial guess of kappa)
    # kappas <<- 10^(seq(log10(kappaEst)-2, log10(kappaEst)+2, l=100))
    # kappaWidths <<- log10(kappas[2]) - log10(kappas[1])
    # margVars <<- sapply(kappas, getMultiMargVar, rho=1, tod=2.5, nx=nx, ny=ny, nu=nu, nLayer=nLayer)[1,]
    # logMargVarSpline <<- splinefun(log(kappas), margVars, "natural")
    
    first.time <<- FALSE
  }
  
  # theta is of the form:
  # c(betas, effectiveCor, sigmaSq, kappa, rho, nu, alphas)
  interpret.theta = function()
  {
    ## internal helper-function to map the parameters from the internal-scale to the
    ## user-scale
    effectiveCor = exp(theta[1])
    sigmaSq = exp(theta[2])
    
    # precomputations: get lattice grid cell width, convert parameters from effective correlation 
    # and marginal variance to kappa and rho.  Use spline to convert from marginal variance to kappa
    latticeWidth = (xRange[2]-xRange[1])/(nx-1)
    kap = 2.3/effectiveCor * latticeWidth
    
    # since we are normalizing the process, rho is just sigmaSq
    rho = sigmaSq
    
    # # kap = 2.3/effectiveCor * latticeWidth
    # 
    # # If we're at a value of kappa outside out current reparameterization range, 
    # # adjust the range by refitting the spline function
    # if(kap > max(kappas)) {
    #   newKappas = 10^(seq(log10(kappas[length(kappas)] + kappaWidths), log10(kap)+.5, by=kappaWidths))
    #   kappas <<- c(kappas, newKappas)
    #   kappas = c(kappas, newKappas)
    #   margVars <<- c(margVars, sapply(newKappas, getMultiMargVar, rho=1, tod=2.5, nx=nx, ny=ny, nu=nu, nLayer=nLayer)[1,])
    #   logMargVarSpline <<- splinefun(log(kappas), margVars, "natural")
    # }
    # else if(kap < min(kappas)) {
    #   newKappas = 10^(seq(log10(kap) - .5, log10(kappas[1] - kappaWidths), by=kappaWidths))
    #   kappas <<- c(newKappas, kappas)
    #   margVars <<- c(sapply(newKappas, getMultiMargVar, rho=1, tod=2.5, nx=nx, ny=ny, nu=nu, nLayer=nLayer)[1,], margVars)
    #   logMargVarSpline <<- splinefun(log(kappas), margVars, "natural")
    # }
    
    # # now find rho using reparameterization:
    # # h(log kappa) = sigma^2/rho
    # # rho = sigma^2/h(log kappa)
    # 
    # # rho = sigmaSq * 4*pi * kappa^2
    # rho = sigmaSq / logMargVarSpline(log(kap))
    
    # compute layer weights, alpha_1, ..., alpha_L
    L = nLayer
    alphas = getAlphas(L, nu)
    
    list(effectiveCor = effectiveCor, 
         sigmaSq = sigmaSq, 
         kappa = kap, 
         rho = rho, 
         alphas = alphas)
  }
  
  if(cmd != "initial") {
    pars = interpret.theta()
    effectiveCor = pars$effectiveCor
    sigmaSq = pars$sigmaSq
    kap = pars$kappa
    rho = pars$rho
    alphas = pars$alphas
  }
  
  # returns matrix just like Q except ony zero and 1 for nonzero elements
  graph = function(){
    graph = makeGraph(nx, ny, nLayer=nLayer, xRange=xRange, yRange=yRange, 
                      xRangeDat=xRangeDat, yRangeDat=yRangeDat, nBuffer=nBuffer, NC=NC)
    
    return(graph)
  }
  
  # compute the precision matrix
  Q = function() {
    ## renormalize basis coefficients
    
    # make mid lattice point
    # xi = ceiling(nx/2)
    # yi = ceiling(ny/2)
    # midPt = matrix(c(seq(xRange[1], xRange[2], l=nx)[xi], seq(yRange[1], yRange[2], l=ny)[yi]), nrow=1)
    
    # now construct relevant row of A for value at midpoint, and Q matrix
    # Ai = makeA(midPt, xRange, nx, yRange, ny, nLayer=nLayer, xRangeDat=xRangeDat, 
    #            yRangeDat=yRangeDat, nBuffer=nBuffer)
    
    Q = makeQ(kap, rho, nx=nx, ny=ny, nLayer=nLayer, alphas=alphas, 
              xRangeDat=xRangeDat, yRangeDat=yRangeDat, nBuffer=nBuffer, normalized=normalize, 
              fastNormalize=fastNormalize, NC=NC)
    return(Q)
  }
  
  # get mean of each latent coefficient
  mu = function() {
    
    # get LatticeKrig grid parameters, number of knots in each layer
    ms = getMs(xRangeDat, yRangeDat, NC, nBuffer, nLayer)
    
    # total number of basis functions should be this number
    # return(rep(0, nx*ny*(4^nLayer - 1)/3))
    rep(0, sum(ms))
  }
  
  log.norm.const = function()
  {
    ## let INLA compute it as -n/2 log(2pi) + 1/2 * log(|Q|)
    return (numeric(0))
  }
  
  log.prior = function() { ##### TODO: is nu fit or fixed?
    require(invgamma)
    
    # get prior (note the Jacobian factors)
    dinvexp(effectiveCor, rate=corRange, log=TRUE) + log(effectiveCor) + 
      dinvgamma(sigmaSq, shape=alpha1, rate=alpha2, log=TRUE) + log(sigmaSq)
  }
  
  initial = function() {
    # # initialize fixed effects and process variances
    # mod = lm(ys ~ X-1)
    # betas = coef(mod)
    # 
    # # remove from observations
    # ysCntr = ys - X %*% betas
    # 
    # # initialize covariance parameters
    # effectiveRangeInit = ((xRange[2] - xRange[1]) + (yRange[2] - yRange[1]))/2
    # xLength = xRange[2] - xRange[1]
    # squareWidth = xLength/(nx-1)
    # sigmaSq = sum((ysCntr - 0)^2)/n
    
    # initialize process variance by estimating spatial process with OLS
    AMat = makeA(datCoords, xNKnot=nx, yNKnot=ny, xRangeDat=xRangeDat, 
                 yRangeDat=yRangeDat, nBuffer=nBuffer, NC=NC)
    mod = lm(ys ~ cbind(X, as.matrix(AMat)) - 1)
    r2 = summary(mod)$r.squared
    s2 = summary(mod)$sigma^2
    # sigmaSq = (r2/(1-r2)) * s2
    sigmaSq = var(ys) * r2
    
    # initialize covariance parameters
    effectiveRangeInit = ((xRangeDat[2] - xRangeDat[1]) + (yRangeDat[2] - yRangeDat[1]))/8
    
    return(c(log(effectiveRangeInit), log(sigmaSq)))
  }
  
  quit = function() {
    return(invisible())
  }
  
  val = do.call(match.arg(cmd), args = list())
  return (val)
}

# test simple regeneric model using data simulated from simple LK model
# buffer: buffer distance between domain edge of basis lattice and domain edge of data.
# n: number of observations
# xRange: range of x coordinates
# yRange: range of y coordinates
# nx: number of basis function lattice points in x directions
# NOTE: ny is determined automatically to match scale of x lattice points
# Xmat: design matrix
# ys: observations
# first.time: is first time evaluating function.  User should always set to FALSE
testLKModelSimple = function(buffer=2.5, kappa=1, rho=1, nu=1.5, seed=1, nLayer=3, nx=20, ny=nx, n=50, sigma2 = .1^2, 
                             nBuffer=5, normalize=TRUE, NC=5) {
  set.seed(seed)
  
  # compute alphas, the variance weights for each layer, depending on nu:
  alphas = getAlphas(nLayer, nu)
  
  # get lattice points, prediction points
  xRangeDat = c(0,1)
  yRangeDat = c(0,1)
  
  # generate lattice and simulate observations
  coords = matrix(runif(2*n), ncol=2)
  
  # get lattice points, prediction points
  xRangeDat = range(coords[,1])
  yRangeDat = range(coords[,2])
  latInfo = makeLatGrid(xRangeDat, yRangeDat, NC, nBuffer)
  
  xRangeBasis = latInfo$xRangeKnots
  yRangeBasis = latInfo$yRangeKnots
  nx = latInfo$nx
  ny = latInfo$ny
  
  AObs = makeA(coords, xNKnot=nx, yNKnot=ny, nLayer=nLayer, xRangeDat=xRangeDat, 
               yRangeDat=yRangeDat, nBuffer=nBuffer, NC=NC)
  Q = makeQ(kappa=kappa, rho=rho, nx=nx, ny=ny, 
            nLayer=nLayer, alphas=alphas, xRangeDat=xRangeDat, yRangeDat=yRangeDat, 
            nBuffer=nBuffer, normalized = normalize, NC=NC) 
  L = as.matrix(t(chol(solve(Q))))
  zsims = matrix(rnorm(nrow(Q)), ncol=1)
  fieldSims = L %*% zsims
  ys = as.numeric(AObs %*% fieldSims) + 1 # add a constant unit mean term to be estimated by INLA
  # ys = 1 + as.numeric(AObs %*% fieldSims) + coords[,1] # x-valued mean term to be estimated by INLA
  errs = rnorm(n, sd=sqrt(sigma2))
  ys = ys + errs
  
  # plot the observations
  pdf(file="Figures/simpleObservations.pdf", width=5, height=5)
  par(mfrow=c(1,1))
  quilt.plot(coords, ys)
  dev.off()
  
  # make prediction coordinates on a grid
  xRange=c(0,1)
  yRange=c(0,1)
  mx = 20
  my = 20
  predPts = make.surface.grid(list(x=seq(xRange[1], xRange[2], l=mx), y=seq(yRange[1], yRange[2], l=my)))
  
  # generate hyperparameters based on median and quantiles of inverse exponential and inverse gamma
  priorPar = getPrior(.1, .1, 10)
  X = matrix(rep(1, n), ncol=1)
  # X = matrix(coords[,1], ncol=1)
  XPred = matrix(rep(1, mx*my), ncol=1)
  
  # define the LK model
  first.time <<- TRUE
  xRange=xRangeBasis
  yRange=yRangeBasis
  corRange=priorPar$corScalePar
  alpha1=priorPar$varPar1
  alpha2=priorPar$varPar2
  rgen = inla.rgeneric.define(model=inla.rgeneric.lk.model.simple, n=n, xRange=xRangeBasis, yRange=yRangeBasis, nx=nx, 
                              ys=ys, corRange=priorPar$corScalePar, alpha1=priorPar$varPar1, 
                              alpha2=priorPar$varPar2, first.time=TRUE, X=X, nu=nu, nLayer=nLayer, datCoords=coords, 
                              nBuffer=nBuffer, xRangeDat=xRangeDat, yRangeDat=yRangeDat, normalize=normalize, 
                              fastNormalize=fastNormalize)
  
  ## generate inla stack:
  # Stacked A matrix (A_s from notation of LR2015 Bayesian Spatial Modelling with R_INLA):
  # (AEst   0  )
  # ( 0   APred)
  # eta_s = (c^T c^T)^T
  # where c is the vector of lattice coefficients
  AEst = makeA(coords, xNKnot=nx, yNKnot=ny, nLayer=nLayer, xRangeDat=xRangeDat, 
               yRangeDat=yRangeDat, nBuffer=nBuffer, NC=NC)
  APred = makeA(predPts, xNKnot=nx, yNKnot=ny, nLayer=nLayer, xRangeDat=xRangeDat, 
                yRangeDat=yRangeDat, nBuffer=nBuffer, NC=NC)
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
  
  # show priors on effective correlation and marginal variance:
  xs1 = seq(.01, 1, l=500)
  pdf(file="Figures/simplePriorEffRange.pdf", width=5, height=5)
  plot(xs1, dinvexp(xs1, rate=priorPar$corScalePar), type="l", col="blue", 
       xlab="Effective Correlation Range", main="Effective Correlation Prior", 
       ylab="Prior Density")
  abline(v=qinvexp(.5, rate=priorPar$corScalePar), col="red")
  dev.off()
  
  xs2 = seq(.01, 10.5, l=500)
  pdf(file="Figures/simplePriorMargVar.pdf", width=5, height=5)
  plot(xs2, dinvgamma(xs2, shape=priorPar$varPar1, rate=priorPar$varPar2), type="l", col="blue", 
       xlab="Marginal Variance", main="Marginal Variance Prior", 
       ylab="Prior Density")
  abline(v=qinvgamma(.1, shape=priorPar$varPar1, rate=priorPar$varPar2), col="red")
  abline(v=qinvgamma(.9, shape=priorPar$varPar1, rate=priorPar$varPar2), col="red")
  dev.off()
  
  # fit the model
  mod = inla(y ~ - 1 + X + f(field, model=rgen), data=dat, 
             control.predictor=list(A=inla.stack.A(stack.full), compute=TRUE), 
             family="normal", verbose=TRUE)
  
  # show a model summary
  print(summary(mod))
  
  # show predictive surface, SD, and data
  index = inla.stack.index(stack.full, "pred")$data
  linpred.mean = mod[["summary.linear.predictor"]]$mean
  linpred.sd = mod[["summary.linear.predictor"]]$sd
  
  pdf(file="Figures/simplePreds.pdf", width=15, height=10)
  par(mfrow=c(2,3))
  obsInds = 1:n
  predInds = (n+1):(n+mx*my)
  coefInds = (n+mx*my+1):(n+mx*my+nx*ny)
  colRangeDat = range(c(ys-errs, linpred.mean[c(obsInds, predInds)]))
  colRangeSD = range(linpred.sd)
  gridPtsL1 = make.surface.grid(list(xs=seq(xRangeBasis[1], xRangeBasis[2], l=nx), 
                                     ys=seq(yRangeBasis[1], yRangeBasis[2], l=ny)))
  quilt.plot(coords, ys-errs, main="True Process", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
  quilt.plot(gridPtsL1[,1], gridPtsL1[,2], linpred.mean[coefInds], main="Basis Coefficient Mean (Layer 1)", xlim=xRangeDat, ylim=yRangeDat)
  quilt.plot(gridPtsL1[,1], gridPtsL1[,2], linpred.sd[coefInds], main="Basis Coefficient SD (Layer 1)", zlim=colRangeSD, xlim=xRangeDat, ylim=yRangeDat)
  
  quilt.plot(coords, linpred.mean[obsInds], main="Observation Mean", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
  quilt.plot(predPts[,1], predPts[,2], linpred.mean[predInds], main="Prediction Mean", zlim=colRangeDat, 
             xlim=xRangeDat, ylim=yRangeDat)
  quilt.plot(predPts[,1], predPts[,2], linpred.sd[predInds], main="Prediction SD", 
             xlim=xRangeDat, ylim=yRangeDat)
  dev.off()
  
  # calculate true effective range and marginal variance:
  latticeWidth = (xRangeBasis[2] - xRangeBasis[1])/(nx-1)
  effRange = 2.3/kappa * latticeWidth
  # marginalVar = rho/(4*pi * kappa^2)
  # marginalVar = getMultiMargVar(kappa, rho, nLayer=nLayer, nu=nu, xRange=xRangeBasis, 
  #                               yRange=yRangeBasis, nx=nx, ny=ny)[1]
  marginalVar = rho
  
  # plot marginals on interpretable scale (effective range, marginal variance)
  effRangeMarg = inla.tmarginal(exp, mod$marginals.hyperpar$`Theta1 for field`)
  varMarg = inla.tmarginal(exp, mod$marginals.hyperpar$`Theta2 for field`)
  sigma2Marg = inla.tmarginal(function(x) {1/x}, mod$marginals.hyperpar$`Precision for the Gaussian observations`)
  XMarginal = inla.tmarginal(function(x) {x}, mod$marginals.fixed$X)
  
  par(mfrow=c(1,1))
  pdf(file="Figures/simpleEffRange.pdf", width=5, height=5)
  plot(effRangeMarg, type="l", main="Marginal for effective range")
  abline(v=effRange, col="green")
  abline(v=inla.qmarginal(c(.025, .975), effRangeMarg), col="purple", lty=2)
  dev.off()
  # plot(mod$marginals.hyperpar$`Theta1 for field`, type="l", main="Marginal for log range")
  pdf(file="Figures/simpleVar.pdf", width=5, height=5)
  plot(varMarg, type="l", main="Marginal for process variance")
  abline(v=marginalVar, col="green")
  abline(v=inla.qmarginal(c(.025, .975), varMarg), col="purple", lty=2)
  dev.off()
  # plot(mod$marginals.hyperpar$`Theta2 for field`, type="l", main="Marginal for log variance")
  pdf(file="Figures/simpleSigma2.pdf", width=5, height=5)
  plot(sigma2Marg, type="l", main="Marginal for error variance")
  abline(v=sigma2, col="green")
  abline(v=inla.qmarginal(c(.025, .975), sigma2Marg), col="purple", lty=2)
  dev.off()
  pdf(file="Figures/simpleX.pdf", width=5, height=5)
  plot(XMarginal, type="l", main="Marginal for fixed effect")
  abline(v=1, col="green")
  abline(v=inla.qmarginal(c(.025, .975), XMarginal), col="purple", lty=2)
  dev.off()
  
  # pdf(file="Figures/simpleRho.pdf", width=5, height=5)
  # plot(sigma2Marg, type="l", main=TeX("Marginal for $\\rho$"), xlab=TeX("$\\rho$"))
  # abline(v=rho, col="green")
  # dev.off()
  
  
  # do the same for kappa, rho
  # in order to get distribution for rho, must sample from joint hyperparameters
  kappaMarg = inla.tmarginal(function(x) {2.3/exp(x) * latticeWidth}, mod$marginals.hyperpar$`Theta1 for field`)
  # thetasToRho = function(xs) {
  #   logCor = xs[2]
  #   logVar = xs[3]
  #   kappa = 2.3/exp(logCor) * latticeWidth
  #   sigma2 = exp(logVar)
  #   sigma2 * 4*pi * kappa^2
  # }
  # samples = inla.hyperpar.sample(50000, mod, TRUE)
  # rhos = apply(samples, 1, thetasToRho)
  
  pdf(file="Figures/simpleKappa.pdf", width=5, height=5)
  plot(kappaMarg, type="l", xlab="kappa", main="Marginal for kappa")
  abline(v=kappa, col="green")
  abline(v=inla.qmarginal(c(.025, .975), kappaMarg), col="purple", lty=2)
  dev.off()
  
  # pdf(file="Figures/simpleRho.pdf", width=5, height=5)
  # hist(rhos, xlab="rho", main="Marginal for Rho", breaks=1000, freq=F, xlim=c(0, quantile(probs=.95, rhos)))
  # abline(v=rho, col="green")
  # dev.off()
  invisible(NULL)
}


#####
#####
#####
#####
#####

# fit simple LK model given data and relevant parameters
# coords: spatial coordinates of the data
# ys: observations at the coordinates
# predPts: places at which to make predictions
# nu: Matern smoothness parameter from "simple" LKrig model
# X: observation design matrix (intercept and covariates at obs locations)
# XPred: prediction design matrix (intercept and covariates at pred locations)
# NC: number of coarse lattice points over the longest dimension of the data
# int.strategy: "auto" or "eb" for empirical Bayes
fitLKINLASimple = function(coords, ys, predPts=coords, nu=1.5, seed=1, nLayer=3, NC=5,
                           nBuffer=5, priorPar=getPrior(.1, .1, 10), 
                           X=cbind(1, coords), XPred=cbind(1, predPts), normalize=TRUE, 
                           int.strategy="auto", strategy="gaussian", fastNormalize=FALSE) {
  set.seed(seed)
  
  # get lattice points, prediction points
  xRangeDat = range(coords[,1])
  yRangeDat = range(coords[,2])
  latInfo = makeLatGrid(xRangeDat, yRangeDat, NC, nBuffer)
  
  xRangeBasis = latInfo$xRangeKnots
  yRangeBasis = latInfo$yRangeKnots
  nx = latInfo$nx
  ny = latInfo$ny
  
  # generate lattice basis matrix
  AObs = makeA(coords, xNKnot=nx, yNKnot=ny, nLayer=nLayer, xRangeDat=xRangeDat, 
               yRangeDat=yRangeDat, nBuffer=nBuffer, NC=NC)
  
  # define the LK model
  first.time <<- TRUE
  xRange=xRangeBasis
  yRange=yRangeBasis
  corRange=priorPar$corScalePar
  alpha1=priorPar$varPar1
  alpha2=priorPar$varPar2
  n = length(ys)
  rgen = inla.rgeneric.define(model=inla.rgeneric.lk.model.simple, n=n, nx=nx, NC=NC, 
                              ys=ys, corRange=priorPar$corScalePar, alpha1=priorPar$varPar1, normalize=normalize, 
                              alpha2=priorPar$varPar2, first.time=TRUE, X=X, nu=nu, nLayer=nLayer, datCoords=coords, 
                              nBuffer=nBuffer, xRangeDat=xRangeDat, yRangeDat=yRangeDat, fastNormalize=fastNormalize)
  # first.time <<- TRUE; n<<-n; NC <<- NC; nx<<-nx; ys<<-ys; normalize<<-normalize; X<<-X; nu<<-nu; nLayer<<-nLayer; datCoords<<-coords;
  # nBuffer<<-nBuffer; xRangeDat<<-xRangeDat; yRangeDat<<-yRangeDat; fastNormalize<<-fastNormalize; alpha2<<-priorPar$varPar2;
  # alpha1<<-priorPar$varPar1; corRange<<-priorPar$corScalePar
  ## generate inla stack:
  # Stacked A matrix (A_s from notation of LR2015 Bayesian Spatial Modelling with R_INLA):
  # (AEst   0  )
  # ( 0   APred)
  # eta_s = (c^T c^T)^T
  # where c is the vector of lattice coefficients
  AEst = makeA(coords, xNKnot=nx, yNKnot=ny, nLayer=nLayer, xRangeDat=xRangeDat, 
               yRangeDat=yRangeDat, nBuffer=nBuffer, NC=NC)
  APred = makeA(predPts, xNKnot=nx, yNKnot=ny, nLayer=nLayer, xRangeDat=xRangeDat, 
                yRangeDat=yRangeDat, nBuffer=nBuffer, NC=NC)
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
  control.inla = list(strategy=strategy, int.strategy=int.strategy) 
  mod = inla(y ~ - 1 + X + f(field, model=rgen), data=dat, 
             control.predictor=list(A=inla.stack.A(stack.full), compute=TRUE), 
             family="normal", verbose=TRUE, control.inla=control.inla)
  
  # get predictive surface, SD, and data
  index = inla.stack.index(stack.full, "pred")$data
  linpred.mean = mod[["summary.linear.predictor"]]$mean
  linpred.sd = mod[["summary.linear.predictor"]]$sd
  obsInds = 1:n
  predInds = (n+1):(n+nrow(predPts))
  preds = linpred.mean[predInds]
  predSDs = linpred.sd[predInds]
  obsPreds = linpred.mean[obsInds]
  obsSDs = linpred.sd[obsInds]
  
  startI = n+nrow(predPts)+1
  endI = n+nrow(predPts)+nx*ny
  coefPreds = list(layer1 = linpred.mean[startI:endI])
  coefSDs = list(layer1 = linpred.mean[startI:endI])
  for(i in 2:nLayer) {
    startI = endI + 1
    endI = startI + nrow(predPts)
    coefPreds = c(coefPreds, list(linpred.mean[startI:endI]))
    coefSDs = c(coefSDs, list(linpred.sd[startI:endI]))
  }
  
  # calculate true effective range and marginal variance:
  latticeWidth = (xRangeBasis[2] - xRangeBasis[1])/(nx-1)
  
  list(mod=mod, preds=preds, SDs=predSDs, latInfo=latInfo, latWidth=latticeWidth, obsPreds=obsPreds, 
       obsSDs=obsSDs, coefPreds=coefPreds, coefSDs=coefSDs)
}

# same as testLKModelSimple, but uses fitLKINLASimple as a subfunction
# buffer: buffer distance between domain edge of basis lattice and domain edge of data.
# n: number of observations
# xRange: range of x coordinates
# yRange: range of y coordinates
# nx: number of basis function lattice points in x directions
# NOTE: ny is determined automatically to match scale of x lattice points
# Xmat: design matrix
# ys: observations
# first.time: is first time evaluating function.  User should always set to FALSE
testLKModelSimple2 = function(buffer=2.5, kappa=1, rho=1, nu=1.5, seed=1, nLayer=3, nx=20, ny=nx, n=50, sigma2 = .1^2, 
                              nBuffer=5, normalize=TRUE, fastNormalize=TRUE, NC=5, testCovs=TRUE) {
  set.seed(seed)
  
  # compute alphas, the variance weights for each layer, depending on nu:
  alphas = getAlphas(nLayer, nu)
  
  # generate lattice and simulate observations
  coords = matrix(runif(2*n), ncol=2)
  xRangeDat = range(coords[,1])
  yRangeDat = range(coords[,2])
  latInfo = makeLatGrid(xRangeDat, yRangeDat, NC, nBuffer)
  
  xRangeBasis = latInfo$xRangeKnots
  yRangeBasis = latInfo$yRangeKnots
  nx = latInfo$nx
  ny = latInfo$ny
  
  AObs = makeA(coords, xNKnot=nx, yNKnot=ny, nLayer=nLayer, xRangeDat=xRangeDat, 
               yRangeDat=yRangeDat, nBuffer=nBuffer, NC=NC)
  Q = makeQ(kappa=kappa, rho=rho, nx=nx, ny=ny, 
            nLayer=nLayer, alphas=alphas, xRangeDat=xRangeDat, yRangeDat=yRangeDat, 
            nBuffer=nBuffer, normalized = normalize, NC=NC, fastNormalize=fastNormalize) 
  L = as.matrix(t(chol(solve(Q))))
  zsims = matrix(rnorm(nrow(Q)), ncol=1)
  fieldSims = L %*% zsims
  ys = as.numeric(AObs %*% fieldSims) + 1 # add a constant unit mean term to be estimated by INLA
  # ys = 1 + as.numeric(AObs %*% fieldSims) + coords[,1] # x-valued mean term to be estimated by INLA
  errs = rnorm(n, sd=sqrt(sigma2))
  ys = ys + errs
  
  # plot the observations
  pdf(file="Figures/simpleObservations.pdf", width=5, height=5)
  par(mfrow=c(1,1))
  quilt.plot(coords, ys)
  dev.off()
  
  # make prediction coordinates on a grid
  xRange=c(0,1)
  yRange=c(0,1)
  mx = 20
  my = 20
  predPts = make.surface.grid(list(x=seq(xRange[1], xRange[2], l=mx), y=seq(yRange[1], yRange[2], l=my)))
  
  # generate hyperparameters based on median and quantiles of inverse exponential and inverse gamma
  priorPar = getPrior(.1, .1, 10)
  X = matrix(rep(1, n), ncol=1)
  # X = matrix(coords[,1], ncol=1)
  XPred = matrix(rep(1, mx*my), ncol=1)
  
  # add linear terms in lat/lon to covariate matrices if requested
  if(testCovs) {
    X = cbind(X, coords)
    XPred = cbind(XPred, predPts)
  }
  
  # show priors on effective correlation and marginal variance:
  xs1 = seq(.01, 1, l=500)
  pdf(file="Figures/simplePriorEffRange.pdf", width=5, height=5)
  plot(xs1, dinvexp(xs1, rate=priorPar$corScalePar), type="l", col="blue", 
       xlab="Effective Correlation Range", main="Effective Correlation Prior", 
       ylab="Prior Density")
  abline(v=qinvexp(.5, rate=priorPar$corScalePar), col="red")
  dev.off()
  
  xs2 = seq(.01, 10.5, l=500)
  pdf(file="Figures/simplePriorMargVar.pdf", width=5, height=5)
  plot(xs2, dinvgamma(xs2, shape=priorPar$varPar1, rate=priorPar$varPar2), type="l", col="blue", 
       xlab="Marginal Variance", main="Marginal Variance Prior", 
       ylab="Prior Density")
  abline(v=qinvgamma(.1, shape=priorPar$varPar1, rate=priorPar$varPar2), col="red")
  abline(v=qinvgamma(.9, shape=priorPar$varPar1, rate=priorPar$varPar2), col="red")
  dev.off()
  
  # fit the model
  out = fitLKINLASimple(coords, ys, predPts=predPts, nu=nu, seed=seed, nLayer=nLayer, NC=NC,
                        nBuffer=nBuffer, priorPar=priorPar, X=X, XPred=XPred, normalize=normalize, 
                        int.strategy="auto", strategy="gaussian", fastNormalize=fastNormalize)
  mod = out$mod
  preds=out$preds
  predSDs=out$SDs
  latInfo=out$latInfo
  latWidth=out$latWidth
  obsPreds=out$obsPreds
  obsSDs=out$obsSDs
  coefPreds = out$coefPreds
  coefSDs = out$coefSDs
  
  # show a model summary
  print(summary(mod))
  
  # show predictive surface, SD, and data
  
  pdf(file="Figures/simplePreds.pdf", width=15, height=10)
  par(mfrow=c(2,3))
  # obsInds = 1:n
  # predInds = (n+1):(n+mx*my)
  # coefInds = (n+mx*my+1):(n+mx*my+nx*ny)
  colRangeDat = range(c(ys-errs, obsPreds, preds, coefPreds))
  colRangeSD = range(c(predSDs, obsSDs, coefSDs))
  gridPtsL1 = make.surface.grid(list(xs=seq(xRangeBasis[1], xRangeBasis[2], l=nx), 
                                     ys=seq(yRangeBasis[1], yRangeBasis[2], l=ny)))
  quilt.plot(coords, ys-errs, main="True Process", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
  quilt.plot(gridPtsL1[,1], gridPtsL1[,2], coefPreds[[1]], main="Basis Coefficient Mean (Layer 1)", xlim=xRangeDat, ylim=yRangeDat)
  quilt.plot(gridPtsL1[,1], gridPtsL1[,2], coefSDs[[1]], main="Basis Coefficient SD (Layer 1)", zlim=colRangeSD, xlim=xRangeDat, ylim=yRangeDat)
  
  quilt.plot(coords, obsPreds, main="Observation Mean", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
  quilt.plot(predPts[,1], predPts[,2], preds, main="Prediction Mean", zlim=colRangeDat, 
             xlim=xRangeDat, ylim=yRangeDat)
  quilt.plot(predPts[,1], predPts[,2], predSDs, main="Prediction SD", 
             xlim=xRangeDat, ylim=yRangeDat)
  dev.off()
  
  # calculate true effective range and marginal variance:
  latticeWidth = (xRangeBasis[2] - xRangeBasis[1])/(nx-1)
  effRange = 2.3/kappa * latticeWidth
  # marginalVar = rho/(4*pi * kappa^2)
  # marginalVar = getMultiMargVar(kappa, rho, nLayer=nLayer, nu=nu, xRange=xRangeBasis, 
  #                               yRange=yRangeBasis, nx=nx, ny=ny)[1]
  marginalVar = rho
  
  # plot marginals on interpretable scale (effective range, marginal variance)
  effRangeMarg = inla.tmarginal(exp, mod$marginals.hyperpar$`Theta1 for field`)
  varMarg = inla.tmarginal(exp, mod$marginals.hyperpar$`Theta2 for field`)
  sigma2Marg = inla.tmarginal(function(x) {1/x}, mod$marginals.hyperpar$`Precision for the Gaussian observations`)
  covNames = names(mod$marginals.fixed)
  XMarginals = list()
  for(i in 1:length(covNames)) {
    XMarginal = inla.tmarginal(function(x) {x}, mod$marginals.fixed[[covNames[i]]])
    XMarginals = c(XMarginals, list(XMarginal))
  }
  
  par(mfrow=c(1,1))
  pdf(file="Figures/simpleEffRange.pdf", width=5, height=5)
  plot(effRangeMarg, type="l", main="Marginal for effective range")
  abline(v=effRange, col="green")
  abline(v=inla.qmarginal(c(.025, .975), effRangeMarg), col="purple", lty=2)
  dev.off()
  # plot(mod$marginals.hyperpar$`Theta1 for field`, type="l", main="Marginal for log range")
  pdf(file="Figures/simpleVar.pdf", width=5, height=5)
  plot(varMarg, type="l", main="Marginal for process variance")
  abline(v=marginalVar, col="green")
  abline(v=inla.qmarginal(c(.025, .975), varMarg), col="purple", lty=2)
  dev.off()
  # plot(mod$marginals.hyperpar$`Theta2 for field`, type="l", main="Marginal for log variance")
  pdf(file="Figures/simpleSigma2.pdf", width=5, height=5)
  plot(sigma2Marg, type="l", main="Marginal for error variance")
  abline(v=sigma2, col="green")
  abline(v=inla.qmarginal(c(.025, .975), sigma2Marg), col="purple", lty=2)
  dev.off()
  for(i in 1:length(covNames)) {
    XMarginal = XMarginals[[i]]
    pdf(file=paste0("Figures/simple", covNames[i], ".pdf"), width=5, height=5)
    plot(XMarginal, type="l", main="Marginal for fixed effect")
    if(i==1)
      abline(v=1, col="green")
    else
      abline(v=0, col="green")
    abline(v=inla.qmarginal(c(.025, .975), XMarginal), col="purple", lty=2)
    dev.off()
  }
  
  # pdf(file="Figures/simpleRho.pdf", width=5, height=5)
  # plot(sigma2Marg, type="l", main=TeX("Marginal for $\\rho$"), xlab=TeX("$\\rho$"))
  # abline(v=rho, col="green")
  # dev.off()
  
  
  # do the same for kappa, rho
  # in order to get distribution for rho, must sample from joint hyperparameters
  kappaMarg = inla.tmarginal(function(x) {2.3/exp(x) * latticeWidth}, mod$marginals.hyperpar$`Theta1 for field`)
  # thetasToRho = function(xs) {
  #   logCor = xs[2]
  #   logVar = xs[3]
  #   kappa = 2.3/exp(logCor) * latticeWidth
  #   sigma2 = exp(logVar)
  #   sigma2 * 4*pi * kappa^2
  # }
  # samples = inla.hyperpar.sample(50000, mod, TRUE)
  # rhos = apply(samples, 1, thetasToRho)
  
  pdf(file="Figures/simpleKappa.pdf", width=5, height=5)
  plot(kappaMarg, type="l", xlab="kappa", main="Marginal for kappa")
  abline(v=kappa, col="green")
  abline(v=inla.qmarginal(c(.025, .975), kappaMarg), col="purple", lty=2)
  dev.off()
  
  # pdf(file="Figures/simpleRho.pdf", width=5, height=5)
  # hist(rhos, xlab="rho", main="Marginal for Rho", breaks=1000, freq=F, xlim=c(0, quantile(probs=.95, rhos)))
  # abline(v=rho, col="green")
  # dev.off()
  invisible(NULL)
}

# converts from LatticeKrig lattice parameters to LKinla's
# xRange: range of x-coords of data
# yRange: range of y-coords of data
# NC: number of basis functions to put on largest dimension of data
# nBuffer: number of lattice points outside of range of data to include as buffer
makeLatGrid = function(xRange, yRange, NC=5, nBuffer=5) {
  xLength = xRange[2] - xRange[1]
  yLength = yRange[2] - yRange[1]
  if(xLength >= yLength) {
    # lattice point xs are on edge of data domain, lattice point ys aren't
    latWidth = xLength/(NC-1)
    Buff = latWidth*nBuffer
    xRangeKnots = c(xRange[1] - Buff, xRange[2] + Buff)
    nx = NC + 2*nBuffer
    ny = floor(yLength/latWidth) + 2*nBuffer + 1
    yMidPt = mean(yRange)
    extraYLength = yLength - floor(yLength/latWidth)*latWidth
    yRangeKnots = c(yRange[1] + 0.5*extraYLength - Buff, yRange[2] - 0.5*extraYLength + Buff)
  }
  else {
    # lattice point ys are on edge of data domain, lattice point xs aren't
    latWidth = yLength/(NC-1)
    Buff = latWidth*nBuffer
    yRangeKnots = c(yRange[1] - Buff, yRange[2] + Buff)
    ny = NC + 2*nBuffer
    nx = floor(xLength/latWidth) + 2*nBuffer + 1
    xMidPt = mean(xRange)
    extraXLength = xLength - floor(xLength/latWidth)*latWidth
    xRangeKnots = c(xRange[1] + 0.5*extraXLength - Buff, xRange[2] - 0.5*extraXLength + Buff)
  }
  
  list(xRangeKnots=xRangeKnots, nx=nx, yRangeKnots=yRangeKnots, ny=ny)
}

testMakeLatGrid = function(xRange=c(0,1), yRange=c(0,1)) {
  test = makeLatGrid(xRange, yRange)
  xs = seq(test$xRangeKnots[1], test$xRangeKnots[2], l=test$nx)
  ys = seq(test$yRangeKnots[1], test$yRangeKnots[2], l=test$ny)
  pts = make.surface.grid((list(xs, ys)))
  border = rbind(c(xRange[1], yRange[1]), 
              c(xRange[1], yRange[2]), 
              c(xRange[2], yRange[2]), 
              c(xRange[2], yRange[1]))
  plot(pts[,1], pts[,2], xlab="x", ylab="y", main="Data domain and buffer")
  polygon(border[,1], border[,2], border="blue")
}
# testMakeLatGrid()
# testMakeLatGrid(c(0,1.2))

# generate lattice points for all layers
# NOTE: if any of the possibly NULL input variables are NULL, they are all set to the default
makeLatGrids = function(xRangeDat=c(0,1), yRangeDat=c(0,1), NC=5, nBuffer=5, nLayer=1) {
  
  # # set lattice parameters to defaults if necessary
  # if(is.null(xRangeKnot) || is.null(yRangeKnot) || is.null(xNKnot) || is.null(yNKnot)) {
  #   latGrid = makeLatGrid(xRangeDat, yRangeDat, NC, nBuffer)
  #   xRangeKnot = latGrid$xRangeKnots
  #   xNKnot=latGrid$nx
  #   yRangeKnot = latGrid$yRangeKnots
  #   yNKnot = latGrid$ny
  # }
  # origXNKnot = xNKnot
  # origYNKnot = yNKnot
  
  latCoords = list()
  nx = c()
  ny = c()
  ms = c()
  for(thisLayer in 1:nLayer) {
    i = thisLayer
    # # adjust basis elements to depend on layer
    # xNKnot = origXNKnot * 2^(thisLayer-1)
    # yNKnot = origYNKnot * 2^(thisLayer-1)
    # 
    # # generate knot lattice locations and filter out locations 
    # # too far outside of the data domain
    # knotXs = seq(xRangeKnot[1], xRangeKnot[2], l=xNKnot)
    # knotYs = seq(yRangeKnot[1], yRangeKnot[2], l=yNKnot)
    # if(sum(knotXs > xRangeDat[2]) > nBuffer)
    #   knotXs = knotXs[1:(length(knotXs) - (sum(knotXs > xRangeDat[2]) - nBuffer))]
    # if(sum(knotXs < xRangeDat[1]) > nBuffer)
    #   knotXs = knotXs[(1 + sum(knotXs < xRangeDat[1]) - nBuffer):length(knotXs)]
    # if(sum(knotYs > yRangeDat[2]) > nBuffer)
    #   knotYs = knotYs[1:(length(knotYs) - (sum(knotYs > yRangeDat[2]) - nBuffer))]
    # if(sum(knotYs < yRangeDat[1]) > nBuffer)
    #   knotYs = knotYs[(1 + sum(knotYs < yRangeDat[1]) - nBuffer):length(knotYs)]
    # knotPts = make.surface.grid(list(x=knotXs, y=knotYs))
    layerNC = (NC-1)*2^(i-1) + 1
    latGrid = makeLatGrid(xRangeDat, yRangeDat, layerNC, nBuffer)
    xRangeKnot = latGrid$xRangeKnots
    xNKnot=latGrid$nx
    yRangeKnot = latGrid$yRangeKnots
    yNKnot = latGrid$ny
    xGrid = seq(xRangeKnot[1], xRangeKnot[2], l=xNKnot)
    yGrid = seq(yRangeKnot[1], yRangeKnot[2], l=yNKnot)
    knotPts = make.surface.grid(list(x=xGrid, y=yGrid))
    
    latCoords = c(latCoords, list(knotPts))
    nx = c(nx, xNKnot)
    ny = c(ny, yNKnot)
    ms = nx*ny
  }
  list(latCoords=latCoords, nx=nx, ny=ny, ms=ms)
}

# compute number of basis elements per layer.  Returns a list with element at index i 
# equal to the number of basis elements in layer i.
getMs = function(xRangeDat=c(0,1), yRangeDat=c(0,1), NC=5, nBuffer=5, nLayer=1) {
  makeLatGrids(xRangeDat=xRangeDat, yRangeDat=yRangeDat, NC=NC, nBuffer=nBuffer, nLayer=nLayer)$ms
}

# convert from explicit grid parameters to LatticeKrig grid parameters
rawGridToLK = function(xRangeKnot=c(0,1), xNKnot=10, yRangeKnot=c(0,1), yNKnot=10, nBuffer=5) {
  # first compute NC, the number of lattice points in the largest dimension
  xNKnotInRange = xNKnot - 2*nBuffer
  yNKnotInRange = yNKnot - 2*nBuffer
  NC = max(xNKnotInRange, yNKnotInRange)
  
  # subtract off the buffer to get the "range of the data" under the LatticeKrig parameterization
  latWidth = diff(xRangeKnot)/(xNKnot-1)
  xGrid = seq(xRangeKnot[1], xRangeKnot[2], l=xNKnot)
  yGrid = seq(yRangeKnot[1], yRangeKnot[2], l=yNKnot)
  xRangeDat = c(xGrid[nBuffer+1], xGrid[xNKnot - nBuffer])
  yRangeDat = c(yGrid[nBuffer+1], yGrid[yNKnot - nBuffer])
  list(NC=NC, nBuffer=nBuffer, xRangeDat=xRangeDat, yRangeDat=yRangeDat)
}


# exiting from: inla.rgeneric.lk.model.simple("Q", thetas)
# Browse[1]> dim(test)
# [1] 533 533
# Browse[1]> 210 + 342
# [1] 552
# Browse[1]> 210 + 342 + 675
# [1] 1227
# Browse[1]> xRange
# [1] -1.209754  2.215051
# Browse[1]> xRangeDat
# [1] 0.01339033 0.99190609
# Browse[1]> xRangeBasis
# [1] -1.209754  2.215051
# Browse[1]> yRange
# [1] -1.080312  2.099864
# Browse[1]> yRangeDat
# [1] 0.05893438 0.96061800