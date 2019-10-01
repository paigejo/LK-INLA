##### make "simple" (alphas follow fixed relationship based on nu) model for multi-layer LK-INLA:

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
  envir = environment(sys.call()[[1]])
  
  # load relevant external functions
  require(Matrix)
  require(spam)
  require(fields)
  require(LatticeKrig)
  require(invgamma)
  require(INLA)
  source("testCasePt2.R")
  
  # get range of data
  xRangeDat = latInfo[[1]]$xRangeDat
  yRangeDat = latInfo[[1]]$yRangeDat
  xRange = latInfo[[1]]$xRange
  yRange = latInfo[[1]]$yRange
  nx = latInfo[[1]]$nx
  ny = latInfo[[1]]$ny
  xLength = latInfo[[1]]$xLength
  yLength = latInfo[[1]]$yLength
  knotPts = latInfo[[1]]$knotPts
  latWidth = latInfo[[1]]$latWidth
  nBuffer = latInfo[[1]]$nBuffer
  NC = latInfo[[1]]$NC
  n = nrow(X)
  nLayer = length(latInfo)
  
  # # convert to theoretical kappa estimate
  # kappaEst = 2.3/effectiveRangeInit * latWidth
  # 
  # # now to the important part: get an initial estimate of marginal variance and a range (go with /100 and *100 of initial guess of kappa)
  # kappas <<- 10^(seq(log10(kappaEst)-2, log10(kappaEst)+2, l=100))
  # kappaWidths <<- log10(kappas[2]) - log10(kappas[1])
  # margVars <<- sapply(kappas, getMultiMargVar, rho=1, tod=2.5, nx=nx, ny=ny, nu=nu, nLayer=nLayer)[1,]
  # logMargVarSpline <<- splinefun(log(kappas), margVars, "natural")
  
  # theta is of the form:
  # c(betas, effectiveCor, sigmaSq, kappa, rho, nu, alphas)
  interpret.theta = function()
  {
    ## internal helper-function to map the parameters from the internal-scale to the
    ## user-scale
    if(!is.null(theta)) {
      effectiveCor = exp(theta[1])
      sigmaSq = exp(theta[2])
    }
    else {
      effectiveCor = NULL
      sigmaSq = NULL
    }
    
    # precomputations: get lattice grid cell width, convert parameters from effective correlation 
    # and marginal variance to kappa and rho.  Use spline to convert from marginal variance to kappa
    latticeWidth = latWidth
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
    L = nLayer = length(latInfo)
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
    makeGraph(latInfo)
  }
  
  # compute the precision matrix
  Q = function() {
    makeQ(kap, rho, latInfo, alphas=alphas, normalized=normalize, 
          fastNormalize=fastNormalize)
  }
  
  # get mean of each latent coefficient
  mu = function() {
    return(numeric(0))
  }
  
  log.norm.const = function()
  {
    ## let INLA compute it as -n/2 log(2pi) + 1/2 * log(|Q|)
    return (numeric(0))
  }
  
  log.prior = function() { ##### TODO: is nu fit or fixed?
    require(invgamma)
    
    # get prior (note the Jacobian factors)
    if(prior$priorType == "orig") {
      dinvexp(effectiveCor, rate=prior$corScalePar, log=TRUE) + log(effectiveCor) + 
        invgamma::dinvgamma(sigmaSq, shape=prior$varPar1, rate=prior$varPar2, log=TRUE) + log(sigmaSq)
    } else if(prior$priorType == "pc") {
      dinvexp(effectiveCor, rate=prior$corScalePar, log=TRUE) + log(effectiveCor) + 
        dexp(sqrt(sigmaSq), rate=prior$sdRate, log=TRUE) - log(2) - log(sqrt(sigmaSq))
    }
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
    AMat = makeA(datCoords, latInfo)
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
  
  if(is.null(theta)) { theta = initial() }
  val = do.call(match.arg(cmd), args = list())
  return (val)
}

## function for making prior on corelation scale and marginal variance
# corScaleMed: median correlation scale
# var10quant: 10th percentile of the marginal variance
# var90quant: 90th percentile of the marginal variance
getPrior = function(corScaleMed, var10quant, var90quant, dirichletConcentration=1, nLayer=3) {
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
  
  ## set the layer weight hyperparameters
  alphaPar = rep(dirichletConcentration / nLayer, nLayer)
  
  list(corScalePar=corScalePar, varPar1=varPar1, varPar2=varPar2, alphaPar=alphaPar, priorType="orig")
}

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
                           intStrategy="auto", strategy="gaussian", fastNormalize=FALSE, 
                           predictionType=c("median", "mean")) {
  set.seed(seed)
  
  # get the type of prediction the user wants
  predictionType = match.arg(predictionType)
  
  # get lattice points, prediction points
  xRangeDat = range(coords[,1])
  yRangeDat = range(coords[,2])
  latInfo = makeLatGrids(xRangeDat, yRangeDat, NC, nBuffer, nLayer)
  
  # generate lattice basis matrix
  AObs = makeA(coords, latInfo)
  
  # define the LK model
  n = length(ys)
  
  rgen = inla.rgeneric.define(model=inla.rgeneric.lk.model.simple, latInfo=latInfo, 
                              ys=ys, prior=priorPar, normalize=normalize, 
                              X=X, nu=nu, datCoords=coords, fastNormalize=fastNormalize)
  # use these global variables for testing calls to inla.rgeneric.lk.model.simple
  # ys<<-ys; normalize<<-normalize; X<<-X; nu<<-nu; datCoords<<-coords
  # fastNormalize<<-fastNormalize; alpha2<<-priorPar$varPar2; prior<<-priorPar
  # alpha1<<-priorPar$varPar1; corRange<<-priorPar$corScalePar; latInfo<<-latInfo
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
  controls = list(strategy=strategy, int.strategy=intStrategy) 
  mod = inla(y ~ - 1 + X + f(field, model=rgen), data=dat, 
             control.predictor=list(A=inla.stack.A(stack.full), compute=TRUE), 
             family="normal", verbose=TRUE, control.inla=controls)
  
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
  endI = n+nrow(predPts)+nrow(latInfo[[1]]$latCoords)
  coefPreds = list(layer1 = linpreds[startI:endI])
  coefSDs = list(layer1 = linpred.sd[startI:endI])
  for(i in 2:nLayer) {
    startI = endI + 1
    endI = startI + nrow(predPts)
    coefPreds = c(coefPreds, list(linpreds[startI:endI]))
    coefSDs = c(coefSDs, list(linpred.sd[startI:endI]))
  }
  
  # get true effective range and marginal variance:
  latticeWidth = latInfo[[1]]$latWidth
  
  list(mod=mod, preds=preds, SDs=predSDs, latInfo=latInfo, latWidth=latticeWidth, obsPreds=obsPreds, 
       obsSDs=obsSDs, coefPreds=coefPreds, coefSDs=coefSDs)
}

# buffer: buffer distance between domain edge of basis lattice and domain edge of data.
# n: number of observations
# xRange: range of x coordinates
# yRange: range of y coordinates
# nx: number of basis function lattice points in x data range for the first layer
# NOTE: ny is determined automatically to match scale of x lattice points
# Xmat: design matrix
# ys: observations
testLKModelSimple = function(buffer=2.5, kappa=1, rho=1, nu=1.5, seed=1, nLayer=3, n=50, sigma2 = .1^2, 
                              nBuffer=5, normalize=TRUE, fastNormalize=TRUE, NC=5, testCovs=TRUE) {
  set.seed(seed)
  
  # compute alphas, the variance weights for each layer, depending on nu:
  alphas = getAlphas(nLayer, nu)
  
  # generate lattice and simulate observations
  coords = matrix(runif(2*n), ncol=2)
  xRangeDat = range(coords[,1])
  yRangeDat = range(coords[,2])
  latInfo = makeLatGrids(xRangeDat, yRangeDat, NC, nBuffer, nLayer)
  
  AObs = makeA(coords, latInfo)
  Q = makeQ(kappa=kappa, rho=rho, latInfo, alphas=alphas, normalized=normalize, fastNormalize=fastNormalize) 
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
  
  # # show priors on effective correlation and marginal variance:
  # xs1 = seq(.01, 1, l=500)
  # pdf(file="Figures/simplePriorEffRange.pdf", width=5, height=5)
  # plot(xs1, dinvexp(xs1, rate=priorPar$corScalePar), type="l", col="blue", 
  #      xlab="Effective Correlation Range", main="Effective Correlation Prior", 
  #      ylab="Prior Density")
  # abline(v=qinvexp(.5, rate=priorPar$corScalePar), col="red")
  # dev.off()
  # 
  # xs2 = seq(.01, 10.5, l=500)
  # pdf(file="Figures/simplePriorMargVar.pdf", width=5, height=5)
  # plot(xs2, invgamma::dinvgamma(xs2, shape=priorPar$varPar1, rate=priorPar$varPar2), type="l", col="blue", 
  #      xlab="Marginal Variance", main="Marginal Variance Prior", 
  #      ylab="Prior Density")
  # abline(v=qinvgamma(.1, shape=priorPar$varPar1, rate=priorPar$varPar2), col="red")
  # abline(v=qinvgamma(.9, shape=priorPar$varPar1, rate=priorPar$varPar2), col="red")
  # dev.off()
  
  # fit the model
  out = fitLKINLASimple(coords, ys, predPts=predPts, nu=nu, seed=seed, nLayer=nLayer, NC=NC,
                        nBuffer=nBuffer, priorPar=priorPar, X=X, XPred=XPred, normalize=normalize, 
                        intStrategy="auto", strategy="gaussian", fastNormalize=fastNormalize)
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
  gridPtsL1 = latInfo[[1]]$latCoords
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
  latticeWidth = latInfo[[1]]$latWidth
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

testLKModelSimple()