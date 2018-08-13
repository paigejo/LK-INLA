library(Matrix)
library(spam)
library(fields)
library(LatticeKrig)
library(invgamma)
library(INLA)
library(latex2exp)
library(xtable)
# setwd("~/git/LK-INLA/")
# source('~/git/M9/exploratoryAnalysisFuns.R')

## code for calculating precision and basis matrices for coding LatticeKrig into INLA:
library(profvis)
# profvis({
  # library(fields)
# computing precision matrix for a single layer or a block diagonal sparse 
# precision matrix for multiple layers
# kappa: scale of Matern covariance with smoothness 1
# xRange, yRange: x and y intervals in space over which basis elements are placed
# nx, ny: number of basis elements when counting along the lattive in 
#         x and y directions respectively (for the first layer)
# rho: sill (marginal variance) for first layer
# nLayer: number of lattice layers
# thisLayer: user should always set to 1 to get full block diagonal precision matrix 
#            for all layers
# alphas: weights on the variances of each layer.  Scales rho for each layer.
# fastNormalize: simple normalization to make marginal variance = rho in center. Basis coefficients may have different variances
# assume there is a single kappa, rho is adjusted with alpha when there are multiple layers
makeQ = function(kappa=1, rho=1, xRange=c(0,1), yRange=c(0,1), nx=10, ny=10, nLayer=1, thisLayer=1, alphas=NULL, nu=NULL, 
                 xRangeDat=xRange, yRangeDat=yRange, nBuffer=5, normalized=FALSE, fastNormalize=FALSE, NC=5, rawKnotRange=FALSE) {
  
  # save base layer input quantities
  origRho = rho
  origNx = nx
  origNy = ny
  
  # make alphas according to nu relation if nu is set
  if(is.null(alphas) && !is.null(nu)) {
    alphas = getAlphas(nLayer, nu)
  }
  
  # # adjusted resolution and process params depending on layer
  # nx = nx*2^(thisLayer-1)
  # ny = ny*2^(thisLayer-1)
  # 
  # # generate knot lattice locations and filter out locations 
  # # too far outside of the data domain
  # knotXs = seq(xRange[1], xRange[2], l=nx)
  # knotYs = seq(yRange[1], yRange[2], l=ny)
  # if(sum(knotXs > xRangeDat[2]) > nBuffer)
  #   knotXs = knotXs[1:(length(knotXs) - (sum(knotXs > xRangeDat[2]) - nBuffer))]
  # if(sum(knotXs < xRangeDat[1]) > nBuffer)
  #   knotXs = knotXs[(1 + sum(knotXs < xRangeDat[1]) - nBuffer):length(knotXs)]
  # if(sum(knotYs > yRangeDat[2]) > nBuffer)
  #   knotYs = knotYs[1:(length(knotYs) - (sum(knotYs > yRangeDat[2]) - nBuffer))]
  # if(sum(knotYs < yRangeDat[1]) > nBuffer)
  #   knotYs = knotYs[(1 + sum(knotYs < yRangeDat[1]) - nBuffer):length(knotYs)]
  # gridPts = make.surface.grid(list(x=knotXs, y=knotYs))
  # nx = length(knotXs)
  # ny = length(knotYs)
  
  # make range of knots based on data range
  if(!rawKnotRange) {
    out = makeLatGrid(xRangeDat, yRangeDat, NC=NC, nBuffer=nBuffer)
    xRange=out$xRangeKnots
    yRange=out$yRangeKnots
  }
  else {
    # make sure the raw grid has equal cell width in x and y directions
    xWidth = diff(xRange)/(nx-1)
    yWidth = diff(yRange)/(ny-1)
    if(xWidth != yWidth)
      stop("raw x and y grid widths not equal")
  }
  
  # convert from raw grid parameterization to LatticeKrig grid parameterization, 
  # and generate the knot points for this layer
  gridPar = makeLatGrids(xRangeDat, yRangeDat, NC, nBuffer, nLayer)
  knotPts = gridPar$latCoords[[thisLayer]]
  nx = gridPar$nx[thisLayer]
  ny = gridPar$ny[thisLayer]
  
  # make (Laplacian) differential operators
  Dnx = bandSparse(nx, k=0:1, diag=list(rep(-2, nx), rep(1, nx-1)), symmetric=TRUE)
  Dny = bandSparse(ny, k=0:1, diag=list(rep(-2, ny), rep(1, ny-1)), symmetric=TRUE)
  
  # generate x and y (Laplacian) differential operators
  Inx = Diagonal(n=nx)
  Iny = Diagonal(n=ny)
  Bx = kronecker(Iny, Dnx)
  By = kronecker(Dny, Inx)
  
  # make B, SAR regression matrix for Bc = e
  B = Diagonal(n=nx*ny, x=kappa^2) - (Bx + By)
  
  # compute precision matrix
  Q = (1/rho) * t(B) %*% B
  
  if(normalized) {
    # make mid lattice point
    # xi = ceiling(nx/2)
    # yi = ceiling(ny/2)
    # midPt = matrix(c(seq(xRange[1], xRange[2], l=nx)[xi], seq(yRange[1], yRange[2], l=ny)[yi]), nrow=1)
    midPt = matrix(c((xRange[1] + xRange[2])/2, (yRange[1] + yRange[2])/2), nrow=1)
    
    # now construct relevant row of A for value at midpoint at this layer, and Q matrix
    Ai = makeA(midPt, xNKnot=origNx, yNKnot=origNy, thisLayer=thisLayer, nLayer=thisLayer, nBuffer=nBuffer, xRangeDat=xRangeDat, 
               yRangeDat=yRangeDat, NC=NC)
    
    # # test
    # sds2 = 1/diag(Q)
    # sdMat = Diagonal(x=sds2)
    # # Qnorm2 = sweep(sweep(Q, 1, sds2, "*"), 2, sds2, "*")
    # Qnorm2 = sdMat %*% Q %*% sdMat
    # 
    # # QnormInv = sweep(sweep(Qinv, 1, 1/sds), 2, 1/sds)
    # procVar2 = as.numeric(Ai %*% inla.qsolve(Qnorm2, t(Ai)))
    # # procVar = as.numeric(Ai %*% QnormInv %*% t(Ai))
    # Q2 = Qnorm2 * (procVar2 / rho)
    # # Q = Q2 # system.time(out <- makeQ(nLayer=3, nx=15, ny=15, nu=1, normalized=TRUE, newnormalize=TRUE)): ~5.8s

    #  test 2
    if(fastNormalize) {
      ctilde = as.numeric(Ai %*% inla.qsolve(Q, t(Ai)))/rho
      Qtilde = ctilde * Q
      Q = Qtilde # system.time(out <- makeQ(nLayer=3, nx=15, ny=15, nu=1, normalized=TRUE, newnormalize=TRUE)): 2.72
    }
    else {
      # renormalize basis coefficients to have constant variance, and the process to have unit variance
      Qinv = inla.qsolve(Q, diag(nrow(Q)))
      sds = sqrt(diag(Qinv))
      sdMat = Diagonal(x=sds)
      # Qnorm = sweep(sweep(Q, 1, sds, "*"), 2, sds, "*")
      Qnorm = sdMat %*% Q %*% sdMat
      # QnormInv = sweep(sweep(Qinv, 1, 1/sds), 2, 1/sds)
      procVar = as.numeric(Ai %*% inla.qsolve(Qnorm, t(Ai)))
      # procVar = as.numeric(Ai %*% QnormInv %*% t(Ai))
      Q = Qnorm * (procVar / rho) # system.time(out <- makeQ(nLayer=3, nx=15, ny=15, nu=1, normalized=TRUE) ~5.87
    }
    # # compute how good an approximation it was
    # hist(diag(solve(Q)))
    # hist(diag(solve(Q2)))
    # hist(diag(solve(Qtilde)))
    # image(Q)
    # image(Q2)
    # image(Qtilde)
  }
  if(nLayer == 1 && is.null(alphas))
    alphas = 1
  Q = Q * (1/alphas[thisLayer])
  
  # return results
  # If multiple layers, return block diagonal sparse matrix
  if(thisLayer == nLayer) {
    return(Q)
  }
  else if((thisLayer == 1) && (nLayer != 1)) {
    Q = bdiag(c(list(Q), makeQ(kappa, origRho, xRange, yRange, origNx, origNy, nLayer, thisLayer+1, alphas, 
                               xRangeDat=xRangeDat, yRangeDat=yRangeDat, nBuffer=nBuffer, normalized=normalized, 
                               fastNormalize=fastNormalize, NC=NC)))
    return(Q)
  }
  else {
    return(c(list(Q), makeQ(kappa, origRho, xRange, yRange, origNx, origNy, nLayer, thisLayer+1, alphas, 
                            xRangeDat=xRangeDat, yRangeDat=yRangeDat, nBuffer=nBuffer, normalized=normalized,
                            fastNormalize=fastNormalize, NC=NC)))
  }
}

# test for the makeQ function.  Same inputs as makeQ
testMakeQ = function(kappa=1, rho=1, xRange=c(0,1), yRange=c(0,1), nx=40, ny=40, nLayer=1, nu=1.5) {
  
  if(nLayer != 1) {
    alphas = getAlphas(nLayer, nu)
  }
  else
    alphas = NULL
  
  # make precision matrix
  Q = makeQ(kappa, rho, xRange, yRange, nx, ny, nLayer, alphas=alphas)
  
  # simulate field
  L = t(chol(inla.qsolve(Q, diag(nrow(Q)))))
  Zs = matrix(rnorm(nrow(Q)), ncol=1)
  Cs = matrix(L %*% Zs, nrow=ny, ncol=nx, byrow = TRUE)
  Cs
}
# image.plot(testMakeQ(kappa=.1, rho=5))

testMakeQBuffer = function(kappa=1, rho=1, xRange=c(-1,2), yRange=xRange, nx=20, ny=20, nu=1.5, 
                           xRangeDat=c(0,1), yRangeDat=xRangeDat, nBuffer=5, savePlot=FALSE) {
  nLayer=3
  
  if(nLayer != 1) {
    alphas = getAlphas(nLayer, nu)
  }
  else
    alphas = NULL
  
  # make precision matrix
  Q = makeQ(kappa, rho, xRange, yRange, nx, ny, nLayer, alphas=alphas, 
            xRangeDat=xRangeDat, yRangeDat=yRangeDat, nBuffer=nBuffer)
  
  # simulate field
  L = t(chol(inla.qsolve(Q, diag(nrow(Q)))))
  Zs = matrix(rnorm(nrow(Q)), ncol=1)
  Cs = L %*% Zs
  
  ## plot the coefficients of each layer
  if(savePlot)
    pdf(file="Figures/exampleSim.pdf", width=5, height=5)
  par(mfrow=c(2,2))
  startI = 1
  for(l in 1:nLayer) {
    thisnx = nx * 2^(l-1)
    thisny = ny * 2^(l-1)
    
    # generate knot lattice locations and filter out locations 
    # too far outside of the data domain
    knotXs = seq(xRange[1], xRange[2], l=thisnx)
    knotYs = seq(yRange[1], yRange[2], l=thisny)
    out = rawGridToLK(xRange, nx, yRange, ny, nBuffer)
    gridPar = makeLatGrids(out$xRangeDat, out$yRangeDat, out$NC, out$nBuffer, nLayer)
    gridPts = gridPar$latCoords[[l]]
    
    # filter out simulated coefficients for the layer
    layerL = nrow(gridPts)
    endI = startI + layerL - 1
    layerCs = Cs[startI:endI]
    
    quilt.plot(gridPts, layerCs, main=paste0("Layer ", l, " coefficients"), 
               xlim=xRange, ylim=yRange, xlab="x", ylab="y")
    
    startI = endI+1
  }
  if(savePlot)
    dev.off()
  
  ## plot the coefficient location in each layer together
  cols = rainbow(nLayer)
  for(l in 1:nLayer) {
    thisnx = nx * 2^(l-1)
    thisny = ny * 2^(l-1)
    
    # generate knot lattice locations and filter out locations 
    # too far outside of the data domain
    knotXs = seq(xRange[1], xRange[2], l=thisnx)
    knotYs = seq(yRange[1], yRange[2], l=thisny)
    out = rawGridToLK(xRange, nx, yRange, ny, nBuffer)
    gridPar = makeLatGrids(out$xRangeDat, out$yRangeDat, out$NC, out$nBuffer, nLayer)
    gridPts = gridPar$latCoords[[l]]
    
    if(l == 1) {
      plot(gridPts[,1], gridPts[,2], main="Lattice points", col=cols[l], pch=19, cex=1.5-l*.45, 
           xlim=xRange, ylim=yRange, xlab="x", ylab="y")
    }
    else {
      points(gridPts[,1], gridPts[,2], col=cols[l], pch=19, cex=1.5-l*.45, 
             xlim=xRange, ylim=yRange)
    }
  }
}

# test normalization of Q for 1 layer
testMakeNormalizedQ = function(buffer=1, kappa=1, rho=1, nu=1.5, seed=123, nx=10, ny=nx, mx=20, my=mx, sigma2 = .1^2) {
  nLayer=1
  
  if(nLayer != 1) {
    alphas = getAlphas(nLayer, nu)
  }
  else
    alphas = NULL
  
  # get lattice points, prediction points
  xRangeBasis = c(0-buffer, 1+buffer)
  yRangeBasis = c(0-buffer, 1+buffer)
  xRangeDat = c(0,1)
  yRangeDat = c(0,1)
  xs = seq(xRangeBasis[1], xRangeBasis[2], l=nx)
  ys = seq(yRangeBasis[1], yRangeBasis[2], l=ny)
  latPts = make.surface.grid(list(x=xs, y=ys))
  xmask =  (latPts[,1] >= xRangeDat[1]) & (latPts[,1] <= xRangeDat[2])
  ymask =  (latPts[,2] >= yRangeDat[1]) & (latPts[,2] <= yRangeDat[2])
  mask = xmask & ymask
  latPtsInDom = latPts[mask,]
  
  # generate lattice and simulate observations
  coords = make.surface.grid(list(x=seq(xRangeDat[1], xRangeDat[2], l=mx), y=seq(yRangeDat[1], yRangeDat[2], l=my)))
  AObs = makeA(coords, xRangeBasis, nx, yRangeBasis, ny, nLayer=nLayer)
  Q = makeQ(kappa=kappa, rho=1, xRange=xRangeBasis, yRange=yRangeBasis, nx=nx, ny=ny, nLayer=nLayer, alphas=alphas)
  Qinv = inla.qsolve(Q, diag(nrow(Q)))
  
  # make mid lattice point
  xi = ceiling(nx/2)
  yi = ceiling(ny/2)
  midPt = matrix(c(seq(xRangeBasis[1], xRangeBasis[2], l=nx)[xi], seq(yRangeBasis[1], yRangeBasis[2], l=ny)[yi]), nrow=1)
  
  # now construct relevant row of A for value at midpoint, and Q matrix
  Ai = makeA(midPt, xRange, nx, yRange, ny, nLayer=nLayer)
  
  # renormalize basis coefficients to have constant variance, and the process to have unit variance
  sds = sqrt(diag(Qinv))
  Qnorm = sweep(sweep(Q, 1, sds, "*"), 2, sds, "*")
  procVar = as.numeric(Ai %*% inla.qsolve(Qnorm, t(Ai)))
  Qfinal = Qnorm * (procVar / rho)
  
  ## plot simulations
  # simulate lattice coefficients
  L = t(chol(Qinv))
  Zs = matrix(rnorm(nrow(Q)), ncol=1)
  Cs = L %*% Zs
  CsMat = matrix(Cs, nrow=ny, ncol=nx, byrow = TRUE)
  image.plot(CsMat, main="Simulated lattice coefficients")
  
  # now simulate the continuous process:
  image.plot(matrix(AObs %*% Cs, nrow=my, ncol=mx, byrow=TRUE), main="Simulated Process")
  
  ## plot normalized and unnormalized coefficient variances
  # unnormalized coefficients
  vars = diag(Qinv)
  quilt.plot(latPtsInDom, vars[mask], main="Marginal variance of unnormalized coefficients", nx=20, ny=20)
  
  # normalized coefficeints
  vars = diag(inla.qinv(Qnorm))
  quilt.plot(latPtsInDom, vars[mask] + jitter(rep(0, sum(mask)), amount=.000001), main="Marginal variance of normalized coefficients", nx=20, ny=20)
  
  # final coefficients
  vars = diag(inla.qinv(Qfinal))
  quilt.plot(latPtsInDom, vars[mask] + jitter(rep(0, sum(mask)), amount=.000001), main="Marginal variance of final coefficients", nx=20, ny=20)
  
  ## plot the normalized and unnormalized process variances
  # unnormalized
  vars = diag(AObs %*% inla.qsolve(Q, t(AObs)))
  quilt.plot(coords, vars, main="Marginal variance of unnormalized process")
  
  # normalized
  vars = diag(AObs %*% inla.qsolve(Qnorm, t(AObs)))
  quilt.plot(coords, vars, main="Marginal variance of normalized process")
  
  # final
  vars = diag(AObs %*% inla.qsolve(Qfinal, t(AObs)))
  quilt.plot(coords, vars, main="Marginal variance of final process")
}

# test normalization of Q for multiple layers
testMakeNormalizedQMulti = function(buffer=1, kappa=1, rho=1, nu=1.5, seed=123, nx=10, ny=nx, mx=20, my=mx, 
                                    sigma2 = .1^2, nLayer=3, nBuffer=5, plotMinMax=FALSE) {
  
  if(nLayer != 1) {
    alphas = getAlphas(nLayer, nu)
  }
  else
    alphas = NULL
  
  # get lattice points, prediction points
  xRangeBasis = c(0-buffer, 1+buffer)
  yRangeBasis = c(0-buffer, 1+buffer)
  xRangeDat = c(0,1)
  yRangeDat = c(0,1)
  xs = seq(xRangeBasis[1], xRangeBasis[2], l=nx)
  ys = seq(yRangeBasis[1], yRangeBasis[2], l=ny)
  latPts = make.surface.grid(list(x=xs, y=ys))
  xmask =  (latPts[,1] >= xRangeDat[1]) & (latPts[,1] <= xRangeDat[2])
  ymask =  (latPts[,2] >= yRangeDat[1]) & (latPts[,2] <= yRangeDat[2])
  mask = xmask & ymask
  latPtsInDom = latPts[mask,]
  
  # generate lattice and simulate observations
  coords = make.surface.grid(list(x=seq(xRangeDat[1], xRangeDat[2], l=mx), y=seq(yRangeDat[1], yRangeDat[2], l=my)))
  AObs = makeA(coords, xRangeBasis, nx, yRangeBasis, ny, nLayer=nLayer, xRangeDat=xRangeDat, yRangeDat=yRangeDat, nBuffer=nBuffer)
  Q = makeQ(kappa=kappa, rho=1, xRange=xRangeBasis, yRange=yRangeBasis, nx=nx, ny=ny, nLayer=nLayer, 
            alphas=alphas, nu=nu, xRangeDat=xRangeDat, yRangeDat=yRangeDat, nBuffer=nBuffer, normalized=FALSE)
  Qfinal = makeQ(kappa=kappa, rho=1, xRange=xRangeBasis, yRange=yRangeBasis, nx=nx, ny=ny, nLayer=nLayer, 
                 alphas=alphas, nu=nu, xRangeDat=xRangeDat, yRangeDat=yRangeDat, nBuffer=nBuffer, normalized=TRUE)
  
  ## plot the normalized and unnormalized process variances
  # unnormalized
  vars = diag(AObs %*% inla.qsolve(Q, t(AObs)))
  pdf(file="Figures/varRatioUnnormalized.pdf", width=5, height=5)
  quilt.plot(coords, vars, main="Marginal variance of unnormalized process")
  maxI = which.max(vars)
  minI = which.min(vars)
  if(plotMinMax) {
    points(coords[maxI,1], coords[maxI, 2], cex=2)
    points(coords[minI,1], coords[minI, 2], cex=2)
  }
  unnormRatio = max(vars)/min(vars)
  dev.off()
  print(paste0("Ratio of largest to smallest variance for unnormalized process: ", unnormRatio))
  
  # normalized
  # vars = diag(AObs %*% inla.qsolve(Qnorm, t(AObs)))
  # quilt.plot(coords, vars, main="Marginal variance of normalized (not final) process")
  # normRatio = max(vars)/min(vars)
  # print(paste0("Ratio of largest to smallest variance for normalized process: ", normRatio))
  # 
  # final
  vars = diag(AObs %*% inla.qsolve(Qfinal, t(AObs)))
  pdf(file="Figures/varRatioFinal.pdf", width=5, height=5)
  quilt.plot(coords, vars, main="Marginal variance of normalized process")
  maxI = which.max(vars)
  minI = which.min(vars)
  if(plotMinMax) {
    points(coords[maxI,1], coords[maxI, 2], cex=2)
    points(coords[minI,1], coords[minI, 2], cex=2)
  }
  dev.off()
  finalRatio = max(vars)/min(vars)
  print(paste0("Ratio of largest to smallest variance for final process: ", finalRatio))
  
  c(unnormRatio=unnormRatio, finalRatio=finalRatio)
}
# testMakeNormalizedQMulti(kappa=.5, buffer=2.5, nx=20, ny=20)

makeNormalizedQMultiTable = function() {
  # testMakeNormalizedQMulti()
  # testMakeNormalizedQMulti(kappa=.01)
  # testMakeNormalizedQMulti(kappa=100)
  # testMakeNormalizedQMulti(nu=.5)
  # testMakeNormalizedQMulti(nu=.5, kappa=.01)
  # testMakeNormalizedQMulti(nu=.5, kappa=100)
  # testMakeNormalizedQMulti(nu=2.5)
  # testMakeNormalizedQMulti(nu=2.5, kappa=.01)
  # testMakeNormalizedQMulti(nu=2.5, kappa=100)
  require(xtable)
  
  tab5 = sapply(c(.01, .1, .5, 1, 5, 10, 100), testMakeNormalizedQMulti, buffer=2.5, nx=20, ny=20, nu=.5)
  tab15 = sapply(c(.01, .1, .5, 1, 5, 10, 100), testMakeNormalizedQMulti, buffer=2.5, nx=20, ny=20)
  tab25 = sapply(c(.01, .1, .5, 1, 5, 10, 100), testMakeNormalizedQMulti, buffer=2.5, nx=20, ny=20, nu=2.5)
  
  fullTab = rbind(tab5, tab15, tab25)
  colnames(fullTab) = c("kappa=0.01", "kappa=0.1", "kappa=0.5", "kappa=1", "kappa=5", "kappa=10", "kappa=100")
  xtable(fullTab, digits=2)
}


# computing basis matrix with Wendland covariance
# theta: 
makeA = function(predPts=NULL, xRangeKnot=c(0,1), xNKnot=10, yRangeKnot=c(0,1), yNKnot=10, theta=NULL, nLayer=1, thisLayer=1, 
                 xRangeDat=xRangeKnot, yRangeDat=yRangeKnot, nBuffer=5, rawKnotRange=FALSE, NC=5) {
  
  # get knot range based on data range
  if(!rawKnotRange) {
    out = makeLatGrid(xRangeDat, yRangeDat, NC=NC, nBuffer=nBuffer)
    xRangeKnot=out$xRangeKnots
    yRangeKnot=out$yRangeKnots
    xNKnot=out$nx
    yNKnot=out$ny
  }
  else {
    # make sure the raw grid has equal cell width in x and y directions
    xWidth = diff(xRange)/(nx-1)
    yWidth = diff(yRange)/(ny-1)
    if(xWidth != yWidth)
      stop("raw x and y grid widths not equal")
  }
  
  # setup theta if null
  if(is.null(theta))
    theta = 2.5*(xRangeKnot[2]-xRangeKnot[1])/(xNKnot-1)
  
  # set up prediction locations if necessary
  # have some sort of default behavior for setting predPts for ease of testing
  if(is.null(predPts)) {
    xRangePred = xRangeKnot
    yRangePred = yRangeKnot
    xNPred = xNKnot*3
    yNPred = yNKnot*3
    
    predPts = make.surface.grid(list(x=seq(xRangePred[1], xRangePred[2], l=xNPred), 
                                     y=seq(yRangePred[1], yRangePred[2], l=yNPred)))
  }
  
  # adjust basis elements to depend on layer
  origXNKnot = xNKnot
  origYNKnot = yNKnot
  origTheta = theta
  # xNKnot = xNKnot * 2^(thisLayer-1)
  # yNKnot = yNKnot * 2^(thisLayer-1)
  theta = theta / 2^(thisLayer-1)
  
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
  
  # generate the knot points for this layer
  knotPts = makeLatGrids(xRangeDat, yRangeDat, NC, nBuffer, nLayer)$latCoords[[thisLayer]]
  
  # NOTE: the wendland.cov function in `fields' does fancy things for only computing 
  # covariances for distances smaller than theta.  Returns a spam matrix
  thisA = as.dgCMatrix.spam(wendland.cov(predPts, knotPts, theta=theta, k=2))
  
  # return results recursively
  if(thisLayer == nLayer) {
    return(thisA)
  }
  else {
    return(cbind(thisA, makeA(predPts, xRangeKnot, origXNKnot, yRangeKnot, origYNKnot, origTheta, nLayer, thisLayer+1, 
                              xRangeDat=xRangeDat, yRangeDat=yRangeDat, nBuffer=nBuffer)))
  }
}
# out = makeQ(nLayer=3, nx=15, ny=15, nu=1, normalized=TRUE, newnormalize= TRUE)})

# makes alpha vector based on theoretical relationship in Nychka 2015 paper Cor 4.2
# NOTE: used to be based on ?LKrig, but it doesn't match the paper or 
#       what the package actually does (in LKrigSetupAlpha.default)
getAlphas = function(nLayer=3, nu=.5) {
  # 
  # alphas = exp(-2*(1:nLayer) * nu)
  # alphas = alphas/sum(alphas)
  # alphas
  thetaL=  2^(-(1:nLayer))
  alphas = thetaL^(2*nu)
  alphas = alphas/sum(alphas)
  alphas
}

# test for makeA and makeQ:
testMakeAQ = function(xRange=c(0,1), yRange=c(0,1), nx=10, ny=10, 
                      kappa=.1, rho=1, resFac=3, nLayer=1, nu=1.5, 
                      xRangeDat=xRange, yRangeDat=yRange, nBuffer=5, 
                      normalize=FALSE, seed=1) {
  set.seed(seed)
  
  if(nLayer != 1) {
    alphas = getAlphas(nLayer, nu)
  }
  else
    alphas = NULL
  
  # first step: make Q matrix
  Qmat = makeQ(kappa, rho, xRange, yRange, nx, ny, nLayer, alphas=alphas, xRangeDat=xRangeDat, 
               yRangeDat=yRangeDat, nBuffer=nBuffer, normalize=normalize)
  
  # second step: make A matrix
  par(mfrow=c(1,1))
  image(Qmat)
  predPts = make.surface.grid(list(x=seq(xRange[1], xRange[2], l=nx*resFac), 
                                   y=seq(yRange[1], yRange[2], l=ny*resFac)))
  Amat = makeA(predPts, xRange, nx, yRange, ny, nLayer=nLayer,  xRangeDat=xRangeDat, 
               yRangeDat=yRangeDat, nBuffer=nBuffer)
  
  # third step: simulate field on knot lattice and prediction points:
  L = t(chol(solve(Qmat)))
  Zs = matrix(rnorm(nrow(Qmat)), ncol=1)
  Cs = matrix(as.numeric(L %*% Zs), ncol=1)
  predSim = matrix(Amat %*% Cs, nrow=ny*resFac, ncol=nx*resFac)
  
  # make color scale for plots
  allVals = c(Cs, predSim)
  valRange = range(allVals)
  
  # make basis coefficient/resulting spatial field plots
  if(nLayer == 1) {
    par(mfrow=c(1,2))
    CsMat = matrix(Cs, nrow=ny, ncol=nx, byrow = FALSE)
    image.plot(CsMat, main="Basis Coefficients", breaks=seq(valRange[1], valRange[2], l=65))
    image.plot(predSim, main="Spatial Field", breaks=seq(valRange[1], valRange[2], l=65))
  }
  else if(nLayer < 4) {
    par(mfrow=c(2,2))
    layerInds = 1:(nx*ny)
    thisNX = nx
    thisNY = ny
    for(lay in 1:nLayer) {
      # plot this basis layer
      CsMat = matrix(Cs[layerInds], nrow=thisNY, ncol=thisNX, byrow = FALSE)
      image.plot(CsMat, main=paste0("Basis Coefficients (Layer ", lay, ")"), 
                 breaks=seq(valRange[1], valRange[2], l=65))
      
      # update layer indices
      maxInd = max(layerInds)
      thisNX = thisNX*2
      thisNY = thisNY*2
      layerInds = (maxInd + 1):(maxInd + thisNX*thisNY)
    }
    image.plot(predSim, main="Spatial Field", breaks=seq(valRange[1], valRange[2], l=65))
  }
  else {
    image.plot(predSim, main="Spatial Field", breaks=seq(valRange[1], valRange[2], l=65))
  }
}
# testMakeAQ(kappa=.001, rho=100, resFac=10)
# testMakeAQ(kappa=.001, rho=100, nLayer=3)
# testMakeAQ(kappa=.001, rho=1, nLayer=3, nu=.5)

# make the graph of the SAR model for LatticeKrig
makeGraph = function(nx=10, ny=10, nLayer=1, thisLayer=1, xRange=c(0,1), yRange=c(0,1), 
                     xRangeDat=xRange, yRangeDat=xRangeDat, nBuffer=5, NC=5) {
  # adjust resolution for layer number
  origNX = nx
  origNY = ny
  # nx = nx * 2^(thisLayer-1)
  # ny = ny * 2^(thisLayer-1)
  
  # # generate knot lattice locations and filter out locations 
  # # too far outside of the data domain
  # knotXs = seq(xRange[1], xRange[2], l=nx)
  # knotYs = seq(yRange[1], yRange[2], l=ny)
  # if(sum(knotXs > xRangeDat[2]) > nBuffer)
  #   knotXs = knotXs[1:(length(knotXs) - (sum(knotXs > xRangeDat[2]) - nBuffer))]
  # if(sum(knotXs < xRangeDat[1]) > nBuffer)
  #   knotXs = knotXs[(1 + sum(knotXs < xRangeDat[1]) - nBuffer):length(knotXs)]
  # if(sum(knotYs > yRangeDat[2]) > nBuffer)
  #   knotYs = knotYs[1:(length(knotYs) - (sum(knotYs > yRangeDat[2]) - nBuffer))]
  # if(sum(knotYs < yRangeDat[1]) > nBuffer)
  #   knotYs = knotYs[(1 + sum(knotYs < yRangeDat[1]) - nBuffer):length(knotYs)]
  # nx = length(knotXs)
  # ny = length(knotYs)
  
  # convert from raw grid parameterization to LatticeKrig grid parameterization, 
  # and generate the knot points for this layer
  gridPar = makeLatGrids(xRangeDat, yRangeDat, NC, nBuffer, nLayer)
  knotPts = gridPar$latCoords[[thisLayer]]
  nx = gridPar$nx[thisLayer]
  ny = gridPar$ny[thisLayer]
  
  # make (Laplacian) differential operators
  Dx = bandSparse(nx, k=0:1, diag=list(rep(1, nx), rep(1, nx-1)), symmetric=TRUE)
  Dy = bandSparse(ny, k=0:1, diag=list(rep(1, ny), rep(1, ny-1)), symmetric=TRUE)
  Dxx = bandSparse(nx, k=0:2, diag=list(rep(1, nx), rep(1, nx-1), rep(1, nx-2)), symmetric=TRUE)
  Dyy = bandSparse(ny, k=0:2, diag=list(rep(1, ny), rep(1, ny-1), rep(1, ny-2)), symmetric=TRUE)
  
  # generate x and y (Laplacian) differential operators
  Ix = Diagonal(n=nx)
  Iy = Diagonal(n=ny)
  Bx = kronecker(Iy, Dxx)
  By = kronecker(Dyy, Ix)
  Bxy = kronecker(Dy, Dx)
  
  
  # based on B, SAR regression matrix for Bc = e
  G = (Bx + By + Bxy) > 0
  # G = Bx + By + Bxy
  
  # if only need one layer, return what we have.  Otherwise, make sparse 
  # block diagonal matrix
  if(nLayer == thisLayer)
    return(G)
  else if((nLayer != 1) && (thisLayer == 1))
    return(bdiag(c(list(G), makeGraph(origNX, origNY, nLayer, thisLayer+1, 
                                      xRange=xRange, yRange=yRange, 
                                      xRangeDat=xRangeDat, yRangeDat=yRangeDat, nBuffer=nBuffer))))
  else
    return(c(list(G), makeGraph(origNX, origNY, nLayer, thisLayer+1, xRange=xRange, yRange=yRange, 
                                xRangeDat=xRangeDat, yRangeDat=yRangeDat, nBuffer=nBuffer)))
}

testMakeGraph1 = function(nx=3, ny=4) {
  image(makeGraph(nx, ny))
}
testMakeGraph2 = function(nx=3, ny=4) {
  Q = makeQ(nx=nx, ny=ny)
  image(Q != 0)
}
# testMakeGraph1()
# testMakeGraph2()

plotExampleDatasets = function() {
  # plot example datasets:
  # ozone2 dataset:
  # 8 hr average surface O3 levels (from 9am-4pm) in ppb
  # 153 sites in midwest
  # June 3 1987 - August 31 1987 (89 days)
  # 153 x 89 = 13,617 observations
  data("ozone2")
  names(ozone2)
  dim(ozone2$y)
  lonLat = ozone2$lon.lat
  y1 = ozone2$y[1,]
  quilt.plot(lonLat, y1, main="ozone2 (t=1 of 89)", xlab="Longitude", ylab="Latitude")
  US(add=TRUE)
  
  # also uses ozone2 in some examples
  
  # monthly min/max temperatures (deg C) and precip (total mm) from 1895 to 1997 (103 years).  
  # 376 stations, MAM= March, April, May
  # 376 x 103 = 38,728 observations
  data(COmonthlyMet)
  ?COmonthlyMet
  quilt.plot( CO.loc,CO.tmax.MAM[103,], main="Recorded MAM max temperatures (t=103 of 103)", 
              xlab="Longitude", ylab="Latitude")
  US( add=TRUE)
  
  # LatticeKrig paper also uses NorthAmericanRainfall dataset
  # 1720 stations precip in JJA=June,July,August based on data from 1950-2010 (61 years).  
  # Also includes elevation.
  # dataset in LK package only includes linear trends and intercepts for each station
  # 1720 x 1 = 1720 observations (or, with full dataset, 1720 x 61 = 104,920 observations)
  data(NorthAmericanRainfall)
  x<- cbind(NorthAmericanRainfall$longitude,  NorthAmericanRainfall$latitude)
  y<- NorthAmericanRainfall$precip
  quilt.plot( x,y/10, main="Mean JJA Precipitation, 1950-2010 (mm)", xlab="Longitude", ylab="Latitude")
  world( add=TRUE)
}



#####
# test out the correlation induced by kappa and rho
testCorrelation = function(nsim=1000, kappa=1, rho=1, N=30, NPred=60) {
  # get lattice points, prediction points, basis matrix APred
  gridPts = make.surface.grid(list(xs=seq(0, 1, l=10), 
                                   ys=seq(0, 1, l=10)))
  
  predPts = make.surface.grid(list(xs=seq(0,1,l=20), 
                                   ys=seq(0,1,l=20)))
  APred = makeA(predPts)
  
  # generate simulations
  Q = makeQ(kappa=kappa, rho=rho)
  L = t(chol(solve(Q)))
  zsims = matrix(rnorm(nsim*100), ncol=nsim)
  fieldSims = L %*% zsims
  pointSims = APred %*% fieldSims
  
  # estimate correlogram
  fullCorGram = NULL
  fullCorGramPred = NULL
  for(i in 1:nsim) {
    corG = myvgram(gridPts, fieldSims[,i], type="correlogram", colMeans=0)
    corGPred = myvgram(predPts, pointSims[,i], type="correlogram", colMeans=0)
    if(is.null(fullCorGram)) {
      fullCorGram = corG
      fullCorGramPred = corGPred
    }
    else {
      fullCorGram$d = c(fullCorGram$d, corG$d)
      fullCorGram$vgram = c(fullCorGram$vgram, corG$vgram)
      
      fullCorGramPred$d = c(fullCorGramPred$d, corGPred$d)
      fullCorGramPred$vgram = c(fullCorGramPred$vgram, corGPred$vgram)
    }
  }
  
  ## plot empirical versus theoretical correlograms for basis coefficients
  corMean = getVGMean(fullCorGram, N=N)
  notnans = !is.nan(corMean$ys)
  corRange = range(c(corMean$ys[notnans], 0, 1))
  plot(corMean$centers[notnans], corMean$ys[notnans], ylim=corRange, xlim=c(0, sqrt(2)), type="o", pch=19, 
       xlab="Distance", ylab="Correlation", main="Lattice Coefficient Correlation")
  xs = seq(0, sqrt(2), l=100)
  # > Matern(2.3)
  # [1] 0.1002588
  # effectiveRange=sqrt(8*1)/kappa # the distance at which Matern correlation is roughly .1
  squareWidth = .1
  effectiveRange=2.687/kappa*squareWidth
  lines(xs, Matern(xs, smoothness=1, range=effectiveRange/2.3))
  
  ## plot empirical versus theoretical correlograms for basis coefficients
  corMean = getVGMean(fullCorGramPred, N=NPred)
  notnans = !is.nan(corMean$ys)
  corRange = range(c(corMean$ys[notnans], 0, 1))
  plot(corMean$centers[notnans], corMean$ys[notnans], ylim=corRange, xlim=c(0, sqrt(2)), type="o", pch=19, 
       xlab="Distance", ylab="Correlation", main="Latent Field Correlation")
  xs = seq(0, sqrt(2), l=100)
  # > Matern(2.3)
  # [1] 0.1002588
  # effectiveRange=sqrt(8*1)/kappa # the distance at which Matern correlation is roughly .1
  squareWidth = .1
  effectiveRange=2.687/kappa*squareWidth
  lines(xs, Matern(xs, smoothness=1, range=effectiveRange/2.3))
}
# control.family(list(hyper=list(prec=list(initial=starting_value_for_log_prec, fixed=TRUE))))
# reduce observation noise by a lot (sigma=.1, sigma2=.01)
# look at marginal variance of the field
# compare results of lattice krig, inla.spde
# simulate a realization from the LK model with params near peak of priors
# fit a model
# test against LK, inla.spde
#####
# test out the correlation induced by kappa and rho
testCorrelation2 = function(kappa=1, rho=1, nx=30, mx=20, maxPlot=300, minEdgeDist=1, buffer=1) {
  # get lattice points, prediction points, basis matrix APred
  gridPts = make.surface.grid(list(xs=seq(0-buffer, 1+buffer, l=nx), 
                                   ys=seq(0-buffer, 1+buffer, l=nx)))
  
  predPts = make.surface.grid(list(xs=seq(0,1,l=mx), 
                                   ys=seq(0,1,l=mx)))
  latWidth = (1+buffer - (0 - buffer))/(nx-1)
  theta = 2.5*latWidth
  APred = makeA(predPts, theta=theta, xNKnot=nx, yNKnot=nx, xRangeKnot=c(0-buffer,1+buffer), yRangeKnot=c(0-buffer,1+buffer))
  
  # compute which points are far from the edge, and convert into the appropriate vectorized matrix using outer
  gridPtsCntr = (gridPts[,1] > 0-buffer + minEdgeDist) & (gridPts[,1] < 1+buffer - minEdgeDist) & 
    (gridPts[,2] > 0-buffer + minEdgeDist) & (gridPts[,2] < 1+buffer - minEdgeDist)
  predPtsCntr = (predPts[,1] > 0-buffer + minEdgeDist) & (predPts[,1] < 1+buffer - minEdgeDist) & 
    (predPts[,2] > 0-buffer + minEdgeDist) & (predPts[,2] < 1+buffer - minEdgeDist)
  gridPtsCntr = c(outer(gridPtsCntr, gridPtsCntr))
  predPtsCntr = c(outer(predPtsCntr, predPtsCntr))
  
  # compute variance and correlation matrices for lattice coefficients (c) and pred points (Ac)
  Q = makeQ(kappa=kappa, rho=rho, nx=nx, ny=nx, xRange=c(0-buffer,1+buffer), yRange=c(0-buffer,1+buffer))
  varC = solve(Q)
  varAC = APred %*% varC %*% t(APred)
  sigmasC = sqrt(diag(varC))
  sigmasAC = sqrt(diag(varAC))
  corC = c(sweep(sweep(as.matrix(varC), 1, 1/sigmasC, "*"), 2, 1/sigmasC, "*"))
  corAC = c(sweep(sweep(as.matrix(varAC), 1, 1/sigmasAC, "*"), 2, 1/sigmasAC, "*"))
  D = c(rdist(gridPts))
  DPred = c(rdist(predPts))
  
  ## plot empirical versus theoretical correlograms for basis coefficients
  # get unique distances and correlations and points far enough from edge
  uniqueDI = !duplicated(corC)
  uniqueDPredI = !duplicated(corAC)
  uniqueDI = uniqueDI & gridPtsCntr
  uniqueDPredI = uniqueDPredI & predPtsCntr
  uniqueD = D[uniqueDI]
  uniqueDPred = DPred[uniqueDPredI]
  uniqueCorC = corC[uniqueDI]
  uniqueCorAC = corAC[uniqueDPredI]
  
  if(length(uniqueD) > maxPlot) {
    inds = sample(1:length(uniqueD), maxPlot)
    uniqueD = uniqueD[inds]
    uniqueCorC = uniqueCorC[inds]
  }
  if(length(uniqueDPred) > maxPlot) {
    inds = sample(1:length(uniqueDPred), maxPlot)
    uniqueDPred = uniqueDPred[inds]
    uniqueCorAC = uniqueCorAC[inds]
  }
  
  # make sure all correlations at zero distance are 1
  zeroDistCors = corC[D == 0]
  if(any(rdist(zeroDistCors, 1) > 1e-6))
    warning(paste0("Bad lattice correlation at dist=0: ", zeroDistCors[match(1, rdist(zeroDistCors, 1) > 1e-6)]))
  zeroDistCors = corAC[DPred == 0]
  if(any(zeroDistCors != 1))
    warning(paste0("Bad latent field correlation at dist=0: ", zeroDistCors[match(1, rdist(zeroDistCors, 1) > 1e-6)]))
  
  # add a point at (0,1)
  uniqueD = c(uniqueD, 0)
  uniqueCorC = c(uniqueCorC, 1)
  uniqueDPred = c(uniqueDPred, 0)
  uniqueCorAC = c(uniqueCorAC, 1)
  
  plot(uniqueD, uniqueCorC, ylim=c(0,1), xlim=c(0, (latWidth*nx - 2*buffer)*sqrt(2)), pch=19, cex=.5, 
       xlab="Distance", ylab="Correlation", main="Lattice Coefficient Correlation")
  xs = seq(0, sqrt(2), l=100)
  # > Matern(2.3)
  # [1] 0.1002588
  # effectiveRange=sqrt(8*1)/kappa # the distance at which Matern correlation is roughly .1
  squareWidth = latWidth
  effectiveRange1=2.3/kappa*squareWidth
  effectiveRange2=sqrt(8)/kappa * latWidth
  lines(xs, Matern(xs, smoothness=1, range=effectiveRange1/2.3), col="blue")
  lines(xs, Matern(xs, smoothness=1, range=effectiveRange2/2.3), col="purple")
  legend("topright", c(TeX("$w/\\kappa$"), bquote(sqrt(8)*w/kappa)), 
         col=c("blue", "purple"), lty=1)
  
  ## plot empirical versus theoretical correlograms for basis coefficients
  plot(uniqueDPred, uniqueCorAC, ylim=c(0,1), xlim=c(0, (latWidth*nx - 2*buffer)*sqrt(2)), pch=19, cex=.3, 
       xlab="Distance", ylab="Correlation", main="Latent Field Correlation")
  xs = seq(0, sqrt(2), l=100)
  # > Matern(2.3)
  # [1] 0.1002588
  lines(xs, Matern(xs, smoothness=1, range=effectiveRange1/2.3), col="blue")
  lines(xs, Matern(xs, smoothness=1, range=effectiveRange2/2.3), col="purple")
  legend("topright", c(TeX("$w/\\kappa$"), bquote(sqrt(8)*w/kappa)), 
         col=c("blue", "purple"), lty=1)
}

# test out the correlation induced by kappa and rho for multi-layer
# model.  Plots saved as:
# margVarTest.pdf
# paste0("Figures/corNx", nx, "Nu", round(nu, digits=2), "L", nLayer, ".pdf")
testCorrelation3 = function(kappa=1, rho=1, nx=30, mx=20, maxPlot=300, minEdgeDist=1, buffer=1, nu=1, nLayer=3) {
  
  # get lattice points, prediction points, basis matrix APred
  gridPts = make.surface.grid(list(xs=seq(0-buffer, 1+buffer, l=nx), 
                                   ys=seq(0-buffer, 1+buffer, l=nx)))
  
  predPts = make.surface.grid(list(xs=seq(0,1,l=mx), 
                                   ys=seq(0,1,l=mx)))
  latWidth = (1+buffer - (0 - buffer))/(nx-1)
  theta = 2.5*latWidth
  APred = makeA(predPts, xNKnot=nx, yNKnot=nx, xRangeKnot=c(0-buffer,1+buffer), yRangeKnot=c(0-buffer,1+buffer), nLayer=nLayer, 
                xRangeDat=c(0,1), yRangeDat=c(0,1))
  
  # compute variance and correlation matrices for lattice coefficients (c) and pred points (Ac)
  Q = makeQ(kappa=kappa, rho=rho, nx=nx, ny=nx, xRange=c(0-buffer,1+buffer), yRange=c(0-buffer,1+buffer), 
            nLayer=nLayer, nu=nu, xRangeDat=c(0,1), yRangeDat=c(0,1))
  ny=nx
  
  # make mid lattice point
  xi = ceiling(nx/2)
  yi = ceiling(ny/2)
  midPt = matrix(c(seq(0-buffer, 1+buffer, l=nx)[xi], seq(0-buffer, 1+buffer, l=ny)[yi]), nrow=1)
  
  # now construct relevant row of A for value at midpoint, and Q matrix
  Ai = makeA(midPt, c(0-buffer,1+buffer), nx, c(0-buffer, 1+buffer), ny, nLayer=nLayer, 
             xRangeDat=c(0,1), yRangeDat=c(0,1))
  
  # renormalize basis coefficients to have constant variance, and the process to have unit variance
  Qinv = inla.qinv(Q)
  sds = sqrt(diag(Qinv))
  # Qnorm = sweep(sweep(Q, 1, sds, "*"), 2, sds, "*")
  Qnorm = Diagonal(x=sds) %*% Q %*% Diagonal(x=sds)
  procVar = as.numeric(Ai %*% inla.qsolve(Qnorm, t(Ai)))
  Q = Qnorm * (procVar / rho)
  
  # varC = solve(Q)
  varAC = APred %*% inla.qsolve(Q, t(APred))
  sigmasC = sqrt(diag(inla.qinv(Q)))
  sigmasCMatInv = Diagonal(x=1/sigmasC)
  sigmasAC = sqrt(diag(varAC))
  sigmasACMatInv = Diagonal(x=1/sigmasAC)
  # fullCorC = sweep(sweep(as.matrix(varC), 1, 1/sigmasC, "*"), 2, 1/sigmasC, "*")
  fullCorC = sigmasCMatInv %*% inla.qsolve(Q, sigmasCMatInv)
  # corAC = c(sweep(sweep(as.matrix(varAC), 1, 1/sigmasAC, "*"), 2, 1/sigmasAC, "*"))
  corAC = as.numeric(sigmasACMatInv %*% varAC %*% sigmasACMatInv)
  
  # compute distance matrix for prediction locations
  DPred = c(rdist(predPts))
  
  # plot correlation approximation and alpha weight for each layer
  if((nLayer >= 2) && (nLayer <= 3)) {
    pdf(file=paste0("Figures/corNx", nx, "Nu", round(nu, digits=2), "L", nLayer, ".pdf"), width=10, height=8)
    par(mfrow=c(2,2))
  }
  else if(nLayer > 3) {
    if(nLayer > 5)
      stop("nLayer > 5")
    pdf(file=paste0("Figures/corNx", nx, "Nu", round(nu, digits=2), "L", nLayer, ".pdf"), width=10, height=12)
    par(mfrow=c(3,2))
  }
  else {
    pdf(file=paste0("Figures/corNx", nx, "Nu", round(nu, digits=2), "L", nLayer, ".pdf"), width=10, height=5)
    par(mfrow=c(1,2))
  }
  startI=1
  alphas=getAlphas(nLayer, nu)
  for(l in 1:nLayer) {
    # get resolution for this lattice layer
    res = 2^(l-1)
    nLatticePts = (res*nx)^2
    layerInds = startI:(startI+nLatticePts-1)
    startI = startI + nLatticePts
    
    # get C submatrix for this layer in vector form
    corC = as.numeric(fullCorC[layerInds, layerInds])
    
    # get lattice points for this layer
    gridPts = make.surface.grid(list(xs=seq(0-buffer, 1+buffer, l=nx*res), 
                                     ys=seq(0-buffer, 1+buffer, l=nx*res)))
    
    # compute which points are far from the edge, and convert into the appropriate vectorized matrix using outer
    gridPtsCntr = (gridPts[,1] > 0-buffer + minEdgeDist) & (gridPts[,1] < 1+buffer - minEdgeDist) & 
      (gridPts[,2] > 0-buffer + minEdgeDist) & (gridPts[,2] < 1+buffer - minEdgeDist)
    predPtsCntr = (predPts[,1] > 0-buffer + minEdgeDist) & (predPts[,1] < 1+buffer - minEdgeDist) & 
      (predPts[,2] > 0-buffer + minEdgeDist) & (predPts[,2] < 1+buffer - minEdgeDist)
    gridPtsCntr = c(outer(gridPtsCntr, gridPtsCntr))
    predPtsCntr = c(outer(predPtsCntr, predPtsCntr))
    
    # get distance matrix for lattice points
    D = c(rdist(gridPts))
    
    ## plot approximate versus theoretical correlograms for basis coefficients
    # get unique distances and correlations and points far enough from edge
    uniqueDI = !duplicated(corC)
    uniqueDI = uniqueDI & gridPtsCntr
    uniqueD = D[uniqueDI]
    uniqueCorC = corC[uniqueDI]
    
    # subsample correlogram points to plot so there aren't too many
    if(length(uniqueD) > maxPlot) {
      inds = sample(1:length(uniqueD), maxPlot)
      uniqueD = uniqueD[inds]
      uniqueCorC = uniqueCorC[inds]
    }
    
    # make sure all correlations at zero distance are 1
    zeroDistCors = corC[D == 0]
    if(any(rdist(zeroDistCors, 1) > 1e-6))
      warning(paste0("Bad lattice correlation at dist=0: ", zeroDistCors[match(1, rdist(zeroDistCors, 1) > 1e-6)]))
    zeroDistCors = corAC[DPred == 0]
    if(any(zeroDistCors != 1))
      warning(paste0("Bad latent field correlation at dist=0: ", zeroDistCors[match(1, rdist(zeroDistCors, 1) > 1e-6)]))
    
    # add a point at (0,1)
    uniqueD = c(uniqueD, 0)
    uniqueCorC = c(uniqueCorC, 1)
    
    plot(uniqueD, uniqueCorC, ylim=c(0,1), xlim=c(0, (latWidth*nx - 2*buffer)*sqrt(2)), pch=19, cex=.5, 
         xlab="Distance", ylab="Correlation", main=paste0("Lattice Coefficient Correlation, Layer ", l, " alpha = ", round(alphas[l], digits=3)))
    xs = seq(0, sqrt(2), l=100)
    # > Matern(2.3)
    # [1] 0.1002588
    # effectiveRange=sqrt(8*1)/kappa # the distance at which Matern correlation is roughly .1
    squareWidth = latWidth/res
    effectiveRange=2.3/kappa*squareWidth
    lines(xs, Matern(xs, smoothness=1, range=effectiveRange/2.3), col="blue")
    
    # plot variance of layer basis coefficients
    # thisVarC = as.matrix(varC[layerInds, layerInds])
    # coefVars = diag(inla.qinv(Q)[layerInds, layerInds])
    # quilt.plot(gridPts, coefVars, main=paste0("Basis coefficient variance (layer ", l, ")"))
    
    # plot variance of layer predictions/process
    # tmp = as.matrix(APred)[,layerInds]
    # predVars = tmp %*% thisVarC %*% t(tmp)
    # quilt.plot(predPts, diag(predVars), main=paste0("Layer ", l, " process variance"))
  }
  
  # get correlations at prediction points
  uniqueDPredI = !duplicated(corAC)
  uniqueDPredI = uniqueDPredI & predPtsCntr
  uniqueCorAC = corAC[uniqueDPredI]
  uniqueDPred = DPred[uniqueDPredI]
  
  # subsample correlogram points to plot so there aren't too many
  if(length(uniqueDPred) > maxPlot) {
    inds = sample(1:length(uniqueDPred), maxPlot)
    uniqueDPred = uniqueDPred[inds]
    uniqueCorAC = uniqueCorAC[inds]
  }
  
  # add a point at (0,1)
  uniqueDPred = c(uniqueDPred, 0)
  uniqueCorAC = c(uniqueCorAC, 1)
  
  ## plot empirical versus theoretical correlograms for basis coefficients
  plot(uniqueDPred, uniqueCorAC, ylim=c(0,1), xlim=c(0, (latWidth*nx - 2*buffer)*sqrt(2)), pch=19, cex=.3, 
       xlab="Distance", ylab="Correlation", main=paste0("Latent Field Correlation, ", nLayer, " Layers"))
  xs = seq(0, sqrt(2), l=100)
  effectiveRange = 2.3/kappa*latWidth
  lines(xs, Matern(xs, smoothness=nu, range=effectiveRange/2.3), col="blue")
  
  # add line representing weighted sum of materns based on alphas
  cors = matrix(nrow=length(xs), ncol=nLayer)
  for(l in 1:nLayer) {
    res = 2^(l-1)
    squareWidth = latWidth/res
    effectiveRange=2.3/kappa*squareWidth
    
    cors[,l] = alphas[l] * Matern(xs, smoothness=1, range=effectiveRange/2.3)
  }
  ys = rowSums(cors)
  lines(xs, ys, col="purple")
  legend("topright", c("single Matern", "weighted Matern sum"), 
         col=c("blue", "purple"), lty=1)
  dev.off()
  
  ##### now for testing marginal variance
  sigma2 = rho/(4*pi*(4+kappa^2-4))
  print(paste0("Theoretical marginal variance is: ", sigma2))
  pdf("Figures/margVarTest_procVar.pdf", width=5, height=5)
  quilt.plot(predPts, sigmasAC^2, main="Process variance")
  dev.off()
  
  tmp = as.matrix(APred)
  quilt.plot(predPts, tmp %*% rep(1, ncol(APred)), main="Sum of basis functions")
  
  quilt.plot(predPts, tmp %*% rep(sigma2, ncol(APred)), main="Sum of basis functions times theoretical variance")
  
  # test the amount of inflation versus kappa
  myFun = function(x, thisRho=rho) { c(getMargVar(x, rho=thisRho)$inflation,  getMargVar(x, rho=thisRho)$theorVar, getMargVar(x, rho=thisRho)$actualVar) }
  kappas = 10^(seq(-4, 1, l=100))
  vars = sapply(kappas, myFun)
  inflates = vars[1,]
  actualVars = vars[3,]
  theorVars = vars[2,]
  options(scipen=5)
  pdf("Figures/margVarTest.pdf", width=5, height=5)
  plot(log10(kappas), theorVars, main="Normalized variance estimates", ylim=c(min(actualVars), max(actualVars)), 
       xlab=TeX("$\\kappa$"), ylab=TeX("$\\sigma^2 / \\rho $"), type="l", col="red", xaxt="n", 
       xlim=c(-2, 1))
  axis(1,at=c(-2, -1, 0, 1),labels=c("0.001", "0.1", "1", "10"))
  lines(log10(kappas), actualVars, col="blue")
  # lines(kappas, inflates, col="purple")
  # legend("topright", c("True Variance", "Asymptotic Variance", "Variance Ratio"), 
  #        lty=1, col=c("blue", "red", "purple"))
  legend("topright", c("True Variance", "Asymptotic Variance"), 
         lty=1, col=c("blue", "red"))
  dev.off()
}
# nu=1; L=7; sqrt(exp(-2*nu*(1:L))/sum(exp(-2*nu*(1:L))))

# get the marginal variance for multi-resolution process
# tod: theta/delta, or theta/latticeWidth
# either nu or alphas must be non-null
getMultiMargVar = function(kappa=1, rho=1, tod=2.5, nLayer=3, nu=NULL, alphas=NULL, nx=NULL, ny=NULL, 
                           xRange=c(0,1), yRange=xRange, xRangeDat=c(-2,1), yRangeDat=xRangeDat, nBuffer=5) {
  # set alphas if nu has been set
  if(!is.null(nu)) {
    alphas = getAlphas(nLayer, nu)
  }
  
  # set nx and ny if necessary and add buffer to avoid edge effects
  if(is.null(nx) || is.null(ny)) {
    maxPt = ceiling(tod)*4 + 1
    nx = maxPt
    ny = maxPt
  }
  
  # generate knot lattice locations and filter out locations 
  # too far outside of the data domain
  origNX = nx
  origNY = ny
  knotXs = seq(xRange[1], xRange[2], l=nx)
  knotYs = seq(yRange[1], yRange[2], l=ny)
  if(sum(knotXs > xRangeDat[2]) > nBuffer)
    knotXs = knotXs[1:(length(knotXs) - (sum(knotXs > xRangeDat[2]) - nBuffer))]
  if(sum(knotXs < xRangeDat[1]) > nBuffer)
    knotXs = knotXs[(1 + sum(knotXs < xRangeDat[1]) - nBuffer):length(knotXs)]
  if(sum(knotYs > yRangeDat[2]) > nBuffer)
    knotYs = knotYs[1:(length(knotYs) - (sum(knotYs > yRangeDat[2]) - nBuffer))]
  if(sum(knotYs < yRangeDat[1]) > nBuffer)
    knotYs = knotYs[(1 + sum(knotYs < yRangeDat[1]) - nBuffer):length(knotYs)]
  nx = length(knotXs)
  ny = length(knotYs)
  
  # sum the variances of each layer weighted by alphas
  totalMargVar = c()
  for(l in 1:nLayer) {
    # get the layer marginal variances
    layerMargVar = as.numeric(getMargVar(kappa, rho, tod, origNX*2^(l-1), origNY*2^(l-1), xRange=xRange, yRange=yRange, 
                                         xRangeDat=xRangeDat, yRangeDat=yRangeDat, nBuffer=nBuffer))
    layerMargVar[1:2] = layerMargVar[1:2]*alphas[l]
    if(l == 1)
      totalMargVar = layerMargVar[1:2]
    else
      totalMargVar = totalMargVar + layerMargVar[1:2]
  }
  
  # add in a variance ratio column
  totalMargVar = c(totalMargVar, totalMargVar[1]/totalMargVar[2])
  names(totalMargVar) = c("actualVar", "theorVar", "inflation")
  
  totalMargVar
}

# compute the marginal variance for a given resolution layer
# tod: theta/delta, or theta/latticeWidth
getMargVar = function(kappa=1, rho=1, tod=2.5, nx=NULL, ny=NULL, xRange=c(-1,2), yRange=xRange, 
                      xRangeDat=c(0,1), yRangeDat=xRangeDat, nBuffer=5) {
  # set nx and ny if necessary and add buffer to avoid edge effects
  if(is.null(nx) || is.null(ny)) {
    maxPt = ceiling(tod)*4 + 1
    nx = maxPt
    ny = maxPt
  }
  
  # generate knot lattice locations and filter out locations 
  # too far outside of the data domain
  knotXs = seq(xRange[1], xRange[2], l=nx)
  knotYs = seq(yRange[1], yRange[2], l=ny)
  delta = knotXs[2]-knotXs[1]
  if(sum(knotXs > xRangeDat[2]) > nBuffer)
    knotXs = knotXs[1:(length(knotXs) - (sum(knotXs > xRangeDat[2]) - nBuffer))]
  if(sum(knotXs < xRangeDat[1]) > nBuffer)
    knotXs = knotXs[(1 + sum(knotXs < xRangeDat[1]) - nBuffer):length(knotXs)]
  if(sum(knotYs > yRangeDat[2]) > nBuffer)
    knotYs = knotYs[1:(length(knotYs) - (sum(knotYs > yRangeDat[2]) - nBuffer))]
  if(sum(knotYs < yRangeDat[1]) > nBuffer)
    knotYs = knotYs[(1 + sum(knotYs < yRangeDat[1]) - nBuffer):length(knotYs)]
  knotPts = make.surface.grid(list(x=knotXs, y=knotYs))
  
  # take the middle knot location
  midPt = matrix(c(knotXs[ceiling(length(knotXs)/2)], 
                   knotYs[ceiling(length(knotYs)/2)]), nrow=1)
  
  # compute variance of process at the middle knot
  A = as.matrix(makeA(midPt, xRange, nx, yRange, ny, tod*delta, 
                xRangeDat=xRangeDat, yRangeDat=yRangeDat, nBuffer=nBuffer))
  Q = makeQ(kappa, rho, xRange, yRange, nx, ny, 
            xRangeDat=xRangeDat, yRangeDat=yRangeDat, nBuffer=nBuffer)
  # VarC = as.matrix(solve(Q))
  varMidPt = A %*% inla.qsolve(Q, t(A))
  
  # compare with theoretical marginal variance
  sigma2 = rho/(4*pi * kappa^2)
  inflation = varMidPt/sigma2
  
  # return results
  list(actualVar=varMidPt, theorVar=sigma2, inflation=inflation)
}


# estimates effective range for a given LatticeKrig model
getEffRange = function(predPts=NULL, xRangeKnot=c(0,1), xNKnot=10, yRangeKnot=c(0,1), yNKnot=10, theta=NULL, nLayer=1, thisLayer=1, 
                       xRangeDat=xRangeKnot, yRangeDat=yRangeKnot, nBuffer=5, mx=20, my=20) {
  
}

##### simulate from a latticeKrig model
# return a function with argument nsim for simulating a number of realizations from a latticeKrig model with no nugget.  
# coords: coordinates at which to simulate
# all other arguments: same as general latticeKrig arguments
LKSimulator = function(coords, NC=5, kappa=1, rho=1, nu=1.5, nBuffer=5, nLayer=3, normalize=TRUE) {
  # first make the grid on which to set the basis functions
  xRangeDat = range(coords[,1])
  yRangeDat = range(coords[,2])
  knotGrid = makeLatGrid(xRange=xRangeDat, yRange=yRangeDat, NC=NC, nBuffer=nBuffer)
  xRangeKnots=knotGrid$xRangeKnots
  nx=knotGrid$nx
  yRangeKnots=knotGrid$yRangeKnots
  ny=knotGrid$ny
  
  # generate layer variance weights
  alphas = getAlphas(nLayer, nu)
  
  # generate basis function and precision matrices
  AObs = makeA(coords, xRangeKnots, nx, yRangeKnots, ny, nLayer=nLayer, xRangeDat=xRangeDat, 
               yRangeDat=yRangeDat, nBuffer=nBuffer)
  Q = makeQ(kappa=kappa, rho=rho, xRange=xRangeBasis, yRange=yRangeBasis, nx=nx, ny=ny, 
            nLayer=nLayer, alphas=alphas, xRangeDat=xRangeDat, yRangeDat=yRangeDat, 
            nBuffer=nBuffer, normalized = normalize)
  L = as.matrix(t(chol(solve(Q))))
  zsim = matrix(rnorm(nrow(Q)), ncol=1)
  fieldSim = L %*% zsim
  
  fieldSims
}

# return a function with argument nsim for simulating a number of realizations from a latticeKrig model with no nugget.  
# same as LKSimulator, but uses marginal variance, effective range parameterization instead of rho, kappa.
# coords: coordinates at which to simulate
# all other arguments: same as general latticeKrig arguments
LKSimulator2 = function(coords, nsim=1, NC=5, effRange=(max(coords[,1])-min(coords[,1]))/3, margVar=1, 
                        nu=1.5, nBuffer=5, nLayer=3, normalize=TRUE) {
  
  # first make the grid on which to set the basis functions
  xRangeDat = range(coords[,1])
  yRangeDat = range(coords[,2])
  knotGrid = makeLatGrid(xRange=xRangeDat, yRange=yRangeDat, NC=NC, nBuffer=nBuffer)
  xRangeKnots=knotGrid$xRangeKnots
  nx=knotGrid$nx
  yRangeKnots=knotGrid$yRangeKnots
  ny=knotGrid$ny
  
  # convert from effRange, margVar to rho, kappa
  latticeWidth = (xRangeKnots[2] - xRangeKnots[1])/(nx-1)
  kappa = 2.3/effRange * latticeWidth
  
  # since we are normalizing the process, rho is just sigmaSq
  rho = margVar
  
  # generate layer variance weights
  alphas = getAlphas(nLayer, nu)
  
  # generate basis function and precision matrices
  AObs = makeA(coords, xRangeKnots, nx, yRangeKnots, ny, nLayer=nLayer, xRangeDat=xRangeDat, 
               yRangeDat=yRangeDat, nBuffer=nBuffer)
  Q = makeQ(kappa=kappa, rho=rho, xRange=xRangeKnots, yRange=yRangeKnots, nx=nx, ny=ny, 
            nLayer=nLayer, alphas=alphas, xRangeDat=xRangeDat, yRangeDat=yRangeDat, 
            nBuffer=nBuffer, normalized = normalize)
  L = as.matrix(t(chol(solve(Q))))
  zsims = matrix(rnorm(nrow(Q)*nsim), ncol=nsim)
  fieldSims = matrix(as.numeric(AObs %*% L %*% zsims), ncol=nsim)
  
  fieldSims
}

