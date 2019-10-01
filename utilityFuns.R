# this script contains some miscellaneous, but useful functions

# Do precomputations for computing precision matrix for a single layer or a block diagonal sparse 
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
# latticeInfo: an object returned by makeLatGrids
precomputationsQ = function(latticeInfo, thisLayer=1) {
  
  # save layer info quantities
  nLayer = length(latticeInfo)
  nx = latticeInfo[[thisLayer]]$nx
  ny = latticeInfo[[thisLayer]]$ny
  xRange=latticeInfo[[thisLayer]]$xRangeKnots
  yRange=latticeInfo[[thisLayer]]$yRangeKnots
  knotPts = latticeInfo[[thisLayer]]$latCoords
  
  # make (Laplacian) differential operators
  Dnx = bandSparse(nx, k=0:1, diag=list(rep(-2, nx), rep(1, nx-1)), symmetric=TRUE)
  Dny = bandSparse(ny, k=0:1, diag=list(rep(-2, ny), rep(1, ny-1)), symmetric=TRUE)
  
  # generate x and y (Laplacian) differential operators
  Inx = Diagonal(n=nx)
  Iny = Diagonal(n=ny)
  Bx = kronecker(Iny, Dnx)
  By = kronecker(Dny, Inx)
  Bxy = Bx + By
  
  ## now construct relevant row of A for value at midpoint at this layer, and Q matrix
  
  # make mid lattice point
  # xi = ceiling(nx/2)
  # yi = ceiling(ny/2)
  # midPt = matrix(c(seq(xRange[1], xRange[2], l=nx)[xi], seq(yRange[1], yRange[2], l=ny)[yi]), nrow=1)
  midPt = matrix(c((xRange[1] + xRange[2])/2, (yRange[1] + yRange[2])/2), nrow=1)
  
  Ai = makeA(midPt, latticeInfo, thisLayer=thisLayer, maxLayer=thisLayer)
  
  # return results
  # If multiple layers, return block diagonal sparse matrix
  if(thisLayer == nLayer) {
    return(list(list(Bxy=Bxy, Ai=Ai)))
  }
  else {
    return(c(list(list(Bxy=Bxy, Ai=Ai)), precomputationsQ(latticeInfo, thisLayer+1)))
  }
}

# Do precomputations for computing precision matrix for a single layer or a block diagonal sparse 
# precision matrix for multiple layers. Return pre- constructed block diagonal matrices instead of 
# individual matrices going into the block diagonal
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
# latticeInfo: an object returned by makeLatGrids
precomputationsQ2 = function(latticeInfo) {
  
  out = precomputationsQ(latticeInfo, 1)
  Bxys = lapply(out, function(x) {x$Bxy})
  Ais = lapply(out, function(x) {x$Ai})
  Bxy = bdiag(Bxys)
  mBxymBxyT = -Bxy - t(Bxy) # Bxy - t(Bxy)
  BxyTBxy = t(Bxy) %*% Bxy # t(Bxy) %*% Bxy
  A = bdiag(Ais)
  At = t(A)
  ms = sapply(out, function(x) {length(x$Ai)})
  
  list(mBxymBxyT=mBxymBxyT, BxyTBxy=BxyTBxy, A=A, At=At, ms=ms)
}

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
# latticeInfo: an object returned by makeLatGrids
makeQ = function(kappa=1, rho=1, latticeInfo, thisLayer=1, alphas=NULL, nu=NULL, 
                 normalized=FALSE, fastNormalize=FALSE, precomputedMatrices=NULL) {
  require(Matrix)
  require(spam)
  require(fields)
  
  # save base layer input quantities
  origRho = rho
  nLayer = length(latticeInfo)
  
  # make alphas according to nu relation if nu is set
  if(is.null(alphas) && !is.null(nu)) {
    alphas = getAlphas(nLayer, nu)
  }
  else if(is.null(nu) && nLayer != 1) {
    warning("Both alphas and nu are NULL. Defaulting to exponential covariance.")
    nu = 0.5
    alphas = getAlphas(nLayer, nu)
  }
  
  nx = latticeInfo[[thisLayer]]$nx
  ny = latticeInfo[[thisLayer]]$ny
  if(is.null(precomputedMatrices)) {
    # make (Laplacian) differential operators
    Dnx = bandSparse(nx, k=0:1, diag=list(rep(-2, nx), rep(1, nx-1)), symmetric=TRUE)
    Dny = bandSparse(ny, k=0:1, diag=list(rep(-2, ny), rep(1, ny-1)), symmetric=TRUE)
    
    # generate x and y (Laplacian) differential operators
    Inx = Diagonal(n=nx)
    Iny = Diagonal(n=ny)
    Bx = kronecker(Iny, Dnx)
    By = kronecker(Dny, Inx)
    Bxy = Bx + By
  }
  else
    Bxy = precomputedMatrices[[thisLayer]]$Bxy
  
  # make B, SAR regression matrix for Bc = e
  B = Diagonal(n=nx*ny, x=kappa^2) - Bxy
  
  # compute precision matrix
  Q = (1/rho) * t(B) %*% B
  
  if(normalized) {
    
    # now construct relevant row of A for value at midpoint at this layer, and Q matrix
    if(is.null(precomputedMatrices)) {
      xRange=latticeInfo[[thisLayer]]$xRangeKnots
      yRange=latticeInfo[[thisLayer]]$yRangeKnots
      
      # make mid lattice point
      # xi = ceiling(nx/2)
      # yi = ceiling(ny/2)
      # midPt = matrix(c(seq(xRange[1], xRange[2], l=nx)[xi], seq(yRange[1], yRange[2], l=ny)[yi]), nrow=1)
      midPt = matrix(c((xRange[1] + xRange[2])/2, (yRange[1] + yRange[2])/2), nrow=1)
      
      Ai = makeA(midPt, latticeInfo, thisLayer=thisLayer, maxLayer=thisLayer)
    }
    else
      Ai = precomputedMatrices[[thisLayer]]$Ai
    
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
    Q = bdiag(c(list(Q), makeQ(kappa, origRho, latticeInfo, thisLayer+1, alphas, normalized=normalized, 
                               fastNormalize=fastNormalize)))
    return(Q)
  }
  else {
    return(c(list(Q), makeQ(kappa, origRho, latticeInfo, thisLayer+1, alphas, normalized=normalized, 
                            fastNormalize=fastNormalize)))
  }
}









