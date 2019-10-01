# load required packages and R scripts
library(Matrix)
library(spam)
library(fields)
library(LatticeKrig)
library(invgamma)
library(INLA)

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
                 normalized=FALSE, fastNormalize=FALSE) {
  
  # save base layer input quantities
  origRho = rho
  nLayer = length(latticeInfo)
  
  # make alphas according to nu relation if nu is set
  if(is.null(alphas) && !is.null(nu)) {
    alphas = getAlphas(nLayer, nu)
  }
  
  xRange=latticeInfo[[thisLayer]]$xRangeKnots
  yRange=latticeInfo[[thisLayer]]$yRangeKnots
  knotPts = latticeInfo[[thisLayer]]$latCoords
  nx = latticeInfo[[thisLayer]]$nx
  ny = latticeInfo[[thisLayer]]$ny
  
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
    Ai = makeA(midPt, latticeInfo, thisLayer=thisLayer, maxLayer=thisLayer)
    
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

# make the graph of the SAR model for LatticeKrig
makeGraph = function(latticeInfo, thisLayer=1) {
  nx = latticeInfo[[thisLayer]]$nx
  ny = latticeInfo[[thisLayer]]$ny
  
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
  # G = (Bx + By + Bxy) > 0
  G = Bx + By + Bxy
  
  # if only need one layer, return what we have.  Otherwise, make sparse 
  # block diagonal matrix
  nLayer = length(latticeInfo)
  if(nLayer == thisLayer)
    return(G)
  else if((nLayer != 1) && (thisLayer == 1))
    return(bdiag(c(list(G), makeGraph(latticeInfo, thisLayer+1))))
  else
    return(c(list(G), makeGraph(latticeInfo, thisLayer+1)))
}

# computing basis matrix with Wendland covariance
# theta: 
# latticeInfo: an object returned by makeLatGrids
makeA = function(predPts=NULL, latticeInfo, thisLayer=1, theta=NULL, maxLayer=length(latticeInfo)) {
  nLayer = length(latticeInfo)
  xRangeKnot = latticeInfo[[thisLayer]]$xRangeKnots
  yRangeKnot = latticeInfo[[thisLayer]]$yRangeKnots
  xNKnot = latticeInfo[[thisLayer]]$nx
  yNKnot = latticeInfo[[thisLayer]]$ny
  xRangeDat = latticeInfo[[thisLayer]]$xRangeDat
  yRangeDat = latticeInfo[[thisLayer]]$yRangeDat
  nBuffer = latticeInfo[[thisLayer]]$nBuffer
  NC = latticeInfo[[thisLayer]]$NC
  knotPts = latticeInfo[[thisLayer]]$latCoords
  
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
  
  # NOTE: the wendland.cov function in `fields' does fancy things for only computing 
  # covariances for distances smaller than theta.  Returns a spam matrix
  thisA = as.dgCMatrix.spam(wendland.cov(predPts, knotPts, theta=theta, k=2))
  
  # return results recursively
  if(thisLayer == min(nLayer, maxLayer)) {
    return(thisA)
  }
  else {
    return(cbind(thisA, makeA(predPts, latticeInfo, thisLayer+1, theta=origTheta)))
  }
}

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

# generate lattice points for all layers
# NOTE: if any of the possibly NULL input variables are NULL, they are all set to the default
makeLatGrids = function(xRangeDat=c(0,1), yRangeDat=c(0,1), NC=5, nBuffer=5, nLayer=1) {
  
  grids = list()
  for(thisLayer in 1:nLayer) {
    i = thisLayer
    
    layerNC = (NC-1)*2^(i-1) + 1
    grids = c(grids, list(makeLatGrid(xRangeDat, yRangeDat, layerNC, nBuffer)))
  }
  grids
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
  xGrid = seq(xRangeKnots[1], xRangeKnots[2], l=nx)
  yGrid = seq(yRangeKnots[1], yRangeKnots[2], l=ny)
  latCoords = make.surface.grid(list(x=xGrid, y=yGrid))
  
  list(xRangeKnots=xRangeKnots, nx=nx, yRangeKnots=yRangeKnots, ny=ny, 
       xRangeDat=xRange, yRangeDat=yRange, NC=NC, nBuffer=nBuffer, latWidth=latWidth, 
       xGrid=xGrid, yGrid=yGrid, latCoords=latCoords)
}

# compute number of basis elements per layer.  Returns a list with element at index i 
# equal to the number of basis elements in layer i.
getMs = function(xRangeDat=c(0,1), yRangeDat=c(0,1), NC=5, nBuffer=5, nLayer=1) {
  # makeLatGrids(xRangeDat=xRangeDat, yRangeDat=yRangeDat, NC=NC, nBuffer=nBuffer, nLayer=nLayer)$ms
  sapply(makeLatGrids(xRangeDat=xRangeDat, yRangeDat=yRangeDat, NC=NC, nBuffer=nBuffer, nLayer=nLayer), function(x) {x$nx*x$ny})
}