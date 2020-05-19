##### the scripture contains functions for testing latticeKrig results

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

testMakeQBuffer = function(kappa=1, rho=1, latticeInfo=makeLatGrids(nLayer=3), nu=1.5, savePlot=FALSE) {
  nLayer=3
  
  if(nLayer != 1) {
    alphas = getAlphas(nLayer, nu)
  }
  else
    alphas = NULL
  
  # make precision matrix
  Q = makeQ(kappa, rho, latticeInfo, alphas=alphas)
  
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
    
    # generate knot lattice locations and filter out locations 
    # too far outside of the data domain
    nx = latticeInfo[[l]]$nx
    ny = latticeInfo[[l]]$ny
    xRange = latticeInfo[[l]]$xRangeKnots
    yRange = latticeInfo[[l]]$yRangeKnots
    nBuffer = latticeInfo[[l]]$nBuffer
    out = rawGridToLK(xRange, nx, yRange, ny, nBuffer)
    gridPts = latticeInfo[[l]]$latCoords
    
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
    gridPts = latticeInfo[[l]]$latCoords
    xRange = latticeInfo[[l]]$xRangeKnots
    yRange = latticeInfo[[l]]$yRangeKnots
    
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


testMakeGraph1 = function(nx=3, ny=4) {
  image(makeGraph(nx, ny))
}
testMakeGraph2 = function(nx=3, ny=4) {
  Q = makeQ(nx=nx, ny=ny)
  image(Q != 0)
}
# testMakeGraph1()
# testMakeGraph2()

plotExampleDatasets = function(savePlots=FALSE) {
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
  if(savePlots)
    pdf("ozoneDataSet.pdf", width=5, height=5)
  quilt.plot(lonLat, y1, main="ozone2 (t=1 of 89)", xlab="Longitude", ylab="Latitude")
  US(add=TRUE)
  if(savePlots)
    dev.off()
  
  # also uses ozone2 in some examples
  
  # monthly min/max temperatures (deg C) and precip (total mm) from 1895 to 1997 (103 years).  
  # 376 stations, MAM= March, April, May
  # 376 x 103 = 38,728 observations
  data(COmonthlyMet)
  ?COmonthlyMet
  if(savePlots)
    pdf("coloradoTempDataSet.pdf", width=5, height=5)
  quilt.plot( CO.loc,CO.tmax.MAM[103,], main="Recorded MAM max temperatures (t=103 of 103)", 
              xlab="Longitude", ylab="Latitude")
  US( add=TRUE)
  if(savePlots)
    dev.off()
  
  # LatticeKrig paper also uses NorthAmericanRainfall dataset
  # 1720 stations precip in JJA=June,July,August based on data from 1950-2010 (61 years).  
  # Also includes elevation.
  # dataset in LK package only includes linear trends and intercepts for each station
  # 1720 x 1 = 1720 observations (or, with full dataset, 1720 x 61 = 104,920 observations)
  data(NorthAmericanRainfall)
  x<- cbind(NorthAmericanRainfall$longitude,  NorthAmericanRainfall$latitude)
  y<- NorthAmericanRainfall$precip
  if(savePlots)
    pdf("precipitationDataSet.pdf", width=5, height=5)
  quilt.plot( x,y/10, main="Mean JJA Precipitation, 1950-2010 (mm)", xlab="Longitude", ylab="Latitude")
  world( add=TRUE)
  if(savePlots)
    dev.off
  
  if(savePlots)
    pdf("precipitationCoordsDataSet.pdf", width=5, height=5)
  plot(x[, 1], x[, 2], main="Mean JJA Precipitation, 1950-2010 (mm)", xlab="Longitude", ylab="Latitude", 
       pch=19, cex=.1)
  world( add=TRUE)
  if(savePlots)
    dev.off()
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
  # > Matern(sqrt(8))
  # [1] 0.1002588
  # effectiveRange=sqrt(8*1)/kappa # the distance at which Matern correlation is roughly .1
  squareWidth = .1
  effectiveRange=2.687/kappa*squareWidth
  lines(xs, Matern(xs, smoothness=1, range=effectiveRange/sqrt(8)))
  
  ## plot empirical versus theoretical correlograms for basis coefficients
  corMean = getVGMean(fullCorGramPred, N=NPred)
  notnans = !is.nan(corMean$ys)
  corRange = range(c(corMean$ys[notnans], 0, 1))
  plot(corMean$centers[notnans], corMean$ys[notnans], ylim=corRange, xlim=c(0, sqrt(2)), type="o", pch=19, 
       xlab="Distance", ylab="Correlation", main="Latent Field Correlation")
  xs = seq(0, sqrt(2), l=100)
  # > Matern(sqrt(8))
  # [1] 0.1002588
  # effectiveRange=sqrt(8*1)/kappa # the distance at which Matern correlation is roughly .1
  squareWidth = .1
  effectiveRange=2.687/kappa*squareWidth
  lines(xs, Matern(xs, smoothness=1, range=effectiveRange/sqrt(8)))
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
  # > Matern(sqrt(8))
  # [1] 0.1002588
  # effectiveRange=sqrt(8*1)/kappa # the distance at which Matern correlation is roughly .1
  squareWidth = latWidth
  effectiveRange1=sqrt(8)/kappa*squareWidth
  effectiveRange2=sqrt(8)/kappa * latWidth
  lines(xs, Matern(xs, smoothness=1, range=effectiveRange1/sqrt(8)), col="blue")
  lines(xs, Matern(xs, smoothness=1, range=effectiveRange2/sqrt(8)), col="purple")
  legend("topright", c(TeX("$w/\\kappa$"), bquote(sqrt(8)*w/kappa)), 
         col=c("blue", "purple"), lty=1)
  
  ## plot empirical versus theoretical correlograms for basis coefficients
  plot(uniqueDPred, uniqueCorAC, ylim=c(0,1), xlim=c(0, (latWidth*nx - 2*buffer)*sqrt(2)), pch=19, cex=.3, 
       xlab="Distance", ylab="Correlation", main="Latent Field Correlation")
  xs = seq(0, sqrt(2), l=100)
  # > Matern(sqrt(8))
  # [1] 0.1002588
  lines(xs, Matern(xs, smoothness=1, range=effectiveRange1/sqrt(8)), col="blue")
  lines(xs, Matern(xs, smoothness=1, range=effectiveRange2/sqrt(8)), col="purple")
  legend("topright", c(TeX("$w/\\kappa$"), bquote(sqrt(8)*w/kappa)), 
         col=c("blue", "purple"), lty=1)
}

# test out the correlation induced by kappa and rho for multi-layer
# model.  Plots saved as:
# margVarTest.pdf
# paste0("Figures/corNx", nx, "Nu", round(nu, digits=2), "L", nLayer, ".pdf")
testCorrelation3 = function(kappa=1, rho=1, nx=30, mx=20, maxPlot=300, minEdgeDist=1, buffer=1, nu=1, nLayer=2) {
  
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
    # > Matern(sqrt(8))
    # [1] 0.1002588
    # effectiveRange=sqrt(8*1)/kappa # the distance at which Matern correlation is roughly .1
    squareWidth = latWidth/res
    effectiveRange=sqrt(8)/kappa*squareWidth
    lines(xs, Matern(xs, smoothness=1, range=effectiveRange/sqrt(8)), col="blue")
    
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
  effectiveRange = sqrt(8)/kappa*latWidth
  lines(xs, Matern(xs, smoothness=nu, range=effectiveRange/sqrt(8)), col="blue")
  
  # add line representing weighted sum of materns based on alphas
  cors = matrix(nrow=length(xs), ncol=nLayer)
  for(l in 1:nLayer) {
    res = 2^(l-1)
    squareWidth = latWidth/res
    effectiveRange=sqrt(8)/kappa*squareWidth
    
    cors[,l] = alphas[l] * Matern(xs, smoothness=1, range=effectiveRange/sqrt(8))
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

testAlphaPrior = function() {
  require(hexbin)
  par(mfrow=c(2,3))
  samples = rdirichlet(1000000, alpha=c(1, 1, 1))
  h <- hexbin(samples[,1:2])
  plot(h)
  h <- hexbin(samples[,c(1,3)])
  plot(h)
  h <- hexbin(samples[,2:3])
  plot(h)
  hist(samples[,1], freq=FALSE)
  hist(samples[,2], freq=FALSE)
  hist(samples[,3], freq=FALSE)
}

# tests the fitLKINLAStandard function using data simulated from the LK model
# buffer: buffer distance between domain edge of basis lattice and domain edge of data.
# n: number of observations
# xRange: range of x coordinates
# yRange: range of y coordinates
# nx: number of basis function lattice points in x directions
# NOTE: ny is determined automatically to match scale of x lattice points
# Xmat: design matrix
# ys: observations
# first.time: is first time evaluating function.  User should always set to FALSE
testLKINLAModelStandard = function(buffer=2.5, kappa=1, rho=1, nu=1.5, seed=1, nLayer=3, nx=20, ny=nx, n=900, sigma2 = .2^2, 
                                   nBuffer=5, normalize=TRUE, fastNormalize=TRUE, NC=5, testCovs=TRUE, alphas=NULL, 
                                   printVerboseTimings=FALSE) {
  set.seed(seed)
  
  # compute alphas, the variance weights for each layer, depending on nu if necessary:
  if(is.null(alphas))
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
  pdf(file="Figures/standardObservations.pdf", width=5, height=5)
  par(mfrow=c(1,1))
  quilt.plot(coords, ys)
  dev.off()
  
  # make prediction coordinates on a grid
  xRange=c(0,1)
  yRange=c(0,1)
  mx = 100
  my = 100
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
  
  # show priors on effective correlation, marginal variance, and error variance:
  xs1 = seq(.01, 1, l=500)
  pdf(file="Figures/standardPriorEffRange.pdf", width=5, height=5)
  plot(xs1, dinvexp(xs1, rate=priorPar$corScalePar), type="l", col="blue", 
       xlab="Effective Correlation Range", main="Effective Correlation Prior", 
       ylab="Prior Density")
  abline(v=qinvexp(.5, rate=priorPar$corScalePar), col="red")
  dev.off()
  
  xs2 = seq(.01, 10.5, l=500)
  pdf(file="Figures/standardPriorMargVar.pdf", width=5, height=5)
  plot(xs2, invgamma::dinvgamma(xs2, shape=priorPar$varPar1, rate=priorPar$varPar2), type="l", col="blue", 
       xlab="Marginal Variance", main="Marginal Variance Prior", 
       ylab="Prior Density")
  abline(v=qinvgamma(.1, shape=priorPar$varPar1, rate=priorPar$varPar2), col="red")
  abline(v=qinvgamma(.9, shape=priorPar$varPar1, rate=priorPar$varPar2), col="red")
  dev.off()
  
  xs2 = seq(.001, invgamma::qinvgamma(.905, shape=0.1, rate=0.1), l=500)
  pdf(file="Figures/standardPriorErrorVar.pdf", width=5, height=5)
  plot(xs2, invgamma::dinvgamma(xs2, shape=0.1, rate=0.1), type="l", col="blue", 
       xlab="Error Variance", main="Error Variance Prior", 
       ylab="Prior Density")
  abline(v=invgamma::qinvgamma(.1, shape=0.1, rate=0.1), col="red")
  abline(v=invgamma::qinvgamma(.9, shape=0.1, rate=0.1), col="red")
  dev.off()
  
  for(l in 1:nLayer) {
    pdf(file=paste0("Figures/standardPriorAlpha", l, ".pdf"), width=5, height=5)
    xs = seq(0, 1, l=500)
    tempYs = dbeta(xs, priorPar$alphaPar[l], sum(priorPar$alphaPar[-l]))
    plot(xs, tempYs, type="l", xlab=TeX(paste0("$\\alpha_", l, "$")), ylab="Density", 
         main=TeX(paste0("Marginal for $\\alpha_", l, "$")), xlim=c(0,1), ylim=c(0, max(tempYs[is.finite(tempYs)])))
    abline(v=alphas[l], col="green")
    abline(v=qbeta(c(0.025, 0.975), priorPar$alphaPar[l], sum(priorPar$alphaPar[-l])), col="purple", lty=2)
    dev.off()
  }
  
  
  # fit the model
  time = system.time(out <- fitLKINLAStandard2(coords, ys, predCoords=predPts, nu=nu, seed=seed, nLayer=nLayer, NC=NC,
                                               nBuffer=nBuffer, priorPar=priorPar, xObs=X, xPred=XPred, normalize=normalize, 
                                               intStrategy="auto", strategy="laplace", fastNormalize=fastNormalize, 
                                               printVerboseTimings=printVerboseTimings))
  mod = out$mod
  preds=out$preds
  predSDs=out$SDs
  latInfo=out$latInfo
  latWidth=out$latWidth
  obsPreds=out$obsPreds
  obsSDs=out$obsSDs
  coefPreds = out$coefPreds
  coefSDs = out$coefSDs
  
  # print out the total time
  print(paste0("Total time: ", time[3]))
  
  # show a model summary
  print(summary(mod))
  
  # function for determining if points are in correct range
  inRange = function(pts, rangeShrink=0) {
    inX = (rangeShrink < pts[,1]) & (pts[,1] < 1-rangeShrink)
    inY = (rangeShrink < pts[,2]) & (pts[,2] < 1-rangeShrink)
    inX & inY
  }
  
  # show predictive surface, SD, and data
  
  pdf(file="Figures/standardPreds.pdf", width=15, height=6)
  if(nLayer==1) {
    par(mfrow=c(2,3))
    
    # obsInds = 1:n
    # predInds = (n+1):(n+mx*my)
    # coefInds = (n+mx*my+1):(n+mx*my+nx*ny)
    colRangeDat = range(c(ys-errs, obsPreds, preds, coefPreds))
    colRangeSD = range(c(range(predSDs[inRange(predPts)]), coefSDs[[1]][inRange(gridPtsL1)], 
                         coefSDs[[2]][inRange(gridPtsL2)], coefSDs[[3]][inRange(gridPtsL3)]))
    gridPtsL1 = latInfo[[1]]$latCoords
    quilt.plot(coords, ys-errs, main="True Process", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL1[,1], gridPtsL1[,2], coefPreds[[1]], main="Basis Coefficient Mean (Layer 1)", xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL1[,1], gridPtsL1[,2], coefSDs[[1]], main="Basis Coefficient SD (Layer 1)", zlim=colRangeSD, xlim=xRangeDat, ylim=yRangeDat)
    
    # quilt.plot(coords, obsPreds, main="Prediction mean", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat, 
    #            zlim=range(predSDs[inRange(predPts)]))
    plot.new()
    quilt.plot(predPts[,1], predPts[,2], preds, main="Prediction Mean", zlim=colRangeDat, 
               xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(predPts[,1], predPts[,2], predSDs, main="Prediction SD", 
               xlim=xRangeDat, ylim=yRangeDat)
  }
  else if(nLayer==2) {
    par(mfrow=c(2,4))
    
    # obsInds = 1:n
    # predInds = (n+1):(n+mx*my)
    # coefInds = (n+mx*my+1):(n+mx*my+nx*ny)
    colRangeDat = range(c(ys-errs, obsPreds, preds, coefPreds))
    colRangeCoef = range(c(coefPreds))
    colRangeSD = range(c(predSDs, obsSDs, coefSDs))
    gridPtsL1 = latInfo[[1]]$latCoords
    gridPtsL2 = latInfo[[2]]$latCoords
    quilt.plot(coords, ys-errs, main="True Process", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(predPts[,1], predPts[,2], preds, main="Prediction Mean", zlim=colRangeDat, 
               xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL1[,1], gridPtsL1[,2], coefPreds[[1]], main="Basis Coefficient Mean (Layer 1)", xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL2[,1], gridPtsL2[,2], coefPreds[[2]], main="Basis Coefficient Mean (Layer 2)", xlim=xRangeDat, ylim=yRangeDat)
    
    quilt.plot(coords, obsPreds, main="Observation Mean", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(predPts[,1], predPts[,2], predSDs, main="Prediction SD", 
               xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL1[,1], gridPtsL1[,2], coefSDs[[1]], main="Basis Coefficient SD (Layer 1)", zlim=colRangeSD, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL2[,1], gridPtsL2[,2], coefSDs[[2]], main="Basis Coefficient SD (Layer 2)", zlim=colRangeSD, xlim=xRangeDat, ylim=yRangeDat)
  }
  else if(nLayer==3) {
    par(mfrow=c(2,5), mar=c(5.1, 4.1, 4.1, 6))
    
    # obsInds = 1:n
    # predInds = (n+1):(n+mx*my)
    # coefInds = (n+mx*my+1):(n+mx*my+nx*ny)
    # colRangeDat = range(c(ys-errs, obsPreds, preds, coefPreds))
    colRangeDat = range(c(ys, obsPreds, preds, coefPreds))
    colRangeCoef = range(c(coefPreds))
    colRangeSD = range(c(predSDs, obsSDs, coefSDs))
    gridPtsL1 = latInfo[[1]]$latCoords
    gridPtsL2 = latInfo[[2]]$latCoords
    gridPtsL3 = latInfo[[3]]$latCoords
    quilt.plot(coords, ys-errs, main="True Process", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(predPts[,1], predPts[,2], preds, main="Prediction Mean", zlim=colRangeDat, 
               xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL1[,1], gridPtsL1[,2], coefPreds[[1]], main="Basis Coefficient Mean (Layer 1)", 
               xlim=xRangeDat, ylim=yRangeDat, zlim=colRangeCoef)
    quilt.plot(gridPtsL2[,1], gridPtsL2[,2], coefPreds[[2]], main="Basis Coefficient Mean (Layer 2)", 
               xlim=xRangeDat, ylim=yRangeDat, zlim=colRangeCoef)
    quilt.plot(gridPtsL3[,1], gridPtsL3[,2], coefPreds[[3]], main="Basis Coefficient Mean (Layer 3)", 
               xlim=xRangeDat, ylim=yRangeDat, zlim=colRangeCoef)
    
    # quilt.plot(coords, obsPreds, main="Observation Mean", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    # plot.new()
    quilt.plot(coords, ys, main="Observations", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(predPts[,1], predPts[,2], predSDs, main="Prediction SD",
               xlim=xRangeDat, ylim=yRangeDat, zlim=range(predSDs[inRange(predPts, rangeShrink=.03)]))
    quilt.plot(gridPtsL1[,1], gridPtsL1[,2], coefSDs[[1]], main="Basis Coefficient SD (Layer 1)", zlim=colRangeSD, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL2[,1], gridPtsL2[,2], coefSDs[[2]], main="Basis Coefficient SD (Layer 2)", zlim=colRangeSD, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL3[,1], gridPtsL3[,2], coefSDs[[3]], main="Basis Coefficient SD (Layer 3)", zlim=colRangeSD, xlim=xRangeDat, ylim=yRangeDat)
  }
  dev.off()
  
  # calculate true effective range and marginal variance:
  latticeWidth = latInfo[[1]]$latWidth
  effRange = sqrt(8)/kappa * latticeWidth
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
  pdf(file="Figures/standardEffRange.pdf", width=5, height=5)
  plot(effRangeMarg, type="l", main="Marginal for effective range")
  abline(v=effRange, col="green")
  abline(v=inla.qmarginal(c(.025, .975), effRangeMarg), col="purple", lty=2)
  dev.off()
  # plot(mod$marginals.hyperpar$`Theta1 for field`, type="l", main="Marginal for log range")
  pdf(file="Figures/standardVar.pdf", width=5, height=5)
  plot(varMarg, type="l", main="Marginal for spatial variance")
  abline(v=marginalVar, col="green")
  abline(v=inla.qmarginal(c(.025, .975), varMarg), col="purple", lty=2)
  dev.off()
  # plot(mod$marginals.hyperpar$`Theta2 for field`, type="l", main="Marginal for log variance")
  pdf(file="Figures/standardSigma2.pdf", width=5, height=5)
  plot(sigma2Marg, type="l", main="Marginal for error variance")
  abline(v=sigma2, col="green")
  abline(v=inla.qmarginal(c(.025, .975), sigma2Marg), col="purple", lty=2)
  dev.off()
  for(i in 1:length(covNames)) {
    XMarginal = XMarginals[[i]]
    pdf(file=paste0("Figures/standard", covNames[i], ".pdf"), width=5, height=5)
    plot(XMarginal, type="l", main="Marginal for fixed effect")
    if(i==1)
      abline(v=1, col="green")
    else
      abline(v=0, col="green")
    abline(v=inla.qmarginal(c(.025, .975), XMarginal), col="purple", lty=2)
    dev.off()
  }
  
  # pdf(file="Figures/standardRho.pdf", width=5, height=5)
  # plot(sigma2Marg, type="l", main=TeX("Marginal for $\\rho$"), xlab=TeX("$\\rho$"))
  # abline(v=rho, col="green")
  # dev.off()
  
  
  # do the same for kappa, rho
  # in order to get distribution for rho, must sample from joint hyperparameters
  kappaMarg = inla.tmarginal(function(x) {sqrt(8)/exp(x) * latticeWidth}, mod$marginals.hyperpar$`Theta1 for field`)
  # thetasToRho = function(xs) {
  #   logCor = xs[2]
  #   logVar = xs[3]
  #   kappa = sqrt(8)/exp(logCor) * latticeWidth
  #   sigma2 = exp(logVar)
  #   sigma2 * 4*pi * kappa^2
  # }
  # samples = inla.hyperpar.sample(50000, mod, TRUE)
  # rhos = apply(samples, 1, thetasToRho)
  
  pdf(file="Figures/standardKappa.pdf", width=5, height=5)
  plot(kappaMarg, type="l", xlab="kappa", main="Marginal for kappa")
  abline(v=kappa, col="green")
  abline(v=inla.qmarginal(c(.025, .975), kappaMarg), col="purple", lty=2)
  dev.off()
  
  # pdf(file="Figures/standardRho.pdf", width=5, height=5)
  # hist(rhos, xlab="rho", main="Marginal for Rho", breaks=1000, freq=F, xlim=c(0, quantile(probs=.95, rhos)))
  # abline(v=rho, col="green")
  # dev.off()
  
  ## Now generate marginals for the alpha parameters. In order to do this, we must generate draws from 
  ## the posterior, and transform them back to the probability scale
  out = inla.hyperpar.sample(10000, mod, improve.marginals=TRUE)
  zSamples = out[,4:(3+nLayer-1)]
  xSamples = apply(zSamples, 1, multivariateExpit)
  xSamples = rbind(xSamples, 1-colSums(xSamples))
  
  for(l in 1:nLayer) {
    pdf(file=paste0("Figures/standardAlpha", l, ".pdf"), width=5, height=5)
    hist(xSamples[l,], xlab=TeX(paste0("$\\alpha_", l, "$")), main=TeX(paste0("Marginal for $\\alpha_", l, "$")), breaks=100, freq=F, xlim=c(0,1))
    abline(v=alphas[l], col="green")
    abline(v=mean(xSamples[l,]), col="purple", lty=1)
    abline(v=quantile(probs=c(.025, .975), xSamples[l,]), col="purple", lty=2)
    dev.off()
  }
  
  invisible(NULL)
}

# tests the fitLKINLAStandard function using data simulated from the LK model with binomial likelihood
# buffer: buffer distance between domain edge of basis lattice and domain edge of data.
# n: number of observations
# xRange: range of x coordinates
# yRange: range of y coordinates
# nx: number of basis function lattice points in x directions
# NOTE: ny is determined automatically to match scale of x lattice points
# Xmat: design matrix
# ys: observations
# first.time: is first time evaluating function.  User should always set to FALSE
testLKINLAModelStandardBinomial = function(buffer=2.5, kappa=1, rho=1, nu=1.5, seed=1, nLayer=3, nx=20, ny=nx, n=900, sigma2 = .2^2, 
                                           nBuffer=5, normalize=TRUE, fastNormalize=TRUE, NC=5, testCovs=TRUE, alphas=NULL, 
                                           printVerboseTimings=FALSE, nObs=rep(25, n), intStrategy="auto", strategy="gaussian") {
  set.seed(seed)
  
  # compute alphas, the variance weights for each layer, depending on nu if necessary:
  if(is.null(alphas))
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
  ps = as.numeric(AObs %*% fieldSims) + 1 # add a constant unit mean term to be estimated by INLA
  # ys = 1 + as.numeric(AObs %*% fieldSims) + coords[,1] # x-valued mean term to be estimated by INLA
  errs = rnorm(n, sd=sqrt(sigma2))
  ps = expit(ps + errs)
  ys = rbinom(n, nObs, ps)
  
  # plot the observations
  pdf(file="Figures/standardObservationsBinom.pdf", width=5, height=5)
  par(mfrow=c(1,1))
  quilt.plot(coords, ys, main="Number observations")
  dev.off()
  
  pdf(file="Figures/standardLatentBinom.pdf", width=5, height=5)
  par(mfrow=c(1,1))
  quilt.plot(coords, ps, main="Latent probabilities")
  dev.off()
  
  # make prediction coordinates on a grid
  xRange=c(0,1)
  yRange=c(0,1)
  mx = 100
  my = 100
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
  
  # show priors on effective correlation, marginal variance, and error variance:
  xs1 = seq(.01, 1, l=500)
  pdf(file="Figures/standardPriorEffRangeBinom.pdf", width=5, height=5)
  plot(xs1, dinvexp(xs1, rate=priorPar$corScalePar), type="l", col="blue", 
       xlab="Effective Correlation Range", main="Effective Correlation Prior", 
       ylab="Prior Density")
  abline(v=qinvexp(.5, rate=priorPar$corScalePar), col="red")
  dev.off()
  
  xs2 = seq(.01, 10.5, l=500)
  pdf(file="Figures/standardPriorMargVarBinom.pdf", width=5, height=5)
  plot(xs2, invgamma::dinvgamma(xs2, shape=priorPar$varPar1, rate=priorPar$varPar2), type="l", col="blue", 
       xlab="Marginal Variance", main="Marginal Variance Prior", 
       ylab="Prior Density")
  abline(v=qinvgamma(.1, shape=priorPar$varPar1, rate=priorPar$varPar2), col="red")
  abline(v=qinvgamma(.9, shape=priorPar$varPar1, rate=priorPar$varPar2), col="red")
  dev.off()
  
  xs2 = seq(.001, invgamma::qinvgamma(.905, shape=0.1, rate=0.1), l=500)
  pdf(file="Figures/standardPriorErrorVarBinom.pdf", width=5, height=5)
  plot(xs2, invgamma::dinvgamma(xs2, shape=0.1, rate=0.1), type="l", col="blue", 
       xlab="Error Variance", main="Error Variance Prior", 
       ylab="Prior Density")
  abline(v=invgamma::qinvgamma(.1, shape=0.1, rate=0.1), col="red")
  abline(v=invgamma::qinvgamma(.9, shape=0.1, rate=0.1), col="red")
  dev.off()
  
  for(l in 1:nLayer) {
    pdf(file=paste0("Figures/standardPriorAlpha", l, "Binom.pdf"), width=5, height=5)
    xs = seq(0, 1, l=500)
    tempYs = dbeta(xs, priorPar$alphaPar[l], sum(priorPar$alphaPar[-l]))
    plot(xs, tempYs, type="l", xlab=TeX(paste0("$\\alpha_", l, "$")), ylab="Density", 
         main=TeX(paste0("Marginal for $\\alpha_", l, "$")), xlim=c(0,1), ylim=c(0, max(tempYs[is.finite(tempYs)])))
    abline(v=alphas[l], col="green")
    abline(v=qbeta(c(0.025, 0.975), priorPar$alphaPar[l], sum(priorPar$alphaPar[-l])), col="purple", lty=2)
    dev.off()
  }
  
  
  # fit the model
  time = system.time(out <- fitLKINLAStandard2(coords, ys, predCoords=predPts, nu=nu, seed=seed, nLayer=nLayer, NC=NC,
                                               nBuffer=nBuffer, priorPar=priorPar, xObs=X, xPred=XPred, normalize=normalize, 
                                               intStrategy=intStrategy, strategy=strategy, fastNormalize=fastNormalize, 
                                               printVerboseTimings=printVerboseTimings, family="binomial", obsNs=nObs))
  mod = out$mod
  preds=out$preds
  predSDs=out$sigmas
  latInfo=out$latInfo
  latWidth=out$latWidth
  obsPreds=out$obsPreds
  obsSDs=out$obsSDs
  coefPreds = out$coefPreds
  coefSDs = out$coefSDs
  
  # print out the total time
  print(paste0("Total time: ", time[3]))
  
  # show a model summary
  print(summary(mod))
  
  # function for determining if points are in correct range
  inRange = function(pts, rangeShrink=0) {
    inX = (rangeShrink < pts[,1]) & (pts[,1] < 1-rangeShrink)
    inY = (rangeShrink < pts[,2]) & (pts[,2] < 1-rangeShrink)
    inX & inY
  }
  
  # show predictive surface, SD, and data
  
  pdf(file="Figures/standardPredsBinom.pdf", width=15, height=6)
  if(nLayer==1) {
    par(mfrow=c(2,3))
    
    # obsInds = 1:n
    # predInds = (n+1):(n+mx*my)
    # coefInds = (n+mx*my+1):(n+mx*my+nx*ny)
    colRangeDat = range(c((ys-errs) / nObs, obsPreds, preds, expit(unlist(coefPreds))))
    colRangeSD = range(c(range(predSDs[inRange(predPts)]), expit(coefSDs[[1]][inRange(gridPtsL1)]), 
                         expit(coefSDs[[2]][inRange(gridPtsL2)]), expit(coefSDs[[3]][inRange(gridPtsL3)])))
    gridPtsL1 = latInfo[[1]]$latCoords
    quilt.plot(coords, (ys-errs) / nObs, main="True Process", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL1[,1], gridPtsL1[,2], expit(coefPreds[[1]]), main="Basis Coefficient Mean (Layer 1)", xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL1[,1], gridPtsL1[,2], expit(coefSDs[[1]]), main="Basis Coefficient SD (Layer 1)", zlim=colRangeSD, xlim=xRangeDat, ylim=yRangeDat)
    
    # quilt.plot(coords, obsPreds, main="Prediction mean", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat, 
    #            zlim=range(predSDs[inRange(predPts)]))
    plot.new()
    quilt.plot(predPts[,1], predPts[,2], preds, main="Prediction Mean", zlim=colRangeDat, 
               xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(predPts[,1], predPts[,2], predSDs, main="Prediction SD", 
               xlim=xRangeDat, ylim=yRangeDat)
  } else if(nLayer==2) {
    par(mfrow=c(2,4))
    
    # obsInds = 1:n
    # predInds = (n+1):(n+mx*my)
    # coefInds = (n+mx*my+1):(n+mx*my+nx*ny)
    colRangeDat = range(c((ys-errs) / nObs, obsPreds, preds, expit(unlist(coefPreds))))
    colRangeCoef = range(c(expit(unlist(coefPreds))))
    colRangeSD = range(c(predSDs, obsSDs, expit(unlist(coefSDs))))
    gridPtsL1 = latInfo[[1]]$latCoords
    gridPtsL2 = latInfo[[2]]$latCoords
    quilt.plot(coords, (ys-errs) / nObs, main="True Process", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(predPts[,1], predPts[,2], preds, main="Prediction Mean", zlim=colRangeDat, 
               xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL1[,1], gridPtsL1[,2], expit(coefPreds[[1]]), main="Basis Coefficient Mean (Layer 1)", xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL2[,1], gridPtsL2[,2], expit(coefPreds[[2]]), main="Basis Coefficient Mean (Layer 2)", xlim=xRangeDat, ylim=yRangeDat)
    
    quilt.plot(coords, obsPreds, main="Observation Mean", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(predPts[,1], predPts[,2], predSDs, main="Prediction SD", 
               xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL1[,1], gridPtsL1[,2], expit(coefSDs[[1]]), main="Basis Coefficient SD (Layer 1)", zlim=colRangeSD, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL2[,1], gridPtsL2[,2], expit(coefSDs[[2]]), main="Basis Coefficient SD (Layer 2)", zlim=colRangeSD, xlim=xRangeDat, ylim=yRangeDat)
  } else if(nLayer==3) {
    par(mfrow=c(2,5), mar=c(5.1, 4.1, 4.1, 6))
    
    # obsInds = 1:n
    # predInds = (n+1):(n+mx*my)
    # coefInds = (n+mx*my+1):(n+mx*my+nx*ny)
    # colRangeDat = range(c(ys-errs, obsPreds, preds, coefPreds))
    colRangeDat = range(c(ys / nObs, obsPreds, preds, expit(unlist(coefPreds))))
    colRangeCoef = range(c(expit(unlist(coefPreds))))
    colRangeSD = range(c(predSDs, obsSDs, expit(unlist(coefSDs))))
    gridPtsL1 = latInfo[[1]]$latCoords
    gridPtsL2 = latInfo[[2]]$latCoords
    gridPtsL3 = latInfo[[3]]$latCoords
    quilt.plot(coords, (ys-errs) / nObs, main="True Process", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(predPts[,1], predPts[,2], preds, main="Prediction Mean", zlim=colRangeDat, 
               xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL1[,1], gridPtsL1[,2], expit(coefPreds[[1]]), main="Basis Coefficient Mean (Layer 1)", 
               xlim=xRangeDat, ylim=yRangeDat, zlim=colRangeCoef)
    quilt.plot(gridPtsL2[,1], gridPtsL2[,2], expit(coefPreds[[2]]), main="Basis Coefficient Mean (Layer 2)", 
               xlim=xRangeDat, ylim=yRangeDat, zlim=colRangeCoef)
    quilt.plot(gridPtsL3[,1], gridPtsL3[,2], expit(coefPreds[[3]]), main="Basis Coefficient Mean (Layer 3)", 
               xlim=xRangeDat, ylim=yRangeDat, zlim=colRangeCoef)
    
    # quilt.plot(coords, obsPreds, main="Observation Mean", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    # plot.new()
    quilt.plot(coords, ys / nObs, main="Empirical Proportions", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(predPts[,1], predPts[,2], predSDs, main="Prediction SD",
               xlim=xRangeDat, ylim=yRangeDat, zlim=range(predSDs[inRange(predPts, rangeShrink=.03)]))
    quilt.plot(gridPtsL1[,1], gridPtsL1[,2], expit(coefSDs[[1]]), main="Basis Coefficient SD (Layer 1)", zlim=colRangeSD, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL2[,1], gridPtsL2[,2], expit(coefSDs[[2]]), main="Basis Coefficient SD (Layer 2)", zlim=colRangeSD, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL3[,1], gridPtsL3[,2], expit(coefSDs[[3]]), main="Basis Coefficient SD (Layer 3)", zlim=colRangeSD, xlim=xRangeDat, ylim=yRangeDat)
  }
  dev.off()
  
  # calculate true effective range and marginal variance:
  latticeWidth = latInfo[[1]]$latWidth
  effRange = sqrt(8)/kappa * latticeWidth
  # marginalVar = rho/(4*pi * kappa^2)
  # marginalVar = getMultiMargVar(kappa, rho, nLayer=nLayer, nu=nu, xRange=xRangeBasis, 
  #                               yRange=yRangeBasis, nx=nx, ny=ny)[1]
  marginalVar = rho
  
  # plot marginals on interpretable scale (effective range, marginal variance)
  effRangeMarg = inla.tmarginal(exp, mod$marginals.hyperpar$`Theta1 for field`)
  varMarg = inla.tmarginal(exp, mod$marginals.hyperpar$`Theta2 for field`)
  sigma2Marg = inla.tmarginal(function(x) {1/x}, mod$marginals.hyperpar$`Precision for clust`)
  covNames = names(mod$marginals.fixed)
  XMarginals = list()
  for(i in 1:length(covNames)) {
    XMarginal = inla.tmarginal(function(x) {x}, mod$marginals.fixed[[covNames[i]]])
    XMarginals = c(XMarginals, list(XMarginal))
  }
  
  par(mfrow=c(1,1))
  pdf(file="Figures/standardEffRangeBinom.pdf", width=5, height=5)
  plot(effRangeMarg, type="l", main="Marginal for effective range")
  abline(v=effRange, col="green")
  abline(v=inla.qmarginal(c(.025, .975), effRangeMarg), col="purple", lty=2)
  dev.off()
  # plot(mod$marginals.hyperpar$`Theta1 for field`, type="l", main="Marginal for log range")
  pdf(file="Figures/standardVarBinom.pdf", width=5, height=5)
  plot(varMarg, type="l", main="Marginal for spatial variance")
  abline(v=marginalVar, col="green")
  abline(v=inla.qmarginal(c(.025, .975), varMarg), col="purple", lty=2)
  dev.off()
  # plot(mod$marginals.hyperpar$`Theta2 for field`, type="l", main="Marginal for log variance")
  pdf(file="Figures/standardSigma2Binom.pdf", width=5, height=5)
  plot(sigma2Marg, type="l", main="Marginal for error variance")
  abline(v=sigma2, col="green")
  abline(v=inla.qmarginal(c(.025, .975), sigma2Marg), col="purple", lty=2)
  dev.off()
  for(i in 1:length(covNames)) {
    XMarginal = XMarginals[[i]]
    pdf(file=paste0("Figures/standard", covNames[i], "Binom.pdf"), width=5, height=5)
    plot(XMarginal, type="l", main="Marginal for fixed effect")
    if(i==1)
      abline(v=1, col="green")
    else
      abline(v=0, col="green")
    abline(v=inla.qmarginal(c(.025, .975), XMarginal), col="purple", lty=2)
    dev.off()
  }
  
  # pdf(file="Figures/standardRhoBinom.pdf", width=5, height=5)
  # plot(sigma2Marg, type="l", main=TeX("Marginal for $\\rho$"), xlab=TeX("$\\rho$"))
  # abline(v=rho, col="green")
  # dev.off()
  
  
  # do the same for kappa, rho
  # in order to get distribution for rho, must sample from joint hyperparameters
  kappaMarg = inla.tmarginal(function(x) {sqrt(8)/exp(x) * latticeWidth}, mod$marginals.hyperpar$`Theta1 for field`)
  # thetasToRho = function(xs) {
  #   logCor = xs[2]
  #   logVar = xs[3]
  #   kappa = sqrt(8)/exp(logCor) * latticeWidth
  #   sigma2 = exp(logVar)
  #   sigma2 * 4*pi * kappa^2
  # }
  # samples = inla.hyperpar.sample(50000, mod, TRUE)
  # rhos = apply(samples, 1, thetasToRho)
  
  pdf(file="Figures/standardKappaBinom.pdf", width=5, height=5)
  plot(kappaMarg, type="l", xlab="kappa", main="Marginal for kappa")
  abline(v=kappa, col="green")
  abline(v=inla.qmarginal(c(.025, .975), kappaMarg), col="purple", lty=2)
  dev.off()
  
  # pdf(file="Figures/standardRhoBinom.pdf", width=5, height=5)
  # hist(rhos, xlab="rho", main="Marginal for Rho", breaks=1000, freq=F, xlim=c(0, quantile(probs=.95, rhos)))
  # abline(v=rho, col="green")
  # dev.off()
  
  ## Now generate marginals for the alpha parameters. In order to do this, we must generate draws from 
  ## the posterior, and transform them back to the probability scale
  test = inla.hyperpar(mod)
  hyperMat = t(inla.hyperpar.sample(10000, mod, improve.marginals=TRUE))
  hyperMat = out$hyperMat
  zSamples = hyperMat[3:(2+nLayer-1),]
  xSamples = apply(zSamples, 2, multivariateExpit)
  xSamples = rbind(xSamples, 1-colSums(xSamples))
  
  for(l in 1:nLayer) {
    pdf(file=paste0("Figures/standardAlpha", l, "Binom.pdf"), width=5, height=5)
    hist(xSamples[l,], xlab=TeX(paste0("$\\alpha_", l, "$")), main=TeX(paste0("Marginal for $\\alpha_", l, "$")), breaks=100, freq=F, xlim=c(0,1))
    abline(v=alphas[l], col="green")
    abline(v=mean(xSamples[l,]), col="purple", lty=1)
    abline(v=quantile(probs=c(.025, .975), xSamples[l,]), col="purple", lty=2)
    dev.off()
  }
  
  invisible(NULL)
}

# tests the fitLK function using data simulated from the LK model
# buffer: buffer distance between domain edge of basis lattice and domain edge of data.
# n: number of observations
# xRange: range of x coordinates
# yRange: range of y coordinates
# nx: number of basis function lattice points in x directions
# NOTE: ny is determined automatically to match scale of x lattice points
# Xmat: design matrix
# ys: observations
# first.time: is first time evaluating function.  User should always set to FALSE
testLKModel = function(buffer=2.5, kappa=1, rho=1, nu=1, seed=1, nLayer=3, nx=20, ny=nx, n=50, sigma2 = .1^2, 
                       nBuffer=5, NC=5, testCovs=TRUE, alphas=NULL, lambdaStart=.1, a.wghtStart=5, maxit=15, 
                       printVerboseTimings=FALSE) {
  set.seed(seed)
  
  # compute alphas, the variance weights for each layer, depending on nu if necessary:
  if(is.null(alphas))
    alphas = getAlphas(nLayer, nu)
  
  # generate lattice and simulate observations
  coords = matrix(runif(2*n), ncol=2)
  xRangeDat = range(coords[,1])
  yRangeDat = range(coords[,2])
  latInfo = makeLatGrids(xRangeDat, yRangeDat, NC, nBuffer, nLayer)
  
  AObs = makeA(coords, latInfo)
  Q = makeQ(kappa=kappa, rho=rho, latInfo, alphas=alphas, normalized=TRUE, fastNormalize=TRUE) 
  L = as.matrix(t(chol(solve(Q))))
  zsims = matrix(rnorm(nrow(Q)), ncol=1)
  fieldSims = L %*% zsims
  ys = as.numeric(AObs %*% fieldSims) + 1 # add a constant unit mean term to be estimated by INLA
  # ys = 1 + as.numeric(AObs %*% fieldSims) + coords[,1] # x-valued mean term to be estimated by INLA
  errs = rnorm(n, sd=sqrt(sigma2))
  ys = ys + errs
  
  # plot the observations
  pdf(file="Figures/standardObservations.pdf", width=5, height=5)
  par(mfrow=c(1,1))
  quilt.plot(coords, ys)
  dev.off()
  
  # make prediction coordinates on a grid
  xRange=c(0,1)
  yRange=c(0,1)
  mx = 20
  my = 20
  predPts = make.surface.grid(list(x=seq(xRange[1], xRange[2], l=mx), y=seq(yRange[1], yRange[2], l=my)))
  
  ## linear term is already included in Lattice Krig model, so not necessary to include in covariates
  # add linear terms in lat/lon to covariate matrices if requested
  # if(testCovs) {
  #   X = cbind(X, coords)
  #   XPred = cbind(XPred, predPts)
  # }
  X = NULL
  XPred = NULL
  
  # fit the model
  time = system.time(out <- fitLKStandard(coords, ys, predCoords=predPts, XObs=X, XPred=XPred, NC=NC, nLayer=nLayer,
                                          nBuffer=nBuffer, lambdaStart=lambdaStart, a.wghtStart=a.wghtStart))
  mod = out$mod
  preds = out$preds
  predSDs = out$sigmas
  LKinfo = out$LKinfo
  lower = out$lower
  upper = out$upper
  
  latInfo=out$latInfo
  latWidth=out$latWidth
  obsPreds=out$obsPreds
  obsSDs=out$obsSDs
  coefPreds = out$coefPreds
  coefSDs = out$coefSDs
  
  # print out the total time
  print(paste0("Total time: ", time[3]))
  
  # show a model summary
  print(summary(mod))
  
  # function for determining if points are in correct range
  inRange = function(pts, rangeShrink=0) {
    inX = (rangeShrink < pts[,1]) & (pts[,1] < 1-rangeShrink)
    inY = (rangeShrink < pts[,2]) & (pts[,2] < 1-rangeShrink)
    inX & inY
  }
  
  # show predictive surface, SD, and data
  
  pdf(file="Figures/standardPreds.pdf", width=15, height=8)
  if(nLayer==1) {
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
  }
  else if(nLayer==2) {
    par(mfrow=c(2,4))
    
    # obsInds = 1:n
    # predInds = (n+1):(n+mx*my)
    # coefInds = (n+mx*my+1):(n+mx*my+nx*ny)
    colRangeDat = range(c(ys-errs, obsPreds, preds, coefPreds))
    colRangeSD = range(c(range(predSDs[inRange(predPts)]), coefSDs[[1]][inRange(gridPtsL1)], 
                         coefSDs[[2]][inRange(gridPtsL2)]))
    gridPtsL1 = latInfo[[1]]$latCoords
    gridPtsL2 = latInfo[[2]]$latCoords
    quilt.plot(coords, ys-errs, main="True Process", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(predPts[,1], predPts[,2], preds, main="Prediction Mean", zlim=colRangeDat, 
               xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL1[,1], gridPtsL1[,2], coefPreds[[1]], main="Basis Coefficient Mean (Layer 1)", xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL2[,1], gridPtsL2[,2], coefPreds[[2]], main="Basis Coefficient Mean (Layer 2)", xlim=xRangeDat, ylim=yRangeDat)
    
    quilt.plot(coords, obsPreds, main="Observation Mean", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(predPts[,1], predPts[,2], predSDs, main="Prediction SD", 
               xlim=xRangeDat, ylim=yRangeDat, zlim=range(predSDs[inRange(predPts)]))
    quilt.plot(gridPtsL1[,1], gridPtsL1[,2], coefSDs[[1]], main="Basis Coefficient SD (Layer 1)", 
               zlim=colRangeSD, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL2[,1], gridPtsL2[,2], coefSDs[[2]], main="Basis Coefficient SD (Layer 2)", 
               zlim=colRangeSD, xlim=xRangeDat, ylim=yRangeDat)
  }
  else if(nLayer==3) {
    par(mfrow=c(2,5))
    
    # obsInds = 1:n
    # predInds = (n+1):(n+mx*my)
    # coefInds = (n+mx*my+1):(n+mx*my+nx*ny)
    colRangeDat = range(c(ys-errs, obsPreds, preds, coefPreds))
    colRangeSD = range(c(range(predSDs[inRange(predPts)]), coefSDs[[1]][inRange(gridPtsL1)], 
                         coefSDs[[2]][inRange(gridPtsL2)], coefSDs[[3]][inRange(gridPtsL3)]))
    gridPtsL1 = latInfo[[1]]$latCoords
    gridPtsL2 = latInfo[[2]]$latCoords
    gridPtsL3 = latInfo[[3]]$latCoords
    quilt.plot(coords, ys-errs, main="True Process", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(predPts[,1], predPts[,2], preds, main="Prediction Mean", zlim=colRangeDat, 
               xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL1[,1], gridPtsL1[,2], coefPreds[[1]], main="Basis Coefficient Mean (Layer 1)", xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL2[,1], gridPtsL2[,2], coefPreds[[2]], main="Basis Coefficient Mean (Layer 2)", xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL3[,1], gridPtsL3[,2], coefPreds[[3]], main="Basis Coefficient Mean (Layer 3)", xlim=xRangeDat, ylim=yRangeDat)
    
    quilt.plot(coords, obsPreds, main="Observation predictions", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(predPts[,1], predPts[,2], predSDs, main="Prediction SD", 
               xlim=xRangeDat, ylim=yRangeDat, zlim=range(predSDs[inRange(predPts)]))
    quilt.plot(gridPtsL1[,1], gridPtsL1[,2], coefSDs[[1]], main="Basis Coefficient SD (Layer 1)", zlim=colRangeSD, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL2[,1], gridPtsL2[,2], coefSDs[[2]], main="Basis Coefficient SD (Layer 2)", zlim=colRangeSD, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL3[,1], gridPtsL3[,2], coefSDs[[3]], main="Basis Coefficient SD (Layer 3)", zlim=colRangeSD, xlim=xRangeDat, ylim=yRangeDat)
  }
  dev.off()
  
  # calculate true effective range and marginal variance:
  latticeWidth = latInfo[[1]]$latWidth
  effRange = sqrt(8)/kappa * latticeWidth
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
  pdf(file="Figures/standardEffRange.pdf", width=5, height=5)
  plot(effRangeMarg, type="l", main="Marginal for effective range")
  abline(v=effRange, col="green")
  abline(v=inla.qmarginal(c(.025, .975), effRangeMarg), col="purple", lty=2)
  dev.off()
  # plot(mod$marginals.hyperpar$`Theta1 for field`, type="l", main="Marginal for log range")
  pdf(file="Figures/standardVar.pdf", width=5, height=5)
  plot(varMarg, type="l", main="Marginal for process variance")
  abline(v=marginalVar, col="green")
  abline(v=inla.qmarginal(c(.025, .975), varMarg), col="purple", lty=2)
  dev.off()
  # plot(mod$marginals.hyperpar$`Theta2 for field`, type="l", main="Marginal for log variance")
  pdf(file="Figures/standardSigma2.pdf", width=5, height=5)
  plot(sigma2Marg, type="l", main="Marginal for error variance")
  abline(v=sigma2, col="green")
  abline(v=inla.qmarginal(c(.025, .975), sigma2Marg), col="purple", lty=2)
  dev.off()
  for(i in 1:length(covNames)) {
    XMarginal = XMarginals[[i]]
    pdf(file=paste0("Figures/standard", covNames[i], ".pdf"), width=5, height=5)
    plot(XMarginal, type="l", main="Marginal for fixed effect")
    if(i==1)
      abline(v=1, col="green")
    else
      abline(v=0, col="green")
    abline(v=inla.qmarginal(c(.025, .975), XMarginal), col="purple", lty=2)
    dev.off()
  }
  
  # pdf(file="Figures/standardRho.pdf", width=5, height=5)
  # plot(sigma2Marg, type="l", main=TeX("Marginal for $\\rho$"), xlab=TeX("$\\rho$"))
  # abline(v=rho, col="green")
  # dev.off()
  
  
  # do the same for kappa, rho
  # in order to get distribution for rho, must sample from joint hyperparameters
  kappaMarg = inla.tmarginal(function(x) {sqrt(8)/exp(x) * latticeWidth}, mod$marginals.hyperpar$`Theta1 for field`)
  # thetasToRho = function(xs) {
  #   logCor = xs[2]
  #   logVar = xs[3]
  #   kappa = sqrt(8)/exp(logCor) * latticeWidth
  #   sigma2 = exp(logVar)
  #   sigma2 * 4*pi * kappa^2
  # }
  # samples = inla.hyperpar.sample(50000, mod, TRUE)
  # rhos = apply(samples, 1, thetasToRho)
  
  pdf(file="Figures/standardKappa.pdf", width=5, height=5)
  plot(kappaMarg, type="l", xlab="kappa", main="Marginal for kappa")
  abline(v=kappa, col="green")
  abline(v=inla.qmarginal(c(.025, .975), kappaMarg), col="purple", lty=2)
  dev.off()
  
  # pdf(file="Figures/standardRho.pdf", width=5, height=5)
  # hist(rhos, xlab="rho", main="Marginal for Rho", breaks=1000, freq=F, xlim=c(0, quantile(probs=.95, rhos)))
  # abline(v=rho, col="green")
  # dev.off()
  
  ## Now generate marginals for the alpha parameters. In order to do this, we must generate draws from 
  ## the posterior, and transform them back to the probability scale
  out = inla.hyperpar.sample(10000, mod, improve.marginals=TRUE)
  zSamples = out[,4:(3+nLayer-1)]
  xSamples = apply(zSamples, 1, multivariateExpit)
  xSamples = rbind(xSamples, 1-colSums(xSamples))
  
  for(l in 1:nLayer) {
    pdf(file=paste0("Figures/standardAlpha", l, ".pdf"), width=5, height=5)
    hist(xSamples[l,], xlab=TeX(paste0("$\\alpha_", l, "$")), main=TeX(paste0("Marginal for $\\alpha_", l, "$")), breaks=100, freq=F, xlim=c(0,1))
    abline(v=alphas[l], col="green")
    abline(v=mean(xSamples[l,]), col="purple", lty=1)
    abline(v=quantile(probs=c(.025, .975), xSamples[l,]), col="purple", lty=2)
    dev.off()
  }
  
  invisible(NULL)
}

# tests the LKINLA model applied to Kenya secondary education prevalence data
# buffer: buffer distance between domain edge of basis lattice and domain edge of data.
# n: number of observations
# xRange: range of x coordinates
# yRange: range of y coordinates
# nx: number of basis function lattice points in x directions
# NOTE: ny is determined automatically to match scale of x lattice points
# Xmat: design matrix
# ys: observations
# first.time: is first time evaluating function.  User should always set to FALSE
testLKINLAKenyaDat = function(seed=1, nLayer=3, NC=13, nBuffer=5, 
                              normalize=TRUE, fastNormalize=TRUE, urbanEffect=TRUE, clusterEffect=TRUE, 
                              printVerboseTimings=FALSE, latInfo=NULL, 
                              plotNameRoot="", effRangeRange=NULL, urbanOverSamplefrac=0, 
                              intStrategy="ccd", strategy="gaussian", 
                              targetPop=c("women", "children"), separateRanges=FALSE, dropUrban=FALSE, 
                              clusterLevel=TRUE, pixelLevel=TRUE, countyLevel=TRUE, regionLevel=TRUE, 
                              family=c("binomial", "betabinomial")) {
  set.seed(seed)
  targetPop = match.arg(targetPop)
  family = match.arg(family)
  
  # make sure input arguments are compatible
  if(dropUrban) {
    if(urbanEffect)
      warning("dropUrban and urbanEffect both set to TRUE. Setting urbanEffect to FALSE...")
    urbanEffect=FALSE
  }
  if(family == "betabinomial") {
    if(clusterEffect)
      stop("cluster effect must not be set to TRUE for betaBinomial model")
  }
  # if(separateRanges) {
    # if(clusterEffect)
    #   warning("separateRanges and clusterEffect both set to TRUE. Setting clusterEffect to FALSE...")
    # clusterEffect=FALSE
  # }
  
  # set default resolution if using separate layer ranges
  if(length(NC) == 1) {
    if(separateRanges)
      NC = c(30, 107) # by default, use two layers with the finest layer having resolution equal to 10km
  }
  if(separateRanges)
    nLayer = length(NC)
  
  # set plotNameRoot
  ncText = ""
  if(length(NC) == 1) {
    if(separateRanges)
      ncText = "_NC30_107"
    else {
      ncText = paste0("_NC", NC)
    }
  } else {
    tempText = do.call("paste0", as.list(c(NC[1], paste0("_", NC[-1]))))
    ncText = paste0("_NC", tempText)
  }
  if(family == "betabinomial")
    familyText = "_BBin"
  else
    familyText = "_Bin"
  if(targetPop == "women") {
    dataType = "ed"
    plotNameRoot = paste0(plotNameRoot, "_Ed", familyText, "_L", nLayer, ncText, "_sepRange", separateRanges, 
                          "_urb", urbanEffect, "_clust", clusterEffect, "_dropUrb", dropUrban)
  }
  else {
    dataType = "mort"
    plotNameRoot = paste0(plotNameRoot, "_Mort", familyText, "_L", nLayer, ncText, "_sepRange", separateRanges, 
                          "_urb", urbanEffect, "_clust", clusterEffect, "_dropUrb", dropUrban)
  }
  
  # load data set
  if(dataType == "mort") {
    out = load("../U5MR/kenyaData.RData")
    dat = mort
  }
  else {
    out = load("../U5MR/kenyaDataEd.RData")
    dat = ed
  }
  if(dropUrban) {
    dat = dat[!dat$urban,]
  }
  coords = cbind(dat$east, dat$north)
  ys = dat$y
  ns = dat$n
  obsUrban = dat$urban
  
  # plot the observations
  pdf(file=paste0("Figures/LKINLAObservations", plotNameRoot, ".pdf"), width=5, height=5)
  par(mfrow=c(1,1))
  quilt.plot(coords, ys/ns, FUN=function(x) {mean(x, na.rm=TRUE)})
  dev.off()
  
  # create default lattice basis
  latInfo = makeLatGridsKenya(nLayer=nLayer, NC=NC, nBuffer=nBuffer)
  
  # generate hyperparameters based on median and quantiles of inverse exponential and inverse gamma
  # priorPar = getPrior(.1, .1, 10)
  # generate hyperparameters for pc priors
  # median effective range is .4 or 200 for kenya data (a fifth of the spatial domain diameter), median spatial variance is 1
  priorPar = getPCPrior(200, .01, 1, nLayer=nLayer, separateRanges=separateRanges, latticeInfo=latInfo)
  
  # show priors on effective correlation, marginal variance, and error variance:
  xs1 = seq(1, 500, l=500)
  if(!separateRanges) {
    pdf(file=paste0("Figures/LKINLAPriorEffRange", plotNameRoot, ".pdf"), width=5, height=5)
    plot(xs1, dinvexp(xs1, rate=priorPar$corScalePar), type="l", col="blue", 
         xlab="Effective Correlation Range", main="Effective Correlation Prior", 
         ylab="Prior Density")
    abline(v=qinvexp(.5, rate=priorPar$corScalePar), col="red")
    dev.off()
  } else {
    for(i in 1:nLayer) {
      if(i == nLayer && identical(NC, c(30, 107)))
        xs1 = seq(1, 200, l=500)
      pdf(file=paste0("Figures/LKINLAPriorEff", i, "Range", plotNameRoot, ".pdf"), width=5, height=5)
      plot(xs1, dinvexp(xs1, rate=priorPar$corScalePar[i]), type="l", col="blue", 
           xlab="Effective Correlation Range", main="Effective Correlation Prior", 
           ylab="Prior Density")
      abline(v=qinvexp(.5, rate=priorPar$corScalePar[i]), col="red")
      dev.off()
    }
  }
  
  if(priorPar$priorType == "orig") {
    xs2 = seq(.01, 10.5, l=500)
    pdf(file=paste0("Figures/LKINLAPriorMargVar", plotNameRoot, ".pdf"), width=5, height=5)
    plot(xs2, invgamma::dinvgamma(xs2, shape=priorPar$varPar1, rate=priorPar$varPar2), type="l", col="blue", 
         xlab="Marginal Variance", main="Marginal Variance Prior", 
         ylab="Prior Density")
    abline(v=qinvgamma(.1, shape=priorPar$varPar1, rate=priorPar$varPar2), col="red")
    abline(v=qinvgamma(.9, shape=priorPar$varPar1, rate=priorPar$varPar2), col="red")
    dev.off()
  } else if(priorPar$priorType == "pc") {
    xs2 = seq(.01, 11.5, l=500)
    pdf(file=paste0("Figures/LKINLAPriorMargVar", plotNameRoot, ".pdf"), width=5, height=5)
    plot(xs2, dpcvar(xs2, alpha=priorPar$alpha, u=priorPar$u), type="l", col="blue", 
         xlab="Marginal Variance", main="Marginal Variance Prior", 
         ylab="Prior Density")
    abline(v=qpcvar(.1, alpha=priorPar$alpha, u=priorPar$u), col="red")
    abline(v=qpcvar(.9, alpha=priorPar$alpha, u=priorPar$u), col="red")
    abline(v=1, col="green")
    dev.off()
  }
  
  # xs2 = seq(.001, invgamma::qinvgamma(.905, shape=0.1, rate=0.1), l=500)
  # pdf(file="Figures/mixtureLKINLAPriorErrorVar.pdf", width=5, height=5)
  # plot(xs2, invgamma::dinvgamma(xs2, shape=0.1, rate=0.1), type="l", col="blue", 
  #      xlab="Error Variance", main="Error Variance Prior", 
  #      ylab="Prior Density")
  # abline(v=invgamma::qinvgamma(.1, shape=0.1, rate=0.1), col="red")
  # abline(v=invgamma::qinvgamma(.9, shape=0.1, rate=0.1), col="red")
  # dev.off()
  
  xs2 = seq(.01, 1, l=500)
  pdf(file=paste0("Figures/LKINLAPriorErrorVar", plotNameRoot, ".pdf"), width=5, height=5)
  plot(xs2, dpcvar(xs2, alpha=.05, u=1), type="l", col="blue", 
       xlab="Marginal Variance", main="Marginal Variance Prior", 
       ylab="Prior Density")
  abline(v=qpcvar(.1, alpha=.05, u=1), col="red")
  abline(v=qpcvar(.9, alpha=.05, u=1), col="red")
  abline(v=sqrt(.1), col="green")
  dev.off()
  
  # browser()
  
  
  for(l in 1:nLayer) {
    pdf(file=paste0("Figures/LKINLAPriorAlpha", l, plotNameRoot, ".pdf"), width=5, height=5)
    xs = seq(0, 1, l=500)
    tempYs = dbeta(xs, priorPar$alphaPar[l], sum(priorPar$alphaPar[-l]))
    plot(xs, tempYs, type="l", xlab=TeX(paste0("$\\alpha_", l, "$")), ylab="Density", 
         main=TeX(paste0("Prior for $\\alpha_", l, "$")), xlim=c(0,1), ylim=c(0, max(tempYs[is.finite(tempYs)])))
    abline(v=qbeta(c(0.025, 0.975), priorPar$alphaPar[l], sum(priorPar$alphaPar[-l])), col="purple", lty=2)
    dev.off()
  }
  
  # prior of overdispersion
  if(family == "betabinomial") {
    # lambda = getLambdapcBeta(U=1, logitU=TRUE, alpha=0.01, p=.5, normalize=TRUE)
    # set median at .04 and upper 97.5th pctile at 0.2
    mu = logit(0.04)
    prec = 1/((logit(.2)-logit(.04))/qnorm(.975))^2
    rhos = seq(0, 1, l=1000)
    # tempYs = dpcBeta2(rhos, lambda=lambda, normalize=TRUE)
    tempYs = dlogitNormal(rhos, mu=mu, prec=prec)
    pdf(file=paste0("Figures/LKINLAPriorRho", plotNameRoot, ".pdf"), width=5, height=5)
    plot(rhos, tempYs, type="l", xlab=TeX(paste0("$\\rho$")), ylab="Density", 
         main=TeX(paste0("Overdispersion Prior")), xlim=c(0,1), ylim=c(0, max(tempYs[is.finite(tempYs)])))
    # abline(v=qpcBeta2(c(0.025, 0.975), lambda=lambda, normalize=TRUE), col="purple", lty=2)
    abline(v=expit(qnorm(c(0.025, 0.975), mu, sqrt(1/prec))), col="purple", lty=2)
    dev.off()
  }
  
  # prior on covariogram
  alphaVals = t(rdirichlet(100, alpha=priorPar$alphaPar))
  rhoVals = rpcvar(100, alpha=priorPar$alpha, u=priorPar$u)
  effectiveRangeVals = t(matrix(rinvexp(100*length(priorPar$corScalePar), rate=priorPar$corScalePar), nrow=length(priorPar$corScalePar)))
  if(separateRanges) {
    latticeWidth = sapply(latInfo, function(x) {x$latWidth})
    kappaVals = t(sweep(sqrt(8)/effectiveRangeVals, 2, latticeWidth, "*"))
  } else {
    latticeWidth = latInfo[[1]]$latWidth
    kappaVals = c(sqrt(8)/exp(effectiveRangeVals) * latticeWidth)
  }
  if(clusterEffect)
    nuggetVarVals = rpcvar(100, alpha=.05, u=1)
  else
    nuggetVarVals = rep(0, 100)
  out = covarianceDistributionLKINLA(latInfo, kappaVals, rhoVals, nuggetVarVals, alphaVals, 
                                     normalize=normalize, fastNormalize=fastNormalize)
  d = out$d
  sortI = sort(d, index.return=TRUE)$ix
  d = d[sortI]
  covMean = out$cov[sortI]
  upperCov=out$upperCov[sortI]
  lowerCov=out$lowerCov[sortI]
  covMat=out$covMat[sortI]
  corMean = out$cor[sortI]
  upperCor=out$upperCor[sortI]
  lowerCor=out$lowerCor[sortI]
  corMat=out$corMat[sortI]
  
  # correlation and covariance functions
  
  # plot the covariance an correlation priors
  yRange = range(c(covMean, lowerCov, upperCov))
  pdf(file=paste0("Figures/LKINLAPriorCov", plotNameRoot, ".pdf"), width=5, height=5)
  plot(d, covMean, type="l", main="Prior of covariance function", xlab="Distance", ylab="Covariance", 
       ylim=yRange)
  lines(d, lowerCov, lty=2)
  lines(d, upperCov, lty=2)
  legend("topright", c("Estimate", "80% CI"), lty=c(1, 2), col=c("black", "black"))
  dev.off()
  
  pdf(file=paste0("Figures/LKINLAPriorCor", plotNameRoot, ".pdf"), width=5, height=5)
  plot(d, corMean, type="l", main="Prior of correlation function", xlab="Distance", ylab="Covariance", 
       ylim=c(0,1))
  lines(d, lowerCor, lty=2)
  lines(d, upperCor, lty=2)
  legend("topright", c("Estimate", "80% CI"), lty=c(1, 2), col=c("black", "black"))
  dev.off()
  
  # fit the model
  time = system.time(fit <- fitLKINLAKenyaDat(dat, dataType, nu=1, seed=seed, nLayer=nLayer, NC=NC, nBuffer=nBuffer, 
                                              priorPar=priorPar, normalize=normalize, fastNormalize=fastNormalize, 
                                              latInfo=latInfo, urbanEffect=urbanEffect, clusterEffect=clusterEffect, 
                                              intStrategy=intStrategy, strategy=strategy, 
                                              printVerboseTimings=printVerboseTimings, separateRanges=separateRanges, 
                                              family=family))
  mod = fit$mod
  preds=fit$preds
  predSDs=fit$sigmas
  latInfo=fit$latInfo
  latWidth=fit$latWidth
  obsPreds=fit$obsPreds
  obsSDs=fit$obsSDs
  coefPreds = fit$coefPreds
  coefSDs = fit$coefSDs
  predPts = fit$predPts
  predsUrban = fit$predsUrban
  
  # print out the total time
  print(paste0("Total time before aggregation: ", time[3]))
  
  # aggregate the model
  aggregationTime = 
    system.time(aggregatedFit <- 
                  aggregateModelResultsKenya(fit, clusterLevel=clusterLevel, pixelLevel=pixelLevel, 
                                             countyLevel=countyLevel, regionLevel=regionLevel, 
                                             targetPop=targetPop)
                )[3]
  allResults = list(fit=fit, aggregatedResults=aggregatedFit)
  
  print(paste0("Time to aggregate: ", aggregationTime))
  
  # show a model summary
  print(summary(mod))
  
  # function for determining if points are in correct range
  inRange = function(pts, rangeShrink=0) {
    rep(TRUE, nrow(pts))
  }
  
  xRangeDat = latInfo[[1]]$xRangeDat
  yRangeDat = latInfo[[1]]$yRangeDat
  
  # show predictive surface, SD, and data
  out = load("../U5MR/adminMapData.RData")
  # kenyaMap = adm0
  kenyaMap = adm1
  if(nLayer==1) {
    pdf(file=paste0("Figures/LKINLAPreds", plotNameRoot, ".pdf"), width=9, height=6)
    par(mfrow=c(2,3))
    
    # obsInds = 1:n
    # predInds = (n+1):(n+mx*my)
    # coefInds = (n+mx*my+1):(n+mx*my+nx*ny)
    colRangeDat = range(c(ys/ns, obsPreds, preds))
    colRangeCoef = range(c(coefPreds))
    colRangeSD = range(c(predSDs, obsSDs))
    colRangeSDCoef = range(c(coefSDs))
    gridPtsL1 = latInfo[[1]]$latCoords
    quilt.plot(coords, ys/ns, main="True Process", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    plotMapDat(mapDat=kenyaMap, lwd=.5, project=TRUE)
    quilt.plot(gridPtsL1[,1], gridPtsL1[,2], coefPreds[[1]], main="Basis Coefficient Mean (Layer 1)", xlim=xRangeDat, ylim=yRangeDat)
    plotMapDat(mapDat=kenyaMap, lwd=.5, project=TRUE)
    quilt.plot(gridPtsL1[,1], gridPtsL1[,2], coefSDs[[1]], main="Basis Coefficient SD (Layer 1)", zlim=colRangeSDCoef, xlim=xRangeDat, ylim=yRangeDat)
    plotMapDat(mapDat=kenyaMap, lwd=.5, project=TRUE)
    
    # quilt.plot(coords, obsPreds, main="Prediction mean", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat, 
    #            zlim=range(predSDs[inRange(predPts)]))
    plot.new()
    quilt.plot(predPts[,1], predPts[,2], preds, main="Prediction Mean", zlim=colRangeDat, 
               xlim=xRangeDat, ylim=yRangeDat)
    plotMapDat(mapDat=kenyaMap, lwd=.5, project=TRUE)
    quilt.plot(predPts[,1], predPts[,2], predSDs, main="Prediction SD", 
               xlim=xRangeDat, ylim=yRangeDat)
    plotMapDat(mapDat=kenyaMap, lwd=.5, project=TRUE)
  }
  else if(nLayer==2) {
    pdf(file=paste0("Figures/LKINLAPreds", plotNameRoot, ".pdf"), width=12, height=6)
    par(mfrow=c(2,4), mar=c(5.1, 4.1, 4.1, 6))
    
    # obsInds = 1:n
    # predInds = (n+1):(n+mx*my)
    # coefInds = (n+mx*my+1):(n+mx*my+nx*ny)
    colRangeDat = range(c(ys/ns, obsPreds, preds))
    colRangeCoef = range(c(coefPreds))
    colRangeSD = range(c(predSDs, obsSDs))
    colRangeSDCoef = range(c(coefSDs))
    gridPtsL1 = latInfo[[1]]$latCoords
    gridPtsL2 = latInfo[[2]]$latCoords
    quilt.plot(coords, ys/ns, main="Empirical Proportions", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    plotMapDat(mapDat=kenyaMap, lwd=.5, project=TRUE)
    quilt.plot(predPts[,1], predPts[,2], preds, main="Prediction Mean", 
               xlim=xRangeDat, ylim=yRangeDat, nx=100, ny=100)
    plotMapDat(mapDat=kenyaMap, lwd=.5, project=TRUE)
    quilt.plot(gridPtsL1[,1], gridPtsL1[,2], coefPreds[[1]], main="Basis Coefficient Mean (Layer 1)", xlim=xRangeDat, ylim=yRangeDat, 
               zlim=colRangeCoef)
    plotMapDat(mapDat=kenyaMap, lwd=.5, project=TRUE)
    quilt.plot(gridPtsL2[,1], gridPtsL2[,2], coefPreds[[2]], main="Basis Coefficient Mean (Layer 2)", xlim=xRangeDat, ylim=yRangeDat, 
               zlim=colRangeCoef, nx=100, ny=100)
    plotMapDat(mapDat=kenyaMap, lwd=.5, project=TRUE)
    
    quilt.plot(coords, obsPreds, main="Observation Predictions", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    plotMapDat(mapDat=kenyaMap, lwd=.5, project=TRUE)
    quilt.plot(predPts[,1], predPts[,2], predSDs, main="Prediction SD", 
               xlim=xRangeDat, ylim=yRangeDat, zlim=colRangeSD, nx=100, ny=100)
    plotMapDat(mapDat=kenyaMap, lwd=.5, project=TRUE)
    quilt.plot(gridPtsL1[,1], gridPtsL1[,2], coefSDs[[1]], main="Basis Coefficient SD (Layer 1)", zlim=colRangeSDCoef, xlim=xRangeDat, ylim=yRangeDat)
    plotMapDat(mapDat=kenyaMap, lwd=.5, project=TRUE)
    quilt.plot(gridPtsL2[,1], gridPtsL2[,2], coefSDs[[2]], main="Basis Coefficient SD (Layer 2)", zlim=colRangeSDCoef, xlim=xRangeDat, ylim=yRangeDat, 
               nx=100, ny=100)
    plotMapDat(mapDat=kenyaMap, lwd=.5, project=TRUE)
  }
  else if(nLayer==3) {
    pdf(file=paste0("Figures/LKINLAPreds", plotNameRoot, ".pdf"), width=15, height=6)
    par(mfrow=c(2,5), mar=c(5.1, 4.1, 4.1, 6))
    
    # obsInds = 1:n
    # predInds = (n+1):(n+mx*my)
    # coefInds = (n+mx*my+1):(n+mx*my+nx*ny)
    # colRangeDat = range(c(ys, obsPreds, preds, coefPreds))
    colRangeDat = range(c(ys/ns, obsPreds, preds))
    colRangeCoef = range(c(coefPreds))
    colRangeSD = range(c(predSDs, obsSDs))
    colRangeSDCoef = range(c(coefSDs))
    gridPtsL1 = latInfo[[1]]$latCoords
    gridPtsL2 = latInfo[[2]]$latCoords
    gridPtsL3 = latInfo[[3]]$latCoords
    quilt.plot(coords, ys/ns, main="True Process", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    plotMapDat(mapDat=kenyaMap, lwd=.5, project=TRUE)
    quilt.plot(predPts[,1], predPts[,2], preds, main="Prediction Mean", zlim=colRangeDat, 
               xlim=xRangeDat, ylim=yRangeDat)
    plotMapDat(mapDat=kenyaMap, lwd=.5, project=TRUE)
    quilt.plot(gridPtsL1[,1], gridPtsL1[,2], coefPreds[[1]], main="Basis Coefficient Mean (Layer 1)", 
               xlim=xRangeDat, ylim=yRangeDat, zlim=colRangeCoef)
    plotMapDat(mapDat=kenyaMap, lwd=.5, project=TRUE)
    quilt.plot(gridPtsL2[,1], gridPtsL2[,2], coefPreds[[2]], main="Basis Coefficient Mean (Layer 2)", 
               xlim=xRangeDat, ylim=yRangeDat, zlim=colRangeCoef)
    plotMapDat(mapDat=kenyaMap, lwd=.5, project=TRUE)
    quilt.plot(gridPtsL3[,1], gridPtsL3[,2], coefPreds[[3]], main="Basis Coefficient Mean (Layer 3)", 
               xlim=xRangeDat, ylim=yRangeDat, zlim=colRangeCoef)
    plotMapDat(mapDat=kenyaMap, lwd=.5, project=TRUE)
    
    # quilt.plot(coords, obsPreds, main="Observation Mean", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    # plot.new()
    quilt.plot(coords, ys/ns, main="Observations", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    plotMapDat(mapDat=kenyaMap, lwd=.5, project=TRUE)
    quilt.plot(predPts[,1], predPts[,2], predSDs, main="Prediction SD",
               xlim=xRangeDat, ylim=yRangeDat, zlim=range(predSDs[inRange(predPts, rangeShrink=.03)]))
    plotMapDat(mapDat=kenyaMap, lwd=.5, project=TRUE)
    quilt.plot(gridPtsL1[,1], gridPtsL1[,2], coefSDs[[1]], main="Basis Coefficient SD (Layer 1)", zlim=colRangeSDCoef, xlim=xRangeDat, ylim=yRangeDat)
    plotMapDat(mapDat=kenyaMap, lwd=.5, project=TRUE)
    quilt.plot(gridPtsL2[,1], gridPtsL2[,2], coefSDs[[2]], main="Basis Coefficient SD (Layer 2)", zlim=colRangeSDCoef, xlim=xRangeDat, ylim=yRangeDat)
    plotMapDat(mapDat=kenyaMap, lwd=.5, project=TRUE)
    quilt.plot(gridPtsL3[,1], gridPtsL3[,2], coefSDs[[3]], main="Basis Coefficient SD (Layer 3)", zlim=colRangeSDCoef, xlim=xRangeDat, ylim=yRangeDat)
    plotMapDat(mapDat=kenyaMap, lwd=.5, project=TRUE)
  }
  else if(nLayer==4) {
    pdf(file=paste0("Figures/LKINLAPreds", plotNameRoot, ".pdf"), width=18, height=6)
    par(mfrow=c(2,6), mar=c(5.1, 4.1, 4.1, 6))
    
    # obsInds = 1:n
    # predInds = (n+1):(n+mx*my)
    # coefInds = (n+mx*my+1):(n+mx*my+nx*ny)
    # colRangeDat = range(c(ys, obsPreds, preds, coefPreds))
    colRangeDat = range(c(ys/ns, obsPreds, preds))
    colRangeCoef = range(c(coefPreds))
    colRangeSD = range(c(predSDs, obsSDs))
    colRangeSDCoef = range(c(coefSDs))
    gridPtsL1 = latInfo[[1]]$latCoords
    gridPtsL2 = latInfo[[2]]$latCoords
    gridPtsL3 = latInfo[[3]]$latCoords
    gridPtsL4 = latInfo[[4]]$latCoords
    quilt.plot(coords, ys/ns, main="True Process", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    plotMapDat(mapDat=kenyaMap, lwd=.5, project=TRUE)
    quilt.plot(predPts[,1], predPts[,2], preds, main="Prediction Mean", zlim=colRangeDat, 
               xlim=xRangeDat, ylim=yRangeDat)
    plotMapDat(mapDat=kenyaMap, lwd=.5, project=TRUE)
    quilt.plot(gridPtsL1[,1], gridPtsL1[,2], coefPreds[[1]], main="Basis Coefficient Mean (Layer 1)", 
               xlim=xRangeDat, ylim=yRangeDat, zlim=colRangeCoef)
    plotMapDat(mapDat=kenyaMap, lwd=.5, project=TRUE)
    quilt.plot(gridPtsL2[,1], gridPtsL2[,2], coefPreds[[2]], main="Basis Coefficient Mean (Layer 2)", 
               xlim=xRangeDat, ylim=yRangeDat, zlim=colRangeCoef)
    plotMapDat(mapDat=kenyaMap, lwd=.5, project=TRUE)
    quilt.plot(gridPtsL3[,1], gridPtsL3[,2], coefPreds[[3]], main="Basis Coefficient Mean (Layer 3)", 
               xlim=xRangeDat, ylim=yRangeDat, zlim=colRangeCoef)
    plotMapDat(mapDat=kenyaMap, lwd=.5, project=TRUE)
    quilt.plot(gridPtsL4[,1], gridPtsL4[,2], coefPreds[[4]], main="Basis Coefficient Mean (Layer 4)", 
               xlim=xRangeDat, ylim=yRangeDat, zlim=colRangeCoef)
    plotMapDat(mapDat=kenyaMap, lwd=.5, project=TRUE)
    
    # quilt.plot(coords, obsPreds, main="Observation Mean", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    # plot.new()
    quilt.plot(coords, ys/ns, main="Observations", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(predPts[,1], predPts[,2], predSDs, main="Prediction SD",
               xlim=xRangeDat, ylim=yRangeDat, zlim=range(predSDs[inRange(predPts, rangeShrink=.03)]))
    plotMapDat(mapDat=kenyaMap, lwd=.5, project=TRUE)
    quilt.plot(gridPtsL1[,1], gridPtsL1[,2], coefSDs[[1]], main="Basis Coefficient SD (Layer 1)", zlim=colRangeSDCoef, xlim=xRangeDat, ylim=yRangeDat)
    plotMapDat(mapDat=kenyaMap, lwd=.5, project=TRUE)
    quilt.plot(gridPtsL2[,1], gridPtsL2[,2], coefSDs[[2]], main="Basis Coefficient SD (Layer 2)", zlim=colRangeSDCoef, xlim=xRangeDat, ylim=yRangeDat)
    plotMapDat(mapDat=kenyaMap, lwd=.5, project=TRUE)
    quilt.plot(gridPtsL3[,1], gridPtsL3[,2], coefSDs[[3]], main="Basis Coefficient SD (Layer 3)", zlim=colRangeSDCoef, xlim=xRangeDat, ylim=yRangeDat)
    plotMapDat(mapDat=kenyaMap, lwd=.5, project=TRUE)
    quilt.plot(gridPtsL4[,1], gridPtsL4[,2], coefSDs[[4]], main="Basis Coefficient SD (Layer 4)", zlim=colRangeSDCoef, xlim=xRangeDat, ylim=yRangeDat)
    plotMapDat(mapDat=kenyaMap, lwd=.5, project=TRUE)
  }
  else if(nLayer==5) {
    pdf(file=paste0("Figures/LKINLAPreds", plotNameRoot, ".pdf"), width=21, height=6)
    par(mfrow=c(2,7), mar=c(5.1, 4.1, 4.1, 6))
    
    # obsInds = 1:n
    # predInds = (n+1):(n+mx*my)
    # coefInds = (n+mx*my+1):(n+mx*my+nx*ny)
    # colRangeDat = range(c(ys, obsPreds, preds, coefPreds))
    colRangeDat = range(c(ys/ns, obsPreds, preds))
    colRangeCoef = range(c(coefPreds))
    colRangeSD = range(c(predSDs, obsSDs))
    colRangeSDCoef = range(c(coefSDs))
    gridPtsL1 = latInfo[[1]]$latCoords
    gridPtsL2 = latInfo[[2]]$latCoords
    gridPtsL3 = latInfo[[3]]$latCoords
    gridPtsL4 = latInfo[[4]]$latCoords
    gridPtsL5 = latInfo[[5]]$latCoords
    quilt.plot(coords, ys/ns, main="True Process", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    plotMapDat(mapDat=kenyaMap, lwd=.5, project=TRUE)
    quilt.plot(predPts[,1], predPts[,2], preds, main="Prediction Mean", zlim=colRangeDat, 
               xlim=xRangeDat, ylim=yRangeDat)
    plotMapDat(mapDat=kenyaMap, lwd=.5, project=TRUE)
    quilt.plot(gridPtsL1[,1], gridPtsL1[,2], coefPreds[[1]], main="Basis Coefficient Mean (Layer 1)", 
               xlim=xRangeDat, ylim=yRangeDat, zlim=colRangeCoef)
    plotMapDat(mapDat=kenyaMap, lwd=.5, project=TRUE)
    quilt.plot(gridPtsL2[,1], gridPtsL2[,2], coefPreds[[2]], main="Basis Coefficient Mean (Layer 2)", 
               xlim=xRangeDat, ylim=yRangeDat, zlim=colRangeCoef)
    plotMapDat(mapDat=kenyaMap, lwd=.5, project=TRUE)
    quilt.plot(gridPtsL3[,1], gridPtsL3[,2], coefPreds[[3]], main="Basis Coefficient Mean (Layer 3)", 
               xlim=xRangeDat, ylim=yRangeDat, zlim=colRangeCoef)
    plotMapDat(mapDat=kenyaMap, lwd=.5, project=TRUE)
    quilt.plot(gridPtsL4[,1], gridPtsL4[,2], coefPreds[[4]], main="Basis Coefficient Mean (Layer 4)", 
               xlim=xRangeDat, ylim=yRangeDat, zlim=colRangeCoef)
    plotMapDat(mapDat=kenyaMap, lwd=.5, project=TRUE)
    quilt.plot(gridPtsL5[,1], gridPtsL5[,2], coefPreds[[5]], main="Basis Coefficient Mean (Layer 5)", 
               xlim=xRangeDat, ylim=yRangeDat, zlim=colRangeCoef)
    plotMapDat(mapDat=kenyaMap, lwd=.5, project=TRUE)
    
    # quilt.plot(coords, obsPreds, main="Observation Mean", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    # plot.new()
    quilt.plot(coords, ys/ns, main="Observations", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    plotMapDat(mapDat=kenyaMap, lwd=.5, project=TRUE)
    quilt.plot(predPts[,1], predPts[,2], predSDs, main="Prediction SD",
               xlim=xRangeDat, ylim=yRangeDat, zlim=range(predSDs[inRange(predPts, rangeShrink=.03)]))
    plotMapDat(mapDat=kenyaMap, lwd=.5, project=TRUE)
    quilt.plot(gridPtsL1[,1], gridPtsL1[,2], coefSDs[[1]], main="Basis Coefficient SD (Layer 1)", zlim=colRangeSDCoef, xlim=xRangeDat, ylim=yRangeDat)
    plotMapDat(mapDat=kenyaMap, lwd=.5, project=TRUE)
    quilt.plot(gridPtsL2[,1], gridPtsL2[,2], coefSDs[[2]], main="Basis Coefficient SD (Layer 2)", zlim=colRangeSDCoef, xlim=xRangeDat, ylim=yRangeDat)
    plotMapDat(mapDat=kenyaMap, lwd=.5, project=TRUE)
    quilt.plot(gridPtsL3[,1], gridPtsL3[,2], coefSDs[[3]], main="Basis Coefficient SD (Layer 3)", zlim=colRangeSDCoef, xlim=xRangeDat, ylim=yRangeDat)
    plotMapDat(mapDat=kenyaMap, lwd=.5, project=TRUE)
    quilt.plot(gridPtsL4[,1], gridPtsL4[,2], coefSDs[[4]], main="Basis Coefficient SD (Layer 4)", zlim=colRangeSDCoef, xlim=xRangeDat, ylim=yRangeDat)
    plotMapDat(mapDat=kenyaMap, lwd=.5, project=TRUE)
    quilt.plot(gridPtsL5[,1], gridPtsL5[,2], coefSDs[[5]], main="Basis Coefficient SD (Layer 5)", zlim=colRangeSDCoef, xlim=xRangeDat, ylim=yRangeDat)
    plotMapDat(mapDat=kenyaMap, lwd=.5, project=TRUE)
  }
  dev.off()
  
  browser()
  plotSingleModelPredictions(dat=dat, allResults, modelName="", targetPop=targetPop, 
                             areaLevels=c("Region", "County", "Pixel", "Cluster"), 
                             plotNameRoot=plotNameRoot, 
                             meanRange=NULL, meanTicks=NULL, meanTickLabels=NULL, widthRange=NULL, widthTicks=NULL, widthTickLabels=NULL, 
                             meanCols=makeRedBlueDivergingColors(64), 
                             widthCols=makeBlueYellowSequentialColors(64), 
                             kenyaLatRange=c(-4.6, 5), kenyaLonRange=c(33.5, 42.0))
  
  pdf(file=paste0("Figures/LKINLAInSampleResidualsLabeled", plotNameRoot, ".pdf"), width=5, height=5)
  ylim = range(ys/ns-obsPreds)
  xlim = range(obsPreds)
  plot(obsPreds, ys/ns-obsPreds, pch=19, cex=.1, main="Residuals versus fitted", 
       ylab="Residuals", xlab="Fitted", xlim=xlim, ylim=ylim, type="n")
  points(obsPreds[!obsUrban], ys[!obsUrban]/ns[!obsUrban]-obsPreds[!obsUrban], pch=19, cex=.1, col="green")
  points(obsPreds[obsUrban], ys[obsUrban]/ns[obsUrban]-obsPreds[obsUrban], pch=19, cex=.1, col="blue")
  abline(h=0, lty=2)
  legend("topright", c("Rural", "Urban"), col=c("green", "blue"), pch=19)
  dev.off()
  
  # calculate true effective range and marginal variance:
  latticeWidth = latInfo[[1]]$latWidth
  if(separateRanges)
    latticeWidth = sapply(latInfo, function(x) {x$latWidth})
  
  # plot marginals on interpretable scale (effective range, marginal variance)
  if(!separateRanges) {
    effRangeMarg = inla.tmarginal(exp, mod$marginals.hyperpar$`Theta1 for field`)
    varMarg = inla.tmarginal(exp, mod$marginals.hyperpar$`Theta2 for field`)
  } else {
    numberThetas = nLayer + 1 + nLayer - 1
    allNames = paste0("Theta", 1:numberThetas, " for field")
    effRangeMargs = list()
    for(i in 1:nLayer) {
      effRangeMargs = c(effRangeMargs, list(inla.tmarginal(exp, mod$marginals.hyperpar[[allNames[i]]])))
    }
    varMarg = inla.tmarginal(exp, mod$marginals.hyperpar[[allNames[nLayer+1]]])
  }
  if(clusterEffect)
    sigma2Marg = inla.tmarginal(function(x) {1/x}, mod$marginals.hyperpar[[length(mod$marginals.hyperpar)]])
  else if(family == "betabinomial")
    overdispersionMarg = mod$marginals.hyperpar$`overdispersion for the betabinomial observations`
  covNames = names(mod$marginals.fixed)
  XMarginals = list()
  for(i in 1:length(covNames)) {
    XMarginal = inla.tmarginal(function(x) {x}, mod$marginals.fixed[[covNames[i]]])
    XMarginals = c(XMarginals, list(XMarginal))
  }
  
  par(mfrow=c(1,1))
  
  if(!separateRanges) {
    pdf(file=paste0("Figures/LKINLAEffRange", plotNameRoot, ".pdf"), width=5, height=5)
    plot(effRangeMarg, type="l", main="Marginal for effective range")
    abline(v=inla.qmarginal(c(.025, .975), effRangeMarg), col="purple", lty=2)
    dev.off()
  } else {
    for(i in 1:nLayer) {
      pdf(file=paste0("Figures/LKINLAEffRange", i, plotNameRoot, ".pdf"), width=5, height=5)
      plot(effRangeMargs[[i]], type="l", main="Marginal for effective range")
      abline(v=inla.qmarginal(c(.025, .975), effRangeMargs[[i]]), col="purple", lty=2)
      dev.off()
    }
  }
  
  # plot(mod$marginals.hyperpar$`Theta1 for field`, type="l", main="Marginal for log range")
  pdf(file=paste0("Figures/LKINLAVar", plotNameRoot, ".pdf"), width=5, height=5)
  plot(varMarg, type="l", main="Marginal for spatial variance")
  abline(v=inla.qmarginal(c(.025, .975), varMarg), col="purple", lty=2)
  dev.off()
  # plot(mod$marginals.hyperpar$`Theta2 for field`, type="l", main="Marginal for log variance")
  if(clusterEffect) {
    pdf(file=paste0("Figures/LKINLASigma2", plotNameRoot, ".pdf"), width=5, height=5)
    plot(sigma2Marg, type="l", main="Marginal for error variance")
    abline(v=inla.qmarginal(c(.025, .975), sigma2Marg), col="purple", lty=2)
    dev.off()
  } else if(family == "betabinomial") {
    pdf(file=paste0("Figures/LKINLARho", plotNameRoot, ".pdf"), width=5, height=5)
    plot(overdispersionMarg, type="l", main="Marginal for overdispersion")
    abline(v=inla.qmarginal(c(.025, .975), overdispersionMarg), col="purple", lty=2)
    dev.off()
  }
  
  for(i in 1:length(covNames)) {
    XMarginal = XMarginals[[i]]
    pdf(file=paste0("Figures/LKINLA", covNames[i], plotNameRoot, ".pdf"), width=5, height=5)
    plot(XMarginal, type="l", main=paste0("Marginal for fixed effect ", i))
    abline(v=inla.qmarginal(c(.025, .975), XMarginal), col="purple", lty=2)
    dev.off()
  }
  
  # pdf(file="Figures/LKINLARho.pdf", width=5, height=5)
  # plot(sigma2Marg, type="l", main=TeX("Marginal for $\\rho$"), xlab=TeX("$\\rho$"))
  # abline(v=rho, col="green")
  # dev.off()
  
  
  # # do the same for kappa, rho
  # # in order to get distribution for rho, must sample from joint hyperparameters
  # kappaMarg = inla.tmarginal(function(x) {sqrt(8)/exp(x) * latticeWidth}, mod$marginals.hyperpar$`Theta1 for field`)
  # # thetasToRho = function(xs) {
  # #   logCor = xs[2]
  # #   logVar = xs[3]
  # #   kappa = sqrt(8)/exp(logCor) * latticeWidth
  # #   sigma2 = exp(logVar)
  # #   sigma2 * 4*pi * kappa^2
  # # }
  # # samples = inla.hyperpar.sample(50000, mod, TRUE)
  # # rhos = apply(samples, 1, thetasToRho)
  # 
  # pdf(file=paste0("Figures/LKINLAKappa", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(kappaMarg, type="l", xlab="kappa", main="Marginal for kappa")
  # abline(v=inla.qmarginal(c(.025, .975), kappaMarg), col="purple", lty=2)
  # dev.off()
  # 
  # pdf(file="Figures/LKINLARho.pdf", width=5, height=5)
  # hist(rhos, xlab="rho", main="Marginal for Rho", breaks=1000, freq=F, xlim=c(0, quantile(probs=.95, rhos)))
  # abline(v=rho, col="green")
  # dev.off()
  
  ## Now generate marginals for the alpha parameters if they exist. In order to do this, we must generate draws from 
  ## the posterior, and transform them back to the probability scale
  out = inla.hyperpar.sample(20000, mod, improve.marginals=TRUE)
  if(nLayer >= 2) {
    if(!separateRanges) {
      zSamples = out[,3:(2+nLayer-1)]
      xSamples = apply(zSamples, 1, multivariateExpit)
      xSamples = rbind(xSamples, 1-colSums(xSamples))
    } else {
      zSamples = matrix(out[,(nLayer+1+1):(nLayer + 1 + 1 + nLayer-2)], ncol=nLayer-1)
      xSamples = matrix(apply(zSamples, 1, multivariateExpit), nrow=nLayer-1)
      xSamples = rbind(xSamples, 1-colSums(xSamples))
    }
    
    for(l in 1:nLayer) {
      thisQuantileRange = abs(diff(quantile(probs=c(.025, .975), xSamples[l,])))
      if(thisQuantileRange <= 0.2) {
        plotRange = range(xSamples[l,])
      } else {
        plotRange = c(0, 1)
      }
      pdf(file=paste0("Figures/LKINLAAlpha", l, plotNameRoot, ".pdf"), width=5, height=5)
      hist(xSamples[l,], xlab=TeX(paste0("$\\alpha_", l, "$")), main=TeX(paste0("Marginal for $\\alpha_", l, "$")), breaks=100, freq=F, xlim=plotRange)
      abline(v=mean(xSamples[l,]), col="purple", lty=1)
      abline(v=quantile(probs=c(.025, .975), xSamples[l,]), col="purple", lty=2)
      dev.off()
    }
  }
  
  ## plot covariance and correlation functions
  
  # first to transform all the hyperparameter samples to their relevant values
  # !separateRanges:
  # 1: log effective range
  # 2: log spatial variance
  # 3:(3 + nLayer - 2): multivariateLogit alpha
  # 3 + nLayer - 1: error precision
  # separateRanges:
  # 1:nLayer: log effective range
  # nLayer + 1: log spatial variance
  # (nLayer + 2):(nLayer + 2 + nLayer - 2): multivariateLogit alpha
  # nLayer + 1 + nLayer - 1 + 1: error precision
  if(clusterEffect)
    nuggetVarVals = 1/out[,ncol(out)]
  else
    nuggetVarVals = rep(0, nrow(out))
  
  if(separateRanges) {
    kappaVals = t(sweep(sqrt(8)/exp(out[,1:nLayer]), 2, latticeWidth, "*"))
    rhoVals = exp(out[,nLayer+1])
  } else {
    kappaVals = sqrt(8)/exp(out[,1]) * latticeWidth
    rhoVals = exp(out[,2])
  }
  alphaMat = xSamples
  
  # compute the covariance function for many different hyperparameter samples
  out = covarianceDistributionLKINLA(latInfo, kappaVals, rhoVals, nuggetVarVals, alphaMat)
  d = out$d
  sortI = sort(d, index.return=TRUE)$ix
  d = d[sortI]
  covMean = out$cov[sortI]
  upperCov=out$upperCov[sortI]
  lowerCov=out$lowerCov[sortI]
  covMat=out$covMat[sortI]
  corMean = out$cor[sortI]
  upperCor=out$upperCor[sortI]
  lowerCor=out$lowerCor[sortI]
  corMat=out$corMat[sortI]
  
  # plot the covariance function
  
  yRange = range(c(covMean, lowerCov, upperCov))
  pdf(file=paste0("Figures/LKINLACov", plotNameRoot, ".pdf"), width=5, height=5)
  plot(d, covMean, type="l", main="Posterior of covariance function", xlab="Distance", ylab="Covariance", 
       ylim=yRange)
  lines(d, lowerCov, lty=2)
  lines(d, upperCov, lty=2)
  legend("topright", c("Estimate", "80% CI"), lty=c(1, 2), col=c("black", "black"))
  dev.off()
  
  pdf(file=paste0("Figures/LKINLACor", plotNameRoot, ".pdf"), width=5, height=5)
  plot(d, corMean, type="l", main="Posterior of correlation function", xlab="Distance", ylab="Covariance", 
       ylim=c(0,1))
  lines(d, lowerCor, lty=2)
  lines(d, upperCor, lty=2)
  legend("topright", c("Estimate", "80% CI"), lty=c(1, 2), col=c("black", "black"))
  dev.off()
  
  # get scoring rules
  truth = ys / ns
  est = obsPreds
  vars = obsSDs^2
  # lower = fit$obsLower
  # upper = fit$obsUpper
  lower = NULL
  upper = NULL # these will be recalculated including binomial variation in the getScores function
  estMat = fit$obsMat
  estMatBinomial = addBinomialVar(estMat, ns)
  browser()
  
  if(!dropUrban) {
    print("Pooled scores:")
    print(data.frame(c(getScores(truth, est, vars, lower, upper, estMatBinomial, doRandomReject=TRUE), Time=time[3])))
    print("Rural scores:")
    print(data.frame(c(getScores(truth[!obsUrban], est[!obsUrban], vars[!obsUrban], lower[!obsUrban], upper[!obsUrban], estMatBinomial[!obsUrban,], doRandomReject=TRUE), Time=time[3])))
    print("Urban scores:")
    print(data.frame(c(getScores(truth[obsUrban], est[obsUrban], vars[obsUrban], lower[obsUrban], upper[obsUrban], estMatBinomial[obsUrban,], doRandomReject=TRUE), Time=time[3])))
  } else {
    print("Rural scores:")
    print(data.frame(c(getScores(truth, est, vars, lower, upper, estMatBinomial, doRandomReject=TRUE), Time=time[3])))
  }
}

# tests the LKINLA model applied to Kenya secondary education prevalence data
# buffer: buffer distance between domain edge of basis lattice and domain edge of data.
# n: number of observations
# xRange: range of x coordinates
# yRange: range of y coordinates
# nx: number of basis function lattice points in x directions
# NOTE: ny is determined automatically to match scale of x lattice points
# Xmat: design matrix
# ys: observations
# first.time: is first time evaluating function.  User should always set to FALSE
testSPDEKenyaDat = function(seed=1, urbanEffect=TRUE, clusterEffect=TRUE, 
                            plotNameRoot="", nPostSamples=1000, significanceCI=.8, 
                            mesh=getSPDEMeshKenya(), prior=getSPDEPrior(mesh), 
                            intStrategy="ccd", strategy="gaussian", 
                            targetPop=c("women", "children"), dropUrban=FALSE, 
                            clusterLevel=TRUE, pixelLevel=TRUE, countyLevel=TRUE, regionLevel=TRUE, 
                            family=c("binomial", "betabinomial"), verbose=TRUE) {
  
  set.seed(seed)
  targetPop = match.arg(targetPop)
  family = match.arg(family)
  
  # make sure input arguments are compatible
  if(dropUrban) {
    if(urbanEffect)
      warning("dropUrban and urbanEffect both set to TRUE. Setting urbanEffect to FALSE...")
    urbanEffect=FALSE
  }
  if(family == "betabinomial") {
    if(clusterEffect)
      stop("cluster effect must not be set to TRUE for betaBinomial model")
  }
  # if(separateRanges) {
  # if(clusterEffect)
  #   warning("separateRanges and clusterEffect both set to TRUE. Setting clusterEffect to FALSE...")
  # clusterEffect=FALSE
  # }
  
  # set default mesh
  
  
  # set plotNameRoot
  if(family == "betabinomial")
    familyText = "_BBin"
  else
    familyText = "_Bin"
  if(targetPop == "women") {
    dataType = "ed"
    plotNameRoot = paste0(plotNameRoot, "_Ed", familyText, "_urb", urbanEffect, "_clust", clusterEffect, "_dropUrb", dropUrban)
  }
  else {
    dataType = "mort"
    plotNameRoot = paste0(plotNameRoot, "_Mort", familyText, "_urb", urbanEffect, "_clust", clusterEffect, "_dropUrb", dropUrban)
  }
  
  # load data set
  if(dataType == "mort") {
    out = load("../U5MR/kenyaData.RData")
    dat = mort
  }
  else {
    out = load("../U5MR/kenyaDataEd.RData")
    dat = ed
  }
  if(dropUrban) {
    dat = dat[!dat$urban,]
  }
  coords = cbind(dat$east, dat$north)
  ys = dat$y
  ns = dat$n
  obsUrban = dat$urban
  
  # plot the observations
  pdf(file=paste0("Figures/SPDEObservations", plotNameRoot, ".pdf"), width=5, height=5)
  par(mfrow=c(1,1))
  quilt.plot(coords, ys/ns, FUN=function(x) {mean(x, na.rm=TRUE)})
  dev.off()
  
  # generate hyperparameters based on median and quantiles of inverse exponential and inverse gamma
  # priorPar = getPrior(.1, .1, 10)
  # generate hyperparameters for pc priors
  # median effective range is .4 or 200 for kenya data (a fifth of the spatial domain diameter), median spatial variance is 1
  # priorPar = getPCPrior(200, .01, 1, nLayer=nLayer, separateRanges=FALSE, latticeInfo=latInfo)
  # # show priors on effective correlation, marginal variance, and error variance:
  # xs1 = seq(1, 500, l=500)
  # if(!separateRanges) {
  #   pdf(file=paste0("Figures/SPDEPriorEffRange", plotNameRoot, ".pdf"), width=5, height=5)
  #   plot(xs1, dinvexp(xs1, rate=priorPar$corScalePar), type="l", col="blue", 
  #        xlab="Effective Correlation Range", main="Effective Correlation Prior", 
  #        ylab="Prior Density")
  #   abline(v=qinvexp(.5, rate=priorPar$corScalePar), col="red")
  #   dev.off()
  # } else {
  #   for(i in 1:nLayer) {
  #     if(i == nLayer && identical(NC, c(30, 107)))
  #       xs1 = seq(1, 200, l=500)
  #     pdf(file=paste0("Figures/SPDEPriorEff", i, "Range", plotNameRoot, ".pdf"), width=5, height=5)
  #     plot(xs1, dinvexp(xs1, rate=priorPar$corScalePar[i]), type="l", col="blue", 
  #          xlab="Effective Correlation Range", main="Effective Correlation Prior", 
  #          ylab="Prior Density")
  #     abline(v=qinvexp(.5, rate=priorPar$corScalePar[i]), col="red")
  #     dev.off()
  #   }
  # }
  # 
  # if(priorPar$priorType == "orig") {
  #   xs2 = seq(.01, 10.5, l=500)
  #   pdf(file=paste0("Figures/SPDEPriorMargVar", plotNameRoot, ".pdf"), width=5, height=5)
  #   plot(xs2, invgamma::dinvgamma(xs2, shape=priorPar$varPar1, rate=priorPar$varPar2), type="l", col="blue", 
  #        xlab="Marginal Variance", main="Marginal Variance Prior", 
  #        ylab="Prior Density")
  #   abline(v=qinvgamma(.1, shape=priorPar$varPar1, rate=priorPar$varPar2), col="red")
  #   abline(v=qinvgamma(.9, shape=priorPar$varPar1, rate=priorPar$varPar2), col="red")
  #   dev.off()
  # } else if(priorPar$priorType == "pc") {
  #   xs2 = seq(.01, 11.5, l=500)
  #   pdf(file=paste0("Figures/SPDEPriorMargVar", plotNameRoot, ".pdf"), width=5, height=5)
  #   plot(xs2, dpcvar(xs2, alpha=priorPar$alpha, u=priorPar$u), type="l", col="blue", 
  #        xlab="Marginal Variance", main="Marginal Variance Prior", 
  #        ylab="Prior Density")
  #   abline(v=qpcvar(.1, alpha=priorPar$alpha, u=priorPar$u), col="red")
  #   abline(v=qpcvar(.9, alpha=priorPar$alpha, u=priorPar$u), col="red")
  #   abline(v=1, col="green")
  #   dev.off()
  # }
  
  # xs2 = seq(.001, invgamma::qinvgamma(.905, shape=0.1, rate=0.1), l=500)
  # pdf(file="Figures/mixtureSPDEPriorErrorVar.pdf", width=5, height=5)
  # plot(xs2, invgamma::dinvgamma(xs2, shape=0.1, rate=0.1), type="l", col="blue", 
  #      xlab="Error Variance", main="Error Variance Prior", 
  #      ylab="Prior Density")
  # abline(v=invgamma::qinvgamma(.1, shape=0.1, rate=0.1), col="red")
  # abline(v=invgamma::qinvgamma(.9, shape=0.1, rate=0.1), col="red")
  # dev.off()
  
  # prior of cluster variance or overdispersion
  if(clusterEffect) {
    xs2 = seq(.01, 1, l=500)
    pdf(file=paste0("Figures/SPDEPriorErrorVar", plotNameRoot, ".pdf"), width=5, height=5)
    plot(xs2, dpcvar(xs2, alpha=.01, u=1), type="l", col="blue", 
         xlab="Marginal Variance", main="Marginal Variance Prior", 
         ylab="Prior Density")
    abline(v=qpcvar(.1, alpha=.01, u=1), col="red")
    abline(v=qpcvar(.9, alpha=.01, u=1), col="red")
    dev.off()
  } else if(family == "betabinomial") {
    # lambda = getLambdapcBeta(U=1, logitU=TRUE, alpha=0.01, p=.5, normalize=TRUE)
    # set median at .04 and upper 97.5th pctile at 0.2
    mu = logit(0.04)
    prec = 1/((logit(.2)-logit(.04))/qnorm(.975))^2
    rhos = seq(0, 1, l=1000)
    # tempYs = dpcBeta2(rhos, lambda=lambda, normalize=TRUE)
    tempYs = dlogitNormal(rhos, mu=mu, prec=prec)
    pdf(file=paste0("Figures/SPDEPriorRho", plotNameRoot, ".pdf"), width=5, height=5)
    plot(rhos, tempYs, type="l", xlab=TeX(paste0("$\\rho$")), ylab="Density", 
         main=TeX(paste0("Overdispersion Prior")), xlim=c(0,1), ylim=c(0, max(tempYs[is.finite(tempYs)])))
    # abline(v=qpcBeta2(c(0.025, 0.975), lambda=lambda, normalize=TRUE), col="purple", lty=2)
    abline(v=expit(qnorm(c(0.025, 0.975), mu, sqrt(1/prec))), col="purple", lty=2)
    dev.off()
  }
  
  # # prior on covariogram
  # alphaVals = t(rdirichlet(100, alpha=priorPar$alphaPar))
  # rhoVals = rpcvar(100, alpha=priorPar$alpha, u=priorPar$u)
  # effectiveRangeVals = t(matrix(rinvexp(100*length(priorPar$corScalePar), rate=priorPar$corScalePar), nrow=length(priorPar$corScalePar)))
  # if(separateRanges) {
  #   latticeWidth = sapply(latInfo, function(x) {x$latWidth})
  #   kappaVals = t(sweep(sqrt(8)/effectiveRangeVals, 2, latticeWidth, "*"))
  # } else {
  #   latticeWidth = latInfo[[1]]$latWidth
  #   kappaVals = c(sqrt(8)/exp(effectiveRangeVals) * latticeWidth)
  # }
  # if(clusterEffect)
  #   nuggetVarVals = rpcvar(100, alpha=.05, u=1)
  # else
  #   nuggetVarVals = rep(0, 100)
  # out = covarianceDistributionSPDE(latInfo, kappaVals, rhoVals, nuggetVarVals, alphaVals, 
  #                                    normalize=normalize, fastNormalize=fastNormalize)
  # d = out$d
  # sortI = sort(d, index.return=TRUE)$ix
  # d = d[sortI]
  # covMean = out$cov[sortI]
  # upperCov=out$upperCov[sortI]
  # lowerCov=out$lowerCov[sortI]
  # covMat=out$covMat[sortI]
  # corMean = out$cor[sortI]
  # upperCor=out$upperCor[sortI]
  # lowerCor=out$lowerCor[sortI]
  # corMat=out$corMat[sortI]
  # 
  # # correlation and covariance functions
  # 
  # # plot the covariance an correlation priors
  # yRange = range(c(covMean, lowerCov, upperCov))
  # pdf(file=paste0("Figures/SPDEPriorCov", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(d, covMean, type="l", main="Prior of covariance function", xlab="Distance", ylab="Covariance", 
  #      ylim=yRange)
  # lines(d, lowerCov, lty=2)
  # lines(d, upperCov, lty=2)
  # legend("topright", c("Estimate", "80% CI"), lty=c(1, 2), col=c("black", "black"))
  # dev.off()
  # 
  # pdf(file=paste0("Figures/SPDEPriorCor", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(d, corMean, type="l", main="Prior of correlation function", xlab="Distance", ylab="Covariance", 
  #      ylim=c(0,1))
  # lines(d, lowerCor, lty=2)
  # lines(d, upperCor, lty=2)
  # legend("topright", c("Estimate", "80% CI"), lty=c(1, 2), col=c("black", "black"))
  # dev.off()
  
  # fit the model
  time = system.time(fit <- fitSPDEKenyaDat(dat, dataType, mesh, prior, significanceCI, intStrategy, strategy, nPostSamples, 
                                            verbose, seed=seed, urbanEffect=urbanEffect, clusterEffect=clusterEffect, 
                                            family=family))
  mod = fit$mod
  preds=fit$preds
  predSDs=fit$sigmas
  obsPreds=fit$obsPreds
  obsSDs=fit$obsSDs
  coefPreds = fit$coefPreds
  coefSDs = fit$coefSDs
  predPts = fit$predPts
  predsUrban = fit$predsUrban
  
  # print out the total time
  print(paste0("Total time before aggregation: ", time[3]))
  
  # aggregate the model
  aggregationTime = 
    system.time(aggregatedFit <- 
                  aggregateModelResultsKenya(fit, clusterLevel=clusterLevel, pixelLevel=pixelLevel, 
                                             countyLevel=countyLevel, regionLevel=regionLevel, 
                                             targetPop=targetPop)
    )[3]
  allResults = list(fit=fit, aggregatedResults=aggregatedFit)
  
  print(paste0("Time to aggregate: ", aggregationTime))
  
  # show a model summary
  print(summary(mod))
  
  # function for determining if points are in correct range
  inRange = function(pts, rangeShrink=0) {
    rep(TRUE, nrow(pts))
  }
  
  # xRangeDat = latInfo[[1]]$xRangeDat
  # yRangeDat = latInfo[[1]]$yRangeDat
  xRangeDat = range(coords[,1])
  yRangeDat = range(coords[,2])
  
  # show predictive surface, SD, and data
  out = load("../U5MR/adminMapData.RData")
  # kenyaMap = adm0
  kenyaMap = adm1
  pdf(file=paste0("Figures/SPDEPreds", plotNameRoot, ".pdf"), width=8, height=8)
  par(mfrow=c(2,2), mar=c(5.1, 4.1, 4.1, 6))
  
  # obsInds = 1:n
  # predInds = (n+1):(n+mx*my)
  # coefInds = (n+mx*my+1):(n+mx*my+nx*ny)
  # colRangeDat = range(c(ys, obsPreds, preds, coefPreds))
  colRangeDat = range(c(ys, obsPreds, preds))
  colRangeSD = range(c(predSDs, obsSDs))
  quilt.plot(coords, ys, main="True Process", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
  quilt.plot(predPts[,1], predPts[,2], preds, main="Prediction Mean", zlim=colRangeDat, 
             xlim=xRangeDat, ylim=yRangeDat)
  
  quilt.plot(coords, ys, main="Observations", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
  quilt.plot(predPts[,1], predPts[,2], predSDs, main="Prediction SD",
             xlim=xRangeDat, ylim=yRangeDat, zlim=range(predSDs[inRange(predPts, rangeShrink=.03)]))
  dev.off()
  
  plotSingleModelPredictions(dat=dat, allResults, modelName="", targetPop=targetPop, 
                             areaLevels=c("Region", "County", "Pixel", "Cluster"), 
                             plotNameRoot=plotNameRoot, 
                             meanRange=NULL, meanTicks=NULL, meanTickLabels=NULL, widthRange=NULL, widthTicks=NULL, widthTickLabels=NULL, 
                             meanCols=makeRedBlueDivergingColors(64), 
                             widthCols=makeBlueYellowSequentialColors(64), 
                             kenyaLatRange=c(-4.6, 5), kenyaLonRange=c(33.5, 42.0))
  
  pdf(file=paste0("Figures/SPDEInSampleResidualsLabeled", plotNameRoot, ".pdf"), width=5, height=5)
  ylim = range(ys/ns-obsPreds)
  xlim = range(obsPreds)
  plot(obsPreds, ys/ns-obsPreds, pch=19, cex=.1, main="Residuals versus fitted", 
       ylab="Residuals", xlab="Fitted", xlim=xlim, ylim=ylim, type="n")
  points(obsPreds[!obsUrban], ys[!obsUrban]/ns[!obsUrban]-obsPreds[!obsUrban], pch=19, cex=.1, col="green")
  points(obsPreds[obsUrban], ys[obsUrban]/ns[obsUrban]-obsPreds[obsUrban], pch=19, cex=.1, col="blue")
  abline(h=0, lty=2)
  legend("topright", c("Rural", "Urban"), col=c("green", "blue"), pch=19)
  dev.off()
  
  # plot marginals on interpretable scale (effective range, marginal variance)
  effRangeMarg = mod$marginals.hyperpar$`Range for field`
  varMarg = inla.tmarginal(function(x) {x^2}, mod$marginals.hyperpar$`Stdev for field`)
  if(clusterEffect)
    sigma2Marg = inla.tmarginal(function(x) {1/x}, mod$marginals.hyperpar[[length(mod$marginals.hyperpar)]])
  else if(family == "betabinomial")
    overdispersionMarg = mod$marginals.hyperpar$`overdispersion for the betabinomial observations`
  covNames = names(mod$marginals.fixed)
  XMarginals = list()
  for(i in 1:length(covNames)) {
    XMarginal = inla.tmarginal(function(x) {x}, mod$marginals.fixed[[covNames[i]]])
    XMarginals = c(XMarginals, list(XMarginal))
  }
  
  par(mfrow=c(1,1))
  
  pdf(file=paste0("Figures/SPDEEffRange", plotNameRoot, ".pdf"), width=5, height=5)
  plot(effRangeMarg, type="l", main="Marginal for effective range")
  abline(v=inla.qmarginal(c(.025, .975), effRangeMarg), col="purple", lty=2)
  dev.off()
  
  # plot(mod$marginals.hyperpar$`Theta1 for field`, type="l", main="Marginal for log range")
  pdf(file=paste0("Figures/SPDEVar", plotNameRoot, ".pdf"), width=5, height=5)
  plot(varMarg, type="l", main="Marginal for spatial variance")
  abline(v=inla.qmarginal(c(.025, .975), varMarg), col="purple", lty=2)
  dev.off()
  # plot(mod$marginals.hyperpar$`Theta2 for field`, type="l", main="Marginal for log variance")
  if(clusterEffect) {
    pdf(file=paste0("Figures/SPDESigma2", plotNameRoot, ".pdf"), width=5, height=5)
    plot(sigma2Marg, type="l", main="Marginal for error variance")
    abline(v=inla.qmarginal(c(.025, .975), sigma2Marg), col="purple", lty=2)
    dev.off()
  } else if(family == "betabinomial") {
    pdf(file=paste0("Figures/SPDERho", plotNameRoot, ".pdf"), width=5, height=5)
    plot(overdispersionMarg, type="l", main="Marginal for overdispersion")
    abline(v=inla.qmarginal(c(.025, .975), overdispersionMarg), col="purple", lty=2)
    dev.off()
  }
  
  for(i in 1:length(covNames)) {
    XMarginal = XMarginals[[i]]
    pdf(file=paste0("Figures/SPDE", covNames[i], plotNameRoot, ".pdf"), width=5, height=5)
    plot(XMarginal, type="l", main=paste0("Marginal for fixed effect ", i))
    abline(v=inla.qmarginal(c(.025, .975), XMarginal), col="purple", lty=2)
    dev.off()
  }
  
  # pdf(file="Figures/SPDERho.pdf", width=5, height=5)
  # plot(sigma2Marg, type="l", main=TeX("Marginal for $\\rho$"), xlab=TeX("$\\rho$"))
  # abline(v=rho, col="green")
  # dev.off()
  
  
  # # do the same for kappa, rho
  # # in order to get distribution for rho, must sample from joint hyperparameters
  # kappaMarg = inla.tmarginal(function(x) {sqrt(8)/exp(x) * latticeWidth}, mod$marginals.hyperpar$`Theta1 for field`)
  # # thetasToRho = function(xs) {
  # #   logCor = xs[2]
  # #   logVar = xs[3]
  # #   kappa = sqrt(8)/exp(logCor) * latticeWidth
  # #   sigma2 = exp(logVar)
  # #   sigma2 * 4*pi * kappa^2
  # # }
  # # samples = inla.hyperpar.sample(50000, mod, TRUE)
  # # rhos = apply(samples, 1, thetasToRho)
  # 
  # pdf(file=paste0("Figures/SPDEKappa", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(kappaMarg, type="l", xlab="kappa", main="Marginal for kappa")
  # abline(v=inla.qmarginal(c(.025, .975), kappaMarg), col="purple", lty=2)
  # dev.off()
  # 
  # pdf(file="Figures/SPDERho.pdf", width=5, height=5)
  # hist(rhos, xlab="rho", main="Marginal for Rho", breaks=1000, freq=F, xlim=c(0, quantile(probs=.95, rhos)))
  # abline(v=rho, col="green")
  # dev.off()
  
  ## Now generate marginals for the alpha parameters if they exist. In order to do this, we must generate draws from 
  ## the posterior, and transform them back to the probability scale
  out = inla.hyperpar.sample(20000, mod, improve.marginals=TRUE)
  
  ## plot covariance and correlation functions
  
  # first to transform all the hyperparameter samples to their relevant values
  # !separateRanges:
  # 1: log effective range
  # 2: log spatial variance
  # 3:(3 + nLayer - 2): multivariateLogit alpha
  # 3 + nLayer - 1: error precision
  # separateRanges:
  # 1:nLayer: log effective range
  # nLayer + 1: log spatial variance
  # (nLayer + 2):(nLayer + 2 + nLayer - 2): multivariateLogit alpha
  # nLayer + 1 + nLayer - 1 + 1: error precision
  if(clusterEffect)
    nuggetVarVals = 1/out[,ncol(out)]
  else
    nuggetVarVals = rep(0, nrow(out))
  effectiveRangeVals = out[,2]
  varVals = out[,3]^2
  
  # compute the covariance function for many different hyperparameter samples
  # NOTE: this is on a logit scale for binomial and beta binomial case, and does not 
  #       account for overdispersion in beta binomial case
  out = covarianceDistributionSPDE(effectiveRangeVals, varVals, nuggetVarVals, mesh, xRangeDat=xRangeDat, yRangeDat=yRangeDat)
  
  d = out$d
  sortI = sort(d, index.return=TRUE)$ix
  d = d[sortI]
  covMean = out$cov[sortI]
  upperCov=out$upperCov[sortI]
  lowerCov=out$lowerCov[sortI]
  covMat=out$covMat[sortI]
  corMean = out$cor[sortI]
  upperCor=out$upperCor[sortI]
  lowerCor=out$lowerCor[sortI]
  corMat=out$corMat[sortI]
  
  # plot the covariance function
  
  yRange = range(c(covMean, lowerCov, upperCov))
  pdf(file=paste0("Figures/SPDECov", plotNameRoot, ".pdf"), width=5, height=5)
  plot(d, covMean, type="l", main="Posterior of covariance function", xlab="Distance", ylab="Covariance", 
       ylim=yRange)
  lines(d, lowerCov, lty=2)
  lines(d, upperCov, lty=2)
  legend("topright", c("Estimate", "80% CI"), lty=c(1, 2), col=c("black", "black"))
  dev.off()
  
  pdf(file=paste0("Figures/SPDECor", plotNameRoot, ".pdf"), width=5, height=5)
  plot(d, corMean, type="l", main="Posterior of correlation function", xlab="Distance", ylab="Covariance", 
       ylim=c(0,1))
  lines(d, lowerCor, lty=2)
  lines(d, upperCor, lty=2)
  legend("topright", c("Estimate", "80% CI"), lty=c(1, 2), col=c("black", "black"))
  dev.off()
  
  # get scoring rules
  truth = ys / ns
  est = obsPreds
  vars = obsSDs^2
  # lower = fit$obsLower
  # upper = fit$obsUpper
  lower = NULL
  upper = NULL # these will be recalculated including binomial variation in the getScores function
  estMat = fit$obsMat
  estMatBinomial = addBinomialVar(estMat, ns)
  
  if(!dropUrban) {
    print("Pooled scores:")
    print(data.frame(c(getScores(truth, est, vars, lower, upper, estMatBinomial, doRandomReject=TRUE), Time=time[3])))
    print("Rural scores:")
    print(data.frame(c(getScores(truth[!obsUrban], est[!obsUrban], vars[!obsUrban], lower[!obsUrban], upper[!obsUrban], estMatBinomial[!obsUrban,], doRandomReject=TRUE), Time=time[3])))
    print("Urban scores:")
    print(data.frame(c(getScores(truth[obsUrban], est[obsUrban], vars[obsUrban], lower[obsUrban], upper[obsUrban], estMatBinomial[obsUrban,], doRandomReject=TRUE), Time=time[3])))
  } else {
    print("Rural scores:")
    print(data.frame(c(getScores(truth, est, vars, lower, upper, estMatBinomial, doRandomReject=TRUE), Time=time[3])))
  }
}

# modelResults must have matrices "predMat" and "hyperMat"
# B: 
aggregateModelPreds = function(modelResults, B, addClusterError=TRUE) {
  
}

##### 
##### 
##### 
##### 
##### 
# testing functions for the covariance of mixture data set

# tests the fitLKINLAStandard function using data simulated from the LK model
# buffer: buffer distance between domain edge of basis lattice and domain edge of data.
# n: number of observations
# xRange: range of x coordinates
# yRange: range of y coordinates
# nx: number of basis function lattice points in x directions
# NOTE: ny is determined automatically to match scale of x lattice points
# Xmat: design matrix
# ys: observations
# first.time: is first time evaluating function.  User should always set to FALSE
# thetas: originally was c(.1, 3), but 3 is too large
testLKINLAModelMixture = function(seed=1, nLayer=3, nx=20, ny=nx, assumeMeanZero=TRUE, nu=1, 
                                  nBuffer=5, normalize=TRUE, fastNormalize=TRUE, NC=14, testCovs=FALSE, 
                                  printVerboseTimings=FALSE, latInfo=NULL, n=900, thetas=NULL, 
                                  testfrac=1/9, plotNameRoot="", sigma2=.1^2, useKenya=FALSE, 
                                  effRangeRange=NULL, urbanOverSamplefrac=0, nHyperSamples=1000, 
                                  intStrategy="ccd", strategy="gaussian", separateRanges=FALSE, 
                                  leaveOutRegion=TRUE, gscratch=FALSE, 
                                  savePrecomputationResults=FALSE, loadPrecomputationResults=FALSE, 
                                  precomputationFileNameRoot="precomputationResults") {
  set.seed(seed)
  clusterEffect=TRUE
  
  startTime = proc.time()[3]
  
  if(useKenya)
    distanceBreaks = seq(0, 300, l=20)
  else
    distanceBreaks = seq(0, 0.5, l=20)
  
  # set plotNameRoot
  if(length(NC) == 1) {
    if(separateRanges)
      NC = c(14, 126) # by default, use two layers with the finest layer having resolution equal to 10km
  }
  if(separateRanges)
    nLayer = length(NC)
  
  # set plotNameRoot
  # > 2/.08 * 5
  # [1] 125
  # > 2/.8 * 5
  # [1] 12.5
  ncText = ""
  if(length(NC) == 1) {
    if(separateRanges)
      ncText = "_NC14_126"
    else {
      ncText = paste0("_NC", NC)
    }
  } else {
    tempText = do.call("paste0", as.list(c(NC[1], paste0("_", NC[-1]))))
    ncText = paste0("_NC", tempText)
  }
  plotNameRoot = paste0(plotNameRoot, "_L", nLayer, ncText, "_sepRange", separateRanges, "_n", n, "_nu", nu, "_nugV", 
                        round(sigma2, 2), "_Kenya", useKenya, "_noInt", assumeMeanZero, "_urbOversamp", round(urbanOverSamplefrac, 4))
  
  # set true parameter values
  if(useKenya) {
    if(is.null(thetas))
      thetas=c(.08, .8) * (1000/2) / sqrt(8)
  } else {
    if(is.null(thetas))
      thetas=c(.08, .8) / sqrt(8)
  }
  rho = 1
  effectiveRange = thetas * sqrt(8)
  
  # load data set if necessary
  if(is.null(n)) {
    out = load("mixtureDataSet.RData")
  } else {
    spatialCorFun = function(x) {0.5 * stationary.cov(x, theta=thetas[1], Covariance="Matern", smoothness=nu) + 
        0.5 * stationary.cov(x, theta=thetas[2], Covariance="Matern", smoothness=nu)}
    mixtureCorFun = function(x) {0.5 * stationary.cov(x, theta=thetas[1], Covariance="Matern", smoothness=nu) + 
        0.5 * stationary.cov(x, theta=thetas[2], Covariance="Matern", smoothness=nu)}
    nTest = round(testfrac * n)
    if(leaveOutRegion) {
      simulationData = getSimulationDataSetsGivenCovarianceTest(mixtureCorFun, nTotal=n, nTest=nTest, marginalVar=rho, errorVar=sigma2, 
                                                                nDataSets=2, plotNameRoot=paste0("(0.5*Matern(", thetas[1], ") + 0.5*Matern(", thetas[2], "))"), fileNameRoot="mix", 
                                                                saveDataSetPlot=FALSE, doPredGrid=TRUE)
    } else {
      simulationData = getSimulationDataSetsGivenCovariance(mixtureCorFun, nTotal=n, nTest=nTest, marginalVar=rho, errorVar=sigma2, 
                                                            nDataSets=2, plotNameRoot=paste0("(0.5*Matern(", thetas[1], ") + 0.5*Matern(", thetas[2], "))"), fileNameRoot="mix", 
                                                            saveDataSetPlot=FALSE, useKenyaLocations=useKenya, urbanOverSamplefrac=urbanOverSamplefrac)
    }
  }
  coords = cbind(simulationData$xTrain[,1], simulationData$yTrain[,1])
  ys = simulationData$zTrain[,1]
  
  # generate lattice and simulate observations
  # coords = matrix(runif(2*n), ncol=2)
  if(useKenya) {
    xRangeDat = simulationData$xRange
    yRangeDat = simulationData$yRange
    # if(is.null(effRangeRange))
    #   effRangeRange=exp(c(-6, 7))
  } else {
    xRangeDat = c(-1, 1)
    yRangeDat = c(-1, 1)
  }
  if(is.null(latInfo))
    latInfo = makeLatGrids(xRangeDat, yRangeDat, NC, nBuffer, nLayer)
  
  AObs = makeA(coords, latInfo)
  # Q = makeQ(kappa=kappa, rho=rho, latInfo, alphas=alphas, normalized=normalize, fastNormalize=fastNormalize) 
  # L = as.matrix(t(chol(solve(Q))))
  # zsims = matrix(rnorm(nrow(Q)), ncol=1)
  # fieldSims = L %*% zsims
  # ys = as.numeric(AObs %*% fieldSims) + 1 # add a constant unit mean term to be estimated by INLA
  # # ys = 1 + as.numeric(AObs %*% fieldSims) + coords[,1] # x-valued mean term to be estimated by INLA
  # errs = rnorm(n, sd=sqrt(sigma2))
  # ys = ys + errs
  
  # plot the observations
  pdf(file=paste0("Figures/mixtureLKINLAObservations", plotNameRoot, ".pdf"), width=5, height=5)
  par(mfrow=c(1,1))
  quilt.plot(coords, ys)
  dev.off()
  
  # make prediction coordinates on a grid, and add testing points
  mx = 100
  my = 100
  predPts = make.surface.grid(list(x=seq(xRangeDat[1], xRangeDat[2], l=mx), y=seq(yRangeDat[1], yRangeDat[2], l=my)))
  if(useKenya) {
    # remove grid points outside of Kenya national boundaries
    load("../U5MR/adminMapData.RData")
    polys = adm0@polygons
    kenyaPoly = polys[[1]]@Polygons[[77]]@coords
    kenyaPolyProj = projKenya(kenyaPoly)
    inKenya = in.poly(predPts, kenyaPolyProj)
    predPts = predPts[inKenya,]
    
    # add other testing locations to matrix of prediction locations and remember which 
    predPts = rbind(predPts, cbind(simulationData$xTest[,1], simulationData$yTest[,1]))
    predPts = rbind(predPts, cbind(simulationData$xTestRural[,1], simulationData$yTestRural[,1]))
    predPts = rbind(predPts, cbind(simulationData$xTestUrban[,1], simulationData$yTestUrban[,1]))
    plotGridI = 1:sum(inKenya)
    overallTestI = simulationData$overallTestI
    ruralTestI = simulationData$ruralTestI
    urbanTestI = simulationData$urbanTestI
    gridTestI = (max(urbanTestI) + 1):(max(urbanTestI) + length(simulationData$xGrid))
  } else {
    predPts = rbind(predPts, cbind(simulationData$xTest[,1], simulationData$yTest[,1]))
  }
  predPts = rbind(predPts, cbind(simulationData$xGrid, simulationData$yGrid))
  ysTest = c(simulationData$zTest[,1], simulationData$zTestRural[,1], simulationData$zTestUrban[,1], simulationData$zGrid[,1])
  
  # generate hyperparameters based on median and quantiles of inverse exponential and inverse gamma
  # priorPar = getPrior(.1, .1, 10)
  # generate hyperparameters for pc priors
  # median effective range is .4 or 200 for kenya data (a fifth of the spatial domain diameter), median spatial variance is 1
  if(!useKenya)
    priorPar = getPCPrior(.4, .5, 1, nLayer=nLayer, separateRanges=separateRanges, latticeInfo=latInfo)
  else
    priorPar = getPCPrior(200, .5, 1, nLayer=nLayer, separateRanges=separateRanges, latticeInfo=latInfo)
  X = matrix(rep(1, nrow(coords)), ncol=1)
  # X = matrix(coords[,1], ncol=1)
  XPred = matrix(rep(1, nrow(predPts)), ncol=1)
  
  # add linear terms in lat/lon to covariate matrices if requested
  if(testCovs) {
    X = cbind(X, coords)
    XPred = cbind(XPred, predPts)
  }
  
  if(assumeMeanZero) {
    X = NULL
    XPred = NULL
  }
  
  # show priors on effective correlation, marginal variance, and error variance:
  if(!useKenya)
    xs1 = seq(.01, 2, l=500)
  else
    xs1 = seq(.01, 1000, l=500)
  if(!separateRanges) {
    pdf(file=paste0("Figures/mixtureLKINLAPriorEffRange", plotNameRoot, ".pdf"), width=5, height=5)
    plot(xs1, dinvexp(xs1, rate=priorPar$corScalePar), type="l", col="blue", 
         xlab="Effective Correlation Range", main="Effective Correlation Prior", 
         ylab="Prior Density")
    abline(v=qinvexp(.5, rate=priorPar$corScalePar), col="red")
    dev.off()
  } else {
    for(i in 1:nLayer) {
      if(!useKenya)
        xs1 = seq(.01, 2, l=500)
      else
        xs1 = seq(.01, 1000, l=500)
      
      if(i == nLayer && useKenya)
        xs1 = seq(1, 200, l=500)
      else if(!useKenya)
        xs1 = seq(.001, .4, l=500)
      pdf(file=paste0("Figures/mixtureLKINLAPriorEff", i, "Range", plotNameRoot, ".pdf"), width=5, height=5)
      plot(xs1, dinvexp(xs1, rate=priorPar$corScalePar[i]), type="l", col="blue", 
           xlab="Effective Correlation Range", main="Effective Correlation Prior", 
           ylab="Prior Density")
      abline(v=qinvexp(.5, rate=priorPar$corScalePar[i]), col="red")
      dev.off()
    }
  }
  
  if(priorPar$priorType == "orig") {
    xs2 = seq(.01, 10.5, l=500)
    pdf(file=paste0("Figures/mixtureLKINLAPriorMargVar", plotNameRoot, ".pdf"), width=5, height=5)
    plot(xs2, invgamma::dinvgamma(xs2, shape=priorPar$varPar1, rate=priorPar$varPar2), type="l", col="blue", 
         xlab="Marginal Variance", main="Marginal Variance Prior", 
         ylab="Prior Density")
    abline(v=qinvgamma(.1, shape=priorPar$varPar1, rate=priorPar$varPar2), col="red")
    abline(v=qinvgamma(.9, shape=priorPar$varPar1, rate=priorPar$varPar2), col="red")
    dev.off()
  } else if(priorPar$priorType == "pc") {
    xs2 = seq(.01, 11.5, l=500)
    pdf(file=paste0("Figures/mixtureLKINLAPriorMargVar", plotNameRoot, ".pdf"), width=5, height=5)
    plot(xs2, dpcvar(xs2, alpha=priorPar$alpha, u=priorPar$u), type="l", col="blue", 
         xlab="Marginal Variance", main="Marginal Variance Prior", 
         ylab="Prior Density")
    abline(v=qpcvar(.1, alpha=priorPar$alpha, u=priorPar$u), col="red")
    abline(v=qpcvar(.9, alpha=priorPar$alpha, u=priorPar$u), col="red")
    abline(v=1, col="green")
    dev.off()
  }
  
  # xs2 = seq(.001, invgamma::qinvgamma(.905, shape=0.1, rate=0.1), l=500)
  # pdf(file="Figures/mixtureLKINLAPriorErrorVar.pdf", width=5, height=5)
  # plot(xs2, invgamma::dinvgamma(xs2, shape=0.1, rate=0.1), type="l", col="blue", 
  #      xlab="Error Variance", main="Error Variance Prior", 
  #      ylab="Prior Density")
  # abline(v=invgamma::qinvgamma(.1, shape=0.1, rate=0.1), col="red")
  # abline(v=invgamma::qinvgamma(.9, shape=0.1, rate=0.1), col="red")
  # dev.off()
  
  xs2 = seq(.01, 1, l=500)
  pdf(file=paste0("Figures/mixtureLKINLAPriorErrorVar", plotNameRoot, ".pdf"), width=5, height=5)
  plot(xs2, dpcvar(xs2, alpha=.05, u=1), type="l", col="blue", 
       xlab="Error Variance", main="Error Variance Prior", 
       ylab="Prior Density")
  abline(v=qpcvar(.1, alpha=.05, u=1), col="red")
  abline(v=qpcvar(.9, alpha=.05, u=1), col="red")
  abline(v=sqrt(.1), col="green")
  dev.off()
  
  # browser()
  
  for(l in 1:nLayer) {
    pdf(file=paste0("Figures/mixtureLKINLAPriorAlpha", l, plotNameRoot, ".pdf"), width=5, height=5)
    xs = seq(0, 1, l=500)
    tempYs = dbeta(xs, priorPar$alphaPar[l], sum(priorPar$alphaPar[-l]))
    plot(xs, tempYs, type="l", xlab=TeX(paste0("$\\alpha_", l, "$")), ylab="Density", 
         main=TeX(paste0("Marginal for $\\alpha_", l, "$")), xlim=c(0,1), ylim=c(0, max(tempYs[is.finite(tempYs)])))
    abline(v=qbeta(c(0.025, 0.975), priorPar$alphaPar[l], sum(priorPar$alphaPar[-l])), col="purple", lty=2)
    dev.off()
  }
  
  # prior on covariogram
  alphaVals = t(rdirichlet(100, alpha=priorPar$alphaPar))
  rhoVals = rpcvar(100, alpha=priorPar$alpha, u=priorPar$u)
  if(!separateRanges) {
    kappaVals = sqrt(8)/rinvexp(100, rate=priorPar$corScalePar) * latInfo[[1]]$latWidth
  } else {
    kappaVals = matrix(sqrt(8)/rinvexp(100*3, rate=priorPar$corScalePar) * sapply(latInfo, function(x) {x$latWidth}), nrow=nLayer)
  }
  nuggetVarVals = rpcvar(100, alpha=.05, u=1)
  if(loadPrecomputationResults)
    loadFilename = precomputationFileNameRoot
  else
    loadFilename = ""
  out = covarianceDistributionLKINLA(latInfo, kappaVals, rhoVals, nuggetVarVals, alphaVals, 
                                     normalize=normalize, fastNormalize=fastNormalize, precomputationsFileNameRoot=loadFilename)
  d = out$d
  sortI = sort(d, index.return=TRUE)$ix
  d = d[sortI]
  covMean = out$cov[sortI]
  upperCov=out$upperCov[sortI]
  lowerCov=out$lowerCov[sortI]
  covMat=out$covMat[sortI]
  corMean = out$cor[sortI]
  upperCor=out$upperCor[sortI]
  lowerCor=out$lowerCor[sortI]
  corMat=out$corMat[sortI]
  
  # true correlation and covariance functions
  spatialCorFun = function(x) {0.5 * stationary.cov(x, theta=thetas[1], Covariance="Matern", distMat=x, smoothness=nu) + 
      0.5 * stationary.cov(x, theta=thetas[2], Covariance="Matern", distMat=x, smoothness=nu)}
  spatialCovFun = spatialCorFun
  mixtureCovFun = function(x) {
    out = spatialCorFun(x)
    out[x == 0] = 1 + sigma2
    out
  }
  mixtureCorFun = function(x) { mixtureCovFun(x) * (1 / (1 + sigma2)) }
  
  # plot the covariance an correlation priors
  yRange = range(c(covMean, lowerCov, upperCov, mixtureCovFun(d)))
  pdf(file=paste0("Figures/mixtureLKINLAPriorCov", plotNameRoot, ".pdf"), width=5, height=5)
  plot(d, covMean, type="l", main="Posterior of covariance function", xlab="Distance", ylab="Covariance", 
       ylim=yRange)
  lines(d, lowerCov, lty=2)
  lines(d, upperCov, lty=2)
  lines(d, mixtureCovFun(d), col="green")
  legend("topright", c("Truth", "Estimate", "80% CI"), lty=c(1, 1, 2), col=c("green", "black", "black"))
  dev.off()
  
  pdf(file=paste0("Figures/mixtureLKINLAPriorCor", plotNameRoot, ".pdf"), width=5, height=5)
  yRange = range(c(corMean, lowerCor, upperCor, mixtureCorFun(d)))
  plot(d, corMean, type="l", main="Posterior of correlation function", xlab="Distance", ylab="Covariance", 
       ylim=c(0,1))
  lines(d, lowerCor, lty=2)
  lines(d, upperCor, lty=2)
  lines(d, mixtureCorFun(d), col="green")
  legend("topright", c("Truth", "Estimate", "80% CI"), lty=c(1, 1, 2), col=c("green", "black", "black"))
  dev.off()
  
  # fit the model
  time = system.time(fit <- fitLKINLAStandard2(coords, ys, predCoords=predPts, seed=seed, nLayer=nLayer, NC=NC,
                                               nBuffer=nBuffer, priorPar=priorPar, xObs=X, xPred=XPred, normalize=normalize, 
                                               intStrategy=intStrategy, strategy=strategy, fastNormalize=fastNormalize, 
                                               printVerboseTimings=printVerboseTimings, latInfo=latInfo, effRangeRange=effRangeRange, 
                                               separateRanges=separateRanges, clusterEffect=clusterEffect, 
                                               savePrecomputationResults=savePrecomputationResults, 
                                               loadPrecomputationResults=loadPrecomputationResults, 
                                               precomputationFileNameRoot=precomputationFileNameRoot))
  mod = fit$mod
  preds=fit$preds
  predSDs=fit$sigmas
  latInfo=fit$latInfo
  latWidth=fit$latWidth
  obsPreds=fit$obsPreds
  obsSDs=fit$obsSDs
  coefPreds = fit$coefPreds
  coefSDs = fit$coefSDs
  
  # print out the total time
  print(paste0("Total time: ", time[3]))
  
  # show a model summary
  print(summary(mod))
  
  # function for determining if points are in correct range
  if(!useKenya) {
    inRange = function(pts, rangeShrink=0) {
      inX = (rangeShrink < pts[,1]) & (pts[,1] < 1-rangeShrink)
      inY = (rangeShrink < pts[,2]) & (pts[,2] < 1-rangeShrink)
      inX & inY
    }
  } else {
    inRange = function(pts, rangeShrink=0) {
      rep(TRUE, nrow(pts))
    }
  }
  
  # show predictive surface, SD, and data
  
  if(nLayer==1) {
    pdf(file=paste0("Figures/mixtureLKINLAPreds", plotNameRoot, ".pdf"), width=9, height=6)
    par(mfrow=c(2,3))
    
    # obsInds = 1:n
    # predInds = (n+1):(n+mx*my)
    # coefInds = (n+mx*my+1):(n+mx*my+nx*ny)
    colRangeDat = range(c(ys, obsPreds, preds, coefPreds))
    colRangeSD = range(c(range(predSDs[inRange(predPts)]), coefSDs[[1]][inRange(gridPtsL1)], 
                         coefSDs[[2]][inRange(gridPtsL2)], coefSDs[[3]][inRange(gridPtsL3)]))
    gridPtsL1 = latInfo[[1]]$latCoords
    quilt.plot(coords, ys, main="True Process", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL1[,1], gridPtsL1[,2], coefPreds[[1]], main="Basis Coefficient Mean (Layer 1)", xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL1[,1], gridPtsL1[,2], coefSDs[[1]], main="Basis Coefficient SD (Layer 1)", zlim=colRangeSD, xlim=xRangeDat, ylim=yRangeDat)
    
    # quilt.plot(coords, obsPreds, main="Prediction mean", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat, 
    #            zlim=range(predSDs[inRange(predPts)]))
    plot.new()
    quilt.plot(predPts[,1], predPts[,2], preds, main="Prediction Mean", zlim=colRangeDat, 
               xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(predPts[,1], predPts[,2], predSDs, main="Prediction SD", 
               xlim=xRangeDat, ylim=yRangeDat)
  }
  else if(nLayer==2) {
    pdf(file=paste0("Figures/mixtureLKINLAPreds", plotNameRoot, ".pdf"), width=12, height=6)
    par(mfrow=c(2,4))
    
    # obsInds = 1:n
    # predInds = (n+1):(n+mx*my)
    # coefInds = (n+mx*my+1):(n+mx*my+nx*ny)
    colRangeDat = range(c(ys, obsPreds, preds, coefPreds))
    colRangeCoef = range(c(coefPreds))
    colRangeSD = range(c(predSDs, obsSDs, coefSDs))
    gridPtsL1 = latInfo[[1]]$latCoords
    gridPtsL2 = latInfo[[2]]$latCoords
    quilt.plot(coords, ys, main="True Process", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(predPts[,1], predPts[,2], preds, main="Prediction Mean", zlim=colRangeDat, 
               xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL1[,1], gridPtsL1[,2], coefPreds[[1]], main="Basis Coefficient Mean (Layer 1)", xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL2[,1], gridPtsL2[,2], coefPreds[[2]], main="Basis Coefficient Mean (Layer 2)", xlim=xRangeDat, ylim=yRangeDat)
    
    quilt.plot(coords, obsPreds, main="Observation Mean", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(predPts[,1], predPts[,2], predSDs, main="Prediction SD", 
               xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL1[,1], gridPtsL1[,2], coefSDs[[1]], main="Basis Coefficient SD (Layer 1)", zlim=colRangeSD, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL2[,1], gridPtsL2[,2], coefSDs[[2]], main="Basis Coefficient SD (Layer 2)", zlim=colRangeSD, xlim=xRangeDat, ylim=yRangeDat)
  }
  else if(nLayer==3) {
    pdf(file=paste0("Figures/mixtureLKINLAPreds", plotNameRoot, ".pdf"), width=15, height=6)
    par(mfrow=c(2,5), mar=c(5.1, 4.1, 4.1, 6))
    
    # obsInds = 1:n
    # predInds = (n+1):(n+mx*my)
    # coefInds = (n+mx*my+1):(n+mx*my+nx*ny)
    # colRangeDat = range(c(ys, obsPreds, preds, coefPreds))
    colRangeDat = range(c(ys, obsPreds, preds, coefPreds))
    colRangeCoef = range(c(coefPreds))
    colRangeSD = range(c(predSDs, obsSDs, coefSDs))
    gridPtsL1 = latInfo[[1]]$latCoords
    gridPtsL2 = latInfo[[2]]$latCoords
    gridPtsL3 = latInfo[[3]]$latCoords
    quilt.plot(coords, ys, main="True Process", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(predPts[,1], predPts[,2], preds, main="Prediction Mean", zlim=colRangeDat, 
               xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL1[,1], gridPtsL1[,2], coefPreds[[1]], main="Basis Coefficient Mean (Layer 1)", 
               xlim=xRangeDat, ylim=yRangeDat, zlim=colRangeCoef)
    quilt.plot(gridPtsL2[,1], gridPtsL2[,2], coefPreds[[2]], main="Basis Coefficient Mean (Layer 2)", 
               xlim=xRangeDat, ylim=yRangeDat, zlim=colRangeCoef)
    quilt.plot(gridPtsL3[,1], gridPtsL3[,2], coefPreds[[3]], main="Basis Coefficient Mean (Layer 3)", 
               xlim=xRangeDat, ylim=yRangeDat, zlim=colRangeCoef)
    
    # quilt.plot(coords, obsPreds, main="Observation Mean", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    # plot.new()
    quilt.plot(coords, ys, main="Observations", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(predPts[,1], predPts[,2], predSDs, main="Prediction SD",
               xlim=xRangeDat, ylim=yRangeDat, zlim=range(predSDs[inRange(predPts, rangeShrink=.03)]))
    quilt.plot(gridPtsL1[,1], gridPtsL1[,2], coefSDs[[1]], main="Basis Coefficient SD (Layer 1)", zlim=colRangeSD, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL2[,1], gridPtsL2[,2], coefSDs[[2]], main="Basis Coefficient SD (Layer 2)", zlim=colRangeSD, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL3[,1], gridPtsL3[,2], coefSDs[[3]], main="Basis Coefficient SD (Layer 3)", zlim=colRangeSD, xlim=xRangeDat, ylim=yRangeDat)
  }
  else if(nLayer==4) {
    pdf(file=paste0("Figures/mixtureLKINLAPreds", plotNameRoot, ".pdf"), width=18, height=6)
    par(mfrow=c(2,6), mar=c(5.1, 4.1, 4.1, 6))
    
    # obsInds = 1:n
    # predInds = (n+1):(n+mx*my)
    # coefInds = (n+mx*my+1):(n+mx*my+nx*ny)
    # colRangeDat = range(c(ys, obsPreds, preds, coefPreds))
    colRangeDat = range(c(ys, obsPreds, preds, coefPreds))
    colRangeCoef = range(c(coefPreds))
    colRangeSD = range(c(predSDs, obsSDs, coefSDs))
    gridPtsL1 = latInfo[[1]]$latCoords
    gridPtsL2 = latInfo[[2]]$latCoords
    gridPtsL3 = latInfo[[3]]$latCoords
    gridPtsL4 = latInfo[[4]]$latCoords
    quilt.plot(coords, ys, main="True Process", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(predPts[,1], predPts[,2], preds, main="Prediction Mean", zlim=colRangeDat, 
               xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL1[,1], gridPtsL1[,2], coefPreds[[1]], main="Basis Coefficient Mean (Layer 1)", 
               xlim=xRangeDat, ylim=yRangeDat, zlim=colRangeCoef)
    quilt.plot(gridPtsL2[,1], gridPtsL2[,2], coefPreds[[2]], main="Basis Coefficient Mean (Layer 2)", 
               xlim=xRangeDat, ylim=yRangeDat, zlim=colRangeCoef)
    quilt.plot(gridPtsL3[,1], gridPtsL3[,2], coefPreds[[3]], main="Basis Coefficient Mean (Layer 3)", 
               xlim=xRangeDat, ylim=yRangeDat, zlim=colRangeCoef)
    quilt.plot(gridPtsL4[,1], gridPtsL4[,2], coefPreds[[4]], main="Basis Coefficient Mean (Layer 4)", 
               xlim=xRangeDat, ylim=yRangeDat, zlim=colRangeCoef)
    
    # quilt.plot(coords, obsPreds, main="Observation Mean", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    # plot.new()
    quilt.plot(coords, ys, main="Observations", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(predPts[,1], predPts[,2], predSDs, main="Prediction SD",
               xlim=xRangeDat, ylim=yRangeDat, zlim=range(predSDs[inRange(predPts, rangeShrink=.03)]))
    quilt.plot(gridPtsL1[,1], gridPtsL1[,2], coefSDs[[1]], main="Basis Coefficient SD (Layer 1)", zlim=colRangeSD, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL2[,1], gridPtsL2[,2], coefSDs[[2]], main="Basis Coefficient SD (Layer 2)", zlim=colRangeSD, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL3[,1], gridPtsL3[,2], coefSDs[[3]], main="Basis Coefficient SD (Layer 3)", zlim=colRangeSD, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL4[,1], gridPtsL4[,2], coefSDs[[4]], main="Basis Coefficient SD (Layer 4)", zlim=colRangeSD, xlim=xRangeDat, ylim=yRangeDat)
  }
  else if(nLayer==5) {
    pdf(file=paste0("Figures/mixtureLKINLAPreds", plotNameRoot, ".pdf"), width=21, height=6)
    par(mfrow=c(2,7), mar=c(5.1, 4.1, 4.1, 6))
    
    # obsInds = 1:n
    # predInds = (n+1):(n+mx*my)
    # coefInds = (n+mx*my+1):(n+mx*my+nx*ny)
    # colRangeDat = range(c(ys, obsPreds, preds, coefPreds))
    colRangeDat = range(c(ys, obsPreds, preds, coefPreds))
    colRangeCoef = range(c(coefPreds))
    colRangeSD = range(c(predSDs, obsSDs, coefSDs))
    gridPtsL1 = latInfo[[1]]$latCoords
    gridPtsL2 = latInfo[[2]]$latCoords
    gridPtsL3 = latInfo[[3]]$latCoords
    gridPtsL4 = latInfo[[4]]$latCoords
    gridPtsL5 = latInfo[[5]]$latCoords
    quilt.plot(coords, ys, main="True Process", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(predPts[,1], predPts[,2], preds, main="Prediction Mean", zlim=colRangeDat, 
               xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL1[,1], gridPtsL1[,2], coefPreds[[1]], main="Basis Coefficient Mean (Layer 1)", 
               xlim=xRangeDat, ylim=yRangeDat, zlim=colRangeCoef)
    quilt.plot(gridPtsL2[,1], gridPtsL2[,2], coefPreds[[2]], main="Basis Coefficient Mean (Layer 2)", 
               xlim=xRangeDat, ylim=yRangeDat, zlim=colRangeCoef)
    quilt.plot(gridPtsL3[,1], gridPtsL3[,2], coefPreds[[3]], main="Basis Coefficient Mean (Layer 3)", 
               xlim=xRangeDat, ylim=yRangeDat, zlim=colRangeCoef)
    quilt.plot(gridPtsL4[,1], gridPtsL4[,2], coefPreds[[4]], main="Basis Coefficient Mean (Layer 4)", 
               xlim=xRangeDat, ylim=yRangeDat, zlim=colRangeCoef)
    quilt.plot(gridPtsL5[,1], gridPtsL5[,2], coefPreds[[5]], main="Basis Coefficient Mean (Layer 5)", 
               xlim=xRangeDat, ylim=yRangeDat, zlim=colRangeCoef)
    
    # quilt.plot(coords, obsPreds, main="Observation Mean", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    # plot.new()
    quilt.plot(coords, ys, main="Observations", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(predPts[,1], predPts[,2], predSDs, main="Prediction SD",
               xlim=xRangeDat, ylim=yRangeDat, zlim=range(predSDs[inRange(predPts, rangeShrink=.03)]))
    quilt.plot(gridPtsL1[,1], gridPtsL1[,2], coefSDs[[1]], main="Basis Coefficient SD (Layer 1)", zlim=colRangeSD, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL2[,1], gridPtsL2[,2], coefSDs[[2]], main="Basis Coefficient SD (Layer 2)", zlim=colRangeSD, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL3[,1], gridPtsL3[,2], coefSDs[[3]], main="Basis Coefficient SD (Layer 3)", zlim=colRangeSD, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL4[,1], gridPtsL4[,2], coefSDs[[4]], main="Basis Coefficient SD (Layer 4)", zlim=colRangeSD, xlim=xRangeDat, ylim=yRangeDat)
    quilt.plot(gridPtsL5[,1], gridPtsL5[,2], coefSDs[[5]], main="Basis Coefficient SD (Layer 5)", zlim=colRangeSD, xlim=xRangeDat, ylim=yRangeDat)
  }
  dev.off()
  
  pdf(file=paste0("Figures/mixtureLKINLALeftOutResiduals", plotNameRoot, ".pdf"), width=5, height=5)
  testIndices = (length(preds) - length(ysTest) + 1):length(preds)
  plot(preds[testIndices], ysTest-preds[testIndices], pch=19, cex=.5, col="blue", main="Residuals versus fitted", 
       ylab="Residuals", xlab="Fitted")
  abline(h=0, lty=2)
  dev.off()
  
  if(useKenya) {
    pdf(file=paste0("Figures/mixtureLKINLALeftOutResidualsLabeled", plotNameRoot, ".pdf"), width=5, height=5)
    testIndices = (length(preds) - length(ysTest) + 1):length(preds)
    gridIndices = 
    ylim = range(ysTest-preds[testIndices])
    xlim = range(preds[testIndices])
    plot(preds[testIndices][overallTestI], ysTest[overallTestI]-preds[testIndices][overallTestI], pch=19, cex=.1, col="black", main="Residuals versus fitted", 
         ylab="Residuals", xlab="Fitted", xlim=xlim, ylim=ylim)
    points(preds[testIndices][ruralTestI], ysTest[ruralTestI]-preds[testIndices][ruralTestI], pch=19, cex=.1, col="green")
    points(preds[testIndices][urbanTestI], ysTest[urbanTestI]-preds[testIndices][urbanTestI], pch=19, cex=.1, col="blue")
    points(preds[testIndices][gridTestI], ysTest[gridTestI]-preds[testIndices][gridTestI], pch=19, cex=.1, col="red")
    abline(h=0, lty=2)
    legend("topright", c("Overall", "Rural", "Urban", "Grid"), col=c("black", "green", "blue", "red"), pch=19)
    dev.off()
  }
  
  # calculate true effective range and marginal variance:
  latticeWidth = latInfo[[1]]$latWidth
  if(separateRanges)
    latticeWidth = sapply(latInfo, function(x) {x$latWidth})#
  # marginalVar = rho/(4*pi * kappa^2)
  # marginalVar = getMultiMargVar(kappa, rho, nLayer=nLayer, nu=nu, xRange=xRangeBasis, 
  #                               yRange=yRangeBasis, nx=nx, ny=ny)[1]
  marginalVar = rho
  
  # # plot marginals on interpretable scale (effective range, marginal variance)
  # effRangeMarg = inla.tmarginal(exp, mod$marginals.hyperpar$`Theta1 for field`)
  # varMarg = inla.tmarginal(exp, mod$marginals.hyperpar$`Theta2 for field`)
  # sigma2Marg = inla.tmarginal(function(x) {1/x}, mod$marginals.hyperpar$`Precision for the Gaussian observations`)
  # covNames = names(mod$marginals.fixed)
  # if(!assumeMeanZero) {
  #   XMarginals = list()
  #   for(i in 1:length(covNames)) {
  #     XMarginal = inla.tmarginal(function(x) {x}, mod$marginals.fixed[[covNames[i]]])
  #     XMarginals = c(XMarginals, list(XMarginal))
  #   }
  # }
  # 
  # par(mfrow=c(1,1))
  # if(!separateRanges) {
  #   pdf(file=paste0("Figures/mixtureLKINLAEffRange", plotNameRoot, ".pdf"), width=5, height=5)
  #   plot(effRangeMarg, type="l", main="Marginal for effective range")
  #   abline(v=effectiveRange, col="green")
  #   abline(v=inla.qmarginal(c(.025, .975), effRangeMarg), col="purple", lty=2)
  #   dev.off()
  # } else {
  #   for(l in 1:nLayer) {
  #     pdf(file=paste0("Figures/mixtureLKINLAEffRange", "Layer", l, plotNameRoot, ".pdf"), width=5, height=5)
  #     plot(effRangeMarg, type="l", main=paste0("Marginal for effective range (Layer ", l, ")"))
  #     abline(v=effectiveRange, col="green")
  #     abline(v=inla.qmarginal(c(.025, .975), effRangeMarg), col="purple", lty=2)
  #     dev.off()
  #   }
  # }
  
  # plot marginals on interpretable scale (effective range, marginal variance)
  if(!separateRanges) {
    effRangeMarg = inla.tmarginal(exp, mod$marginals.hyperpar$`Theta1 for field`)
    varMarg = inla.tmarginal(exp, mod$marginals.hyperpar$`Theta2 for field`)
  } else {
    numberThetas = nLayer + 1 + nLayer - 1
    allNames = paste0("Theta", 1:numberThetas, " for field")
    effRangeMargs = list()
    for(i in 1:nLayer) {
      effRangeMargs = c(effRangeMargs, list(inla.tmarginal(exp, mod$marginals.hyperpar[[allNames[i]]])))
    }
    varMarg = inla.tmarginal(exp, mod$marginals.hyperpar[[allNames[nLayer+1]]])
  }
  sigma2Marg = inla.tmarginal(function(x) {1/x}, mod$marginals.hyperpar[[1]])
  covNames = names(mod$marginals.fixed)
  XMarginals = list()
  if(length(covNames) != 0) {
    for(i in 1:length(covNames)) {
      XMarginal = inla.tmarginal(function(x) {x}, mod$marginals.fixed[[covNames[i]]])
      XMarginals = c(XMarginals, list(XMarginal))
    }
  }
  
  par(mfrow=c(1,1))
  
  if(!separateRanges) {
    pdf(file=paste0("Figures/mixtureLKINLAEffRange", plotNameRoot, ".pdf"), width=5, height=5)
    plot(effRangeMarg, type="l", main="Marginal for effective range")
    abline(v=inla.qmarginal(c(.025, .975), effRangeMarg), col="purple", lty=2)
    dev.off()
  } else {
    for(i in 1:nLayer) {
      pdf(file=paste0("Figures/mixtureLKINLAEffRange", i, plotNameRoot, ".pdf"), width=5, height=5)
      plot(effRangeMargs[[i]], type="l", main="Marginal for effective range")
      abline(v=inla.qmarginal(c(.025, .975), effRangeMargs[[i]]), col="purple", lty=2)
      dev.off()
    }
  }
  
  # plot(mod$marginals.hyperpar$`Theta1 for field`, type="l", main="Marginal for log range")
  pdf(file=paste0("Figures/mixtureLKINLAVar", plotNameRoot, ".pdf"), width=5, height=5)
  plot(varMarg, type="l", main="Marginal for spatial variance")
  abline(v=marginalVar, col="green")
  abline(v=inla.qmarginal(c(.025, .975), varMarg), col="purple", lty=2)
  dev.off()
  # plot(mod$marginals.hyperpar$`Theta2 for field`, type="l", main="Marginal for log variance")
  pdf(file=paste0("Figures/mixtureLKINLASigma2", plotNameRoot, ".pdf"), width=5, height=5)
  plot(sigma2Marg, type="l", main="Marginal for error variance")
  abline(v=sigma2, col="green")
  abline(v=inla.qmarginal(c(.025, .975), sigma2Marg), col="purple", lty=2)
  dev.off()
  
  if(!assumeMeanZero) {
    for(i in 1:length(covNames)) {
      XMarginal = XMarginals[[i]]
      pdf(file=paste0("Figures/mixtureLKINLA", covNames[i], plotNameRoot, ".pdf"), width=5, height=5)
      plot(XMarginal, type="l", main="Marginal for fixed effect")
      abline(v=0, col="green")
      abline(v=inla.qmarginal(c(.025, .975), XMarginal), col="purple", lty=2)
      dev.off()
    }
  }
  
  # pdf(file="Figures/mixtureLKINLARho.pdf", width=5, height=5)
  # plot(sigma2Marg, type="l", main=TeX("Marginal for $\\rho$"), xlab=TeX("$\\rho$"))
  # abline(v=rho, col="green")
  # dev.off()
  
  
  # do the same for kappa, rho
  # in order to get distribution for rho, must sample from joint hyperparameters
  kappaMarg = inla.tmarginal(function(x) {sqrt(8)/exp(x) * latticeWidth}, mod$marginals.hyperpar$`Theta1 for field`)
  # thetasToRho = function(xs) {
  #   logCor = xs[2]
  #   logVar = xs[3]
  #   kappa = sqrt(8)/exp(logCor) * latticeWidth
  #   sigma2 = exp(logVar)
  #   sigma2 * 4*pi * kappa^2
  # }
  # samples = inla.hyperpar.sample(50000, mod, TRUE)
  # rhos = apply(samples, 1, thetasToRho)
  
  # pdf(file=paste0("Figures/mixtureLKINLAKappa", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(kappaMarg, type="l", xlab="kappa", main="Marginal for kappa")
  # abline(v=inla.qmarginal(c(.025, .975), kappaMarg), col="purple", lty=2)
  # dev.off()
  
  # pdf(file="Figures/mixtureLKINLARho.pdf", width=5, height=5)
  # hist(rhos, xlab="rho", main="Marginal for Rho", breaks=1000, freq=F, xlim=c(0, quantile(probs=.95, rhos)))
  # abline(v=rho, col="green")
  # dev.off()
  
  ## Now generate marginals for the alpha parameters. In order to do this, we must generate draws from 
  ## the posterior, and transform them back to the probability scale
  out = inla.hyperpar.sample(nHyperSamples, mod, improve.marginals=TRUE)
  if(separateRanges)
    alphaI = (1 + nLayer+1 + 1):(1 + nLayer+1 + nLayer-1)
  else
    alphaI = 4:(3+nLayer-1)
  zSamples = matrix(out[,alphaI], ncol=length(alphaI))
  xSamples = matrix(apply(zSamples, 1, multivariateExpit), nrow=length(alphaI))
  xSamples = rbind(xSamples, 1-colSums(xSamples))
  
  for(l in 1:nLayer) {
    pdf(file=paste0("Figures/mixtureLKINLAAlpha", l, plotNameRoot, ".pdf"), width=5, height=5)
    hist(xSamples[l,], xlab=TeX(paste0("$\\alpha_", l, "$")), main=TeX(paste0("Marginal for $\\alpha_", l, "$")), breaks=100, freq=F, xlim=c(0,1))
    abline(v=mean(xSamples[l,]), col="purple", lty=1)
    abline(v=quantile(probs=c(.025, .975), xSamples[l,]), col="purple", lty=2)
    dev.off()
  }
  
  ## plot covariance and correlation functions
  # first get the true covariance an correlation functions
  # spatialCovFun = function(x) {0.4 * Exp.cov(x, theta=0.1) + 0.6 * Exp.cov(x, theta=3)}
  # mixtureCovFun = function(x) {
  #   out = spatialCovFun(x)[1,]
  #   out[x == 0] = 1 + sqrt(.1)
  #   out
  # }
  # mixtureCorFun = function(x) { mixtureCovFun(x) * (1 / (1 + sqrt(.1))) }
  # spatialCorFun = function(x) {0.4 * Exp.cov(x, theta=thetas[1]) + 0.6 * Exp.cov(x, theta=thetas[2])}
  spatialCorFun = function(x) {0.5 * stationary.cov(x, theta=thetas[1], Covariance="Matern", distMat=x, smoothness=nu) + 
      0.5 * stationary.cov(x, theta=thetas[2], Covariance="Matern", distMat=x, smoothness=nu)}
  spatialCovFun = spatialCorFun
  mixtureCovFun = function(x) {
    out = spatialCorFun(x)
    out[x == 0] = 1 + sigma2
    out
  }
  mixtureCorFun = function(x) { mixtureCovFun(x) * (1 / (1 + sigma2)) }
  
  # first to transform all the hyperparameter samples to their relevant values
  # 1: error precision
  # 2: log effective range
  # 3: log spatial variance
  # 4-(3 + nLayer - 1): multivariateLogit alpha
  # nuggetVarVals = 1 / out[,1]
  nuggetVarVals = rep(0, ncol(xSamples))
  if(separateRanges) {
    kappaVals = t(sweep(sqrt(8)/exp(out[,2:(nLayer+1)]), 2, sapply(latInfo, function(x) {x$latWidth}), "*"))
    rhoVals = exp(out[,nLayer+2])
  } else {
    kappaVals = sqrt(8)/exp(out[,2]) * latticeWidth
    rhoVals = exp(out[,3])
  }
  alphaMat = xSamples
  
  # compute the covariance function for many different hyperparameter samples
  out = covarianceDistributionLKINLA(latInfo, kappaVals, rhoVals, nuggetVarVals, alphaMat, precomputationsFileNameRoot=loadFilename)
  d = out$d
  sortI = sort(d, index.return=TRUE)$ix
  d = d[sortI]
  covMean = out$cov[sortI]
  upperCov=out$upperCov[sortI]
  lowerCov=out$lowerCov[sortI]
  covMat=out$covMat[sortI,]
  corMean = out$cor[sortI]
  upperCor=out$upperCor[sortI]
  lowerCor=out$lowerCor[sortI]
  corMat=out$corMat[sortI,]
  covInfo = list(d=d, covMean=covMean, upperCov=upperCov, lowerCov=lowerCov, covMat=covMat, 
                 corMean=corMean, upperCor=upperCor, lowerCor=lowerCor, corMat=corMat)
  
  # plot the covariance function
  pdf(file=paste0("Figures/mixtureLKINLACov", plotNameRoot, ".pdf"), width=5, height=5)
  plot(d, covMean, type="l", main="Posterior of covariance function", xlab="Distance", ylab="Covariance")
  lines(d, lowerCov, lty=2)
  lines(d, upperCov, lty=2)
  lines(d, mixtureCovFun(d), col="green")
  dev.off()
  
  pdf(file=paste0("Figures/mixtureLKINLACor", plotNameRoot, ".pdf"), width=5, height=5)
  plot(d, corMean, type="l", main="Posterior of correlation function", xlab="Distance", ylab="Correlation")
  lines(d, lowerCor, lty=2)
  lines(d, upperCor, lty=2)
  lines(d, mixtureCorFun(d), col="green")
  dev.off()
  
  yRange = range(c(covMean, lowerCov, upperCov, mixtureCovFun(d)))
  pdf(file=paste0("Figures/mixtureLKINLACov", plotNameRoot, ".pdf"), width=5, height=5)
  plot(d, covMean, type="l", main="Posterior of covariance function", xlab="Distance", ylab="Covariance", 
       ylim=yRange)
  lines(d, lowerCov, lty=2)
  lines(d, upperCov, lty=2)
  lines(d, mixtureCovFun(d), col="green")
  legend("topright", c("Truth", "Estimate", "80% CI"), lty=c(1, 1, 2), col=c("green", "black", "black"))
  dev.off()
  
  pdf(file=paste0("Figures/mixtureLKINLACor", plotNameRoot, ".pdf"), width=5, height=5)
  yRange = range(c(corMean, lowerCor, upperCor, mixtureCorFun(d)))
  plot(d, corMean, type="l", main="Posterior of correlation function", xlab="Distance", ylab="Correlation", 
       ylim=c(0,1))
  lines(d, lowerCor, lty=2)
  lines(d, upperCor, lty=2)
  lines(d, mixtureCorFun(d), col="green")
  legend("topright", c("Truth", "Estimate", "80% CI"), lty=c(1, 1, 2), col=c("green", "black", "black"))
  dev.off()
  
  # get scoring rules
  testIndices = (length(preds) - length(ysTest) + 1):length(preds)
  leftOutIndices = (length(preds) - length(ysTest) + 1):(length(preds) - length(ysTest) + length(simulationData$zTest[,1]))
  gridIndices = (length(preds) - length(ysTest) + length(simulationData$zTest[,1]) + 1):length(preds)
  leftOutIndicesTest = match(leftOutIndices, testIndices)
  gridIndicesTest = match(gridIndices, testIndices)
  
  # first calculate scoring rules at grid points
  truth = ysTest[gridIndicesTest]
  est = preds[gridIndices]
  vars = predSDs[gridIndices]^2
  lower = fit$lower[gridIndices]
  upper = fit$upper[gridIndices]
  
  # compute nearest neighbor distances and scores as a function of them
  gridPts = predPts[gridIndices,]
  distMat = rdist(coords, gridPts)
  nndists = apply(distMat, 2, function(x) {min(x[x != 0])})
  print("Binned grid scores:")
  gridScoringRules = getScores(truth, est, vars, lower, upper, distances=nndists, breaks=distanceBreaks)
  print(gridScoringRules$pooledResults)
  print(gridScoringRules$binnedResults)
  
  # now calculate square rules at left out points
  truth = ysTest[leftOutIndicesTest]
  est = preds[leftOutIndices]
  vars = predSDs[leftOutIndices]^2
  lower = fit$lower[leftOutIndices]
  upper = fit$upper[leftOutIndices]
  leftOutScoringRules = getScores(truth, est, vars, lower, upper)
  leftOutScoringRules = data.frame(c(leftOutScoringRules, Time=time[3]))
  print("Binned left out scores:")
  print(leftOutScoringRules)
  
  scoringRules = list(gridScoringRules=gridScoringRules, leftOutScoringRules=leftOutScoringRules)
  # 
  # # get scoring rules
  # truth = ysTest
  # est = preds[testIndices]
  # vars = predSDs[testIndices]^2
  # lower = fit$lower[testIndices]
  # upper = fit$upper[testIndices]
  # 
  # # compute nearest neighbor distances and scores as a function of them
  # testPts = predPts[testIndices,]
  # distMat = rdist(coords, testPts)
  # nndists = apply(distMat, 2, function(x) {min(x[x != 0])})
  # print("Binned scores:")
  # scoringRules = getScores(truth, est, vars, lower, upper, distances=nndists, breaks=distanceBreaks)
  # scoringRules$pooledResults = data.frame(c(scoringRules$pooledResults, Time=time[3]))
  # print(scoringRules$binnedResults)
  
  
  
  # print("Grid scores:")
  # print(getScores(truth[gridTestI], est[gridTestI], vars[gridTestI], lower[gridTestI], upper[gridTestI], 
  #                 distances=nndists[gridTestI], breaks=distanceBreaks)$binnedResults)
  # 
  # # plot scores as a function of distance
  # distanceScores = getScores(truth[gridTestI], est[gridTestI], vars[gridTestI], lower[gridTestI], upper[gridTestI], 
  #                            distances=nndists[gridTestI], breaks=distanceBreaks)$binnedResults
  
  # pdf(file=paste0("Figures/mixtureSPDEScoreBias", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(distanceScores$NNDist, distanceScores$Bias, pch=19, col="blue", main="Bias", ylab="Bias", xlab="Nearest neighbor distance (km)")
  # abline(h=0, lty=2)
  # dev.off()
  # 
  # pdf(file=paste0("Figures/mixtureSPDEScoreVar", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(distanceScores$NNDist, distanceScores$Var, pch=19, col="blue", main="Variance", ylab="Variance", xlab="Nearest neighbor distance (km)", 
  #      ylim=c(0, max(distanceScores$Var)))
  # abline(h=0, lty=2)
  # dev.off()
  # 
  # pdf(file=paste0("Figures/mixtureSPDEScoreMSE", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(distanceScores$NNDist, distanceScores$MSE, pch=19, col="blue", main="MSE", ylab="MSE", xlab="Nearest neighbor distance (km)", 
  #      ylim=c(0, max(distanceScores$MSE)))
  # abline(h=0, lty=2)
  # dev.off()
  # 
  # pdf(file=paste0("Figures/mixtureSPDEScoreRMSE", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(distanceScores$NNDist, distanceScores$RMSE, pch=19, col="blue", main="RMSE", ylab="RMSE", xlab="Nearest neighbor distance (km)", 
  #      ylim=c(0, max(distanceScores$RMSE)))
  # abline(h=0, lty=2)
  # dev.off()
  # 
  # pdf(file=paste0("Figures/mixtureSPDEScoreCRPS", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(distanceScores$NNDist, distanceScores$CRPS, pch=19, col="blue", main="CRPS", ylab="CRPS", xlab="Nearest neighbor distance (km)", 
  #      ylim=c(0, max(distanceScores$CRPS)))
  # abline(h=0, lty=2)
  # dev.off()
  # 
  # pdf(file=paste0("Figures/mixtureSPDEScoreCvg", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(distanceScores$NNDist, distanceScores$Coverage, pch=19, col="blue", main="80% Coverage", ylab="80% Coverage", xlab="Nearest neighbor distance (km)", 
  #      ylim=c(0, max(distanceScores$Coverage)))
  # abline(h=0, lty=2)
  # dev.off()
  # 
  # pdf(file=paste0("Figures/mixtureSPDEScoreWidth", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(distanceScores$NNDist, distanceScores$Width, pch=19, col="blue", main="Width", ylab="Width", xlab="Nearest neighbor distance (km)", 
  #      ylim=c(0, max(distanceScores$Width)))
  # abline(h=0, lty=2)
  # dev.off()
  
  if(!useKenya) {
    # print("Pooled scores:")
    # print(data.frame(c(getScores(truth, est, vars, lower, upper), Time=time[3])))
  } else {
    print("Pooled scores:")
    print(data.frame(c(getScores(truth, est, vars, lower, upper), Time=time[3])))
    print("Overall scores:")
    print(data.frame(c(getScores(truth[overallTestI], est[overallTestI], vars[overallTestI], lower[overallTestI], upper[overallTestI]), Time=time[3])))
    print("Rural scores:")
    print(data.frame(c(getScores(truth[ruralTestI], est[ruralTestI], vars[ruralTestI], lower[ruralTestI], upper[ruralTestI]), Time=time[3])))
    print("Urban scores:")
    print(data.frame(c(getScores(truth[urbanTestI], est[urbanTestI], vars[urbanTestI], lower[urbanTestI], upper[urbanTestI]), Time=time[3])))
    # print("Grid scores:")
    # print(data.frame(c(getScores(truth[gridTestI], est[gridTestI], vars[gridTestI], lower[gridTestI], upper[gridTestI]), Time=time[3])))
  }
  
  # get aggregated predictions
  # A = t(sapply(1:(mx*my), getARow))
  # A = sweep(A, 1, rowSums(A), "/")
  # mx = 100
  # my = 100
  # predPts = make.surface.grid(list(x=seq(xRangeDat[1], xRangeDat[2], l=mx), y=seq(yRangeDat[1], yRangeDat[2], l=my)))
  gridCoords = cbind(simulationData$xGrid, simulationData$yGrid)
  # testIndices = (length(preds) - length(ysTest) + 1):length(preds)
  # gridIndices = (length(preds) - length(simulationData$xGrid) + 1):length(preds)
  # gridIndicesEst = (length(est) - length(simulationData$xGrid) + 1):length(est)
  
  A = makeNumericalIntegralMat(gridCoords, mx=3, my=3)
  aggregatedPreds = A %*% preds[gridIndices]
  
  truth = A %*% ysTest[gridIndicesTest]
  est = A %*% preds[gridIndices]
  aggregatedPredMat = A %*% fit$predMat[gridIndices,]
  vars = apply(aggregatedPredMat, 1, var)
  sds = apply(aggregatedPredMat, 1, sd)
  lower = apply(aggregatedPredMat, 1, function(x) {quantile(x, probs=.1)})
  upper = apply(aggregatedPredMat, 1, function(x) {quantile(x, probs=.9)})
  predictionMatrix = data.frame(Truth=truth, Est=est, SDs=sds, Lower=lower, Upper=upper)
  print("Aggregated prediction summary statistics:")
  print(predictionMatrix)
  
  print("Pooled aggregated scores:")
  pooledAggregatedScores = getScores(truth, est, vars, lower, upper)
  print(pooledAggregatedScores)
  print("Left out region aggregated scores:")
  leftOutAggregatedScores = getScores(truth[5], est[5], vars[5], lower[5], upper[5])
  leftOutAggregatedScores$Var = sds[5]
  names(leftOutAggregatedScores)[2] = "Predictive.SD"
  print(leftOutAggregatedScores)
  print("Included regions aggregated scores:")
  includedAggregatedScores = getScores(truth[-5], est[-5], vars[-5], lower[-5], upper[-5])
  print(includedAggregatedScores)
  
  aggregatedScoringRules = list(pooledAggregatedScores=pooledAggregatedScores, leftOutAggregatedScores=leftOutAggregatedScores, 
                                includedAggregatedScores=includedAggregatedScores)
  
  analysisTime = proc.time()[3] - startTime
  
  fit$mod = NULL
  if(!gscratch)
    save(scoringRules, fit, covInfo, predictionMatrix, aggregatedScoringRules, analysisTime, file=paste0("savedOutput/simulations/mixtureLKINLA", plotNameRoot, ".RData"))
  else
    save(scoringRules, fit, covInfo, predictionMatrix, aggregatedScoringRules, analysisTime, file=paste0("/work/johnpai/mixtureLKINLA", plotNameRoot, ".RData"))
  
  invisible(list(scoringRules, fit, covInfo, predictionMatrix, aggregatedScoringRules, analysisTime))
}

# tests the fitLKINLAStandard function using data simulated from the LK model
# buffer: buffer distance between domain edge of basis lattice and domain edge of data.
# n: number of observations
# xRange: range of x coordinates
# yRange: range of y coordinates
# nx: number of basis function lattice points in x directions
# NOTE: ny is determined automatically to match scale of x lattice points
# Xmat: design matrix
# ys: observations
# first.time: is first time evaluating function.  User should always set to FALSE
# thetas: originally was c(.1, 3), but 3 is too large
# testLKINLAModelMixtureMultiple = function(seed=1, nSamples=10, nLayer=3, nx=20, ny=nx, assumeMeanZero=TRUE, nu=1, 
#                                   nBuffer=5, normalize=TRUE, fastNormalize=TRUE, NC=14, testCovs=FALSE, 
#                                   printVerboseTimings=FALSE, latInfo=NULL, n=900, thetas=NULL, 
#                                   testfrac=.1, plotNameRoot="", sigma2=.1^2, useKenya=FALSE, 
#                                   effRangeRange=NULL, urbanOverSamplefrac=0, 
#                                   intStrategy="ccd", strategy="gaussian", separateRanges=FALSE, 
#                                   leaveOutRegion=TRUE) {
testLKINLAModelMixtureMultiple = function(seed=1, nSamples=100, startI=1, endI=nSamples, loadResults=FALSE, NC=14, nLayer=3, separateRanges=FALSE, n=900, nu=1, sigma2=0.1^2, 
                                          useKenya=FALSE, assumeMeanZero=TRUE, urbanOverSamplefrac=0, gscratch=TRUE, ...) {
  # set random seeds for each simulation
  set.seed(seed)
  allSeeds = sample(1:1000000, nSamples, replace = FALSE)
  
  # load in the results
  print("Loading in simulation results")
  # set plotNameRoot
  if(length(NC) == 1) {
    if(separateRanges)
      NC = c(14, 126) # by default, use two layers with the finest layer having resolution equal to 10km
  }
  if(separateRanges)
    nLayer = length(NC)
  
  # set plotNameRoot
  # > 2/.08 * 5
  # [1] 125
  # > 2/.8 * 5
  # [1] 12.5
  ncText = ""
  if(length(NC) == 1) {
    if(separateRanges)
      ncText = "_NC14_126"
    else {
      ncText = paste0("_NC", NC)
    }
  } else {
    tempText = do.call("paste0", as.list(c(NC[1], paste0("_", NC[-1]))))
    ncText = paste0("_NC", tempText)
  }
  plotNameRoot = paste0("_L", nLayer, ncText, "_sepRange", separateRanges, "_n", n, "_nu", nu, "_nugV", 
                        round(sigma2, 2), "_Kenya", useKenya, "_noInt", assumeMeanZero, "_urbOversamp", round(urbanOverSamplefrac, 4))
  
  
  # call testLKINLAModelMixture for each simulation requested
  precomputationFileNameRoot = paste0("precomputationMixLKINLA", plotNameRoot)
  temp = function(i) {
    print(paste0("Beginning simulation ", i, "/", nSamples))
    thisPlotNameRoot = paste0("sim", i)
    
    # make sure to save the precomputed results if this is the first run
    savePrecomputationResults = i == 1
    loadPrecomputationResults = i != 1
    do.call("testLKINLAModelMixture", c(list(seed = allSeeds[i], NC=NC, nLayer=nLayer, separateRanges=separateRanges, n=n, nu=nu, sigma2=sigma2, 
                                             useKenya=useKenya, urbanOverSamplefrac=urbanOverSamplefrac, assumeMeanZero=assumeMeanZero, 
                                             plotNameRoot=thisPlotNameRoot, gscratch=gscratch, 
                                             savePrecomputationResults=savePrecomputationResults, 
                                             loadPrecomputationResults=loadPrecomputationResults, 
                                             precomputationFileNameRoot=precomputationFileNameRoot), list(...)))
  }
  if(!loadResults)
    sapply(startI:endI, temp)
  
  if(startI != 1 || endI != nSamples)
    return(invisible(NULL))
  
  # save(scoringRules, fit, covInfo, predictionMatrix, aggregatedScoringRules, file=paste0("savedOutput/simulations/mixtureLKINLA", plotNameRoot, ".RData"))
  allScoringRulesGrid = list()
  allScoringRulesLeftOut = list()
  allFits = list()
  allCovInfo = list()
  allPredictionMatrices = list()
  allAggregatedScoringRules = list()
  allAnalysisTimes = c()
  for(i in 1:nSamples) {
    if(!gscratch)
      out = load(paste0("savedOutput/simulations/mixtureLKINLAsim", i, plotNameRoot, ".RData"))
    else
      out = load(paste0("/work/johnpai/mixtureLKINLAsim", i, plotNameRoot, ".RData"))
    allScoringRulesGrid = c(allScoringRulesGrid, list(scoringRules$gridScoringRules))
    allScoringRulesLeftOut = c(allScoringRulesLeftOut, list(scoringRules$leftOutScoringRules))
    allFits = c(allFits, list(fit))
    allCovInfo = c(allCovInfo, list(covInfo))
    allPredictionMatrices = c(allPredictionMatrices, list(predictionMatrix))
    allAggregatedScoringRules = c(allAggregatedScoringRules, list(aggregatedScoringRules))
    allAnalysisTimes = c(allAnalysisTimes, analysisTime)
  }
  
  ##### average results from each simulation
  # pointwise scoring rules
  allPooledScoringRulesGrid = do.call("rbind", lapply(allScoringRulesGrid, function(x) {x$pooledResults}))
  allBinnedScoringRulesGrid = lapply(allScoringRulesGrid, function(x) {x$binnedResults})
  binnedScoringRulesGrid = averageBinnedScores(allBinnedScoringRulesGrid)
  ns = binnedScoringRulesGrid[,2]
  pooledScoringRulesGrid = apply(binnedScoringRulesGrid, 2, function(x) {sum(x * (ns / sum(ns)))})
  pooledScoringRulesGrid = as.data.frame(matrix(pooledScoringRulesGrid, nrow=1))
  names(pooledScoringRulesGrid) = names(binnedScoringRulesGrid)
  
  fullPooledScoringRulesLeftOut = do.call("rbind", allScoringRulesLeftOut)
  pooledScoringRulesLeftOut = colMeans(fullPooledScoringRulesLeftOut)
  
  # covInfo
  # covInfo = list(d=d, covMean=covMean, upperCov=upperCov, lowerCov=lowerCov, covMat=covMat, 
  #                corMean=corMean, upperCor=upperCor, lowerCor=lowerCor, corMat=corMat)
  d = allCovInfo[[1]]$d
  covMean = rowMeans(do.call("cbind", lapply(allCovInfo, function(x) {x$covMean})))
  upperCov = rowMeans(do.call("cbind", lapply(allCovInfo, function(x) {x$upperCov})))
  lowerCov = rowMeans(do.call("cbind", lapply(allCovInfo, function(x) {x$lowerCov})))
  corMean = rowMeans(do.call("cbind", lapply(allCovInfo, function(x) {x$corMean})))
  upperCor = rowMeans(do.call("cbind", lapply(allCovInfo, function(x) {x$upperCor})))
  lowerCor = rowMeans(do.call("cbind", lapply(allCovInfo, function(x) {x$lowerCor})))
  
  # aggregated scoring rules
  # aggregatedScoringRules = list(pooledAggregatedScores=pooledAggregatedScores, leftOutAggregatedScores=leftOutAggregatedScores, 
  #                               includedAggregatedScores=includedAggregatedScores)
  # predictionMatrix = data.frame(Truth=truth, Est=est, SDs=sds, Lower=lower, Upper=upper)
  fullPredictionMatrix = do.call("rbind", allPredictionMatrices)
  leftOutIndices = seq(5, nrow(fullPredictionMatrix), by=9)
  leftOutPredictionMatrix = fullPredictionMatrix[leftOutIndices, ]
  leftInPredictionMatrix = fullPredictionMatrix[-leftOutIndices, ]
  leftOutScores = getScores(leftOutPredictionMatrix[,1], leftOutPredictionMatrix[,2], leftOutPredictionMatrix[,3]^2)
  leftInScores = getScores(leftInPredictionMatrix[,1], leftInPredictionMatrix[,2], leftInPredictionMatrix[,3]^2)
  aggregatedScores = getScores(fullPredictionMatrix[,1], fullPredictionMatrix[,2], fullPredictionMatrix[,3]^2)
  
  ##### Save results
  if(!gscratch) {
    save(#allScoringRules, 
         allFits, 
         allCovInfo, 
         allPredictionMatrices, 
         allAggregatedScoringRules, 
         allAnalysisTimes, 
         binnedScoringRulesGrid, 
         pooledScoringRulesGrid, 
         fullPooledScoringRulesLeftOut, 
         pooledScoringRulesLeftOut, 
         covMean, 
         upperCov, 
         lowerCov, 
         corMean, 
         upperCor, 
         lowerCor, 
         fullPredictionMatrix, 
         leftOutPredictionMatrix, 
         leftInPredictionMatrix, 
         leftOutScores, 
         leftInScores, 
         aggregatedScores, 
         file=paste0("savedOutput/simulations/mixtureLKINLAAll_nsim", nSamples, plotNameRoot, ".RData"))
  } else {
    save(#allScoringRules, 
         allFits, 
         allCovInfo, 
         allPredictionMatrices, 
         allAggregatedScoringRules, 
         binnedScoringRulesGrid, 
         pooledScoringRulesGrid, 
         fullPooledScoringRulesLeftOut, 
         pooledScoringRulesLeftOut, 
         covMean, 
         upperCov, 
         lowerCov, 
         corMean, 
         upperCor, 
         lowerCor, 
         fullPredictionMatrix, 
         leftOutPredictionMatrix, 
         leftInPredictionMatrix, 
         leftOutScores, 
         leftInScores, 
         aggregatedScores, 
         file=paste0("/work/johnpai/simulations/mixtureLKINLAAll_nsim", nSamples, plotNameRoot, ".RData"))
  }
}

# tests the fitLKStandard function using data simulated from the LK model
# buffer: buffer distance between domain edge of basis lattice and domain edge of data.
# n: number of observations
# xRange: range of x coordinates
# yRange: range of y coordinates
# nx: number of basis function lattice points in x directions
# NOTE: ny is determined automatically to match scale of x lattice points
# Xmat: design matrix
# ys: observations
# first.time: is first time evaluating function.  User should always set to FALSE
testLKModelMixture = function(seed=548676, nLayer=3, nx=20, ny=nx, nu=1, assumeMeanZero=TRUE, 
                              nBuffer=5, normalize=TRUE, NC=14, testCovs=TRUE, 
                              printVerboseTimings=FALSE, n=900, separatea.wght=FALSE, 
                              plotNameRoot="", doMatern=FALSE, fixNu=FALSE, thetas=c(.08, .8) / sqrt(8), 
                              testfrac=1/9, leaveOutRegion=TRUE, sigma2 = 0.1^2, extraPlotName=plotNameRoot, 
                              gscratch=FALSE) {
  set.seed(seed)
  
  startTime = proc.time()[3]
  
  # if(useKenya)
  #   distanceBreaks = seq(0, 300, l=20)
  # else
  distanceBreaks = seq(0, 0.5, l=20)
  
  # set true parameter values
  rho = 1
  effectiveRange = thetas * sqrt(8)
  
  # load data set if necessary
  spatialCorFun = function(x) {0.5 * stationary.cov(x, theta=thetas[1], Covariance="Matern", smoothness=nu) + 
      0.5 * stationary.cov(x, theta=thetas[2], Covariance="Matern", smoothness=nu)}
  if(is.null(n)) {
    out = load("mixtureDataSet.RData")
  } else {
    mixtureCorFun = function(x) {0.5 * stationary.cov(x, theta=thetas[1], Covariance="Matern", smoothness=nu) + 
        0.5 * stationary.cov(x, theta=thetas[2], Covariance="Matern", smoothness=nu)}
    nTest = round(testfrac * n)
    if(leaveOutRegion) {
      simulationData = getSimulationDataSetsGivenCovarianceTest(mixtureCorFun, nTotal=n, nTest=nTest, marginalVar=rho, errorVar=sigma2, 
                                                                nDataSets=2, plotNameRoot=paste0("(0.5*Matern(", thetas[1], ") + 0.5*Matern(", thetas[2], "))"), fileNameRoot="mix", 
                                                                saveDataSetPlot=FALSE, doPredGrid=TRUE)
    } else {
      simulationData = getSimulationDataSetsGivenCovariance(mixtureCorFun, nTotal=n, nTest=nTest, marginalVar=rho, errorVar=sigma2, 
                                                            nDataSets=2, plotNameRoot=paste0("(0.5*Matern(", thetas[1], ") + 0.5*Matern(", thetas[2], "))"), fileNameRoot="mix", 
                                                            saveDataSetPlot=FALSE, useKenyaLocations=useKenya, urbanOverSamplefrac=urbanOverSamplefrac)
    }
  }
  coords = cbind(simulationData$xTrain[,1], simulationData$yTrain[,1])
  ys = simulationData$zTrain[,1]
  
  # generate lattice and simulate observations
  # coords = matrix(runif(2*n), ncol=2)
  xRangeDat = c(-1, 1)
  yRangeDat = c(-1, 1)
  
  # plot the observations
  pdf(file=paste0("Figures/mixtureObservations", extraPlotName, ".pdf"), width=5, height=5)
  par(mfrow=c(1,1))
  quilt.plot(coords, ys)
  dev.off()
  
  # make prediction coordinates on a grid, and add testing points
  # urban and rural points are NULL, but they are left in the code in case kenya data used later
  xRange=c(-1,1)
  yRange=c(-1,1)
  mx = 100
  my = 100
  predPts = make.surface.grid(list(x=seq(xRange[1], xRange[2], l=mx), y=seq(yRange[1], yRange[2], l=my)))
  predPts = rbind(predPts, cbind(simulationData$xTest[,1], simulationData$yTest[,1]))
  # ysTest = simulationData$zTest[,1]
  predPts = rbind(predPts, cbind(simulationData$xGrid, simulationData$yGrid))
  ysTest = c(simulationData$zTest[,1], simulationData$zTestRural[,1], simulationData$zTestUrban[,1], simulationData$zGrid[,1])
  
  # fit the model
  if(assumeMeanZero)
    m=0
  else
    m=1
  time = system.time(fit <- fitLKStandard(coords, ys, predCoords=predPts, nLayer=nLayer, NC=NC,
                                          nBuffer=nBuffer, normalize=normalize, fixedFunctionArgs=list(m=m), 
                                          xRangeDat=xRange, yRangeDat=yRange, separatea.wght=separatea.wght, 
                                          doMatern=doMatern, fixNu=fixNu))
  mod = fit$mod
  preds = fit$preds
  predSDs = fit$sigmas
  parameterSummaryTable = fit$parameterSummaryTable
  parSim = fit$parSim
  
  # print out the total time
  print(paste0("Total time: ", time[3]))
  
  # show a model summary
  print(summary(mod))
  
  # function for determining if points are in correct range
  inRange = function(pts, rangeShrink=0) {
    inX = (rangeShrink < pts[,1]) & (pts[,1] < 1-rangeShrink)
    inY = (rangeShrink < pts[,2]) & (pts[,2] < 1-rangeShrink)
    inX & inY
  }
  
  # show predictive surface, SD, and data
  
  pdf(file=paste0("Figures/mixtureLKPreds", extraPlotName, ".pdf"), width=8, height=8)
  par(mfrow=c(2,2))
  
  # obsInds = 1:n
  # predInds = (n+1):(n+mx*my)
  # coefInds = (n+mx*my+1):(n+mx*my+nx*ny)
  colRangeDat = range(c(ys, preds))
  colRangeSD = range(range(predSDs[inRange(predPts)]))
  
  quilt.plot(coords, ys, main="Observations", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
  quilt.plot(predPts[,1], predPts[,2], preds, main="Prediction Mean", zlim=colRangeDat, 
             xlim=xRangeDat, ylim=yRangeDat)
  
  plot.new()
  quilt.plot(predPts[,1], predPts[,2], predSDs, main="Prediction SD", 
             xlim=xRangeDat, ylim=yRangeDat)
  dev.off()
  
  testIndices = (length(preds) - length(ysTest) + 1):length(preds)
  leftOutIndices = (length(preds) - length(ysTest) + 1):(length(preds) - length(ysTest) + length(ys))
  gridIndices = (length(preds) - length(ysTest) + length(ys)):length(preds)
  leftOutIndicesTest = match(leftOutIndices, testIndices)
  gridIndicesTest = match(gridIndices, testIndices)
  
  pdf(file=paste0("Figures/mixtureLKLeftOutResidualsGrid", extraPlotName, ".pdf"), width=5, height=5)
  plot(preds[gridIndices], ysTest[gridIndicesTest]-preds[gridIndices], pch=19, cex=.5, col="blue", main="Residuals versus fitted (Grid)", 
       ylab="Residuals", xlab="Fitted")
  abline(h=0, lty=2)
  dev.off()
  
  pdf(file=paste0("Figures/mixtureLKLeftOutResidualsLeftOut", extraPlotName, ".pdf"), width=5, height=5)
  plot(preds[leftOutIndices], ysTest[leftOutIndicesTest]-preds[leftOutIndices], pch=19, cex=.5, col="blue", main="Residuals versus fitted (Left out)", 
       ylab="Residuals", xlab="Fitted")
  abline(h=0, lty=2)
  dev.off()
  
  # calculate true effective range and marginal variance:
  # marginalVar = rho/(4*pi * kappa^2)
  # marginalVar = getMultiMargVar(kappa, rho, nLayer=nLayer, nu=nu, xRange=xRangeBasis, 
  #                               yRange=yRangeBasis, nx=nx, ny=ny)[1]
  
  # plot estimates on interpretable scale (effective range, marginal variance)
  
  # transform all the hyperparameter samples to their relevant values
  totalVariance = fit$totalVariance
  lambdaVals = fit$lambdaVals
  rhoVals = fit$rhoVals
  nuggetVarVals = fit$nuggetVarVals
  a.wghtVals = fit$a.wghtVals
  alphaVals = fit$alphaVals
  nuVals = fit$nuVals
  
  par(mfrow=c(1,1))
  # plot(mod$marginals.hyperpar$`Theta1 for field`, type="l", main="Marginal for log range")
  pdf(file=paste0("Figures/mixtureLKVar", extraPlotName, ".pdf"), width=5, height=5)
  hist(rhoVals, main="Estimate of spatial variance", freq=FALSE, breaks=20)
  abline(v=1 + sqrt(.1), col="green")
  abline(v=quantile(probs=c(.025, .975), rhoVals), col="purple", lty=2)
  dev.off()
  # plot(mod$marginals.hyperpar$`Theta2 for field`, type="l", main="Marginal for log variance")
  pdf(file=paste0("Figures/mixtureLKSigma2", extraPlotName, ".pdf"), width=5, height=5)
  hist(nuggetVarVals, main="Estimate of error variance", freq=FALSE, breaks=20)
  abline(v=sigma2, col="green")
  abline(v=quantile(probs=c(.025, .975), nuggetVarVals), col="purple", lty=2)
  dev.off()
  # for(i in 1:length(covNames)) {
  #   XMarginal = XMarginals[[i]]
  #   pdf(file=paste0("Figures/mixtureLK", covNames[i], ".pdf"), width=5, height=5)
  #   plot(XMarginal, type="l", main="Marginal for fixed effect")
  #   abline(v=0, col="green")
  #   abline(v=inla.qmarginal(c(.025, .975), XMarginal), col="purple", lty=2)
  #   dev.off()
  # }
  
  pdf(file=paste0("Figures/mixtureLKRho", extraPlotName, ".pdf"), width=5, height=5)
  hist(rhoVals, main=TeX("Estimate of $\\rho$"), xlab=TeX("$\\rho$"), breaks=20, freq=FALSE)
  abline(v=quantile(probs=c(.025, .975), rhoVals), col="purple", lty=2)
  abline(v=rho, col="green")
  dev.off()
  
  pdf(file=paste0("Figures/mixtureLKa.wght", extraPlotName, ".pdf"), width=5, height=5)
  hist(a.wghtVals, xlab="a.wght", main="Marginal for a.wght", breaks=20, freq=FALSE)
  abline(v=quantile(probs=c(.025, .975), a.wghtVals), col="purple", lty=2)
  dev.off()
  
  # pdf(file=paste0("Figures/mixtureLKRho", extraPlotName, ".pdf"), width=5, height=5)
  # hist(rhos, xlab="rho", main="Marginal for Rho", breaks=1000, freq=F, xlim=c(0, quantile(probs=.95, rhos)))
  # abline(v=rho, col="green")
  # dev.off()
  
  ## Now get distribution for the alpha parameters
  xSamples = alphaVals
  for(l in 1:nLayer) {
    pdf(file=paste0("Figures/mixtureLKAlpha", l, extraPlotName, ".pdf"), width=5, height=5)
    hist(xSamples[l,], xlab=TeX(paste0("$\\alpha_", l, "$")), main=TeX(paste0("Marginal for $\\alpha_", l, "$")), breaks=20, freq=F, xlim=c(0,1))
    abline(v=mean(xSamples[l,]), col="purple", lty=1)
    abline(v=quantile(probs=c(.025, .975), xSamples[l,]), col="purple", lty=2)
    dev.off()
  }
  
  ## plot covariance and correlation functions
  # first get the true covariance an correlation functions
  spatialCovFun = spatialCorFun
  mixtureCovFun = function(x) {
    out = spatialCorFun(x)[1,]
    out[x == 0] = 1 + sigma2
    out
  }
  mixtureCorFun = function(x) { mixtureCovFun(x) * (1 / (1 + sigma2)) }
  
  # compute the covariance function for many different hyperparameter samples
  out = covarianceDistributionLK(mod$LKinfo, alphaVals, lambdaVals, a.wghtVals, rhoVals)
  d = out$d
  sortI = sort(d, index.return=TRUE)$ix
  d = d[sortI]
  covMean = out$cov[sortI]
  upperCov=out$upperCov[sortI]
  lowerCov=out$lowerCov[sortI]
  covMat=out$covMat[sortI,]
  corMean = out$cor[sortI]
  upperCor=out$upperCor[sortI]
  lowerCor=out$lowerCor[sortI]
  corMat=out$corMat[sortI,]
  corMatNoNugget=out$corMatNoNugget[sortI,]
  covInfo = list(d=d, covMean=covMean, upperCov=upperCov, lowerCov=lowerCov, covMat=covMat, 
                 corMean=corMean, upperCor=upperCor, lowerCor=lowerCor, corMat=corMat)
  
  getEffectiveRange = function(ds, cors) {
    firstI = match(TRUE, cors<.1)
    if(is.na(firstI)) {
      warning("effective range larger than spatial domain diameter. Returning spatial domain diameter as effective range")
      firstI = length(ds)
    } else if(firstI == 1)
      warning("effective range is extremely small and more resolution required to accurately determine it")
    ds[firstI]
  }
  effectiveRanges = apply(corMatNoNugget, 2, getEffectiveRange, ds=d)
  
  pdf(file=paste0("Figures/mixtureLKEffRange", extraPlotName, ".pdf"), width=5, height=5)
  hist(effectiveRanges, main="Estimated effective range", freq=FALSE, breaks=20)
  abline(v=effectiveRange, col="green")
  abline(v=quantile(probs=c(.025, .975), effectiveRanges), col="purple", lty=2)
  dev.off()
  
  # plot the covariance function
  yRange = range(c(covMean, lowerCov, upperCov, mixtureCovFun(d)))
  pdf(file=paste0("Figures/mixtureLKCov", extraPlotName, ".pdf"), width=5, height=5)
  plot(d, covMean, type="l", main="Estimated covariance function", xlab="Distance", ylab="Covariance", 
       ylim=yRange)
  lines(d, lowerCov, lty=2)
  lines(d, upperCov, lty=2)
  lines(d, mixtureCovFun(d), col="green")
  legend("topright", c("Truth", "Estimate", "80% CI"), lty=c(1, 1, 2), col=c("green", "black", "black"))
  dev.off()
  
  pdf(file=paste0("Figures/mixtureLKCor", extraPlotName, ".pdf"), width=5, height=5)
  plot(d, corMean, type="l", main="Estimated correlation function", xlab="Distance", ylab="Covariance", 
       ylim=c(0, 1))
  lines(d, lowerCor, lty=2)
  lines(d, upperCor, lty=2)
  lines(d, mixtureCorFun(d), col="green")
  legend("topright", c("Truth", "Estimate", "80% CI"), lty=c(1, 1, 2), col=c("green", "black", "black"))
  dev.off()
  
  # get scoring rules
  testIndices = (length(preds) - length(ysTest) + 1):length(preds)
  leftOutIndices = (length(preds) - length(ysTest) + 1):(length(preds) - length(ysTest) + length(simulationData$zTest[,1]))
  gridIndices = (length(preds) - length(ysTest) + length(simulationData$zTest[,1]) + 1):length(preds)
  leftOutIndicesTest = match(leftOutIndices, testIndices)
  gridIndicesTest = match(gridIndices, testIndices)
  
  # first calculate scoring rules at grid points
  truth = ysTest[gridIndicesTest]
  est = preds[gridIndices]
  vars = predSDs[gridIndices]^2
  lower = mod$lower[gridIndices]
  upper = mod$upper[gridIndices]
  
  # compute nearest neighbor distances and scores as a function of them
  gridPts = predPts[gridIndices,]
  distMat = rdist(coords, gridPts)
  nndists = apply(distMat, 2, function(x) {min(x[x != 0])})
  print("Binned grid scores:")
  gridScoringRules = getScores(truth, est, vars, lower, upper, distances=nndists, breaks=distanceBreaks)
  print(gridScoringRules$pooledResults)
  print(gridScoringRules$binnedResults)
  
  # now calculate square rules at left out points
  truth = ysTest[leftOutIndicesTest]
  est = preds[leftOutIndices]
  vars = predSDs[leftOutIndices]^2
  lower = mod$lower[leftOutIndices]
  upper = mod$upper[leftOutIndices]
  leftOutScoringRules = getScores(truth, est, vars, lower, upper)
  leftOutScoringRules = data.frame(c(leftOutScoringRules, Time=time[3]))
  print("Binned left out scores:")
  print(leftOutScoringRules)
  
  scoringRules = list(gridScoringRules=gridScoringRules, leftOutScoringRules=leftOutScoringRules)
  # 
  # # get scoring rules
  # truth = ysTest
  # est = preds[testIndices]
  # vars = predSDs[testIndices]^2
  # lower = mod$lower[testIndices]
  # upper = mod$upper[testIndices]
  # 
  # # compute nearest neighbor distances and scores as a function of them
  # testPts = predPts[testIndices,]
  # distMat = rdist(coords, testPts)
  # nndists = apply(distMat, 2, function(x) {min(x[x != 0])})
  # print("Binned scores:")
  # scoringRules = getScores(truth, est, vars, lower, upper, distances=nndists, breaks=distanceBreaks)
  # scoringRules$pooledResults = data.frame(c(scoringRules$pooledResults, Time=time[3]))
  # print(scoringRules$binnedResults)
  
  
  
  # print("Grid scores:")
  # print(getScores(truth[gridTestI], est[gridTestI], vars[gridTestI], lower[gridTestI], upper[gridTestI], 
  #                 distances=nndists[gridTestI], breaks=distanceBreaks)$binnedResults)
  # 
  # # plot scores as a function of distance
  # distanceScores = getScores(truth[gridTestI], est[gridTestI], vars[gridTestI], lower[gridTestI], upper[gridTestI], 
  #                            distances=nndists[gridTestI], breaks=distanceBreaks)$binnedResults
  
  # pdf(file=paste0("Figures/mixtureSPDEScoreBias", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(distanceScores$NNDist, distanceScores$Bias, pch=19, col="blue", main="Bias", ylab="Bias", xlab="Nearest neighbor distance (km)")
  # abline(h=0, lty=2)
  # dev.off()
  # 
  # pdf(file=paste0("Figures/mixtureSPDEScoreVar", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(distanceScores$NNDist, distanceScores$Var, pch=19, col="blue", main="Variance", ylab="Variance", xlab="Nearest neighbor distance (km)", 
  #      ylim=c(0, max(distanceScores$Var)))
  # abline(h=0, lty=2)
  # dev.off()
  # 
  # pdf(file=paste0("Figures/mixtureSPDEScoreMSE", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(distanceScores$NNDist, distanceScores$MSE, pch=19, col="blue", main="MSE", ylab="MSE", xlab="Nearest neighbor distance (km)", 
  #      ylim=c(0, max(distanceScores$MSE)))
  # abline(h=0, lty=2)
  # dev.off()
  # 
  # pdf(file=paste0("Figures/mixtureSPDEScoreRMSE", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(distanceScores$NNDist, distanceScores$RMSE, pch=19, col="blue", main="RMSE", ylab="RMSE", xlab="Nearest neighbor distance (km)", 
  #      ylim=c(0, max(distanceScores$RMSE)))
  # abline(h=0, lty=2)
  # dev.off()
  # 
  # pdf(file=paste0("Figures/mixtureSPDEScoreCRPS", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(distanceScores$NNDist, distanceScores$CRPS, pch=19, col="blue", main="CRPS", ylab="CRPS", xlab="Nearest neighbor distance (km)", 
  #      ylim=c(0, max(distanceScores$CRPS)))
  # abline(h=0, lty=2)
  # dev.off()
  # 
  # pdf(file=paste0("Figures/mixtureSPDEScoreCvg", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(distanceScores$NNDist, distanceScores$Coverage, pch=19, col="blue", main="80% Coverage", ylab="80% Coverage", xlab="Nearest neighbor distance (km)", 
  #      ylim=c(0, max(distanceScores$Coverage)))
  # abline(h=0, lty=2)
  # dev.off()
  # 
  # pdf(file=paste0("Figures/mixtureSPDEScoreWidth", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(distanceScores$NNDist, distanceScores$Width, pch=19, col="blue", main="Width", ylab="Width", xlab="Nearest neighbor distance (km)", 
  #      ylim=c(0, max(distanceScores$Width)))
  # abline(h=0, lty=2)
  # dev.off()
  
  # get aggregated predictions
  # A = t(sapply(1:(mx*my), getARow))
  # A = sweep(A, 1, rowSums(A), "/")
  # mx = 100
  # my = 100
  # predPts = make.surface.grid(list(x=seq(xRangeDat[1], xRangeDat[2], l=mx), y=seq(yRangeDat[1], yRangeDat[2], l=my)))
  gridCoords = cbind(simulationData$xGrid, simulationData$yGrid)
  # testIndices = (length(preds) - length(ysTest) + 1):length(preds)
  # gridIndices = (length(preds) - length(simulationData$xGrid) + 1):length(preds)
  # gridIndicesEst = (length(est) - length(simulationData$xGrid) + 1):length(est)
  
  A = makeNumericalIntegralMat(gridCoords, mx=3, my=3)
  aggregatedPreds = A %*% preds[gridIndices]
  
  truth = A %*% ysTest[gridIndicesTest]
  est = A %*% preds[gridIndices]
  aggregatedPredMat = A %*% fit$predMat[gridIndices,]
  vars = apply(aggregatedPredMat, 1, var)
  sds = apply(aggregatedPredMat, 1, sd)
  lower = apply(aggregatedPredMat, 1, function(x) {quantile(x, probs=.1)})
  upper = apply(aggregatedPredMat, 1, function(x) {quantile(x, probs=.9)})
  predictionMatrix = data.frame(Truth=truth, Est=est, SDs=sds, Lower=lower, Upper=upper)
  print("Aggregated prediction summary statistics:")
  print(predictionMatrix)
  
  print("Pooled aggregated scores:")
  pooledAggregatedScores = getScores(truth, est, vars, lower, upper)
  print(pooledAggregatedScores)
  print("Left out region aggregated scores:")
  leftOutAggregatedScores = getScores(truth[5], est[5], vars[5], lower[5], upper[5])
  leftOutAggregatedScores$Var = sds[5]
  names(leftOutAggregatedScores)[2] = "Predictive.SD"
  print(leftOutAggregatedScores)
  print("Included regions aggregated scores:")
  includedAggregatedScores = getScores(truth[-5], est[-5], vars[-5], lower[-5], upper[-5])
  print(includedAggregatedScores)
  
  aggregatedScoringRules = list(pooledAggregatedScores=pooledAggregatedScores, leftOutAggregatedScores=leftOutAggregatedScores, 
                                includedAggregatedScores=includedAggregatedScores)
  
  analysisTime = proc.time()[3] - startTime
  
  fit$mod = NULL
  if(!gscratch)
    save(scoringRules, fit, covInfo, predictionMatrix, aggregatedScoringRules, analysisTime, file=paste0("savedOutput/simulations/mixtureLK", plotNameRoot, ".RData"))
  else
    save(scoringRules, fit, covInfo, predictionMatrix, aggregatedScoringRules, analysisTime, file=paste0("/work/johnpai/mixtureLK", plotNameRoot, ".RData"))
  
  invisible(list(scoringRules, fit, covInfo, predictionMatrix, aggregatedScoringRules, analysisTime))
}

# seed=1, nLayer=3, nx=20, ny=nx, nu=1, assumeMeanZero=TRUE, 
# nBuffer=5, normalize=TRUE, NC=14, testCovs=TRUE, 
# printVerboseTimings=FALSE, n=900, separatea.wght=FALSE, 
# plotNameRoot="", doMatern=FALSE, fixNu=FALSE, thetas=c(.08, .8) / sqrt(8), 
# testfrac=.1, leaveOutRegion=TRUE, sigma2 = 0.1^2, extraPlotName=plotNameRoot
testLKModelMixtureMultiple = function(seed=1, nSamples=100, gscratch=TRUE, loadResults=FALSE, startI=1, endI=nSamples, ...) {
  # set random seeds for each simulation
  set.seed(seed)
  allSeeds = sample(1:1000000, nSamples, replace = FALSE)
  
  # call testLKINLAModelMixture for each simulation requested
  temp = function(i) {
    print(paste0("Beginning simulation ", i, "/", nSamples))
    thisPlotNameRoot = paste0("sim", i)
    do.call("testLKModelMixture", c(list(seed = allSeeds[i], plotNameRoot=thisPlotNameRoot, gscratch=gscratch), list(...)))
  }
  if(!loadResults)
    sapply(startI:endI, temp)
  
  if(startI != 1 || endI != nSamples)
    return(invisible(NULL))
  
  # load in the results
  allScoringRulesGrid = list()
  allScoringRulesLeftOut = list()
  allFits = list()
  allCovInfo = list()
  allPredictionMatrices = list()
  allAggregatedScoringRules = list()
  allAnalysisTimes = c()
  for(i in 1:nSamples) {
    if(!gscratch)
      out = load(paste0("savedOutput/simulations/mixtureLKsim", i, ".RData"))
    else
      out = load(paste0("/work/johnpai/mixtureLKsim", i, ".RData"))
    allScoringRulesGrid = c(allScoringRulesGrid, list(scoringRules$gridScoringRules))
    allScoringRulesLeftOut = c(allScoringRulesLeftOut, list(scoringRules$leftOutScoringRules))
    allFits = c(allFits, list(fit))
    allCovInfo = c(allCovInfo, list(covInfo))
    allPredictionMatrices = c(allPredictionMatrices, list(predictionMatrix))
    allAggregatedScoringRules = c(allAggregatedScoringRules, list(aggregatedScoringRules))
    allAnalysisTimes = c(allAnalysisTimes, analysisTime)
  }
  
  ##### average results from each simulation
  # pointwise scoring rules
  allPooledScoringRulesGrid = do.call("rbind", lapply(allScoringRulesGrid, function(x) {x$pooledResults}))
  allBinnedScoringRulesGrid = lapply(allScoringRulesGrid, function(x) {x$binnedResults})
  binnedScoringRulesGrid = averageBinnedScores(allBinnedScoringRulesGrid)
  ns = binnedScoringRulesGrid[,2]
  pooledScoringRulesGrid = apply(binnedScoringRulesGrid, 2, function(x) {sum(x * (ns / sum(ns)))})
  pooledScoringRulesGrid = as.data.frame(matrix(pooledScoringRulesGrid, nrow=1))
  names(pooledScoringRulesGrid) = names(binnedScoringRulesGrid)
  
  fullPooledScoringRulesLeftOut = do.call("rbind", allScoringRulesLeftOut)
  pooledScoringRulesLeftOut = colMeans(fullPooledScoringRulesLeftOut)
  
  # covInfo
  # covInfo = list(d=d, covMean=covMean, upperCov=upperCov, lowerCov=lowerCov, covMat=covMat, 
  #                corMean=corMean, upperCor=upperCor, lowerCor=lowerCor, corMat=corMat)
  d = allCovInfo[[1]]$d
  covMean = rowMeans(do.call("cbind", lapply(allCovInfo, function(x) {x$covMean})))
  upperCov = rowMeans(do.call("cbind", lapply(allCovInfo, function(x) {x$upperCov})))
  lowerCov = rowMeans(do.call("cbind", lapply(allCovInfo, function(x) {x$lowerCov})))
  corMean = rowMeans(do.call("cbind", lapply(allCovInfo, function(x) {x$corMean})))
  upperCor = rowMeans(do.call("cbind", lapply(allCovInfo, function(x) {x$upperCor})))
  lowerCor = rowMeans(do.call("cbind", lapply(allCovInfo, function(x) {x$lowerCor})))
  
  # aggregated scoring rules
  # aggregatedScoringRules = list(pooledAggregatedScores=pooledAggregatedScores, leftOutAggregatedScores=leftOutAggregatedScores, 
  #                               includedAggregatedScores=includedAggregatedScores)
  # predictionMatrix = data.frame(Truth=truth, Est=est, SDs=sds, Lower=lower, Upper=upper)
  fullPredictionMatrix = do.call("rbind", allPredictionMatrices)
  leftOutIndices = seq(5, nrow(fullPredictionMatrix), by=9)
  leftOutPredictionMatrix = fullPredictionMatrix[leftOutIndices, ]
  leftInPredictionMatrix = fullPredictionMatrix[-leftOutIndices, ]
  leftOutScores = getScores(leftOutPredictionMatrix[,1], leftOutPredictionMatrix[,2], leftOutPredictionMatrix[,3]^2)
  leftInScores = getScores(leftInPredictionMatrix[,1], leftInPredictionMatrix[,2], leftInPredictionMatrix[,3]^2)
  aggregatedScores = getScores(fullPredictionMatrix[,1], fullPredictionMatrix[,2], fullPredictionMatrix[,3]^2)
  
  ##### Save results
  if(!gscratch) {
    save(#allScoringRules, 
         allFits, 
         allCovInfo, 
         allPredictionMatrices, 
         allAggregatedScoringRules, 
         allAnalysisTimes, 
         binnedScoringRulesGrid, 
         pooledScoringRulesGrid, 
         fullPooledScoringRulesLeftOut, 
         pooledScoringRulesLeftOut, 
         covMean, 
         upperCov, 
         lowerCov, 
         corMean, 
         upperCor, 
         lowerCor, 
         fullPredictionMatrix, 
         leftOutPredictionMatrix, 
         leftInPredictionMatrix, 
         leftOutScores, 
         leftInScores, 
         aggregatedScores, 
         file=paste0("savedOutput/simulations/mixtureLKAll_nsim", nSamples, ".RData"))
  } else {
    save(#allScoringRules, 
         allFits, 
         allCovInfo, 
         allPredictionMatrices, 
         allAggregatedScoringRules, 
         binnedScoringRulesGrid, 
         pooledScoringRulesGrid, 
         fullPooledScoringRulesLeftOut, 
         pooledScoringRulesLeftOut, 
         covMean, 
         upperCov, 
         lowerCov, 
         corMean, 
         upperCor, 
         lowerCor, 
         fullPredictionMatrix, 
         leftOutPredictionMatrix, 
         leftInPredictionMatrix, 
         leftOutScores, 
         leftInScores, 
         aggregatedScores, 
         file=paste0("/work/johnpai/simulations/mixtureLKAll_nsim", nSamples, ".RData"))
  }
}

# recalculate covariances for the LatticeKrig model
fixLKModelMixtureCovariance = function(seed=1, nSamples=100, gscratch=TRUE, startI=1, endI=nSamples, ...) {
  # set random seeds for each simulation
  set.seed(seed)
  allSeeds = sample(1:1000000, nSamples, replace = FALSE)
  
  # load in the results
  for(i in startI:endI) {
    set.seed(allSeeds[i])
    print(paste0("Recalculating covariance for simulation ", i, "/", nSamples))
    if(!gscratch)
      out = load(paste0("savedOutput/simulations/mixtureLKsim", i, ".RData"))
    else
      out = load(paste0("/work/johnpai/mixtureLKsim", i, ".RData"))
    
    # calculate covariance
    out = covarianceDistributionLK(fit$LKinfo, fit$alphaVals, fit$lambdaVals, fit$a.wghtVals, fit$rhoVals)
    d = out$d
    sortI = sort(d, index.return=TRUE)$ix
    d = d[sortI]
    covMean = out$cov[sortI]
    upperCov=out$upperCov[sortI]
    lowerCov=out$lowerCov[sortI]
    covMat=out$covMat[sortI,]
    corMean = out$cor[sortI]
    upperCor=out$upperCor[sortI]
    lowerCor=out$lowerCor[sortI]
    corMat=out$corMat[sortI,]
    corMatNoNugget=out$corMatNoNugget[sortI,]
    covInfoNew = list(d=d, covMean=covMean, upperCov=upperCov, lowerCov=lowerCov, covMat=covMat, 
                   corMean=corMean, upperCor=upperCor, lowerCor=lowerCor, corMat=corMat)
    
    # save results
    covInfoOld = covInfo
    covInfo = covInfoNew
    if(!gscratch)
      save(scoringRules, fit, covInfo, covInfoOld, predictionMatrix, aggregatedScoringRules, file=paste0("savedOutput/simulations/mixtureLK", i, ".RData"))
    else
      save(scoringRules, fit, covInfo, covInfoOld, predictionMatrix, aggregatedScoringRules, file=paste0("/work/johnpai/mixtureLK", i, ".RData"))
  }
}

# tests the fitSPDE function using data simulated from the LK model
# buffer: buffer distance between domain edge of basis lattice and domain edge of data.
# n: number of observations
# xRange: range of x coordinates
# yRange: range of y coordinates
# nx: number of basis function lattice points in x directions
# NOTE: ny is determined automatically to match scale of x lattice points
# Xmat: design matrix
# ys: observations
# first.time: is first time evaluating function.  User should always set to FALSE
# thetas: originally was c(.1, 3), but 3 is too large
testSPDEModelMixture = function(seed=1, nx=20, ny=nx, assumeMeanZero=TRUE, 
                                testCovs=FALSE, n=900, thetas=NULL, 
                                int.strategy="auto", strategy="gaussian", 
                                nPostSamples=1000, mesh=NULL, 
                                prior=NULL, testfrac=1/9, nu=1, nHyperSamples=1000, 
                                plotNameRoot="", sigma2 = 0.1^2, useKenya=FALSE, 
                                urbanOverSamplefrac=0, leaveOutRegion=TRUE, gscratch=FALSE) {
  set.seed(seed)
  startTime = proc.time()[3]
  
  if(useKenya)
    distanceBreaks = seq(0, 300, l=20)
  else
    distanceBreaks = seq(0, 0.5, l=20)
  
  # set plotNameRoot
  plotNameRoot = paste0(plotNameRoot, "_n", n, "_nu", nu, "_nugV", round(sigma2, 2), "_Kenya", useKenya, 
                        "_noInt", assumeMeanZero, "_urbOversamp", round(urbanOverSamplefrac, 4))
  
  # set true parameter values
  if(useKenya) {
    if(is.null(thetas))
      thetas=c(.08, .8) * (1000/2) / sqrt(8)
  } else {
    if(is.null(thetas))
      thetas=c(.08, .8) / sqrt(8)
  }
  rho = 1
  effectiveRange = thetas * sqrt(8)
  
  # set the SPDE mesh if necessary
  if(is.null(mesh)) {
    if(useKenya)
      mesh = getSPDEMeshKenya()
    else
      mesh = getSPDEMesh()
  }
  
  # load data set if necessary
  if(is.null(n)) {
    out = load("mixtureDataSet.RData")
  } else {
    mixtureCorFun = function(x) {0.5 * stationary.cov(x, theta=thetas[1], Covariance="Matern", smoothness=nu) + 
        0.5 * stationary.cov(x, theta=thetas[2], Covariance="Matern", smoothness=nu)}
    nTest = round(testfrac * n)
    if(leaveOutRegion) {
      simulationData = getSimulationDataSetsGivenCovarianceTest(mixtureCorFun, nTotal=n, nTest=nTest, marginalVar=rho, errorVar=sigma2, 
                                                                nDataSets=2, plotNameRoot=paste0("(0.5*Matern(", thetas[1], ") + 0.5*Matern(", thetas[2], "))"), fileNameRoot="mix", 
                                                                saveDataSetPlot=FALSE, doPredGrid=TRUE)
    } else {
      simulationData = getSimulationDataSetsGivenCovariance(mixtureCorFun, nTotal=n, nTest=nTest, marginalVar=rho, errorVar=sigma2, 
                                                            nDataSets=2, plotNameRoot=paste0("(0.5*Matern(", thetas[1], ") + 0.5*Matern(", thetas[2], "))"), fileNameRoot="mix", 
                                                            saveDataSetPlot=FALSE, useKenyaLocations=useKenya, urbanOverSamplefrac=urbanOverSamplefrac)
    }
    
  }
  coords = cbind(simulationData$xTrain[,1], simulationData$yTrain[,1])
  ys = simulationData$zTrain[,1]
  
  # generate lattice and simulate observations
  # coords = matrix(runif(2*n), ncol=2)
  if(useKenya) {
    xRangeDat = simulationData$xRange
    yRangeDat = simulationData$yRange
  } else {
    xRangeDat = c(-1, 1)
    yRangeDat = c(-1, 1)
  }
  
  AObs = inla.spde.make.A(mesh, loc = coords)
  # Q = makeQ(kappa=kappa, rho=rho, latInfo, alphas=alphas, normalized=normalize, fastNormalize=fastNormalize) 
  # L = as.matrix(t(chol(solve(Q))))
  # zsims = matrix(rnorm(nrow(Q)), ncol=1)
  # fieldSims = L %*% zsims
  # ys = as.numeric(AObs %*% fieldSims) + 1 # add a constant unit mean term to be estimated by INLA
  # # ys = 1 + as.numeric(AObs %*% fieldSims) + coords[,1] # x-valued mean term to be estimated by INLA
  # errs = rnorm(n, sd=sqrt(sigma2))
  # ys = ys + errs
  
  # plot the observations
  pdf(file=paste0("Figures/mixtureSPDEObservations", plotNameRoot, ".pdf"), width=5, height=5)
  par(mfrow=c(1,1))
  quilt.plot(coords, ys)
  dev.off()
  
  # make prediction coordinates on a grid, and add testing points
  mx = 100
  my = 100
  predPts = make.surface.grid(list(x=seq(xRangeDat[1], xRangeDat[2], l=mx), y=seq(yRangeDat[1], yRangeDat[2], l=my)))
  if(useKenya) {
    # remove grid points outside of Kenya national boundaries
    load("../U5MR/adminMapData.RData")
    polys = adm0@polygons
    kenyaPoly = polys[[1]]@Polygons[[77]]@coords
    kenyaPolyProj = projKenya(kenyaPoly)
    inKenya = in.poly(predPts, kenyaPolyProj)
    predPts = predPts[inKenya,]
    
    # add other testing locations to matrix of prediction locations and remember which 
    predPts = rbind(predPts, cbind(simulationData$xTest[,1], simulationData$yTest[,1]))
    predPts = rbind(predPts, cbind(simulationData$xTestRural[,1], simulationData$yTestRural[,1]))
    predPts = rbind(predPts, cbind(simulationData$xTestUrban[,1], simulationData$yTestUrban[,1]))
    plotGridI = 1:sum(inKenya)
    overallTestI = simulationData$overallTestI
    ruralTestI = simulationData$ruralTestI
    urbanTestI = simulationData$urbanTestI
    gridTestI = (max(urbanTestI) + 1):(max(urbanTestI) + length(simulationData$xGrid))
  } else {
    predPts = rbind(predPts, cbind(simulationData$xTest[,1], simulationData$yTest[,1]))
  }
  predPts = rbind(predPts, cbind(simulationData$xGrid, simulationData$yGrid))
  ysTest = c(simulationData$zTest[,1], simulationData$zTestRural[,1], simulationData$zTestUrban[,1], simulationData$zGrid[,1])
  
  # generate hyperparameters based on median and quantiles of inverse exponential and inverse gamma
  # priorPar = getPrior(.1, .1, 10)
  # generate hyperparameters for pc priors
  # median effective range is .4 or 200 for Kenya data (a fifth of the spatial domain diameter), median spatial variance is 1
  # priorPar = getPCPrior(.4, .5, 1) 
  if(is.null(prior)) {
    if(!useKenya) {
      prior = getSPDEPrior(mesh, U=1, alpha=.5, medianRange=.4)
    } else {
      prior = getSPDEPrior(mesh, U=1, alpha=.5, medianRange=.4 * 1000 / 2)
    }
  }
  
  X = matrix(rep(1, nrow(coords)), ncol=1)
  # X = matrix(coords[,1], ncol=1)
  XPred = matrix(rep(1, nrow(predPts)), ncol=1)
  
  # add linear terms in lat/lon to covariate matrices if requested
  if(testCovs) {
    X = cbind(X, coords)
    XPred = cbind(XPred, predPts)
  }
  
  if(assumeMeanZero) {
    X = NULL
    XPred = NULL
  }
  
  # show priors on effective correlation, marginal variance, and error variance:
  if(!useKenya)
    xs1 = seq(.01, 1, l=500)
  else
    xs1 = seq(1, 1000, l=500)
  # pdf(file=paste0("Figures/mixtureSPDEPriorEffRange", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(xs1, dinvexp(xs1, rate=priorPar$corScalePar), type="l", col="blue", 
  #      xlab="Effective Correlation Range", main="Effective Correlation Prior", 
  #      ylab="Prior Density")
  # abline(v=qinvexp(.5, rate=priorPar$corScalePar), col="red")
  # dev.off()
  
  # if(priorPar$priorType == "orig") {
  #   xs2 = seq(.01, 10.5, l=500)
  #   pdf(file=paste0("Figures/mixtureSPDEPriorMargVar", plotNameRoot, ".pdf"), width=5, height=5)
  #   plot(xs2, invgamma::dinvgamma(xs2, shape=priorPar$varPar1, rate=priorPar$varPar2), type="l", col="blue", 
  #        xlab="Marginal Variance", main="Marginal Variance Prior", 
  #        ylab="Prior Density")
  #   abline(v=qinvgamma(.1, shape=priorPar$varPar1, rate=priorPar$varPar2), col="red")
  #   abline(v=qinvgamma(.9, shape=priorPar$varPar1, rate=priorPar$varPar2), col="red")
  #   dev.off()
  # } else if(priorPar$priorType == "pc") {
  #   xs2 = seq(.01, 11.5, l=500)
  #   pdf(file=paste0("Figures/mixtureSPDEPriorMargVar", plotNameRoot, ".pdf"), width=5, height=5)
  #   plot(xs2, dpcvar(xs2, alpha=priorPar$alpha, u=priorPar$u), type="l", col="blue", 
  #        xlab="Marginal Variance", main="Marginal Variance Prior", 
  #        ylab="Prior Density")
  #   abline(v=qpcvar(.1, alpha=priorPar$alpha, u=priorPar$u), col="red")
  #   abline(v=qpcvar(.9, alpha=priorPar$alpha, u=priorPar$u), col="red")
  #   abline(v=1, col="green")
  #   dev.off()
  # }
  
  # xs2 = seq(.001, invgamma::qinvgamma(.905, shape=0.1, rate=0.1), l=500)
  # pdf(file="Figures/mixtureSPDEPriorErrorVar.pdf", width=5, height=5)
  # plot(xs2, invgamma::dinvgamma(xs2, shape=0.1, rate=0.1), type="l", col="blue", 
  #      xlab="Error Variance", main="Error Variance Prior", 
  #      ylab="Prior Density")
  # abline(v=invgamma::qinvgamma(.1, shape=0.1, rate=0.1), col="red")
  # abline(v=invgamma::qinvgamma(.9, shape=0.1, rate=0.1), col="red")
  # dev.off()
  
  xs2 = seq(.01, 1, l=500)
  pdf(file=paste0("Figures/mixtureSPDEPriorErrorVar", plotNameRoot, ".pdf"), width=5, height=5)
  plot(xs2, dpcvar(xs2, alpha=.05, u=1), type="l", col="blue", 
       xlab="Marginal Variance", main="Marginal Variance Prior", 
       ylab="Prior Density")
  abline(v=qpcvar(.1, alpha=.05, u=1), col="red")
  abline(v=qpcvar(.9, alpha=.05, u=1), col="red")
  abline(v=sqrt(.1), col="green")
  dev.off()
  
  # browser()
  
  # fit the model
  time = system.time(fit <- fitSPDE(coords, ys, predCoords=predPts, seed=seed, prior=prior, 
                                    xObs=X, xPred=XPred, int.strategy=int.strategy, strategy=strategy, 
                                    mesh=mesh, nPostSamples=nPostSamples))
  mod = fit$mod
  preds=fit$preds
  predSDs=fit$sigmas
  latInfo=fit$latInfo
  latWidth=fit$latWidth
  obsPreds=fit$obsPreds
  obsSDs=fit$obsSDs
  
  # print out the total time
  print(paste0("Total time: ", time[3]))
  
  # show a model summary
  print(summary(mod))
  
  # function for determining if points are in correct range
  if(!useKenya) {
    inRange = function(pts, rangeShrink=0) {
      inX = (rangeShrink < pts[,1]) & (pts[,1] < 1-rangeShrink)
      inY = (rangeShrink < pts[,2]) & (pts[,2] < 1-rangeShrink)
      inX & inY
    }
  } else {
    inRange = function(pts, rangeShrink=0) {
      rep(TRUE, nrow(pts))
    }
  }
  
  # show predictive surface, SD, and data
  
  pdf(file=paste0("Figures/mixtureSPDEPreds", plotNameRoot, ".pdf"), width=8, height=8)
  par(mfrow=c(2,2), mar=c(5.1, 4.1, 4.1, 6))
  
  # obsInds = 1:n
  # predInds = (n+1):(n+mx*my)
  # coefInds = (n+mx*my+1):(n+mx*my+nx*ny)
  # colRangeDat = range(c(ys, obsPreds, preds, coefPreds))
  colRangeDat = range(c(ys, obsPreds, preds))
  colRangeSD = range(c(predSDs, obsSDs))
  quilt.plot(coords, ys, main="True Process", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
  quilt.plot(predPts[,1], predPts[,2], preds, main="Prediction Mean", zlim=colRangeDat, 
             xlim=xRangeDat, ylim=yRangeDat)
  
  quilt.plot(coords, ys, main="Observations", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
  quilt.plot(predPts[,1], predPts[,2], predSDs, main="Prediction SD",
             xlim=xRangeDat, ylim=yRangeDat, zlim=range(predSDs[inRange(predPts, rangeShrink=.03)]))
  dev.off()
  
  pdf(file=paste0("Figures/mixtureSPDELeftOutResiduals", plotNameRoot, ".pdf"), width=5, height=5)
  testIndices = (length(preds) - length(ysTest) + 1):length(preds)
  plot(preds[testIndices], ysTest-preds[testIndices], pch=19, cex=.5, col="blue", main="Residuals versus fitted", 
       ylab="Residuals", xlab="Fitted")
  abline(h=0, lty=2)
  dev.off()
  
  if(useKenya) {
    pdf(file=paste0("Figures/mixtureSPDELeftOutResidualsLabeled", plotNameRoot, ".pdf"), width=5, height=5)
    testIndices = (length(preds) - length(ysTest) + 1):length(preds)
    gridIndices = 
      ylim = range(ysTest-preds[testIndices])
    xlim = range(preds[testIndices])
    plot(preds[testIndices][overallTestI], ysTest[overallTestI]-preds[testIndices][overallTestI], pch=19, cex=.1, col="black", main="Residuals versus fitted", 
         ylab="Residuals", xlab="Fitted", xlim=xlim, ylim=ylim)
    points(preds[testIndices][ruralTestI], ysTest[ruralTestI]-preds[testIndices][ruralTestI], pch=19, cex=.1, col="green")
    points(preds[testIndices][urbanTestI], ysTest[urbanTestI]-preds[testIndices][urbanTestI], pch=19, cex=.1, col="blue")
    points(preds[testIndices][gridTestI], ysTest[gridTestI]-preds[testIndices][gridTestI], pch=19, cex=.1, col="red")
    abline(h=0, lty=2)
    legend("topright", c("Overall", "Rural", "Urban", "Grid"), col=c("black", "green", "blue", "red"), pch=19)
    dev.off()
  }
  
  # calculate true effective range and marginal variance:
  marginalVar = rho
  
  # plot marginals on interpretable scale (effective range, marginal variance)
  effRangeMarg = mod$marginals.hyperpar$`Range for field`
  varMarg = inla.tmarginal(function(x) {x^2}, mod$marginals.hyperpar$`Stdev for field`)
  sigma2Marg = inla.tmarginal(function(x) {1/x}, mod$marginals.hyperpar$`Precision for the Gaussian observations`)
  covNames = names(mod$marginals.fixed)
  
  if(!assumeMeanZero) {
    XMarginals = list()
    for(i in 1:length(covNames)) {
      XMarginal = inla.tmarginal(function(x) {x}, mod$marginals.fixed[[covNames[i]]])
      XMarginals = c(XMarginals, list(XMarginal))
    }
  }
  
  par(mfrow=c(1,1))
  pdf(file=paste0("Figures/mixtureSPDEEffRange", plotNameRoot, ".pdf"), width=5, height=5)
  plot(effRangeMarg, type="l", main="Marginal for effective range")
  abline(v=effectiveRange, col="green")
  abline(v=inla.qmarginal(c(.025, .975), effRangeMarg), col="purple", lty=2)
  dev.off()
  # plot(mod$marginals.hyperpar$`Theta1 for field`, type="l", main="Marginal for log range")
  pdf(file=paste0("Figures/mixtureSPDEVar", plotNameRoot, ".pdf"), width=5, height=5)
  plot(varMarg, type="l", main="Marginal for spatial variance")
  abline(v=marginalVar, col="green")
  abline(v=inla.qmarginal(c(.025, .975), varMarg), col="purple", lty=2)
  dev.off()
  # plot(mod$marginals.hyperpar$`Theta2 for field`, type="l", main="Marginal for log variance")
  pdf(file=paste0("Figures/mixtureSPDESigma2", plotNameRoot, ".pdf"), width=5, height=5)
  plot(sigma2Marg, type="l", main="Marginal for error variance")
  abline(v=sigma2, col="green")
  abline(v=inla.qmarginal(c(.025, .975), sigma2Marg), col="purple", lty=2)
  dev.off()
  
  if(!assumeMeanZero) {
    for(i in 1:length(covNames)) {
      XMarginal = XMarginals[[i]]
      pdf(file=paste0("Figures/mixtureSPDE", covNames[i], plotNameRoot, ".pdf"), width=5, height=5)
      plot(XMarginal, type="l", main="Marginal for fixed effect")
      abline(v=0, col="green")
      abline(v=inla.qmarginal(c(.025, .975), XMarginal), col="purple", lty=2)
      dev.off()
    }
  }
  
  # pdf(file="Figures/mixtureSPDERho.pdf", width=5, height=5)
  # plot(sigma2Marg, type="l", main=TeX("Marginal for $\\rho$"), xlab=TeX("$\\rho$"))
  # abline(v=rho, col="green")
  # dev.off()
  
  ## Now generate marginals for the alpha parameters. In order to do this, we must generate draws from 
  ## the posterior, and transform them back to the probability scale
  out = inla.hyperpar.sample(nHyperSamples, mod, improve.marginals=TRUE)
  
  ## plot covariance and correlation functions
  # first get the true covariance an correlation functions
  spatialCovFun = function(x) {0.5 * stationary.cov(x, theta=thetas[1], Covariance="Matern", distMat=x, smoothness=nu) + 
      0.5 * stationary.cov(x, theta=thetas[2], Covariance="Matern", smoothness=nu, distMat=x)}
  mixtureCovFun = function(x) {
    out = spatialCovFun(x)
    out[x == 0] = 1 + sigma2
    out
  }
  mixtureCorFun = function(x) { mixtureCovFun(x) * (1 / (1 + sigma2)) }
  
  # first to transform all the hyperparameter samples to their relevant values
  # 1: error precision
  # 2: log effective range
  # 3: log spatial variance
  # 4-(3 + nLayer - 1): multivariateLogit alpha
  nuggetVarVals = 1/out[,1]
  effectiveRangeVals = out[,2]
  varVals = out[,3]^2
  
  # compute the covariance function for many different hyperparameter samples
  out = covarianceDistributionSPDE(effectiveRangeVals, varVals, nuggetVarVals, mesh, xRangeDat=xRangeDat, yRangeDat=yRangeDat)
  d = out$d
  sortI = sort(d, index.return=TRUE)$ix
  d = d[sortI]
  covMean = out$cov[sortI]
  upperCov=out$upperCov[,sortI]
  lowerCov=out$lowerCov[,sortI]
  covMat=out$covMat[sortI,]
  corMean = out$cor[sortI]
  upperCor=out$upperCor[,sortI]
  lowerCor=out$lowerCor[,sortI]
  corMat=out$corMat[sortI,]
  covInfo = list(d=d, covMean=covMean, upperCov=upperCov, lowerCov=lowerCov, covMat=covMat, 
                 corMean=corMean, upperCor=upperCor, lowerCor=lowerCor, corMat=corMat)
  
  # plot the covariance function
  yRange = range(c(covMean, lowerCov, upperCov, mixtureCovFun(d)))
  pdf(file=paste0("Figures/mixtureSPDECov", plotNameRoot, ".pdf"), width=5, height=5)
  plot(d, covMean, type="l", main="Posterior of covariance function", xlab="Distance", ylab="Covariance", 
       ylim=yRange)
  lines(d, lowerCov[1,], lty=2)
  lines(d, upperCov[1,], lty=2)
  lines(d, lowerCov[2,], lty=3)
  lines(d, upperCov[2,], lty=3)
  lines(d, mixtureCovFun(d), col="green")
  legend("topright", c("Truth", "Estimate", "80% CI", "95% CI"), lty=c(1, 1, 2, 3), col=c("green", "black", "black", "black"))
  dev.off()
  
  pdf(file=paste0("Figures/mixtureSPDECor", plotNameRoot, ".pdf"), width=5, height=5)
  yRange = range(c(corMean, lowerCor, upperCor, mixtureCorFun(d)))
  plot(d, corMean, type="l", main="Posterior of correlation function", xlab="Distance", ylab="Covariance", 
       ylim=yRange)
  lines(d, lowerCor[1,], lty=2)
  lines(d, upperCor[1,], lty=2)
  lines(d, lowerCor[2,], lty=3)
  lines(d, upperCor[2,], lty=3)
  lines(d, mixtureCorFun(d), col="green")
  legend("topright", c("Truth", "Estimate", "80% CI", "95% CI"), lty=c(1, 1, 2, 3), col=c("green", "black", "black", "black"))
  dev.off()
  
  # get scoring rules
  testIndices = (length(preds) - length(ysTest) + 1):length(preds)
  leftOutIndices = (length(preds) - length(ysTest) + 1):(length(preds) - length(ysTest) + length(simulationData$zTest[,1]))
  gridIndices = (length(preds) - length(ysTest) + length(simulationData$zTest[,1]) + 1):length(preds)
  leftOutIndicesTest = match(leftOutIndices, testIndices)
  gridIndicesTest = match(gridIndices, testIndices)
  
  # first calculate scoring rules at grid points
  truth = ysTest[gridIndicesTest]
  est = preds[gridIndices]
  vars = predSDs[gridIndices]^2
  lower = fit$lower[gridIndices]
  upper = fit$upper[gridIndices]
  
  # compute nearest neighbor distances and scores as a function of them
  gridPts = predPts[gridIndices,]
  distMat = rdist(coords, gridPts)
  nndists = apply(distMat, 2, function(x) {min(x[x != 0])})
  print("Binned grid scores:")
  gridScoringRules = getScores(truth, est, vars, lower, upper, distances=nndists, breaks=distanceBreaks)
  print(gridScoringRules$pooledResults)
  print(gridScoringRules$binnedResults)
  
  # now calculate square rules at left out points
  truth = ysTest[leftOutIndicesTest]
  est = preds[leftOutIndices]
  vars = predSDs[leftOutIndices]^2
  lower = fit$lower[leftOutIndices]
  upper = fit$upper[leftOutIndices]
  leftOutScoringRules = getScores(truth, est, vars, lower, upper)
  leftOutScoringRules = data.frame(c(leftOutScoringRules, Time=time[3]))
  print("Binned left out scores:")
  print(leftOutScoringRules)
  
  scoringRules = list(gridScoringRules=gridScoringRules, leftOutScoringRules=leftOutScoringRules)
  # 
  # # get scoring rules
  # truth = ysTest
  # est = preds[testIndices]
  # vars = predSDs[testIndices]^2
  # lower = fit$lower[testIndices]
  # upper = fit$upper[testIndices]
  # 
  # # compute nearest neighbor distances and scores as a function of them
  # testPts = predPts[testIndices,]
  # distMat = rdist(coords, testPts)
  # nndists = apply(distMat, 2, function(x) {min(x[x != 0])})
  # print("Binned scores:")
  # scoringRules = getScores(truth, est, vars, lower, upper, distances=nndists, breaks=distanceBreaks)
  # scoringRules$pooledResults = data.frame(c(scoringRules$pooledResults, Time=time[3]))
  # print(scoringRules$binnedResults)
  
  
  
  # print("Grid scores:")
  # print(getScores(truth[gridTestI], est[gridTestI], vars[gridTestI], lower[gridTestI], upper[gridTestI], 
  #                 distances=nndists[gridTestI], breaks=distanceBreaks)$binnedResults)
  # 
  # # plot scores as a function of distance
  # distanceScores = getScores(truth[gridTestI], est[gridTestI], vars[gridTestI], lower[gridTestI], upper[gridTestI], 
  #                            distances=nndists[gridTestI], breaks=distanceBreaks)$binnedResults
  
  # pdf(file=paste0("Figures/mixtureSPDEScoreBias", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(distanceScores$NNDist, distanceScores$Bias, pch=19, col="blue", main="Bias", ylab="Bias", xlab="Nearest neighbor distance (km)")
  # abline(h=0, lty=2)
  # dev.off()
  # 
  # pdf(file=paste0("Figures/mixtureSPDEScoreVar", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(distanceScores$NNDist, distanceScores$Var, pch=19, col="blue", main="Variance", ylab="Variance", xlab="Nearest neighbor distance (km)", 
  #      ylim=c(0, max(distanceScores$Var)))
  # abline(h=0, lty=2)
  # dev.off()
  # 
  # pdf(file=paste0("Figures/mixtureSPDEScoreMSE", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(distanceScores$NNDist, distanceScores$MSE, pch=19, col="blue", main="MSE", ylab="MSE", xlab="Nearest neighbor distance (km)", 
  #      ylim=c(0, max(distanceScores$MSE)))
  # abline(h=0, lty=2)
  # dev.off()
  # 
  # pdf(file=paste0("Figures/mixtureSPDEScoreRMSE", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(distanceScores$NNDist, distanceScores$RMSE, pch=19, col="blue", main="RMSE", ylab="RMSE", xlab="Nearest neighbor distance (km)", 
  #      ylim=c(0, max(distanceScores$RMSE)))
  # abline(h=0, lty=2)
  # dev.off()
  # 
  # pdf(file=paste0("Figures/mixtureSPDEScoreCRPS", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(distanceScores$NNDist, distanceScores$CRPS, pch=19, col="blue", main="CRPS", ylab="CRPS", xlab="Nearest neighbor distance (km)", 
  #      ylim=c(0, max(distanceScores$CRPS)))
  # abline(h=0, lty=2)
  # dev.off()
  # 
  # pdf(file=paste0("Figures/mixtureSPDEScoreCvg", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(distanceScores$NNDist, distanceScores$Coverage, pch=19, col="blue", main="80% Coverage", ylab="80% Coverage", xlab="Nearest neighbor distance (km)", 
  #      ylim=c(0, max(distanceScores$Coverage)))
  # abline(h=0, lty=2)
  # dev.off()
  # 
  # pdf(file=paste0("Figures/mixtureSPDEScoreWidth", plotNameRoot, ".pdf"), width=5, height=5)
  # plot(distanceScores$NNDist, distanceScores$Width, pch=19, col="blue", main="Width", ylab="Width", xlab="Nearest neighbor distance (km)", 
  #      ylim=c(0, max(distanceScores$Width)))
  # abline(h=0, lty=2)
  # dev.off()
  
  if(!useKenya) {
    # print("Pooled scores:")
    # print(data.frame(c(getScores(truth, est, vars, lower, upper), Time=time[3])))
  } else {
    print("Pooled scores:")
    print(data.frame(c(getScores(truth, est, vars, lower, upper), Time=time[3])))
    print("Overall scores:")
    print(data.frame(c(getScores(truth[overallTestI], est[overallTestI], vars[overallTestI], lower[overallTestI], upper[overallTestI]), Time=time[3])))
    print("Rural scores:")
    print(data.frame(c(getScores(truth[ruralTestI], est[ruralTestI], vars[ruralTestI], lower[ruralTestI], upper[ruralTestI]), Time=time[3])))
    print("Urban scores:")
    print(data.frame(c(getScores(truth[urbanTestI], est[urbanTestI], vars[urbanTestI], lower[urbanTestI], upper[urbanTestI]), Time=time[3])))
    # print("Grid scores:")
    # print(data.frame(c(getScores(truth[gridTestI], est[gridTestI], vars[gridTestI], lower[gridTestI], upper[gridTestI]), Time=time[3])))
  }
  
  # get aggregated predictions
  # A = t(sapply(1:(mx*my), getARow))
  # A = sweep(A, 1, rowSums(A), "/")
  # mx = 100
  # my = 100
  # predPts = make.surface.grid(list(x=seq(xRangeDat[1], xRangeDat[2], l=mx), y=seq(yRangeDat[1], yRangeDat[2], l=my)))
  gridCoords = cbind(simulationData$xGrid, simulationData$yGrid)
  # testIndices = (length(preds) - length(ysTest) + 1):length(preds)
  # gridIndices = (length(preds) - length(simulationData$xGrid) + 1):length(preds)
  # gridIndicesEst = (length(est) - length(simulationData$xGrid) + 1):length(est)
  
  A = makeNumericalIntegralMat(gridCoords, mx=3, my=3)
  aggregatedPreds = A %*% preds[gridIndices]
  
  truth = A %*% ysTest[gridIndicesTest]
  est = A %*% preds[gridIndices]
  aggregatedPredMat = A %*% fit$predMat[gridIndices,]
  vars = apply(aggregatedPredMat, 1, var)
  sds = apply(aggregatedPredMat, 1, sd)
  lower = apply(aggregatedPredMat, 1, function(x) {quantile(x, probs=.1)})
  upper = apply(aggregatedPredMat, 1, function(x) {quantile(x, probs=.9)})
  predictionMatrix = data.frame(Truth=truth, Est=est, SDs=sds, Lower=lower, Upper=upper)
  print("Aggregated prediction summary statistics:")
  print(predictionMatrix)
  
  print("Pooled aggregated scores:")
  pooledAggregatedScores = getScores(truth, est, vars, lower, upper)
  print(pooledAggregatedScores)
  print("Left out region aggregated scores:")
  leftOutAggregatedScores = getScores(truth[5], est[5], vars[5], lower[5], upper[5])
  leftOutAggregatedScores$Var = sds[5]
  names(leftOutAggregatedScores)[2] = "Predictive.SD"
  print(leftOutAggregatedScores)
  print("Included regions aggregated scores:")
  includedAggregatedScores = getScores(truth[-5], est[-5], vars[-5], lower[-5], upper[-5])
  print(includedAggregatedScores)
  
  aggregatedScoringRules = list(pooledAggregatedScores=pooledAggregatedScores, leftOutAggregatedScores=leftOutAggregatedScores, 
                                includedAggregatedScores=includedAggregatedScores)
  
  analysisTime = proc.time()[3] - startTime
  
  fit$mod = NULL
  if(!gscratch)
    save(scoringRules, fit, covInfo, predictionMatrix, aggregatedScoringRules, analysisTime, file=paste0("savedOutput/simulations/mixtureSPDE", plotNameRoot, ".RData"))
  else
    save(scoringRules, fit, covInfo, predictionMatrix, aggregatedScoringRules, analysisTime, file=paste0("/work/johnpai/mixtureSPDE", plotNameRoot, ".RData"))
  
  invisible(list(scoringRules, fit, covInfo, predictionMatrix, aggregatedScoringRules, analysisTime))
}

# runs the testSPDEModelMixture function for multiple realizations, saves results
testSPDEModelMixtureMultiple = function(seed=1, nSamples=100, n=900, nu=1, sigma2=0.1^2, 
                                        useKenya=FALSE, assumeMeanZero=TRUE, urbanOverSamplefrac=0, 
                                        gscratch=TRUE, loadResults=FALSE, ...) {
  # set random seeds for each simulation
  set.seed(seed)
  allSeeds = sample(1:1000000, nSamples, replace = FALSE)
  
  # call testSPDEModelMixture for each simulation requested
  temp = function(i) {
    print(paste0("Beginning simulation ", i, "/", nSamples))
    thisPlotNameRoot = paste0("sim", i)
    do.call("testSPDEModelMixture", c(list(seed = allSeeds[i], n=n, nu=nu, sigma2=sigma2, gscratch=gscratch, 
                                             useKenya=useKenya, urbanOverSamplefrac=urbanOverSamplefrac, assumeMeanZero=assumeMeanZero, plotNameRoot=thisPlotNameRoot), list(...)))
  }
  if(!loadResults)
    sapply(1:nSamples, temp)
  
  # set plotNameRoot
  plotNameRoot = paste0("_n", n, "_nu", nu, "_nugV", round(sigma2, 2), "_Kenya", useKenya, 
                        "_noInt", assumeMeanZero, "_urbOversamp", round(urbanOverSamplefrac, 4))
  
  # save(scoringRules, fit, covInfo, predictionMatrix, aggregatedScoringRules, file=paste0("savedOutput/simulations/mixtureSPDE", plotNameRoot, ".RData"))
  allScoringRulesGrid = list()
  allScoringRulesLeftOut = list()
  allFits = list()
  allCovInfo = list()
  allPredictionMatrices = list()
  allAggregatedScoringRules = list()
  allAnalysisTimes = c()
  for(i in 1:nSamples) {
    if(!gscratch)
      out = load(paste0("savedOutput/simulations/mixtureSPDEsim", i, plotNameRoot, ".RData"))
    else
      out = load(paste0("/work/johnpai/mixtureSPDEsim", i, plotNameRoot, ".RData"))
    allScoringRulesGrid = c(allScoringRulesGrid, list(scoringRules$gridScoringRules))
    allScoringRulesLeftOut = c(allScoringRulesLeftOut, list(scoringRules$leftOutScoringRules))
    allFits = c(allFits, list(fit))
    allCovInfo = c(allCovInfo, list(covInfo))
    allPredictionMatrices = c(allPredictionMatrices, list(predictionMatrix))
    allAggregatedScoringRules = c(allAggregatedScoringRules, list(aggregatedScoringRules))
    allAnalysisTimes = c(allAnalysisTimes, analysisTime)
  }
  
  ##### average results from each simulation
  # pointwise scoring rules
  allPooledScoringRulesGrid = do.call("rbind", lapply(allScoringRulesGrid, function(x) {x$pooledResults}))
  allBinnedScoringRulesGrid = lapply(allScoringRulesGrid, function(x) {x$binnedResults})
  binnedScoringRulesGrid = averageBinnedScores(allBinnedScoringRulesGrid)
  ns = binnedScoringRulesGrid[,2]
  pooledScoringRulesGrid = apply(binnedScoringRulesGrid, 2, function(x) {sum(x * (ns / sum(ns)))})
  pooledScoringRulesGrid = as.data.frame(matrix(pooledScoringRulesGrid, nrow=1))
  names(pooledScoringRulesGrid) = names(binnedScoringRulesGrid)
  
  fullPooledScoringRulesLeftOut = do.call("rbind", allScoringRulesLeftOut)
  pooledScoringRulesLeftOut = colMeans(fullPooledScoringRulesLeftOut)
  
  # covInfo
  # covInfo = list(d=d, covMean=covMean, upperCov=upperCov, lowerCov=lowerCov, covMat=covMat, 
  #                corMean=corMean, upperCor=upperCor, lowerCor=lowerCor, corMat=corMat)
  d = allCovInfo[[1]]$d
  covMean = rowMeans(do.call("cbind", lapply(allCovInfo, function(x) {x$covMean})))
  upperCov = rowMeans(do.call("cbind", lapply(allCovInfo, function(x) {x$upperCov[1,]})))
  lowerCov = rowMeans(do.call("cbind", lapply(allCovInfo, function(x) {x$lowerCov[1,]})))
  corMean = rowMeans(do.call("cbind", lapply(allCovInfo, function(x) {x$corMean})))
  upperCor = rowMeans(do.call("cbind", lapply(allCovInfo, function(x) {x$upperCor[1,]})))
  lowerCor = rowMeans(do.call("cbind", lapply(allCovInfo, function(x) {x$lowerCor[1,]})))
  
  # aggregated scoring rules
  # aggregatedScoringRules = list(pooledAggregatedScores=pooledAggregatedScores, leftOutAggregatedScores=leftOutAggregatedScores, 
  #                               includedAggregatedScores=includedAggregatedScores)
  # predictionMatrix = data.frame(Truth=truth, Est=est, SDs=sds, Lower=lower, Upper=upper)
  fullPredictionMatrix = do.call("rbind", allPredictionMatrices)
  leftOutIndices = seq(5, nrow(fullPredictionMatrix), by=9)
  leftOutPredictionMatrix = fullPredictionMatrix[leftOutIndices, ]
  leftInPredictionMatrix = fullPredictionMatrix[-leftOutIndices, ]
  leftOutScores = getScores(leftOutPredictionMatrix[,1], leftOutPredictionMatrix[,2], leftOutPredictionMatrix[,3]^2)
  leftInScores = getScores(leftInPredictionMatrix[,1], leftInPredictionMatrix[,2], leftInPredictionMatrix[,3]^2)
  aggregatedScores = getScores(fullPredictionMatrix[,1], fullPredictionMatrix[,2], fullPredictionMatrix[,3]^2)
  
  ##### Save results
  if(!gscratch) {
    save(#allScoringRules, 
         allFits, 
         allCovInfo, 
         allPredictionMatrices, 
         allAggregatedScoringRules, 
         allAnalysisTimes, 
         binnedScoringRulesGrid, 
         pooledScoringRulesGrid, 
         fullPooledScoringRulesLeftOut, 
         pooledScoringRulesLeftOut, 
         covMean, 
         upperCov, 
         lowerCov, 
         corMean, 
         upperCor, 
         lowerCor, 
         fullPredictionMatrix, 
         leftOutPredictionMatrix, 
         leftInPredictionMatrix, 
         leftOutScores, 
         leftInScores, 
         aggregatedScores, 
         file=paste0("savedOutput/simulations/mixtureSPDEAll_nsim", nSamples, plotNameRoot, ".RData"))
  } else {
    save(#allScoringRules, 
         allFits, 
         allCovInfo, 
         allPredictionMatrices, 
         allAggregatedScoringRules, 
         binnedScoringRulesGrid, 
         pooledScoringRulesGrid, 
         fullPooledScoringRulesLeftOut, 
         pooledScoringRulesLeftOut, 
         covMean, 
         upperCov, 
         lowerCov, 
         corMean, 
         upperCor, 
         lowerCor, 
         fullPredictionMatrix, 
         leftOutPredictionMatrix, 
         leftInPredictionMatrix, 
         leftOutScores, 
         leftInScores, 
         aggregatedScores, 
         file=paste0("/work/johnpai/simulations/mixtureSPDEAll_nsim", nSamples, plotNameRoot, ".RData"))
  }
}

# test how close we can get to the spatial correlation function:
testLKINLACorrelationApproximation = function(seed=1, nLayer=3, NP=200, 
                                        nBuffer=5, normalize=TRUE, fastNormalize=TRUE, NC=13, 
                                        latInfo=NULL, thetas=c(.1, .4), 
                                        initialEffectiveRange=1, initialAlphas=rep(1/nLayer, nLayer-1), 
                                        nu=.5) {
  set.seed(seed)
  
  # set true parameter values
  rho = 1
  effectiveRange = thetas * sqrt(8)
  sigma2 = sqrt(.1)
  # spatialCorFun = function(x) {0.4 * Exp.cov(x, theta=thetas[1], distMat=x) + 0.6 * Exp.cov(x, theta=thetas[2], distMat=x)}
  spatialCorFun = function(x) {0.5 * stationary.cov(x, theta=thetas[1], Covariance="Matern", distMat=x, smoothness=nu) + 
      0.5 * stationary.cov(x, theta=thetas[2], Covariance="Matern", distMat=x, smoothness=nu)}
  
  # construct the lattice
  xRangeDat = c(-1, 1)
  yRangeDat = c(-1, 1)
  if(is.null(latInfo))
    latInfo = makeLatGrids(xRangeDat, yRangeDat, NC, nBuffer, nLayer)
  
  # set initial parameter values
  initialKappa = (sqrt(8) * latInfo[[1]]$latWidth /initialEffectiveRange)
  
  xlim <- latInfo[[1]]$xRangeDat
  ux <- seq(xlim[1], xlim[2], , NP)
  ylim <- latInfo[[1]]$yRangeDat
  uy <- seq(ylim[1], ylim[2], , NP)
  center <- rbind(c(ux[NP/2], uy[NP/2]))
  
  # precompute relevant matrices
  Qprecomputations = precomputationsQ2(latInfo)
  Acenter = makeA(center, latInfo)
  Ax = makeA(cbind(ux, uy[NP/2]), latInfo)
  Ay = makeA(cbind(ux[NP/2], uy), latInfo)
  
  testFun = function(parameters) {
    kappa = exp(parameters[1])
    logitAlphas = parameters[2:nLayer]
    alphas = multivariateExpit(logitAlphas)
    alphas = c(alphas, 1 - sum(alphas))
    test = getLKInlaCovarianceFun(kappa, 1, 0, alphas, latticeInfo = latInfo, 
                                  precomputedMatrices=Qprecomputations, precomputedAcenter=Acenter, precomputedAx=Ax, precomputedAy=Ay)
    
    ds = test[-1,1]
    mean((spatialCorFun(ds) - test[-1,2])^2 * (1/ds)) # inverse distance weighted mean squared error of the correlation function fit
    # mean((spatialCorFun(ds) - test[-1,2])^2)
  }
  
  print("Beginning optimization...")
  temp = optim(c(log(initialKappa), multivariateLogit(initialAlphas)), testFun)
  outPar = temp$par
  outkappa = exp(outPar[1])
  outlogitAlphas = outPar[2:nLayer]
  outalphas = multivariateExpit(outlogitAlphas)
  outalphas = c(outalphas, 1 - sum(outalphas))
  outEffectiveRange = (1/outkappa) * sqrt(8) * latInfo[[1]]$latWidth
  
  trueEffectiveRangeOverall = getTrueLKEffectiveRange(nLayer, NP, sigma2=0, rho=1, nBuffer, normalize, fastNormalize, 
                                                      NC, latInfo, outEffectiveRange, outalphas)
  trueEffectiveRanges = rep(0, nLayer)
  for(i in 1:nLayer) {
    theseAlphas = rep(0, nLayer)
    theseAlphas[i] = 1
    trueEffectiveRanges[i] = getTrueLKEffectiveRange(nLayer, NP, sigma2=0, rho=1, nBuffer, normalize, fastNormalize, 
                                                     NC, latInfo, outEffectiveRange, theseAlphas)
    
    print(paste0("Layer ", i, " weight: ", outalphas[i], ", true effective range: ", trueEffectiveRanges[i]))
  }
  
  test = getLKInlaCovarianceFun(outkappa, 1, 0, outalphas, latticeInfo = latInfo)
  ds = test[,1]
  pdf(paste0("Figures/approxCorrelation_nLayer", nLayer, "_NC", NC, "_nBuffer", nBuffer, "_nu", nu, ".pdf"), width=5, height=5)
  plot(ds, spatialCorFun(ds), type="l", col="green", ylim=c(0,1), ylab="Correlation", main="Correlation", xlab="Distance")
  lines(ds, test[,2], col="blue")
  legend("topright", c("Truth", "Approximate"), col=c("green", "blue"), lty=1)
  dev.off()
}

# get correlation function from the Matern family with smoothness given by nu that is best approximation to
# the exponential mixture correlation function
getCorrelationApproximation = function(thetas=c(.1, .4), nu=0.5) {
  NP = 200
  # set true parameter values
  rho = 1
  effectiveRange = thetas * sqrt(8)
  sigma2 = sqrt(.1)
  # spatialCorFun = function(x) {0.4 * Exp.cov(x, theta=thetas[1], distMat=x) + 0.6 * Exp.cov(x, theta=thetas[2], distMat=x)}
  spatialCorFun = function(x) {0.5 * stationary.cov(x, theta=thetas[1], Covariance="Matern", distMat=x, smoothness=nu) + 
      0.5 * stationary.cov(x, theta=thetas[2], Covariance="Matern", distMat=x, smoothness=nu)}
  
  # construct the lattice
  xRangeDat = c(-1, 1)
  yRangeDat = c(-1, 1)
  ds = seq(0, 1, l=NP)[-1]
  
  testFun = function(parameters) {
    theta = exp(parameters[1])
    
    mean((spatialCorFun(ds) - Matern(ds, range=theta, smoothness=1))^2 / ds) # inverse distance weighted mean squared error of the correlation function fit
    # mean((spatialCorFun(ds) - test[-1,2])^2)
  }
  
  print("Beginning optimization...")
  temp = optim(log(.2), testFun)
  theta = exp(temp$par)
  
  pdf(paste0("Figures/approxCorrelation_MaternNu1.pdf"), width=5, height=5)
  plot(ds, spatialCorFun(ds), type="l", col="green", ylim=c(0,1), ylab="Correlation", main="Correlation", xlab="Distance")
  lines(ds, Matern(ds, range=theta, smoothness=1), col="blue")
  legend("topright", c("Truth", "Approximate"), col=c("green", "blue"), lty=1)
  dev.off()
  theta # 0.1256181
}

maternMixture = function(x, thetas=c(.1, .4), nu=0.5, alphas=c(0.5, 0.5)) {
  alphas[1] * Matern(x, range=thetas[1], smoothness=nu) + 
    alphas[2] * Matern(x, range=thetas[2], smoothness=nu)
}
maternMixtureCor = function(x1, x2=NULL, thetas=c(.1, .4), nu=0.5, alphas=c(0.5, 0.5), distMat=NA) {
  tempFun = function(d) {maternMixture(d, thetas, nu, alphas)}
  stationary.cov(x1, x2, Covariance=tempFun, distMat=distMat)
}
maternApproxCor = function(x1, x2=NULL, theta=0.1256181, smoothness=1, distMat=NA) {
  approxCorrelation=function(d) {Matern(d, range=theta, smoothness=smoothness)}
  stationary.cov(x1, x2, Covariance=approxCorrelation, distMat=distMat)
}

# nx, ny: resolution of prediction grid in the x and y directions
# mx, my: the number of aggregation regions (rectangles) in the x and y directions
testCorrelationApproximation = function(trueCorrelation = maternMixtureCor, 
                                        approxCorrelation = maternApproxCor, 
                                        nx=75, ny=75, seed=1, n=900, 
                                        mx=3, my=3) {
  set.seed(seed)
  
  # set true parameter values
  rho = 1
  # sigma2 = sqrt(.1)
  sigma2 = 0
  
  # generate data set
  nTest = 100
  simulationData = getSimulationDataSetsGivenCovarianceTest(trueCorrelation, nTotal=n, nTest=nTest, marginalVar=rho, errorVar=0, 
                                                        nDataSets=2, plotNameRoot=paste0("(0.5*Exp(.1) + 0.5*Exp(.4))"), fileNameRoot="mix", 
                                                        saveDataSetPlot=FALSE, doPredGrid=TRUE)
  coords = cbind(simulationData$xTrain[,1], simulationData$yTrain[,1])
  testCoords = cbind(simulationData$xTest[,1], simulationData$yTest[,1])
  ys = simulationData$zTrain[,1]
  gridCoords = cbind(simulationData$xGrid, simulationData$yGrid)
  gridValues = simulationData$zGrid[,1]
  
  # generate lattice and simulate observations
  # coords = matrix(runif(2*n), ncol=2)
  xRangeDat = c(-1, 1)
  yRangeDat = c(-1, 1)
  
  # plot the observations
  pdf(file="Figures/mixtureMaternObservations.pdf", width=5, height=5)
  par(mfrow=c(1,1))
  quilt.plot(coords, ys)
  dev.off()
  
  # make prediction coordinates on a grid, and add testing points
  xRange=c(-1,1)
  yRange=c(-1,1)
  predPts = rbind(gridCoords, testCoords)
  ysGrid = gridValues
  ysTest = simulationData$zTest[,1]
  
  # construct aggregation matrix for predictions by testing which prediction locations 
  # are in which aggregation regions
  xRegionGrid = seq(-1, 1, l=mx + 1)[-1]
  yRegionGrid = seq(-1, 1, l=my + 1)[-1]
  xRegion = function(x) {
    match(TRUE, x <= xRegionGrid)
  }
  yRegion = function(y) {
    match(TRUE, y <= yRegionGrid)
  }
  xRegionI = sapply(gridCoords[,1], xRegion)
  yRegionI = sapply(gridCoords[,2], yRegion)
  regionI = (yRegionI-1)*mx + xRegionI
  getARow = function(ai) {regionI == ai}
  
  # A = t(sapply(1:(mx*my), getARow))
  # A = sweep(A, 1, rowSums(A), "/")
  A = makeNumericalIntegralMat(gridCoords, mx=mx, my=my)
  
  # generate the true aggregated values of the field
  # NOTE: this includes any nugget effect. Is that a good idea??
  gridValuesAggregated = A %*% gridValues
  
  # fit the model true and approximate models
  time = system.time(fit <- GPpreds(coords, ys, predCoords=predPts, 
                                    xObs=NULL, xPred=NULL, cov.fun=trueCorrelation, A=A))
  predsGrid=fit$preds[1:ncol(A)]
  predSDsGrid=fit$sigmas[1:ncol(A)]
  predsTest=fit$preds[-(1:ncol(A))]
  predSDsTest=fit$sigmas[-(1:ncol(A))]
  predsAggregated=fit$predsAggregated
  predSDsAggregated=fit$sigmasAggregated
  
  time = system.time(fit <- GPpreds(coords, ys, predCoords=predPts, 
                                    xObs=NULL, xPred=NULL, cov.fun=approxCorrelation, A=A))
  predsApproxGrid=fit$preds[1:ncol(A)]
  predSDsApproxGrid=fit$sigmas[1:ncol(A)]
  predsApproxTest=fit$preds[-(1:ncol(A))]
  predSDsApproxTest=fit$sigmas[-(1:ncol(A))]
  predsApproxAggregated=fit$predsAggregated
  predSDsApproxAggregated=fit$sigmasAggregated
  
  # print results
  print(data.frame(list(truth=gridValuesAggregated, truePredsAggregated=predsAggregated, approxPredsAggregated=predsApproxAggregated, trueSDsAggregated=predSDsAggregated, approxSDsAggregated=predSDsApproxAggregated)))
  print(data.frame(list(trueMean=mean(ysTest), truePredsMean=mean(predsTest), approxPredsMean=mean(predsApproxTest), 
                        trueSDsMean=mean(predSDsTest), approxSDsMean=mean(predSDsApproxTest), 
                        trueMSE=mean((ysTest-predsTest)^2), approxMSE=mean((ysTest-predsApproxTest)^2))))
  
  # Calculate all scoring rules and print them
  aggregatedScoresTrue = getScores(gridValuesAggregated, predsAggregated, predSDsAggregated)
  aggregatedScoresApprox = getScores(gridValuesAggregated, predsApproxAggregated, predSDsApproxAggregated)
  scoresTrue = getScores(ysTest, predsTest, predSDsTest)
  scoresApprox = getScores(ysTest, predsApproxTest, predSDsApproxTest)
  
  print(aggregatedScoresTrue)
  print(aggregatedScoresApprox)
  print(scoresTrue)
  print(scoresApprox)
  
  # function for determining if points are in correct range
  inRange = function(pts, rangeShrink=0) {
    inX = (rangeShrink < pts[,1]) & (pts[,1] < 1-rangeShrink)
    inY = (rangeShrink < pts[,2]) & (pts[,2] < 1-rangeShrink)
    inX & inY
  }
  
  # show predictive surface, SD, and data
  browser()
  
  pdf(file="Figures/mixtureMaternPreds.pdf", width=12, height=8)
  par(mfrow=c(2,3), mar=c(5.1, 4.1, 4.1, 6))
  
  colRangeDat = range(c(ys, predsGrid, predsApproxGrid, gridValues))
  colRangeSD = range(c(predSDsGrid[inRange(gridCoords, rangeShrink=.01)], predSDsApproxGrid[inRange(gridCoords, rangeShrink=.01)]))
  quilt.plot(gridCoords[,1], gridCoords[,2], gridValues, main="True Process", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
  quilt.plot(gridCoords[,1], gridCoords[,2], predsGrid, main="Prediction Mean (True Covariance)", zlim=colRangeDat, 
             xlim=xRangeDat, ylim=yRangeDat)
  quilt.plot(gridCoords[,1], gridCoords[,2], predsApproxGrid, main="Prediction Mean (Approx. Covariance)", zlim=colRangeDat, 
             xlim=xRangeDat, ylim=yRangeDat)
  
  quilt.plot(coords, ys, main="Observations", zlim=colRangeDat, xlim=xRangeDat, ylim=yRangeDat)
  quilt.plot(gridCoords[,1], gridCoords[,2], predSDsGrid, main="Prediction SD (True Covariance)",
             xlim=xRangeDat, ylim=yRangeDat, zlim=colRangeSD)
  quilt.plot(gridCoords[,1], gridCoords[,2], predSDsApproxGrid, main="Prediction SD (Approx. Covariance)",
             xlim=xRangeDat, ylim=yRangeDat, zlim=colRangeSD)
  dev.off()
  
  my_line <- function(x,y,...){
    abline(a = 0,b = 1,...)
    points(x,y,..., col="blue")
  }
  
  pdf(file="Figures/mixtureMaternLeftOutPairGrid.pdf", width=6, height=6)
  thisRange = range(c(cbind(Truth=gridValues, BLUP=predsGrid, ApproxCov=predsApproxGrid)))
  pairs(cbind(Truth=gridValues, BLUP=predsGrid, ApproxCov=predsApproxGrid), pch=19, cex=.1, xlim=thisRange, ylim=thisRange, 
        upper.panel=my_line, lower.panel=my_line, main="Full Domain Grid Pair Plot")
  dev.off()
  
  pdf(file="Figures/mixtureMaternLeftOutPairTest.pdf", width=6, height=6)
  thisRange = range(c(cbind(Truth=ysTest, BLUP=predsTest, ApproxCov=predsApproxTest)))
  pairs(cbind(Truth=ysTest, BLUP=predsTest, ApproxCov=predsApproxTest), pch=19, cex=.4, xlim=thisRange, ylim=thisRange, 
        upper.panel=my_line, lower.panel=my_line, main="Left Out Area Pair Plot")
  dev.off()
  
  # pdf(file="Figures/mixtureMaternLeftOutResidualsGrid.pdf", width=5, height=5)
  # testIndices = (length(predsGrid) - length(ysGrid) + 1):length(predsGrid)
  # plot(predsGrid[testIndices], ysGrid-predsGrid[testIndices], pch=19, cex=.5, col="blue", main="Residuals versus fitted", 
  #      ylab="Residuals", xlab="Fitted")
  # abline(h=0, lty=2)
  # dev.off()
  # 
  # pdf(file="Figures/mixtureMaternLeftOutResidualsTest.pdf", width=5, height=5)
  # testIndices = (length(predsGrid) - length(ysGrid) + 1):length(predsGrid)
  # plot(predsGrid[testIndices], ysGrid-predsGrid[testIndices], pch=19, cex=.5, col="blue", main="Residuals versus fitted", 
  #      ylab="Residuals", xlab="Fitted")
  # abline(h=0, lty=2)
  # dev.off()
}

# The same as getSimulationDataSetsGivenCovariance, but test observations all come from a left out region in the center
getSimulationDataSetsGivenCovarianceTest = function(corFun, nTotal=900, nTest=round(nTotal / 9), marginalVar=1, errorVar=sqrt(.1), nDataSets=100, 
                                                    printEvery=10, saveDataSetPlot=TRUE, fileNameRoot="", plotNameRoot="", 
                                                    doPredGrid=FALSE, nTestGrid=70^2) {
  
  # set the spatial domain
  xRange = c(-1, 1)
  yRange = c(-1, 1)
  
  # set the region grid
  mx=3
  my=3
  xRegionGrid = seq(-1, 1, l=mx + 1)[-1]
  yRegionGrid = seq(-1, 1, l=my + 1)[-1]
  centerxLeftI = floor((mx - 1) / 2)
  
  # if generating prediction grid, make that grid here
  if(doPredGrid) {
    nx = round(sqrt(nTestGrid))
    ny = nx
    if(nx * ny != nTestGrid)
      stop("nTestGrid must be a square number if we are constructing a prediction grid")
    xValuesGrid = seq(-1, 1, l=nx)
    yValuesGrid = xValuesGrid
    glist = make.surface.grid(list(x=xValuesGrid, y=yValuesGrid))
    xValuesGrid = glist[,1]
    yValuesGrid = glist[,2]
  } else {
    nx = 0
    ny = 0
    xValuesGrid = c()
    yValuesGrid = c()
  }
  
  # simulate observation spatial locations
  # xValues = matrix(runif(nTotal * nDataSets, xRange[1], xRange[2]), ncol=nDataSets)
  # yValues = matrix(runif(nTotal * nDataSets, yRange[1], yRange[2]), ncol=nDataSets)
  nTrain = nTotal - nTest
  coordsTrain = runifsqMissingRectangle(nTrain * nDataSets)
  coordsTest = runifsq(nTest * nDataSets, c(-1/3, 1/3), c(-1/3, 1/3))
  xValuesTrain = matrix(coordsTrain[,1], ncol=nDataSets)
  yValuesTrain = matrix(coordsTrain[,2], ncol=nDataSets)
  xValuesTest = matrix(coordsTest[,1], ncol=nDataSets)
  yValuesTest = matrix(coordsTest[,2], ncol=nDataSets)
  xValues = rbind(xValuesTrain, xValuesTest)
  yValues = rbind(yValuesTrain, yValuesTest)
  
  # preallocate observation matrix, and pregenerate standard normal draws
  observations = matrix(nrow=nTotal+nx*ny, ncol=nDataSets)
  zsims = matrix(rnorm((nTotal + nx*ny) * nDataSets), ncol=nDataSets)
  
  # generate spatial component of observation values
  for(i in 1:nDataSets) {
    if(i %% printEvery == 0 || i == 1)
      print(paste0("Simulating data set ", i, "/", nDataSets))
    
    thisx = xValues[,i]
    thisy = yValues[,i]
    L = t(chol(corFun(cbind(c(thisx, xValuesGrid), c(thisy, yValuesGrid)))))
    observations[,i] = L %*% zsims[,i]
  }
  
  # scale by marginal standard deviation and add in error variance
  observations = observations * sqrt(marginalVar) + matrix(rnorm((nTotal + nx*ny) * nDataSets, sd=sqrt(errorVar)), ncol=nDataSets)
  
  # separate out test and train results
  trainI = 1:(nTotal - nTest)
  testI = (nTotal - nTest + 1):nTotal
  gridI = (nTotal + 1):(nTotal + nx * ny)
  if(nTotal - nTest != 0) {
    xTrain = xValues[trainI,]
    yTrain = yValues[trainI,]
    zTrain = observations[trainI,]
  } else {
    xTrain = c()
    yTrain = c()
    zTrain = c()
  }
  if(nTest != 0) {
    xTest = xValues[testI,]
    yTest = yValues[testI,]
    zTest = observations[testI,]
  }
  else {
    xTest = c()
    yTest = c()
    zTest = c()
  }
  # get grid results if necessary
  if(doPredGrid) {
    xGrid = xValuesGrid
    yGrid = yValuesGrid
    zGrid = observations[gridI,]
  } else {
    xGrid = NULL
    yGrid = NULL
    zGrid = NULL
  }
  
  # put relevant values into a list
  out = list(xTrain=xTrain, yTrain=yTrain, zTrain=zTrain, xTest=xTest, yTest=yTest, zTest=zTest, 
             xGrid=xGrid, yGrid=yGrid, zGrid=zGrid, xValuesGrid=xValuesGrid, yValuesGrid=yValuesGrid, nx=nx, ny=ny, 
             corFun=corFun, marginalVar=marginalVar, errorVar=errorVar, xRange=xRange, yRange=yRange)
  
  # plot the results
  plotExampleDataSets(out, saveDataSetPlot=saveDataSetPlot, plotNameRoot=plotNameRoot, fileNameRoot=fileNameRoot)
  
  # return the results
  out
}

# try out several bases, and determine the number of bases elements and their resolution
testBasisResolution = function() {
  latticeInfo = makeLatGridsKenya(nLayer=3, NC=28, nBuffer=5)
  print(paste0("NC=", 28, ", nLayer=", 3))
  print(paste0("lattice width=", min(sapply(latticeInfo, function(x) {x$latWidth}))))
  print(paste0("total basis functions=", sum(sapply(latticeInfo, function(x) {x$nx * x$ny}))))
  
  latticeInfo = makeLatGridsKenya(nLayer=3, NC=14, nBuffer=5)
  print(paste0("NC=", 14, ", nLayer=", 3))
  print(paste0("lattice width=", min(sapply(latticeInfo, function(x) {x$latWidth}))))
  print(paste0("total basis functions=", sum(sapply(latticeInfo, function(x) {x$nx * x$ny}))))
  
  latticeInfo = makeLatGridsKenya(nLayer=2, NC=54, nBuffer=5)
  print(paste0("NC=", 54, ", nLayer=", 2))
  print(paste0("lattice width=", min(sapply(latticeInfo, function(x) {x$latWidth}))))
  print(paste0("total basis functions=", sum(sapply(latticeInfo, function(x) {x$nx * x$ny}))))
  
  latticeInfo = makeLatGridsKenya(nLayer=1, NC=107, nBuffer=5)
  print(paste0("NC=", 107, ", nLayer=", 1))
  print(paste0("lattice width=", min(sapply(latticeInfo, function(x) {x$latWidth}))))
  print(paste0("total basis functions=", sum(sapply(latticeInfo, function(x) {x$nx * x$ny}))))
  
  latticeInfo2 = makeLatGridsKenya(nLayer=1, NC=30, nBuffer=5)
  print(paste0("NC=", 30, ", nLayer=", 1))
  print(paste0("lattice width=", min(sapply(latticeInfo2, function(x) {x$latWidth}))))
  print(paste0("total basis functions=", sum(sapply(latticeInfo2, function(x) {x$nx * x$ny}))))
  
  latticeInfo2 = makeLatGridsKenya(nLayer=1, NC=14, nBuffer=5)
  print(paste0("NC=", 14, ", nLayer=", 1))
  print(paste0("lattice width=", min(sapply(latticeInfo2, function(x) {x$latWidth}))))
  print(paste0("total basis functions=", sum(sapply(latticeInfo2, function(x) {x$nx * x$ny}))))
  
  test = c(latticeInfo2, latticeInfo)
  
  # NC=c(30, 107), 13725
  
  invisible(NULL)
}

testBasisResolution2 = function(nLayer, NC, nBuffer=5) {
  latticeInfo = makeLatGridsKenya(nLayer=nLayer, NC=NC, nBuffer=nBuffer)
  print(paste0("NC=", NC, ", nLayer=", nLayer))
  print(paste0("lattice width=", min(sapply(latticeInfo, function(x) {x$latWidth}))))
  print(paste0("total basis functions=", sum(sapply(latticeInfo, function(x) {x$nx * x$ny}))))
}

testBasisResolution3 = function(nLayer, NC, nBuffer=5) {
  latticeInfo = makeLatGrids(c(-1,1), c(-1,1), nLayer=nLayer, NC=NC, nBuffer=nBuffer)
  print(paste0("NC=", NC, ", nLayer=", nLayer))
  print(paste0("lattice width=", min(sapply(latticeInfo, function(x) {x$latWidth}))))
  print(paste0("total basis functions=", sum(sapply(latticeInfo, function(x) {x$nx * x$ny}))))
}










