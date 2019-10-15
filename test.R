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
                                   printVerboseTimings=FALSE, nObs=rep(25, n)) {
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
                                               intStrategy="auto", strategy="gaussian", fastNormalize=fastNormalize, 
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
  effRange = 2.3/kappa * latticeWidth
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
testLKModel = function(buffer=2.5, kappa=1, rho=1, nu=1.5, seed=1, nLayer=3, nx=20, ny=nx, n=50, sigma2 = .1^2, 
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
testLKINLAEd = function(buffer=2.5, kappa=1, rho=1, nu=1.5, seed=1, nLayer=3, nx=20, ny=nx, n=900, 
                        nBuffer=5, normalize=TRUE, fastNormalize=TRUE, NC=5, spatialFixedEffect=TRUE, 
                        printVerboseTimings=FALSE, nObs=rep(25, n), urbanEffect=TRUE, popGrid=NULL, 
                        kmres=5) {
  set.seed(seed)
  
  # load secondary education prevalence data
  out= load("~/git/U5MR/kenyaDataEd.RData")
  
  # load population density data (adjusted for target population)
  if(is.null(popGrid)) {
    if(kmres == 5)
      load("~/git/U5MR/popGridAdjustedWomen.RData")
    else
      stop("kmres other than 5km not supported")
  }
  predPts = cbind(popGrid$east, popGrid$north)
  
  # generate lattice and simulate observations
  coords = cbind(ed$east, ed$north) # TODO: lon,lat?
  xRangeDat = range(coords[,1])
  yRangeDat = range(coords[,2])
  latInfo = makeLatGrids(xRangeDat, yRangeDat, NC, nBuffer, nLayer)
  
  # set observations
  nObs = ed$n
  ys = ed$y
  
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
  
  ## set fixed effects
  # intercept
  X = matrix(rep(1, n), ncol=1)
  # X = matrix(coords[,1], ncol=1)
  XPred = matrix(rep(1, mx*my), ncol=1)
  
  # add linear terms in lat/lon to covariate matrices if requested
  if(spatialFixedEffect) {
    X = cbind(X, coords)
    XPred = cbind(XPred, predPts)
  }
  
  # add urban effect
  if(urbanEffect) {
    X = cbind(X, ed$urban)
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
  time = system.time(out <- fitLKINLAStandard2(coords, ys, predCoords=predPts, seed=seed, nLayer=nLayer, NC=NC,
                                               nBuffer=nBuffer, priorPar=priorPar, xObs=X, xPred=XPred, normalize=normalize, 
                                               intStrategy="auto", strategy="laplace", fastNormalize=fastNormalize, 
                                               printVerboseTimings=printVerboseTimings, family="binomial", obsNs=nObs))
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
  
  pdf(file="Figures/standardPredsBinom.pdf", width=15, height=6)
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
  out = inla.hyperpar.sample(10000, mod, improve.marginals=TRUE)
  zSamples = out[,4:(3+nLayer-1)]
  xSamples = apply(zSamples, 1, multivariateExpit)
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

# modelResults must have matrices "predMat" and "hyperMat"
# B: 
aggregateModelPreds = function(modelResults, B, addClusterError=TRUE) {
  
}







