# script for comparing models
library(LatticeKrig)
library(geosphere)
library(viridis)
library(sp)
library(raster)
# Outline of comparison:
# Generate n simulations from a Matern Gaussian process
# Fit each simulation with:
#   Lattice Krig
#   Lattice Krig INLA
#   Matern GP (covariance par fitted)
#   Matern GP (true covariance par assumed)?
#  For each simulation calculate:
#   CRPS (energy score?) 
#   MSE
#   MAE
#   95% production intervals and whether the observations are within
# 

# generate simulations on rectangular domain (e.g. unit square) from a 
# Matern Gaussian Process with zero mean
genSimsMatern = function(xRange=c(0,1), yRange=c(0,1), n=100, nsim=100, 
                         beta=(xRange[2]-xRange[1])/10, nu=1.5, sigmaSq=1, tauSq=sigmaSq/10) {
  # generate locations of data
  xs = matrix(runif(n*nsim)*(xRange[2] - xRange[1]) + xRange[1], ncol=nsim)
  ys = matrix(runif(n*nsim)*(yRange[2] - yRange[1]) + yRange[1], ncol=nsim)
  
  # simulate from standard normal
  zSims = matrix(rnorm(n*nsim), ncol=nsim)
  errs = matrix(rnorm(n*nsim, sd=sqrt(tauSq)), ncol=nsim)
  
  # generate a simulation
  genSim = function(i) {
    # print progress
    if(i %% 10 == 0)
      print(paste0("iteration ", i, "/", nsim))
    
    # use Cholesky decomposition of covariance matrix to simulate
    Sigma = stationary.cov(cbind(xs[,i], ys[,i]), Covariance="MaternLR", theta=beta, phi=sigmaSq, nu=nu)
    L = t(chol(Sigma))
    L %*% zSims[,i]
  }
  trueMat = sapply(1:nsim, genSim)
  obsMat = trueMat + errs
  
  list(xs=xs, ys=ys, trueMat=trueMat, obsMat=obsMat)
}

# modified version of fields packages Matern function to code LR2011
# parameterization of the Matern covariance
# range = beta (this is the effective range for 99% correlation, roughly)
# all other parameters are the same as in ?Matern from fields package
MaternLR = function (d, range = 1, beta=range, alpha = 1/beta, smoothness = 0.5, nu = smoothness, 
                     phi = 1) {
  Matern(sqrt(8*nu)*d*alpha, range, alpha, smoothness, nu, phi)
}

# fit North American rainfall data set
# LatticeKrig paper also uses NorthAmericanRainfall dataset
# 1720 stations precip in JJA=June,July,August based on data from 1950-2010 (61 years).  
# Also includes elevation.
# dataset in LK package only includes linear trends and intercepts for each station
# 1720 x 1 = 1720 observations (or, with full dataset, 1720 x 61 = 104,920 observations)
# int.strategy: integration strategy of INLA. either "auto" or "eb" for empirical Bayes.
fitNAmData = function(normalize=TRUE, NC=5, nLayer=3, nBuffer=5, nu=1.5, seed=1, maxit=25, 
                      int.strategy="auto", fastNormalize=FALSE, doCovs=TRUE) {
  # load it
  data(NorthAmericanRainfall)
  x<- cbind(NorthAmericanRainfall$longitude,  NorthAmericanRainfall$latitude)
  y<- log(NorthAmericanRainfall$precip/10) # to mm units (take log first?)
  x = mercator(x)/1000 # m to km units
  X = matrix(NorthAmericanRainfall$elevation, ncol=1) # include elevation as covariate if desired
  
  # plot it
  par(mfrow=c(1,1))
  quilt.plot( x,y, main="Mean JJA Precipitation, 1950-2010 (log mm)", xlab="Longitude", ylab="Latitude")
  world( add=TRUE)
  
  # set up prediction locations for comparison with LatticeKrig
  latLim = c(31.75, 47.25)
  lonLim = c(-115, -101)
  predLon = seq(lonLim[1], lonLim[2], length = 50)
  predLat = seq(latLim[1], latLim[2], length = 50)
  predCoords = make.surface.grid(list(lon = predLon, lat = predLat))
  predCoords = mercator(predCoords)/1000
  
  # fit it
  if(!doCovs) {
    LKITime = system.time(LKIout <- fitLKINLASimple(x, y, predCoords, nu, seed, nLayer, NC, nBuffer, normalize=normalize,
                                                    priorPar=getPrior(2500, (700/10)^2, (700/2)^2), int.strategy=int.strategy, 
                                                    fastNormalize=fastNormalize))
  }
  else {
    stop("covariates knot yet implemented for LatticeKrigINLA")
  }
  # normalize: 394.113/361.118  (basically same for "eb"?)
  # fast normalize: 247.252/244.319
  if(!doCovs) {
    LKTime = system.time(LKout <- fitLK(x, y, predCoords, NC=NC, nLayer=nLayer, simpleMod=TRUE, normalize=normalize, 
                                        nBuffer=nBuffer, nu=nu, maxit=maxit))
  }
  else {
    LKTime = system.time(LKout <- fitLK(x, y, predCoords, NC=NC, nLayer=nLayer, simpleMod=TRUE, normalize=normalize, 
                                        nBuffer=nBuffer, nu=nu, maxit=maxit, XObs=X))
  }
  # normalize: 105.155, 241.348/219.526/226.204 ?
  
  # show predictions just for LatticeKrig
  zlimPreds = range(LKout$preds)
  zlimSEs = range(LKout$SEs)
  pdf(file="Figures/NAmTemperaturePredsLK.pdf", width=7, height=7)
  par(mfrow=c(1,2))
  LKMeans = exp(LKout$preds + LKout$SEs^2/2)
  LKVars = (exp(LKout$SEs^2)-1) * LKMeans^2
  quilt.plot(predCoords, LKMeans, main="LatticeKrig predictions (in mm)", zlim=zlimPreds, 
             col=viridis(64))
  points(x, pch=19, cex=.4)
  quilt.plot(predCoords, sqrt(LKVars), main="LatticeKrig SEs", zlim=zlimSEs)
  points(x, pch=19, cex=.4)
  US(add=TRUE)
  dev.off()
  
  # show it all
  zlimPreds = range(c(LKout$preds, LKIout$preds))
  zlimSEs = range(c(LKout$SEs, LKIout$SDs))
  pdf(file="Figures/NAmTemperaturePreds.pdf", width=7, height=7)
  par(mfrow=c(2,2))
  quilt.plot(predCoords, LKout$preds, main="LatticeKrig predictions", zlim=zlimPreds)
  quilt.plot(predCoords, LKout$SEs, main="LatticeKrig SEs", zlim=zlimSEs)
  quilt.plot(predCoords, LKIout$preds, main="LatticeKrig-INLA predictions", zlim=zlimPreds)
  quilt.plot(predCoords, LKIout$SDs, main="LatticeKrig-INLA SDs", zlim=zlimSEs)
  US(add=TRUE)
  dev.off()
  
  pdf(file="Figures/NAmTemperaturePredsVs.pdf", width=10, height=5)
  par(mfrow=c(1,2))
  plot(LKout$preds, LKIout$preds, pch=".", main="LK-INLA vs LatticeKrig predictions", 
       xlab="LatticeKrig", ylab="LK-INLA")
  abline(0,1, col="green")
  plot(LKout$SEs, LKIout$SDs, pch=".", log="xy", xlab="LatticeKrig SEs", ylab="LK-INLA SDs", 
       main="LK-INLA vs LatticeKrig predictive uncertainty")
  abline(0,1, col="green")
  dev.off()
  
  # produce figures like in the LatticeKrig paper of central predictions and standard errors
  
  
  list(LKout=LKout, LKIout=LKIout, LKTime=LKTime, LKITime=LKITime)
  # how good are the INLA results?  Compare with MCMC?
}

# 
doSimStudy = function() {
  3
}

# get elevation at the given set of lat/lon coordinates
getElevation = function(lon, lat) {
  # ogrListLayers("gtopo30/gtopo30.shp")
  shape = readOGR("gtopo30/")
  shape <- shapefile("gtopo30/gtopo30.shp")
  datapol <- data.frame(shape)
  pointtoplot <- data.frame(x=-20, y=40)
  coordinates(pointtoplot) <- ~ x + y 
  proj4string(pointtoplot) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  test <- data.frame(xx=over(shape, pointtoplot))
  combine <- cbind(test, datapol)
  combine <- na.omit(combine)
}


