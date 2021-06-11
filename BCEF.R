library(spNNGP)
makeAllPlots = TRUE
blueGreenCols = rev(makeGreenBlueSequentialColors(64))
yellowBlueCols = makeBlueGreenYellowSequentialColors(64)
# data(MI_TSCA)
# 
# head(MI_TSCA)
# fromproj4 = "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
# coordsOrig = SpatialPoints(cbind(MI_TSCA$long, MI_TSCA$lat), proj4string=CRS(fromproj4))
# lonLatCoords = spTransform(coordsOrig, CRS("+init=epsg:4326"))
# lonLatCoords = attr(lonLatCoords, "coords")
# head(lonLatCoords)
# range(lonLatCoords[,1])
# range(lonLatCoords[,2])
# 
# plot(lonLatCoords, pch=".", col="blue")
# US(add=TRUE)
# 
# 
# 
# 
# 
# data(nor2k)
# summary(nor2k)
# nor2kNGO <- spTransform(nor2k, CRS("+init=epsg:4273"))
# summary(nor2kNGO)
# all.equal(coordinates(nor2k)[,3], coordinates(nor2kNGO)[,3])
# 
# 
# #added after posting by Don MacQueen 
# crds <- cbind(c(-121.524764291826, -121.523480804667), c(37.6600366036405, 37.6543604613483))
# ref <- cbind(c(1703671.30566227, 1704020.20113366), c(424014.398045834, 421943.708664294))
# crs.step1.cf <- CRS(paste("+proj=lcc +lat_1=38.43333333333333",
#                           "+lat_2=37.06666666666667 +lat_0=36.5 +lon_0=-120.5",
#                           "+x_0=2000000.0 +y_0=500000.0 +ellps=GRS80 +units=us-ft +no_defs",
#                           "+towgs84=-0.991,1.9072,0.5129,0.025789908,0.0096501,0.0116599,0.0"))
# locs.step1.cf <- spTransform(SpatialPoints(crds,
#                                            proj4string=CRS("+proj=longlat +datum=WGS84")), crs.step1.cf)
# suppressWarnings(proj4string(locs.step1.cf) <- CRS(paste("+proj=lcc",
#                                                          "+lat_1=38.43333333333333 +lat_2=37.06666666666667 +lat_0=36.5",
#                                                          "+lon_0=-120.5 +x_0=2000000.0 +y_0=500000.0 +ellps=GRS80 +units=us-ft",
#                                                          "+no_defs +nadgrids=@null")))
# locs.step2.cfb <- spTransform(locs.step1.cf, CRS("+init=epsg:26743"))
# coordinates(locs.step2.cfb) - ref
# all.equal(unname(coordinates(locs.step2.cfb)), ref)
# 
# 
# 
# 
# 
# install.packages("spNNGP")
# library(spNNGP)
data(BCEF)

# 2000 percent tree cover from:
# https://storage.googleapis.com/earthenginepartners-hansen/GFC-2020-v1.8/download.html
library(tiff)
library(raster)
# out = readTIFF("~/git/LK-INLA/PTC_2000.tif")
# test = inla.nonconvex.hull(cbind(BCEF$x, BCEF$y))

if(FALSE) {
  plot(cbind(BCEF$x, BCEF$y), pch=".")
}

# leave out central points before making convex hull ----
xLows = c(260, 265.5, 268, 262, 275, 270)
xHighs = c(271, 277.5, 278, 268, 280, 272.5)
yLows = c(1644.5, 1649.6, 1653, 1643, 1657, 1646)
yHighs = c(1649.5, 1653.5, 1657, 1645, 1659.5, 1650)
allPoints = cbind(BCEF$x, BCEF$y)
leaveOutI = rep(FALSE, nrow(allPoints))
for(i in 1:length(xLows)) {
  leaveOutI = leaveOutI | ((BCEF$x >= xLows[i]) & (BCEF$x <= xHighs[i]) & 
                             (BCEF$y >= yLows[i]) & (BCEF$y <= yHighs[i]))
}
keepPoints = allPoints[!leaveOutI,]
leaveOutPoints = allPoints[leaveOutI,]
xlim = range(allPoints[,1])
ylim = range(allPoints[,2])
if(FALSE) {
  plot(keepPoints, pch=".", col="green", xlim=xlim, ylim=ylim)
  points(leaveOutPoints, pch=".", col="red")
}

nrow(allPoints)
sum(!leaveOutI)

# get the convex hull ----
chullI = chull(keepPoints[,1], keepPoints[,2])
chullPoints = keepPoints[chullI,]
if(FALSE) {
  plot(keepPoints, pch=".", col="green", xlim=xlim, ylim=ylim)
  points(leaveOutPoints, pch=".", col="red")
  points(chullPoints, pch=19, cex=.3, col="purple")
  polygon(chullPoints[,1], chullPoints[,2])
}

##### construct prediction points ----
xRange = range(allPoints[,1])
yRange = range(allPoints[,2])
xGrid = seq(xRange[1], xRange[2], by=.1)
yGrid = seq(yRange[1], yRange[2], by=.1)
predPoints = expand.grid(list(x=xGrid, y=yGrid))
inds = in.poly(predPoints, chullPoints)
sum(inds)
predPoints = predPoints[inds,]
fromproj4 = "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=km +no_defs"
coordsOrig = SpatialPoints(predPoints, proj4string=CRS(fromproj4))
predPointsLonLat = spTransform(coordsOrig, CRS("+init=epsg:4326"))
predPointsLonLat = attr(predPointsLonLat, "coords")

# get percent tree cover ----
ptc = raster("~/git/LK-INLA/PTC_2000.tif", values= TRUE)
predPTC = extract(ptc, SpatialPoints(predPoints, proj4string=CRS("+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=km +no_defs")),method="bilinear")
range(predPTC, na.rm=TRUE)

# plot percent tree cover
if(FALSE) {
  png("Figures/BCEF/PTC.png", width=500, height=500)
  quilt.plot(predPoints, predPTC, col=blueGreenCols, 
             nx=length(xGrid)-1, ny=length(yGrid)-1, 
             xlab="Easting (km)", ylab="Northing (km)")
  dev.off()
  
  png("Figures/BCEF/FCH.png", width=500, height=500)
  plotWithColor(BCEF$x, BCEF$y, BCEF$FCH, xlab="Easting (km)", ylab="Northing (km)", 
                colScale=yellowBlueCols, pch=19, cex=.1)
  dev.off()
  
  # A few exploratory plots
  pdf("Figures/BCEF/PTCvFCH.pdf", width=5, height=5)
  plot(BCEF$PTC, BCEF$FCH, pch=".", 
       xlab="Percent Tree Cover", ylab="Forest Canopy Height (m)")
  dev.off()
  
  set.seed(123)
  BCEFSubset = BCEF[sample(1:nrow(BCEF), 50000, replace=F),]
  out = loess(FCH~PTC, data=BCEFSubset)
  pdf("Figures/BCEF/PTCvFCH_loess.pdf", width=5, height=5)
  vals = predict(out, seq(0, 100))
  plot(out, pch=".", 
       xlab="Percent Tree Cover", ylab="Forest Canopy Height (m)")
  lines(seq(0, 100), vals, col="red")
  dev.off()
  
  cuts = cut(BCEF$PTC, breaks=seq(0, 100, by=5), include.lowest=TRUE)
  xs = seq(2.5, 97.5, by=5)
  breakInd = as.integer(cuts)
  xVals = xs[breakInd]
  meanVals = aggregate(BCEF$FCH, by=list(PTC=xVals), FUN=mean)$x
  
  pdf("Figures/BCEF/PTCvFCH_mean.pdf", width=5, height=5)
  plot(BCEF$PTC, BCEF$FCH, pch=".", 
       xlab="Percent Tree Cover", ylab="Forest Canopy Height (m)")
  points(xs, meanVals, col="red", pch=19, cex=.5)
  dev.off()
  
  meanVals = aggregate(BCEF$FCH, by=list(PTC=xVals), FUN=mean)$x
  
  pdf("Figures/BCEF/PTCvLogFCH_mean.pdf", width=5, height=5)
  plot(BCEF$PTC, BCEF$FCH, pch=".", 
       xlab="Percent Tree Cover", ylab="Forest Canopy Height (m)", 
       log="y")
  points(xs, meanVals, col="red", pch=19, cex=.5)
  dev.off()
}

##### construct rotated grid of lattice points ----
# first get range of shifted and rotated data domain for lattice
allPoints = cbind(BCEF$x, BCEF$y)
rotationAngle = 49.5 * (pi/180)
center = colMeans(allPoints)
centered = sweep(allPoints, 2, center, "-")
rotationMat = rbind(c(cos(rotationAngle), -sin(rotationAngle)), 
                    c(sin(rotationAngle), cos(rotationAngle)))
rotatedCentered = t(rotationMat %*% t(centered))
xRange = range(rotatedCentered[,1])
yRange = range(rotatedCentered[,2])

# construct basis functions on the rotated centered coordinate system
NCs = c(13, 50) # 2km and .5km resolution (this results in memory explosion...)
separateRanges = TRUE
NCs = c(25, 50, 99) # 1km, .5km, and .25km resolution (matrix factorization explosion...)
separateRanges = FALSE
NCs = c(25, 100)
separateRanges = TRUE # 1km and .25km (killed on the cluster, bad sampling on laptop)
NCs = c(25, 50) # 1km and .5km (killed on cluster)
separateRanges = TRUE
NCs = c(13, 25, 50) # 2km, 1km and .5km
separateRanges = FALSE

NCsText = paste(c("NC", NCs), collapse="_")

latInfo = makeLatGrids(xRange, yRange, NC=NCs, nLayer=length(NCs))
Ms = getMs(latInfo=latInfo)
Ms
sapply(latInfo, function(x) {x$latWidth})

# transform back to the original coordinate system
for(i in 1:length(latInfo)) {
  latInfo[[i]]$latCoords = t(t(rotationMat) %*% t(latInfo[[i]]$latCoords))
  latInfo[[i]]$latCoords = sweep(latInfo[[i]]$latCoords, 2, center, "+")
}

# plot observations and basis knots
xlim = range(c(latInfo[[1]]$latCoords[,1], allPoints[,1]))
ylim = range(c(latInfo[[1]]$latCoords[,2], allPoints[,2]))

# plot the lattice
if(makeAllPlots) {
  pdf(paste0("Figures/BCEF/lattice_", NCsText, ".pdf"), width=5, height=5)
  plot(latInfo[[1]]$latCoords, type="n", 
       xlim=xlim, ylim=ylim, asp=1, 
       xlab="Easting (km)", ylab="Northing (km)")
  cols = c("blue", "purple", "red")
  sizes = c(.3, .75, .05)
  for(i in 1:length(latInfo)) {
    plotPlusAtAngle(latInfo[[i]]$latCoords, col=cols[i], size=sizes[i], angle=-rotationAngle)
  }
  
  polygon(chullPoints[,1], chullPoints[,2])
  points(allPoints, pch=".", col="green")
  dev.off()
  
  png(paste0("Figures/BCEF/latticeFCH_", NCsText, ".png"), width=500, height=500)
  par(mar=c(5.1, 4.1, 4.1, 7))
  plot(latInfo[[1]]$latCoords, type="n", 
       xlim=xlim, ylim=ylim, asp=1, 
       xlab="Easting (km)", ylab="Northing (km)")
  cols = c("blue", "purple", "red")
  sizes = c(.3, .75, .05)
  for(i in 1:length(latInfo)) {
    plotPlusAtAngle(latInfo[[i]]$latCoords, col=cols[i], size=sizes[i], angle=-rotationAngle)
  }
  polygon(chullPoints[,1], chullPoints[,2])
  plotWithColor(BCEF$x, BCEF$y, BCEF$FCH, xlab="Easting (km)", ylab="Northing (km)", 
                colScale=yellowBlueCols, pch=19, cex=.1, new=FALSE)
  dev.off()
}

# Run the ELK analysis ----
Ns = c(1000, 5000, 25000)
Ns = c(1000, 5000, 25000, sum(BCEF$holdout == 0))
fitModels=TRUE
lastMod = NULL
for(i in 1:length(Ns)) {
  sampleN = ceiling(nrow(BCEF)/2)
  sampleN = sum(BCEF$holdout == 0)
  sampleN = Ns[i]
  
  set.seed(123)
  if(sampleN == sum(BCEF$holdout == 0)) {
    sampleDataI = BCEF$holdout == 0
  } else {
    sampleDataI = sample(1:nrow(BCEF), sampleN, replace=FALSE)
  }
  
  BCEFSubset = BCEF[sampleDataI,]
  ptcSubset = extract(ptc, SpatialPoints(cbind(BCEFSubset$x, BCEFSubset$y), proj4string=CRS("+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=km +no_defs")),method="bilinear")
  
  # check to make sure PTC from downloaded data agrees with BCEF dataset
  # mean(abs(ptcSubset - BCEFSubset$PTC) / BCEFSubset$PTC)
  pdf(paste0("Figures/BCEF/PTCcheck_N", sampleN, ".pdf"), width=5, height=5)
  plot(BCEFSubset$PTC, ptcSubset, pch=19, cex=.1, 
       xlab="Dataset PTC", ylab="Downloaded raster PTC")
  abline(0, 1)
  dev.off()
  
  # fit ELK model ----
  # priorPar are the spatial parameters
  # priorPar = getPCPrior(diff(range(BCEF$x))/5, .01, 1, nLayer=2, separateRanges=TRUE, 
  #                       latticeInfo=latInfo, useUrbanPrior=FALSE)
  priorPar = getPCPrior(diff(latInfo[[1]]$yRangeDat)/5, .01, 1, nLayer=length(NCs), separateRanges=separateRanges, 
                        latticeInfo=latInfo, useUrbanPrior=FALSE)
  if(fitModels) {
    precomputationFileNameRoot=paste(c("BCEFprecomputations", NCsText, "sepR", separateRanges), collapse="_")
    savePrecomputationResults=FALSE
    startTime = proc.time()
    if(sampleN > 25000) {
      bcefELK <- modBCEF(BCEFSubset, predPoints, predPTC, latInfo=latInfo, 
                                                 seed=1, rwModel="rw1", nNonlinearBasis=30, 
                                                 normalize=TRUE, fastNormalize=TRUE, 
                                                 intStrategy="eb", strategy="gaussian", 
                                                 printVerboseTimings=FALSE, priorPar=priorPar, 
                                                 loadPrecomputationResults=!savePrecomputationResults, separateRanges=separateRanges, 
                                                 savePrecomputationResults=savePrecomputationResults, 
                                                 precomputationFileNameRoot=precomputationFileNameRoot, 
                                                 previousFit=lastMod, diagonal=1000)
      lastMod = bcefELK$mod
    }
    bcefELK <- modBCEF(BCEFSubset, predPoints, predPTC, latInfo=latInfo, 
                       seed=1, rwModel="rw1", nNonlinearBasis=30, 
                       normalize=TRUE, fastNormalize=TRUE, 
                       intStrategy="ccd", strategy="gaussian", 
                       printVerboseTimings=FALSE, priorPar=priorPar, 
                       loadPrecomputationResults=!savePrecomputationResults, separateRanges=separateRanges, 
                       savePrecomputationResults=savePrecomputationResults, 
                       precomputationFileNameRoot=precomputationFileNameRoot, 
                       previousFit=lastMod)
    totalTime = proc.time() - startTime
    lastMod = bcefELK$mod
    bcefELK$mod = NULL
    save(bcefELK, totalTime, file=paste0("savedOutput/BCEF/bcefELK_", NCsText, "_N", sampleN, "_grid.RData"))
  } else {
    out = load(paste0("savedOutput/BCEF/bcefELK_", NCsText, "_N", sampleN, "_grid.RData"))
  }
  print(paste0("N: ", sampleN, ", total time (min): ", totalTime[3]/60))
  
  
  # (N, minutes) for laptop NC: 25, 100
  # (1000, 10), (5000, 17), (25000, 46), (100000, 368 or 6:08)
  
  # (N, minutes) for cluster NC: 25, 50, 99
  # (1000, 16), 
  
  # plot predictions ----
  preds = exp(bcefELK$preds)
  predSDs = bcefELK$sigmas
  predCIWidths = exp(bcefELK$upper) - exp(bcefELK$lower)
  
  png(paste0("Figures/BCEF/preds_", NCsText, "_N", sampleN, ".png"), width=500, height=500)
  quilt.plot(predPoints, preds, col=yellowBlueCols, 
             nx=length(xGrid)-1, ny=length(yGrid)-1, 
             xlab="Easting (km)", ylab="Northing (km)")
  dev.off()
  
  require(viridis)
  magmaCols = magma(64)
  png(paste0("Figures/BCEF/CIWidth_", NCsText, "_N", sampleN, ".png"), width=500, height=500)
  quilt.plot(predPoints, predCIWidths, col=magmaCols, 
             nx=length(xGrid)-1, ny=length(yGrid)-1, 
             xlab="Easting (km)", ylab="Northing (km)")
  dev.off()
  
  # plot effect of PTC
  rwKnots = bcefELK$rwKnots
  rwSummary = bcefELK$rwSummary
  fixedSummary = bcefELK$fixedEffectSummary
  fixedMat = bcefELK$fixedMat
  ptcMat = exp(bcefELK$rwMat + outer(rwKnots, fixedMat[2,]))
  # rwLower = exp(rwSummary[,5] + fixedSummary[2,4]*rwKnots)
  # rwUpper = exp(rwSummary[,6] + fixedSummary[2,4]*rwKnots)
  # rwEst = exp(rwSummary[,3] + fixedSummary[2,4]*rwKnots)
  rwLower = 100*(apply(ptcMat, 1, quantile, prob=.1)-1)
  rwUpper = 100*(apply(ptcMat, 1, quantile, prob=.9)-1)
  rwEst = 100*(rowMeans(ptcMat)-1)
  ylim=range(c(rwEst, rwLower, rwUpper))
  pdf(paste0("Figures/BCEF/PTCNonlinear_", NCsText, "_N", sampleN, ".pdf"), width=5, height=5)
  plot(rwKnots, rwEst, ylim=ylim, type="l", 
       xlab="Percent Tree Cover", ylab="Percent Change Forest Canopy Height")
  lines(rwKnots, rwLower, lty=2)
  lines(rwKnots, rwUpper, lty=2)
  dev.off()
  
  # plot individual spatial layers
  for(i in 1:length(latInfo)) {
    Amati = makeA(predPoints, latInfo, maxLayer=i, thisLayer=i)
    predsi = as.numeric(Amati %*% bcefELK$coefPreds[[i]])
    
    png(paste0("Figures/BCEF/logPredsLayer", i, "_", NCsText, "_N", sampleN, ".png"), width=500, height=500)
    quilt.plot(predPoints, predsi, col=yellowBlueCols, 
               nx=length(xGrid)-1, ny=length(yGrid)-1, 
               xlab="Easting (km)", ylab="Northing (km)")
    dev.off()
  }
}



# Construct SPDE mesh ----
mesh = 







