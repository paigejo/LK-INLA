library(spNNGP)
data(MI_TSCA)

head(MI_TSCA)
fromproj4 = "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
coordsOrig = SpatialPoints(cbind(MI_TSCA$long, MI_TSCA$lat), proj4string=CRS(fromproj4))
lonLatCoords = spTransform(coordsOrig, CRS("+init=epsg:4326"))
lonLatCoords = attr(lonLatCoords, "coords")
head(lonLatCoords)
range(lonLatCoords[,1])
range(lonLatCoords[,2])

plot(lonLatCoords, pch=".", col="blue")
US(add=TRUE)





data(nor2k)
summary(nor2k)
nor2kNGO <- spTransform(nor2k, CRS("+init=epsg:4273"))
summary(nor2kNGO)
all.equal(coordinates(nor2k)[,3], coordinates(nor2kNGO)[,3])


#added after posting by Don MacQueen 
crds <- cbind(c(-121.524764291826, -121.523480804667), c(37.6600366036405, 37.6543604613483))
ref <- cbind(c(1703671.30566227, 1704020.20113366), c(424014.398045834, 421943.708664294))
crs.step1.cf <- CRS(paste("+proj=lcc +lat_1=38.43333333333333",
                          "+lat_2=37.06666666666667 +lat_0=36.5 +lon_0=-120.5",
                          "+x_0=2000000.0 +y_0=500000.0 +ellps=GRS80 +units=us-ft +no_defs",
                          "+towgs84=-0.991,1.9072,0.5129,0.025789908,0.0096501,0.0116599,0.0"))
locs.step1.cf <- spTransform(SpatialPoints(crds,
                                           proj4string=CRS("+proj=longlat +datum=WGS84")), crs.step1.cf)
suppressWarnings(proj4string(locs.step1.cf) <- CRS(paste("+proj=lcc",
                                                         "+lat_1=38.43333333333333 +lat_2=37.06666666666667 +lat_0=36.5",
                                                         "+lon_0=-120.5 +x_0=2000000.0 +y_0=500000.0 +ellps=GRS80 +units=us-ft",
                                                         "+no_defs +nadgrids=@null")))
locs.step2.cfb <- spTransform(locs.step1.cf, CRS("+init=epsg:26743"))
coordinates(locs.step2.cfb) - ref
all.equal(unname(coordinates(locs.step2.cfb)), ref)





install.packages("spNNGP")
library(spNNGP)
data(BCEF)
head(BCEF)
summary(BCEF)
fromproj4 = "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=km +no_defs"
coordsOrig = SpatialPoints(cbind(BCEF$x, BCEF$y), proj4string=CRS(fromproj4))
lonLatCoords = spTransform(coordsOrig, CRS("+init=epsg:4326"))
lonLatCoords = attr(lonLatCoords, "coords")

# the result is not in Alaska
head(lonLatCoords)
range(lonLatCoords[,1])
# [1] -148.5642 -148.0947
range(lonLatCoords[,2])
# [1] 64.64402 64.79030

plot(lonLatCoords, pch=".", col="blue")
US(add=TRUE)

# 2000 percent tree cover from:
# https://storage.googleapis.com/earthenginepartners-hansen/GFC-2020-v1.8/download.html
library(tiff)
library(raster)
# out = readTIFF("~/git/LK-INLA/PCT_2000.tif")
test = inla.nonconvex.hull(cbind(BCEF$x, BCEF$y))

plot(cbind(BCEF$x, BCEF$y), pch=".")
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
plot(keepPoints, pch=".", col="green", xlim=xlim, ylim=ylim)
points(leaveOutPoints, pch=".", col="red")
nrow(allPoints)
sum(!leaveOutI)

scaleFun = function(pts, factor) {
  center = colMeans(pts)
  centered = sweep(pts, 2, center, "-")
  scaledCentered = centered * factor
  scaled = sweep(scaledCentered, 2, center, "+")
  scaled
}

scaleFunAnisotropic = function(pts, factors, rotationAngle=0, shift=c(0, 0)) {
  center = colMeans(pts)
  centered = sweep(pts, 2, center, "-")
  rotationMat = rbind(c(cos(rotationAngle), -sin(rotationAngle)), 
                      c(sin(rotationAngle), cos(rotationAngle)))
  rotatedCentered = t(rotationMat %*% t(centered))
  scaledRotatedCentered = sweep(rotatedCentered, 2, factors, "*")
  scaledCentered = t(t(rotationMat) %*% t(scaledRotatedCentered))
  scaled = sweep(scaledCentered, 2, center + shift, "+")
  scaled
}

system.time(test <- inla.nonconvex.hull(keepPoints, convex=2.5))

factors = c(1, 1)
factors = c(.78, .85)
rotationAngle = 49.5 * (pi/180)
shift = c(0, 0)
shift = c(.2, 0)

# for testing purposes
centered = sweep(allPoints, 2, colMeans(allPoints), "-")
rotationMat = rbind(c(cos(rotationAngle), -sin(rotationAngle)), 
                    c(sin(rotationAngle), cos(rotationAngle)))
plot(t(rotationMat %*% t(centered)), pch=".", col="green")

# scale the hull and plot to make sure it's a good fit
hull = scaleFunAnisotropic(test$loc, factors, rotationAngle, c(.2, 0))
xlim = range(c(allPoints[,1], hull[,1]))
ylim = range(c(allPoints[,2], hull[,2]))
plot(keepPoints, pch=".", col="green", xlim=xlim, ylim=ylim)
points(leaveOutPoints, pch=".", col="red")
points(hull, pch=19, cex=.3, col="purple")

# compare to the convex hull
chullI = chull(keepPoints[,1], keepPoints[,2])
chullPoints = keepPoints[chullI,]
plot(keepPoints, pch=".", col="green", xlim=xlim, ylim=ylim)
points(leaveOutPoints, pch=".", col="red")
points(chullPoints, pch=19, cex=.3, col="purple")
polygon(chullPoints[,1], chullPoints[,2])

spatialDomain <- Polygon(chullPoints)

##### construct rotated grid of lattice points ----
# first rotate points
rotationAngle = 49.5 * (pi/180)
center = colMeans(allPoints)
centered = sweep(allPoints, 2, center, "-")
rotationMat = rbind(c(cos(rotationAngle), -sin(rotationAngle)), 
                    c(sin(rotationAngle), cos(rotationAngle)))
rotatedCentered = t(rotationMat %*% t(centered))
plot(rotatedCentered, pch=".", col="green")

# construct basis functions on the rotated centered coordinate system
xRange = range(rotatedCentered[,1])
yRange = range(rotatedCentered[,2])
latInfo = makeLatGrids(xRange, yRange, NC=c(25, 100), nLayer=2)
Ms = getMs(latInfo=latInfo)
Ms

# transform back to the original coordinate system
for(i in 1:length(latInfo)) {
  latInfo[[i]]$latCoords = t(t(rotationMat) %*% t(latInfo[[i]]$latCoords))
  latInfo[[i]]$latCoords = sweep(latInfo[[i]]$latCoords, 2, center, "+")
}

# plot observations and basis knots
xlim = range(c(latInfo[[1]]$latCoords[,1], allPoints[,1]))
ylim = range(c(latInfo[[1]]$latCoords[,2], allPoints[,2]))

# plot the lattice
pdf("Figures/BCEF/lattice.pdf", width=5, height=5)
plot(latInfo[[1]]$latCoords, type="n", 
       xlim=xlim, ylim=ylim, asp=1, 
       xlab="Easting (km)", ylab="Northing (km)")
plotPlusAtAngle(latInfo[[1]]$latCoords, col="blue", size=.3, angle=-rotationAngle)
plotPlusAtAngle(latInfo[[2]]$latCoords, col="purple", size=.075, angle=-rotationAngle)
polygon(chullPoints[,1], chullPoints[,2])
points(allPoints, pch=".", col="green")
dev.off()

png("Figures/BCEF/latticeFCH.png", width=500, height=500)
par(mar=c(5.1, 4.1, 4.1, 7))
plot(latInfo[[1]]$latCoords, type="n", 
     xlim=xlim, ylim=ylim, asp=1, 
     xlab="Easting (km)", ylab="Northing (km)")
plotPlusAtAngle(latInfo[[1]]$latCoords, col="blue", size=.3, angle=-rotationAngle)
plotPlusAtAngle(latInfo[[2]]$latCoords, col="purple", size=.075, angle=-rotationAngle)
polygon(chullPoints[,1], chullPoints[,2])
plotWithColor(BCEF$x, BCEF$y, BCEF$FCH, xlab="Easting (km)", ylab="Northing (km)", 
              colScale=yellowBlueCols, pch=19, cex=.1, new=FALSE)
dev.off()

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
pct = raster("~/git/LK-INLA/PCT_2000.tif", values= TRUE)

predPCT = extract(pct, SpatialPoints(predPoints, proj4string=CRS("+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=km +no_defs")),method="bilinear")
range(predPCT, na.rm=TRUE)

# plot percent tree cover
blueGreenCols = rev(makeGreenBlueSequentialColors(64))
yellowBlueCols = makeBlueGreenYellowSequentialColors(64)
png("Figures/BCEF/PTC.png", width=500, height=500)
quilt.plot(predPoints, pcts, col=blueGreenCols, 
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
meanVals = aggregate(BCEF$FCH, by=list(PCT=xVals), FUN=mean)$x

pdf("Figures/BCEF/PTCvFCH_mean.pdf", width=5, height=5)
plot(BCEF$PTC, BCEF$FCH, pch=".", 
     xlab="Percent Tree Cover", ylab="Forest Canopy Height (m)")
points(xs, meanVals, col="red", pch=19, cex=.5)
dev.off()

meanVals = aggregate(BCEF$FCH, by=list(PCT=xVals), FUN=mean)$x

pdf("Figures/BCEF/PTCvLogFCH_mean.pdf", width=5, height=5)
plot(BCEF$PTC, BCEF$FCH, pch=".", 
     xlab="Percent Tree Cover", ylab="Forest Canopy Height (m)", 
     log="y")
points(xs, meanVals, col="red", pch=19, cex=.5)
dev.off()

# Run the analysis ----
sampleN = 25000
set.seed(123)
sampleDataI = sample(1:nrow(BCEF), sampleN, replace=FALSE)
BCEFSubset = BCEF[sampleDataI,]
pctSubset = extract(pct, SpatialPoints(cbind(BCEFSubset$x, BCEFSubset$y), proj4string=CRS("+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=km +no_defs")),method="bilinear")
head(pctSubset)
head(BCEFSubset$PTC)

# check to make sure PTC from downloaded data agrees with BCEF dataset
mean(abs(pctSubset - BCEFSubset$PTC) / BCEFSubset$PTC)
pdf(paste0("Figures/BCEF/PTCcheck_N", sampleN, ".pdf"), width=5, height=5)
plot(BCEFSubset$PTC, pctSubset, pch=19, cex=.1, 
     xlab="Dataset PTC", ylab="Downloaded raster PTC")
abline(0, 1)
dev.off()


# fit model ----
system.time(bcefELK <- modBCEF(BCEFSubset, predPoints, predPCT, latInfo=latInfo, 
                  seed=1, rwModel="rw1", nNonlinearBasis=20, 
                  normalize=TRUE, fastNormalize=TRUE, 
                  intStrategy="ccd", strategy="gaussian", 
                  printVerboseTimings=TRUE, 
                  loadPrecomputationResults=TRUE))
bcefELK$timings$totalTime/60

# (N, minutes)
# (1000, ~9.5), (5000, 24)


# plot predictions ----
preds = exp(bcefELK$preds)
predSDs = bcefELK$sigmas
predCIWidths = exp(bcefELK$upper) - exp(bcefELK$lower)

png(paste0("Figures/BCEF/predsNC_25_100_N", sampleN, ".png"), width=500, height=500)
quilt.plot(predPoints, preds, col=yellowBlueCols, 
           nx=length(xGrid)-1, ny=length(yGrid)-1, 
           xlab="Easting (km)", ylab="Northing (km)")
dev.off()

library(viridis)
magmaCols = magma(64)
png(paste0("Figures/BCEF/CIWidthNC_25_100_N", sampleN, ".png"), width=500, height=500)
quilt.plot(predPoints, predCIWidths, col=magmaCols, 
           nx=length(xGrid)-1, ny=length(yGrid)-1, 
           xlab="Easting (km)", ylab="Northing (km)")
dev.off()

# plot effect of PTC
rwKnots = bcefELK$rwKnots
rwSummary = bcefELK$rwSummary
rwLower = exp(rwSummary[,5])
rwUpper = exp(rwSummary[,6])
rwEst = exp(rwSummary[,3])
ylim=range(c(rwEst, rwLower, rwUpper))
pdf(paste0("PTCNonlinear_NC25_100_N", sampleN, ".pdf"), width=5, height=5)
plot(rwKnots, rwEst, ylim=ylim, lty=1, 
     xlab="Percent Tree Cover", ylab="Forest Canopy Height (m)")
lines(rwKnots, rwLower, lty=2)
lines(rwKnots, rwUpper, lty=2)
dev.off()



