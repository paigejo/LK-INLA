# the script for making the plots for Jon's grant
setwd('~/Google Drive/UW/Wakefield/WakefieldShared/LK-INLA/code/')
source('~/Google Drive/UW/Wakefield/WakefieldShared/LK-INLA/code/LKinla.R')
source('~/Google Drive/UW/Wakefield/WakefieldShared/LK-INLA/code/LKinla_rgeneric.R')

##### unincluded plot, of the data locations
library(LatticeKrig)
data(NorthAmericanRainfall)

# project the locations stereographically
library(mapproj)
xStereo<- mapproject( NorthAmericanRainfall$lon,NorthAmericanRainfall$lat, projection="stereographic")
NorthAmericanRainfall$x.s<- cbind( xStereo$x, xStereo$y)
NorthAmericanRainfall$projection<- .Last.projection

# x<- cbind(NorthAmericanRainfall$longitude,  NorthAmericanRainfall$latitude)
# y<- NorthAmericanRainfall$precip

x = NorthAmericanRainfall$x.s
y = log10(NorthAmericanRainfall$precip/10) #in log10(mm) now

pdf("Figures/precipitationCoordsDataSet.pdf", width=5, height=5)
plot(x[, 1], x[, 2], main="Mean JJA Precipitation, 1950-2010 (mm)", xlab="Longitude", ylab="Latitude", 
     pch=19, cex=.1)
world( add=TRUE, projection=NorthAmericanRainfall$projection()$projection, 
       orientation=NorthAmericanRainfall$projection()$orientation)
dev.off()

##### first plot, of the spatial domain and the basis function locations
coords = x
latLim = range(NorthAmericanRainfall$latitude)
lonLim = range(NorthAmericanRainfall$longitude)
NC=13
nBuffer=5

# get lattice points, prediction points
xRangeDat = range(coords[,1])
yRangeDat = range(coords[,2])
latInfo = makeLatGrids(xRangeDat, yRangeDat, NC, nBuffer, nLayer=3)

# get first, coarsest layer
pts = latInfo$latCoords[[1]]
border = rbind(c(xRangeDat[1], yRangeDat[1]), 
               c(xRangeDat[1], yRangeDat[2]), 
               c(xRangeDat[2], yRangeDat[2]), 
               c(xRangeDat[2], yRangeDat[1]))
pdf("Figures/grantBasisKnots.pdf", width = 5, height = 5)
plot(pts[,1], pts[,2], type="n", xlab="Easting", ylab="Northing", main="Data domain and buffer")
polygon(border[,1], border[,2], col=rgb(.5, .5, .5, .5), border=rgb(0, 0, 0, 0))
points(pts[,1], pts[,2], pch="+", cex=.6)
world( add=TRUE, projection=NorthAmericanRainfall$projection()$projection, 
       orientation=NorthAmericanRainfall$projection()$orientation, xlim=lonLim, 
       ylim=latLim)

# plot second layer
pts = latInfo$latCoords[[2]]
points(pts[,1], pts[,2], pch="+", cex=.3, col="darkGreen")

# plot third layer
pts = latInfo$latCoords[[3]]
points(pts[,1], pts[,2], pch="+", cex=.2, col="red")
dev.off()

##### second plot, of the posterior median
# first generate the prediction locations (predict)
predCoords = latInfo$latCoords[[3]]
predCoords = predCoords[in.poly(predCoords, border),] # 1920 x 2

# fit the model
# y = NorthAmericanRainfall$precip
# out = fitLKINLASimple(coords, y, predPts=predCoords, nu=1.5, seed=1, nLayer=3, NC=13,
#                       nBuffer=5, priorPar=getPrior(.1, .1, 10), 
#                       X=cbind(1, coords), XPred=cbind(1, predCoords), normalize=TRUE, 
#                       intStrategy="eb", strategy="gaussian", fastNormalize=TRUE)
# quilt.plot(predCoords, out$preds) # there still seems to be a problem with the fast normalization
# world( add=TRUE, projection=NorthAmericanRainfall$projection()$projection, 
#        orientation=NorthAmericanRainfall$projection()$orientation, xlim=lonLim, 
#        ylim=latLim)

testI = 1:100
testI = 1:nrow(coords)
testPredI = 1:100
testPredI = 1:nrow(predCoords)
prior = getPrior(1/4, diff(range(y)) / 100, diff(range(y)) / 5)
prior = getPrior2(1 / 4, 3)
prior = getPrior2(1 / 4, 6)
xs = seq(qexp(.05, rate = prior$sdRate), qexp(.95, rate = prior$sdRate), length=100)
plot(xs, dexp(xs, rate = prior$sdRate), type = "l")
inla.setOption(keep = FALSE)
outSlow = fitLKINLASimple(coords[testI,], y[testI], predPts=predCoords[testPredI,], nu=1.5, seed=1, nLayer=3, NC=13,
                      nBuffer=5, priorPar=prior, 
                      X=cbind(1, coords[testI,]), XPred=cbind(1, predCoords[testPredI,]), normalize=TRUE, 
                      intStrategy="eb", strategy="gaussian", fastNormalize=FALSE)

# make sure to delete temporary files in /var/folders/6y/zncx66tn3z13crx_82g9kdnh0000gp/T if necessary

png("Figures/predictiveField.png", width=560, height=560)
quilt.plot(predCoords, 10^(outSlow$preds), nx=latInfo$nx[3] -11, ny=latInfo$ny[3]-10, 
           main="Median Precipitation, 1950-2010 (mm)", xlab="Easting", ylab="Northing")
world( add=TRUE, projection=NorthAmericanRainfall$projection()$projection, 
       orientation=NorthAmericanRainfall$projection()$orientation, xlim=lonLim, 
       ylim=latLim)
dev.off()

##### generate both plots together:
png("Figures/fullPlot6.png", width=1240, height=560)
par(mfrow=c(1, 2))
coords = x
latLim = range(NorthAmericanRainfall$latitude)
lonLim = range(NorthAmericanRainfall$longitude)
NC=13
nBuffer=5

# get lattice points, prediction points
xRangeDat = range(coords[,1])
yRangeDat = range(coords[,2])
latInfo = makeLatGrids(xRangeDat, yRangeDat, NC, nBuffer, nLayer=3)

# get first, coarsest layer
pts = latInfo$latCoords[[1]]
border = rbind(c(xRangeDat[1], yRangeDat[1]), 
               c(xRangeDat[1], yRangeDat[2]), 
               c(xRangeDat[2], yRangeDat[2]), 
               c(xRangeDat[2], yRangeDat[1]))
plot(pts[,1], pts[,2], type="n", xlab="Easting", ylab="Northing", main="Basis Knot Locations")
polygon(border[,1], border[,2], col=rgb(.5, .5, .5, .5), border=rgb(0, 0, 0, 0))
points(pts[,1], pts[,2], pch="+", cex=1)
world( add=TRUE, projection=NorthAmericanRainfall$projection()$projection, 
       orientation=NorthAmericanRainfall$projection()$orientation, xlim=lonLim, 
       ylim=latLim)

# plot second layer
pts = latInfo$latCoords[[2]]
points(pts[,1], pts[,2], pch="+", cex=.75, col="darkGreen")

# plot third layer
pts = latInfo$latCoords[[3]]
points(pts[,1], pts[,2], pch="+", cex=.4, col="red")

# quilt.plot(predCoords, outSlow$preds / 10, nx=latInfo$nx[3] -11, ny=latInfo$ny[3]-10, 
#            main="Mean JJA Precipitation, 1950-2010 (mm)", xlab="Easting", ylab="Northing")
quilt.plot(predCoords, 10^(outSlow$preds), nx=latInfo$nx[3] -11, ny=latInfo$ny[3]-10, 
           main="Median Precipitation (mm)", xlab="Easting", ylab="Northing", 
           legend.mar=10)
world( add=TRUE, projection=NorthAmericanRainfall$projection()$projection, 
       orientation=NorthAmericanRainfall$projection()$orientation, xlim=lonLim, 
       ylim=latLim)

dev.off()


##### subset predictions over land based on example 7 from:
# http://www.geog.uoregon.edu/GeogR/examples/maps_examples02.htm
library(maptools)

outlines = map("world", c("usa"), projection=NorthAmericanRainfall$projection()$projection, 
               orientation=NorthAmericanRainfall$projection()$orientation, 
               fill=TRUE, col="red", xlim=lonLim, ylim = c(15, 50))
pruned = pruneMap(outlines)
usa.sp = map2SpatialPolygons(pruned, 1:52, checkHoles = TRUE)
plot(usa.sp)

outlines = map("world", c("canada"), projection=NorthAmericanRainfall$projection()$projection, 
               orientation=NorthAmericanRainfall$projection()$orientation, 
               fill=TRUE, col="red", xlim=lonLim, ylim = c(40, 90))
pruned = pruneMap(outlines)
canada.sp = map2SpatialPolygons(pruned, 1:(sum(is.na(pruned$x)) + 1), checkHoles = TRUE)
plot(canada.sp)

outlines = map("world", c("mexico"), projection=NorthAmericanRainfall$projection()$projection, 
               orientation=NorthAmericanRainfall$projection()$orientation, 
               fill=TRUE, col="red", xlim=lonLim, ylim = c(0, 40))
pruned = pruneMap(outlines)
mexico.sp = map2SpatialPolygons(pruned, 1:(sum(is.na(pruned$x)) + 1), checkHoles = TRUE)
plot(mexico.sp)

predPoints = SpatialPoints(predCoords)
inUsa = !is.na(over(predPoints, usa.sp))
inCanada = !is.na(over(predPoints, canada.sp))
inMexico = !is.na(over(predPoints, mexico.sp))
overLand = inUsa | inCanada | inMexico

png("Figures/fullPlotRestricted6.png", width=1240, height=560)
par(mfrow=c(1, 2))
coords = x
latLim = range(NorthAmericanRainfall$latitude)
lonLim = range(NorthAmericanRainfall$longitude)
NC=13
nBuffer=5

# get lattice points, prediction points
xRangeDat = range(coords[,1])
yRangeDat = range(coords[,2])
latInfo = makeLatGrids(xRangeDat, yRangeDat, NC, nBuffer, nLayer=3)

# get first, coarsest layer
pts = latInfo$latCoords[[1]]
border = rbind(c(xRangeDat[1], yRangeDat[1]), 
               c(xRangeDat[1], yRangeDat[2]), 
               c(xRangeDat[2], yRangeDat[2]), 
               c(xRangeDat[2], yRangeDat[1]))
plot(pts[,1], pts[,2], type="n", xlab="Easting", ylab="Northing", main="Basis Knot Locations")
polygon(border[,1], border[,2], col=rgb(.5, .5, .5, .5), border=rgb(0, 0, 0, 0))
points(pts[,1], pts[,2], pch="+", cex=1)
world( add=TRUE, projection=NorthAmericanRainfall$projection()$projection, 
       orientation=NorthAmericanRainfall$projection()$orientation, xlim=lonLim, 
       ylim=latLim)

# plot second layer
pts = latInfo$latCoords[[2]]
points(pts[,1], pts[,2], pch="+", cex=.75, col="darkGreen")

# plot third layer
pts = latInfo$latCoords[[3]]
points(pts[,1], pts[,2], pch="+", cex=.4, col="red")

# plot predictive field
quilt.plot(predCoords[overLand,], 10^(outSlow$preds[overLand]), nx=latInfo$nx[3] -11, ny=latInfo$ny[3]-10, 
           main="Median Precipitation (mm)", xlab="Easting", ylab="Northing", 
           legend.mar=10, xlim=range(latInfo$latCoords[[3]][,1]), ylim=range(latInfo$latCoords[[3]][,2]))
world( add=TRUE, projection=NorthAmericanRainfall$projection()$projection, 
       orientation=NorthAmericanRainfall$projection()$orientation, xlim=lonLim, 
       ylim=latLim)
dev.off()

# analyze the data with the latticeKrig package
lkFit = fitLK(coords, y, predCoords, NULL, NULL, 13, 3, TRUE, TRUE, 5, 1.5, TRUE, doSEs = FALSE)
lkFit$preds
png("Figures/predictionsRestrictedLK.png", width=560, height=560)
quilt.plot(predCoords[overLand,], 10^(lkFit$preds[overLand]), nx=latInfo$nx[3] -11, ny=latInfo$ny[3]-10, 
           main="Median Precipitation (mm)", xlab="Easting", ylab="Northing", 
           legend.mar=10, xlim=range(latInfo$latCoords[[3]][,1]), ylim=range(latInfo$latCoords[[3]][,2]))
world( add=TRUE, projection=NorthAmericanRainfall$projection()$projection, 
       orientation=NorthAmericanRainfall$projection()$orientation, xlim=lonLim, 
       ylim=latLim)
dev.off()