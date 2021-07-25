# get color schemes ----
blueGreenCols = rev(makeGreenBlueSequentialColors(64))
yellowBlueCols = makeBlueGreenYellowSequentialColors(64)
infernoCols = inferno(64)
magmaCols = magma(64)
plasmaCols = plasma(64)

makePlots = FALSE

# get datasets ----
# analysis of the heaton et al. dataset
# https://link.springer.com/article/10.1007/s13253-018-00348-w
# https://github.com/finnlindgren/heatoncomparison
out = load("heaton/AllSatelliteTemps.RData")
names(all.sat.temps)
out = load("heaton/SatelliteTemps.RData")
names(sat.temps)
out = load("heaton/SmallTestData.RData")
names(code.test)

# project onto a reasonable coordinate system
# see https://gis.stackexchange.com/questions/141580/which-projection-is-best-for-mapping-the-contiguous-united-states 
# for coordinate systems for contiguous USA. Choose: EPSG 102004
# from lon/lat coords to easting/northing
lonLatCoords = SpatialPoints(cbind(all.sat.temps$Lon, all.sat.temps$Lat), proj4string=CRS("+proj=longlat"))
coordsUTM = spTransform(lonLatCoords, CRS("+proj=lcc +lat_1=33 +lat_2=45 +lat_0=39 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"))
eastNorth = attr(coordsUTM, "coords")
all.sat.temps$east = eastNorth[,1]
all.sat.temps$north = eastNorth[,2]

# # get sun's azimuth angle
# library(suncalc)
# centerEast = mean(eastNorth[,1])
# centerNorth = mean(eastNorth[,2])
# 
# # azimuth is angle from south to west (e.g. 0 is south, pi/2 is west)
# sunPos = getSunlightPosition(date = "2016-08-04", lat = centerLat, lon = centerLon, data = NULL)
# # getSunlightPosition(date = "2016-08-05 01:00:00", lat = centerLat, lon = centerLon, data = NULL)
# centeredEast = eastNorth[,1] - centerEast
# centeredNorth = eastNorth[,2] - centerNorth
# rotationAngle = -(2*pi - (sunPos$azimuth + pi/2)) # shift 0 to eastward and reverse direction from clockwise to counter clockwise
# rotationMat = rbind(c(cos(rotationAngle), -sin(rotationAngle)), 
#                     c(sin(rotationAngle), cos(rotationAngle)))
# sunCoords = t(rotationMat %*% t(cbind(centeredEast, centeredNorth)))
# sunward = sunCoords[,1]

# test plot
if(makePlots) {
  png("Figures/applicationHeaton/sunward.png", width=700, height=700)
  quilt.plot(all.sat.temps$Lon, all.sat.temps$Lat, sunward, 
             col=inferno(64), nx=500, ny=300, 
             xlab="Longitude", ylab="Latitude")
  dev.off()
  
  png("Figures/applicationHeaton/maskedTemp.png", width=700, height=700)
  quilt.plot(all.sat.temps$Lon, all.sat.temps$Lat, all.sat.temps$MaskTemp, 
             col=inferno(64), nx=500, ny=300, 
             xlab="Longitude", ylab="Latitude")
  dev.off()
  
  png("Figures/applicationHeaton/trueTemp.png", width=700, height=700)
  quilt.plot(all.sat.temps$Lon, all.sat.temps$Lat, all.sat.temps$TrueTemp, 
             col=inferno(64), nx=500, ny=300, 
             xlab="Longitude", ylab="Latitude")
  dev.off()
}

# get elevation ----
# temp = raster("elevation.nc")
# elev = extract(temp, SpatialPoints(cbind(all.sat.temps$Lon, all.sat.temps$Lat), proj4string=CRS("+proj=longlat")),method="bilinear")
# GMTED 2010 elevation data downloaded from https://topotools.cr.usgs.gov/gmted_viewer/viewer.htm
elevRaster = raster("gmted2010_30arcsec_US_lon_-90_-120_lat_30_50.tif")
elev = extract(elevRaster, SpatialPoints(cbind(all.sat.temps$Lon, all.sat.temps$Lat), proj4string=CRS("+proj=longlat")),method="bilinear")
# temp2 = raster("gmted2010_30arcsec_US_lon_-90_-120_lat_30_50_breaklineEmphasis.tif")
# breakline = extract(temp2, SpatialPoints(cbind(all.sat.temps$Lon, all.sat.temps$Lat), proj4string=CRS("+proj=longlat")),method="bilinear")
all.sat.temps$elev = elev

if(makePlots) {
  png("Figures/applicationHeaton/elevation.png", width=700, height=700)
  quilt.plot(all.sat.temps$Lon, all.sat.temps$Lat, elev, 
             col=viridis(64), nx=500, ny=300, 
             xlab="Longitude", ylab="Latitude")
  dev.off()
}

# png("Figures/applicationHeaton/breakline.png", width=700, height=700)
# quilt.plot(all.sat.temps$Lon, all.sat.temps$Lat, breakline, 
#            col=viridis(64), nx=500, ny=300, 
#            xlab="Longitude", ylab="Latitude")
# dev.off()

# get ndvi ----
# https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=1299
ndviRaster = raster("MCD13.A2014.unaccum.nc4")
ndvi = extract(ndviRaster, SpatialPoints(cbind(all.sat.temps$Lon, all.sat.temps$Lat), proj4string=CRS("+proj=longlat")),method="bilinear")
all.sat.temps$ndvi = ndvi

if(makePlots) {
  png("Figures/applicationHeaton/ndvi.png", width=700, height=700)
  quilt.plot(all.sat.temps$Lon, all.sat.temps$Lat, ndvi, 
             col=viridis(64), nx=500, ny=300, 
             xlab="Longitude", ylab="Latitude")
  dev.off()
}

# daysInMonth = c(31, 29, 31, 30, 31, 30, 31)
# sum(daysInMonth) # 213, beginning of August
# sum(daysInMonth) + 29 # 242, last day of August
# temp = raster("MOD13Q1.006__250m_16_days_NDVI_doy2016225_aid0001.tif")
# ndvi2 = extract(temp, SpatialPoints(cbind(all.sat.temps$Lon, all.sat.temps$Lat), proj4string=CRS("+proj=longlat")),method="bilinear")
# 
# png("Figures/applicationHeaton/ndvi2.png", width=700, height=700)
# quilt.plot(all.sat.temps$Lon, all.sat.temps$Lat, ndvi2, 
#            col=viridis(64), nx=500, ny=300, 
#            xlab="Longitude", ylab="Latitude")
# dev.off()
# 
# temp = raster("MOD13Q1.006__250m_16_days_EVI_doy2016225_aid0001.tif")
# dvi = extract(temp, SpatialPoints(cbind(all.sat.temps$Lon, all.sat.temps$Lat), proj4string=CRS("+proj=longlat")),method="bilinear")
# 
# png("Figures/applicationHeaton/dvi.png", width=700, height=700)
# quilt.plot(all.sat.temps$Lon, all.sat.temps$Lat, dvi, 
#            col=viridis(64), nx=500, ny=300, 
#            xlab="Longitude", ylab="Latitude")
# dev.off()

# get distance to water ----
library(ncdf4)
# https://catalogue.ceda.ac.uk/uuid/06cef537c5b14a2e871a333b9bc0b482
# nc_data <- nc_open("globolakes-static_distance_to_water_Map-300m-P5Y-2005-ESACCI_WB-fv1.0.nc")
nc_data <- nc_open("distanceFromWater.nc")
# names(nc_data$var)
distanceToWater = ncvar_get(nc_data, "distance_to_water") # 
lonGrid <- nc_data$dim$lon$vals
latGrid <- nc_data$dim$lat$vals
lon = rep(lonGrid, length(latGrid))
lat = rep(latGrid, each=length(lonGrid))
dataExtent = extent(c(min(lonGrid), max(lonGrid), min(latGrid), max(latGrid)))
# temp = rasterFromXYZ(data.frame(x=lon, y=lat, z=c(distanceToWater)))
# temp <- raster(nrows=length(lonGrid), ncols=length(latGrid),
#                xmn=min(lonGrid), xmx=max(lonGrid),
#                ymn=min(latGrid), ymx=max(latGrid),
#                vals=c(distanceToWater), # latitude is backwards, so we must reverse matrix to compensate
#                ext = dataExtent, crs=CRS("+proj=longlat"))
temp <- raster(dataExtent, nrows=length(latGrid), ncols=length(lonGrid))
distToWaterRaster <- rasterize(cbind(lon, lat), temp, distanceToWater, fun=mean)
# plot(distToWaterRaster)

distToWater = extract(distToWaterRaster, SpatialPoints(cbind(all.sat.temps$Lon, all.sat.temps$Lat), proj4string=CRS("+proj=longlat")),method="bilinear")
all.sat.temps$distToWater = distToWater

if(makePlots) {
  png("Figures/applicationHeaton/distanceToWater.png", width=700, height=700)
  quilt.plot(all.sat.temps$Lon, all.sat.temps$Lat, distToWater, 
             col=viridis(64), nx=500, ny=300, 
             xlab="Longitude", ylab="Latitude")
  dev.off()
  
  png("Figures/applicationHeaton/isWater.png", width=700, height=700)
  quilt.plot(all.sat.temps$Lon, all.sat.temps$Lat, distToWater<1, 
             col=rev(viridis(2)), nx=500, ny=300, 
             xlab="Longitude", ylab="Latitude")
  dev.off()
  
  isWater = distanceToWater < 1
  png("Figures/applicationHeaton/isWaterImageTest.png", width=700, height=700)
  image(isWater)
  dev.off()
  
  png("Figures/applicationHeaton/isWaterImageTest2.png", width=700, height=700)
  image(t(isWater))
  dev.off()
  
  png("Figures/applicationHeaton/isWaterImageTest3.png", width=700, height=700)
  image(t(apply(isWater, 1, rev)))
  dev.off()
}

# fit linear models ----
save(all.sat.temps, file="savedOutput/heaton/all.sat.temps-full.rda")
load("savedOutput/heaton/all.sat.temps-full.rda")
out = lm(TrueTemp ~ Lon + Lat + elev + ndvi + distToWater, data=all.sat.temps)
summary(out)
out = lm(TrueTemp ~ Lon + Lat + ndvi + elev, data=all.sat.temps)
summary(out)
out = lm(TrueTemp ~ Lon + Lat + ndvi*elev, data=all.sat.temps)
summary(out)
out = lm(TrueTemp ~ Lon + Lat + elev, data=all.sat.temps)
summary(out)
out = lm(TrueTemp ~ Lon + Lat, data=all.sat.temps)
summary(out)
out = lm(TrueTemp ~ Lon + Lat + ndvi, data=all.sat.temps)
summary(out)
out = lm(TrueTemp ~ ndvi + elev + distToWater, data=all.sat.temps)
summary(out)
out = lm(TrueTemp ~ ndvi*elev, data=all.sat.temps)
summary(out)
out = lm(TrueTemp ~ ndvi*elev + distToWater, data=all.sat.temps)
summary(out)
out = lm(TrueTemp ~ ndvi*elev*distToWater, data=all.sat.temps)
summary(out)

if(makePlots) {
  pdf("Figures/applicationHeaton/pairPlot.pdf", width=7, height=7)
  pairs(data.frame(all.sat.temps$TrueTemp, elev, ndvi, distToWater), 
        pch=".", labels=c("Temperature", "Elevation", "NDVI", "Distance to water"))
  dev.off()
  
  # get residuals of linear model
  masked = is.na(all.sat.temps$MaskTemp)
  predFrame = all.sat.temps[masked,]
  predFrame$ndvi = ndvi[masked]
  predFrame$elev = elev[masked]
  predFrame$distToWater = distToWater[masked]
  fullPredFrame = all.sat.temps
  fullPredFrame$ndvi = ndvi
  fullPredFrame$elev = elev
  fullPredFrame$distToWater = distToWater
  
  out = lm(TrueTemp ~ Lon + Lat, data=all.sat.temps)
  preds = predict(out, fullPredFrame)
  resids = fullPredFrame$TrueTemp - preds
  
  png("Figures/applicationHeaton/linearResidsLonLat.png", width=700, height=700)
  quilt.plot(fullPredFrame$Lon, fullPredFrame$Lat, resids, 
             col=viridis(64), nx=500, ny=300, 
             xlab="Longitude", ylab="Latitude")
  dev.off()
  
  out = lm(TrueTemp ~ ndvi + elev + distToWater, data=all.sat.temps)
  preds = predict(out, fullPredFrame)
  resids = fullPredFrame$TrueTemp - preds
  png("Figures/applicationHeaton/linearResidsNoLonLatAllCovs.png", width=700, height=700)
  quilt.plot(fullPredFrame$Lon, fullPredFrame$Lat, resids, 
             col=viridis(64), nx=500, ny=300, 
             xlab="Longitude", ylab="Latitude")
  dev.off()
  
  out = lm(TrueTemp ~ ndvi*elev*distToWater, data=all.sat.temps)
  preds = predict(out, fullPredFrame)
  resids = fullPredFrame$TrueTemp - preds
  png("Figures/applicationHeaton/linearResidsNoLonLatAllCovsAndInteractions.png", width=700, height=700)
  quilt.plot(fullPredFrame$Lon, fullPredFrame$Lat, resids, 
             col=viridis(64), nx=500, ny=300, 
             xlab="Longitude", ylab="Latitude")
  dev.off()
  
  # plot interaction between ndvi and elev for linear model
  out = lm(TrueTemp ~ ndvi*elev, data=all.sat.temps)
  coefs = coef(out)
  tempPredFrame = make.surface.grid(list(ndvi=seq(min(ndvi), max(ndvi), l=100), elev=seq(min(elev), max(elev), l=100)))
  tempPredFrame = cbind(int=1, tempPredFrame, ndviElev=tempPredFrame[,1]*tempPredFrame[,2])
  preds = tempPredFrame %*% coefs
  png("Figures/applicationHeaton/linearModElevNDVIpredsCovSpaceInteraction.png", width=700, height=700)
  quilt.plot(tempPredFrame[,2:3], preds, 
             col=viridis(64), nx=100, ny=100, 
             xlab="NDVI", ylab="Elevation (m)")
  dev.off()
  
  tempPredFrame = cbind(1, all.sat.temps$ndvi, all.sat.temps$elev, all.sat.temps$ndvi*all.sat.temps$elev)
  preds = tempPredFrame %*% coefs
  png("Figures/applicationHeaton/linearModElevNDVIpredsLonLatInteraction.png", width=700, height=700)
  quilt.plot(all.sat.temps$Lon, all.sat.temps$Lat, preds, 
             col=viridis(64), nx=500, ny=300, 
             xlab="Longitude", ylab="Latitude")
  dev.off()
  
  # plot interaction (when it's not included) between ndvi and elev for linear model
  out = lm(TrueTemp ~ ndvi+elev, data=all.sat.temps)
  coefs = coef(out)
  tempPredFrame = make.surface.grid(list(ndvi=seq(min(ndvi), max(ndvi), l=100), elev=seq(min(elev), max(elev), l=100)))
  tempPredFrame = cbind(int=1, tempPredFrame)
  preds = tempPredFrame %*% coefs
  png("Figures/applicationHeaton/linearModElevNDVIpredsCovSpace.png", width=700, height=700)
  quilt.plot(tempPredFrame[,2:3], preds, 
             col=viridis(64), nx=100, ny=100, 
             xlab="NDVI", ylab="Elevation (m)")
  dev.off()
  
  tempPredFrame = cbind(1, all.sat.temps$ndvi, all.sat.temps$elev)
  preds = tempPredFrame %*% coefs
  png("Figures/applicationHeaton/linearModElevNDVIpredsLonLat.png", width=700, height=700)
  quilt.plot(all.sat.temps$Lon, all.sat.temps$Lat, preds, 
             col=viridis(64), nx=500, ny=300, 
             xlab="Longitude", ylab="Latitude")
  dev.off()
}

##### try fitting LK ----
# * Setup ----
dataCase = "full"
dataCase = "test"
dataCase = "full"
dataCaseText = ifelse(dataCase == "test", "_test", "")
if(dataCase=="test"){
  # load("heaton/heatoncomparison-master/Data/SmallTestData.RData")
  # lon <- code.test$Lon
  # lat <- code.test$Lat
  # temps <- code.test$MaskedData
  # original testing parameters
  # NC<- 16
  # nlevel<- 2
  # a.wght<- 4.4
  # nu<- .5
  NC<- 20
  nlevel<- 2
  a.wght<-  10.25
  nu<- .1
  separateRanges=FALSE
  sampleN = 100000
} else {
  out = load("heaton/AllSatelliteTemps.RData")
  names(all.sat.temps)
  # original parameters:
  # NC<- 40
  # nlevel<- 4
  # a.wght<-  10.25
  # nu<- .1
  NC<- 40
  nlevel<- 3
  a.wght<-  10.25
  nu<- .1
  separateRanges=FALSE
  sampleN = sum(!is.na(all.sat.temps$MaskTemp))
}
kappa = sqrt(a.wght - 4)

# create locations and responses (y)
source("heaton/heatoncomparison-master/Code/LatticeKrig/functions.R")
if(dataCase == "test") {
  dataObject<- makeData(sat.temps$Lon,sat.temps$Lat, sat.temps$Temp)
  sampleI = sample(1:nrow(dataObject$x), sampleN)
  thisDataObject = dataObject
  thisDataObject$x = thisDataObject$x[sampleI,]
  thisDataObject$y = thisDataObject$y[sampleI]
  thisDataObjectI = which(!is.na(all.sat.temps$TrueTemp))[sampleI]
  thisDataObjectMissingI = which(is.na(all.sat.temps$MaskTemp))
} else {
  dataObject<- makeData(sat.temps$Lon,sat.temps$Lat, sat.temps$Temp)
  dataObject$yMissing = all.sat.temps$TrueTemp[is.na(all.sat.temps$MaskTemp)]
  thisDataObject = dataObject
  thisDataObjectI = which(!is.na(all.sat.temps$TrueTemp))
  thisDataObjectMissingI = which(is.na(all.sat.temps$MaskTemp))
  separateRanges=FALSE
}
thisDataObject$TrueTemp = all.sat.temps$TrueTemp[thisDataObjectI]
thisDataObject$TrueTempMissing = all.sat.temps$TrueTemp[thisDataObjectMissingI]
thisDataObject$TrueTempMissingMasked = all.sat.temps$MaskTemp[thisDataObjectMissingI]

elev = extract(elevRaster, SpatialPoints(thisDataObject$x, proj4string=CRS("+proj=longlat")),method="bilinear")
elevPred = extract(elevRaster, SpatialPoints(thisDataObject$xMissing, proj4string=CRS("+proj=longlat")),method="bilinear")
thisDataObject$elev = elev
thisDataObject$elevPred = elevPred

# setup LKrig object
LKinfoFinal<- LKrigSetup( thisDataObject$x,
                          NC = NC,
                          nlevel = nlevel,
                          a.wght = a.wght,
                          nu = nu)

latInfo = LKinfo2ELKBasis(LKinfoFinal)

NCtext = do.call("paste", c(list("_NC"), lapply(latInfo, function(x){x$NC}), sep="_"))
thisLatticeNameRoot = NCtext
thisFileNameRoot = paste0("_N", sampleN, NCtext, "_sepRanges", separateRanges)

# * Fit model ----
require(spam64)
comp.timeLKfit <- system.time(fitLK<- LatticeKrig(thisDataObject$x, thisDataObject$y, 
                                                     LKinfo=LKinfoFinal, Z=thisDataObject$elev, 
                                                     verbose=TRUE))

# * Get SEs ----
set.seed(234)

M <- 100  # Number of conditional draws (used for standard errors) [Final model: 100]
comp.timeLKse = 
  system.time(outputSim<- 
                LKrig.sim.conditional(fitLK,
                                      x.grid = thisDataObject$xMissing, 
                                      Z.grid = matrix(thisDataObject$elevPred, ncol=1), 
                                      M = M))

# * Get CIs ----
comp.timeLKci <- system.time(standardError<- sqrt( 
  apply( outputSim$g.draw, 1, "var") +
    fitLK$sigma.MLE^2))


yHat<- outputSim$ghat
CI95Lower<- yHat + qnorm(.025) * standardError
CI95Upper<- yHat + qnorm(.975) * standardError

comp.timeLK = comp.timeLKfit + comp.timeLKse + comp.timeLKci

finalResults<- list(x=thisDataObject$xMissing,
                    yHat=yHat, 
                    standError=standardError)

# * get scoring rules ----
scoresLK = getScores(datPred$TrueTemp, est=finalResults$yHat, var=finalResults$standError^2, 
                     lower=NULL, upper=NULL, estMat=NULL, significance=.95, 
                     distances=NULL, breaks=30, doFuzzyReject=FALSE, getAverage=FALSE)

# calculate distance to nearest observation
if(nrow(thisDataObject$x) > 50000) {
  mean.neighbor = nrow(thisDataObject$x)/(5*2.5) *pi * .8
  dists = fields.rdist.near(cbind(datPred$Lon, datPred$Lat), thisDataObject$x, 1, mean.neighbor=mean.neighbor)
} else {
  dists = rdist(cbind(datPred$Lon, datPred$Lat), thisDataObject$x)
}
nndists = apply(dists, 1, min)
scoresLK = cbind(NNDist=nndists, scoresLK)

# * Save results ----

save( finalResults, fitLK, comp.timeLK, comp.timeLKfit, comp.timeLKse, comp.timeLKci, scoresLK, 
      file=paste0("savedOutput/heaton/resultsLK", thisFileNameRoot, ".rda") )
load(paste0("savedOutput/heaton/resultsLK", thisFileNameRoot, ".rda")) 

# * Plot results ----
xrange = range(all.sat.temps$Lon)
yrange = range(all.sat.temps$Lat)
pdf(paste0("Figures/applicationHeaton/LKbasicPreds", thisFileNameRoot, ".pdf"), width=5, height=5)
quilt.plot(cbind(thisDataObject$xMissing[,1], thisDataObject$xMissing[,2]), finalResults$yHat, nx=500, ny=300, 
           xlim=xrange, ylim=yrange, main="LK predictions (linear covariate model)", 
           xlab="Longitude", ylab="Latitude")
dev.off()

pdf(paste0("Figures/applicationHeaton/LKbasicSDs", thisFileNameRoot, ".pdf"), width=5, height=5)
quilt.plot(cbind(thisDataObject$xMissing[,1], thisDataObject$xMissing[,2]), finalResults$standError, nx=500, ny=300, 
           xlim=xrange, ylim=yrange, main="LK SD (linear covariate model)", 
           xlab="Longitude", ylab="Latitude")
dev.off()

pdf(paste0("Figures/applicationHeaton/LKbasicResids", thisFileNameRoot, ".pdf"), width=5, height=5)
plot(datPred$elev, finalResults$yHat-thisDataObject$TrueTempMissing, 
     xlab="Elevation (m)", ylab="LK residuals", 
     main="LK residuals Vs. elevation", 
     pch=19, cex=.1, col="blue")
dev.off()

##### try fitting ELK ----

# * Setup ----
latInfo = LKinfo2ELKBasis(LKinfoFinal)
Ms = getMs(latInfo=latInfo)
Ms
sapply(latInfo, function(x) {x$latWidth})

priorPar = getPCPrior(max(c(latInfo[[1]]$xRangeDat, latInfo[[1]]$yRangeDat))/5, .01, 1, nLayer=nlevel, separateRanges=FALSE, latticeInfo=latInfo, useUrbanPrior=FALSE) # 37.06811/5

# * Do precomputations ----

precomputationFileNameRoot = paste0("precomputationResultsHeatonComparison", NCtext)
if(!file.exists(paste0("savedOutput/precomputations/", precomputationFileNameRoot, ".RData"))) {
  comp.timeELKprecomputation = system.time({
    precomputedMatrices = precomputationsQ2(latInfo)
    precomputedNormalizationFun = precomputeNormalization(saveResults=FALSE, latticeInfo=latInfo, effRangeRange=NULL, 
                                                          plotNormalizationSplines=FALSE)
  })
  
  save(precomputedMatrices, precomputedNormalizationFun, comp.timeELKprecomputation, 
       file=paste0("savedOutput/precomputations/", precomputationFileNameRoot, ".RData"))
} else {
  load(paste0("savedOutput/precomputations/", precomputationFileNameRoot, ".RData"))
}

# * Fit model ----
dat = data.frame(Lon=thisDataObject$x[,1], Lat=thisDataObject$x[,2], MaskedTemp=thisDataObject$y, 
                 elev=extract(elevRaster, SpatialPoints(cbind(thisDataObject$x[,1], thisDataObject$x[,2]), proj4string=CRS("+proj=longlat")),method="bilinear"), 
                 distToWater=extract(distToWaterRaster, SpatialPoints(cbind(thisDataObject$x[,1], thisDataObject$x[,2]), proj4string=CRS("+proj=longlat")),method="bilinear"), 
                 ndvi=extract(ndviRaster, SpatialPoints(cbind(thisDataObject$x[,1], thisDataObject$x[,2]), proj4string=CRS("+proj=longlat")),method="bilinear"))
datPred = data.frame(Lon=thisDataObject$xMissing[,1], Lat=thisDataObject$xMissing[,2], MaskedTemp=thisDataObject$TrueTempMissingMasked, TrueTemp=thisDataObject$TrueTempMissing, 
                 elev=extract(elevRaster, SpatialPoints(cbind(thisDataObject$xMissing[,1], thisDataObject$xMissing[,2]), proj4string=CRS("+proj=longlat")),method="bilinear"), 
                 distToWater=extract(distToWaterRaster, SpatialPoints(cbind(thisDataObject$xMissing[,1], thisDataObject$xMissing[,2]), proj4string=CRS("+proj=longlat")),method="bilinear"), 
                 ndvi=extract(ndviRaster, SpatialPoints(cbind(thisDataObject$xMissing[,1], thisDataObject$xMissing[,2]), proj4string=CRS("+proj=longlat")),method="bilinear"))
dat$Covar = dat$elev
datPred$Covar = datPred$elev
thisDiagonal = 0.0
if(sampleN >= 100000) {
  thisDiagonal = c(10.0, 1.0, 0.0)
}
comp.timeELKfit = system.time(fitELK <- fitLKINLAStandard2(thisDataObject$x, thisDataObject$y, 
                                                    predCoords=cbind(datPred$Lon, datPred$Lat), 
                                                    xObs=cbind(1, dat$Lon, dat$Lat, dat$elev, dat$ndvi, dat$elev*dat$ndvi), 
                                                    xPred=cbind(1, datPred$Lon, datPred$Lat, datPred$elev, datPred$ndvi, datPred$elev*datPred$ndvi), 
                                                    nu=nu, seed=234, nLayer=nlevel, NC=NC,
                                                    nBuffer=nBuffer, priorPar=priorPar, normalize=TRUE, 
                                                    intStrategy="eb", strategy="gaussian", fastNormalize=TRUE, 
                                                    predictionType=c("mean", "median"), significanceCI=0.8, 
                                                    printVerboseTimings=FALSE, nPostSamples=1000, family="normal",
                                                    clusterEffect=TRUE, latInfo=latInfo, 
                                                    initialEffectiveRange=3, 
                                                    verbose=TRUE, separateRanges=separateRanges, 
                                                    loadPrecomputationResults=TRUE, 
                                                    precomputationFileNameRoot=precomputationFileNameRoot, 
                                                    diagonal=thisDiagonal))

comp.timeELKall = comp.timeELKfit + comp.timeELKprecomputation

# * get scoring rules ----
scoresELK = getScores(datPred$TrueTemp, est=fitELK$preds, var=fitELK$sigmas^2, lower=NULL, upper=NULL, estMat=NULL, significance=.95, 
                      distances=NULL, breaks=30, doFuzzyReject=FALSE, getAverage=FALSE)

# calculate distance to nearest observation
if(!exists(nndists)) {
  if(nrow(thisDataObject$x) > 50000) {
    mean.neighbor = nrow(thisDataObject$x)/(5*2.5) *pi * .8
    dists = fields.rdist.near(cbind(datPred$Lon, datPred$Lat), thisDataObject$x, 1, mean.neighbor=mean.neighbor)
  } else {
    dists = rdist(cbind(datPred$Lon, datPred$Lat), thisDataObject$x)
  }
  nndists = apply(dists, 1, min)
}

scoresELK = cbind(NNDist=nndists, scoresELK)

# * save results ----
fitELK$mod$.args = NULL
fitELKfinalInt$mod$all.hyper = NULL
save(fitELK, comp.timeELKall, comp.timeELKfit, comp.timeELKprecomputation, scoresELK, 
     file=paste0("savedOutput/heaton/resultsELK", thisFileNameRoot, ".rda"))
out = load(paste0("savedOutput/heaton/resultsELK", thisFileNameRoot, ".rda"))

# * Plot results ----
# replot LK results with same scales
xrange = range(all.sat.temps$Lon)
yrange = range(all.sat.temps$Lat)

# LK plots
pdf(paste0("Figures/applicationHeaton/LKbasicPreds", thisFileNameRoot, ".pdf"), width=5, height=5)
zlim = range(c(finalResults$yHat, fitELK$preds))
quilt.plot(cbind(thisDataObject$xMissing[,1], thisDataObject$xMissing[,2]), finalResults$yHat, nx=500, ny=300, 
           xlim=xrange, ylim=yrange, main="LK predictions (linear covariate model)", 
           xlab="Longitude", ylab="Latitude", zlim=zlim)
dev.off()

pdf(paste0("Figures/applicationHeaton/LKbasicSDs", thisFileNameRoot, ".pdf"), width=5, height=5)
sds = rowMeans(outer(fitELK$sigmasNoNugget^2, fitELK$clusterVars, function(x, y) {sqrt(x + y)}))
zlim = range(c(finalResults$standError, sds))
quilt.plot(cbind(thisDataObject$xMissing[,1], thisDataObject$xMissing[,2]), finalResults$standError, nx=500, ny=300, 
           xlim=xrange, ylim=yrange, main="LK SD (linear covariate model)", 
           xlab="Longitude", ylab="Latitude", zlim=zlim)
dev.off()

pdf(paste0("Figures/applicationHeaton/LKbasicResids", thisFileNameRoot, ".pdf"), width=5, height=5)
ylim = range(c(finalResults$yHat-thisDataObject$TrueTempMissing, fitELK$preds-datPred$TrueTemp), na.rm=TRUE)
plot(datPred$elev, finalResults$yHat-thisDataObject$TrueTempMissing, 
     xlab="Elevation (m)", ylab="LK residuals", 
     main="LK residuals Vs. elevation", 
     pch=19, cex=.1, col="blue", ylim=ylim)
dev.off()

# ELK plots
pdf(paste0("Figures/applicationHeaton/ELKbasicPredsNDVIxELEV", thisFileNameRoot, ".pdf"), width=5, height=5)
zlim = range(c(finalResults$yHat, fitELK$preds))
quilt.plot(cbind(datPred$Lon, datPred$Lat), fitELK$preds, nx=500, ny=300, 
           xlim=xrange, ylim=yrange, main="ELK predictions (linear covariate model)", 
           xlab="Longitude", ylab="Latitude", zlim=zlim)
dev.off()

pdf(paste0("Figures/applicationHeaton/ELKbasicSDsNDVIxELEV", thisFileNameRoot, ".pdf"), width=5, height=5)
sds = rowMeans(outer(fitELK$sigmasNoNugget^2, fitELK$clusterVars, function(x, y) {sqrt(x + y)}))
zlim = range(c(finalResults$standError, sds))
quilt.plot(cbind(datPred$Lon, datPred$Lat), sds, nx=500, ny=300, 
           xlim=xrange, ylim=yrange, main="ELK SD (linear covariate model)", 
           xlab="Longitude", ylab="Latitude", zlim=zlim)
dev.off()

pdf(paste0("Figures/applicationHeaton/ELKbasicResidsNDVIxELEV", thisFileNameRoot, ".pdf"), width=5, height=5)
ylim = range(c(finalResults$yHat-thisDataObject$TrueTempMissing, fitELK$preds-datPred$TrueTemp), na.rm=TRUE)
plot(datPred$elev, fitELK$preds-datPred$TrueTemp, 
     xlab="Elevation (m)", ylab="LK residuals", 
     main="ELK residuals Vs. elevation", 
     pch=19, cex=.1, col="blue", ylim=ylim)
dev.off()

# calculate distance to nearest observation
if(!exists(nndists)) {
  if(nrow(thisDataObject$x) > 50000) {
    mean.neighbor = nrow(thisDataObject$x)/(5*2.5) *pi * .8
    dists = fields.rdist.near(cbind(datPred$Lon, datPred$Lat), thisDataObject$x, 1, mean.neighbor=mean.neighbor)
  } else {
    dists = rdist(cbind(datPred$Lon, datPred$Lat), thisDataObject$x)
  }
  nndists = apply(dists, 1, min)
}

pdf(paste0("Figures/applicationHeaton/LK-ELKbasicPredsComparison", thisFileNameRoot, ".pdf"), width=5, height=5)
plotWithColor(finalResults$yHat, fitELK$preds, nndists, colScale=yellowBlueCols, 
              xlab="LatticeKrig", ylab="ELK", 
              main="Prediction comparison (linear covariate models)", 
              pch=19, cex=.2, ordering="increasing")
dev.off()

pdf(paste0("Figures/applicationHeaton/LK-ELKbasicSDsComparison", thisFileNameRoot, ".pdf"), width=5, height=5)
plotWithColor(finalResults$standError, sds, nndists, colScale=yellowBlueCols, 
              xlab="LatticeKrig", ylab="ELK", 
              main="SD comparison (linear covariate models)", 
              pch=19, cex=.2, ordering="increasing")
abline(0, 1, lty=2)
dev.off()

##### try fitting ELK with nonlinear effects ----

# * Fit model ----
# 
# comp.timeELKnonlinearFit = system.time(fitELKnonlinear <- modHeaton(dat, 
#                                                            predCoords=cbind(datPred$Lon, datPred$Lat), 
#                                                            predCovar=datPred$ndvi, 
#                                                            nu=nu, seed=234, nLayer=nlevel, NC=NC,
#                                                            nBuffer=nBuffer, priorPar=priorPar, normalize=TRUE, 
#                                                            intStrategy="eb", strategy="gaussian", fastNormalize=TRUE, 
#                                                            predictionType=c("mean", "median"), significanceCI=0.8, 
#                                                            printVerboseTimings=FALSE, nPostSamples=1000, family="normal",
#                                                            clusterEffect=TRUE, latInfo=latInfo, 
#                                                            initialEffectiveRange=3, 
#                                                            verbose=TRUE, separateRanges=FALSE, 
#                                                            loadPrecomputationResults=TRUE, 
#                                                            precomputationFileNameRoot=precomputationFileNameRoot, 
#                                                            diagonal=c(0.0), rwConstr=TRUE))
# 
# comp.timeELKnonlinearAll = comp.timeELKnonlinearFit + comp.timeELKprecomputation
# 
# # * save results ----
# fitELKnonlinear$mod$.args = NULL
# save(fitELKnonlinear, comp.timeELKnonlinearAll, comp.timeELKnonlinearFit, comp.timeELKprecomputation,
#      file=paste0("savedOutput/heaton/resultsELKnonlinear", thisFileNameRoot, ".rda"))
# out = load(paste0("savedOutput/heaton/resultsELKnonlinear", thisFileNameRoot, ".rda"))
# 
# # * Plot results ----
# # replot LK results with same scales
# xrange = range(all.sat.temps$Lon)
# yrange = range(all.sat.temps$Lat)
# 
# # LK plots
# pdf(paste0("Figures/applicationHeaton/LKbasicPreds", thisFileNameRoot, ".pdf"), width=5, height=5)
# zlim = range(c(finalResults$yHat, fitELK$preds, fitELKnonlinear$preds))
# quilt.plot(cbind(thisDataObject$xMissing[,1], thisDataObject$xMissing[,2]), finalResults$yHat, nx=500, ny=300, 
#            xlim=xrange, ylim=yrange, main="LK predictions (linear covariate model)", 
#            xlab="Longitude", ylab="Latitude", zlim=zlim)
# dev.off()
# 
# pdf(paste0("Figures/applicationHeaton/LKbasicSDs", thisFileNameRoot, ".pdf"), width=5, height=5)
# sds = rowMeans(outer(fitELK$sigmasNoNugget^2, fitELK$clusterVars, function(x, y) {sqrt(x + y)}))
# sdsNonlinear = rowMeans(outer(fitELKnonlinear$sigmasNoNugget^2, fitELKnonlinear$clusterVars, function(x, y) {sqrt(x + y)}))
# zlim = range(c(finalResults$standError, sds, sdsNonlinear))
# quilt.plot(cbind(thisDataObject$xMissing[,1], thisDataObject$xMissing[,2]), finalResults$standError, nx=500, ny=300, 
#            xlim=xrange, ylim=yrange, main="LK SD (linear covariate model)", 
#            xlab="Longitude", ylab="Latitude", zlim=zlim)
# dev.off()
# 
# pdf(paste0("Figures/applicationHeaton/LKbasicResids", thisFileNameRoot, ".pdf"), width=5, height=5)
# ylim = range(c(finalResults$yHat-thisDataObject$TrueTempMissing, fitELK$preds-datPred$TrueTemp, fitELKnonlinear$preds-datPred$TrueTemp), na.rm=TRUE)
# plot(datPred$elev, 
#      finalResults$yHat-thisDataObject$TrueTempMissing, 
#      xlab="Elevation (m)", ylab="LK residuals", 
#      main="LK residuals Vs. elevation", 
#      pch=19, cex=.1, col="blue", ylim=ylim)
# dev.off()
# 
# # ELK plots
# pdf(paste0("Figures/applicationHeaton/ELKbasicPreds", thisFileNameRoot, ".pdf"), width=5, height=5)
# zlim = range(c(finalResults$yHat, fitELK$preds, fitELKnonlinear$preds))
# quilt.plot(cbind(datPred$Lon, datPred$Lat), fitELK$preds, nx=500, ny=300, 
#            xlim=xrange, ylim=yrange, main="ELK predictions (linear covariate model)", 
#            xlab="Longitude", ylab="Latitude", zlim=zlim)
# dev.off()
# 
# pdf(paste0("Figures/applicationHeaton/ELKbasicSDs", thisFileNameRoot, ".pdf"), width=5, height=5)
# sds = rowMeans(outer(fitELK$sigmasNoNugget^2, fitELK$clusterVars, function(x, y) {sqrt(x + y)}))
# sdsNonlinear = rowMeans(outer(fitELKnonlinear$sigmasNoNugget^2, fitELKnonlinear$clusterVars, function(x, y) {sqrt(x + y)}))
# zlim = range(c(finalResults$standError, sds, sdsNonlinear))
# quilt.plot(cbind(datPred$Lon, datPred$Lat), sds, nx=500, ny=300, 
#            xlim=xrange, ylim=yrange, main="ELK SD (linear covariate model)", 
#            xlab="Longitude", ylab="Latitude", zlim=zlim)
# dev.off()
# 
# pdf(paste0("Figures/applicationHeaton/ELKbasicResids", thisFileNameRoot, ".pdf"), width=5, height=5)
# ylim = range(c(finalResults$yHat-thisDataObject$TrueTempMissing, fitELK$preds-datPred$TrueTemp, fitELKnonlinear$preds-datPred$TrueTemp), na.rm=TRUE)
# plot(datPred$elev, fitELK$preds-datPred$TrueTemp, 
#      xlab="Elevation (m)", ylab="ELK residuals", 
#      main="ELK residuals Vs. elevation (linear covariate)", 
#      pch=19, cex=.1, col="blue", ylim=ylim)
# dev.off()
# 
# # ELK nonlinear plots
# pdf(paste0("Figures/applicationHeaton/ELKnonlinearPreds", thisFileNameRoot, ".pdf"), width=5, height=5)
# zlim = range(c(finalResults$yHat, fitELK$preds, fitELKnonlinear$preds))
# quilt.plot(cbind(datPred$Lon, datPred$Lat), fitELKnonlinear$preds, nx=500, ny=300, 
#            xlim=xrange, ylim=yrange, main="ELK predictions (nonlinear covariate model)", 
#            xlab="Longitude", ylab="Latitude", zlim=zlim)
# dev.off()
# 
# pdf(paste0("Figures/applicationHeaton/ELKnonlinearSDs", thisFileNameRoot, ".pdf"), width=5, height=5)
# sds = rowMeans(outer(fitELK$sigmasNoNugget^2, fitELK$clusterVars, function(x, y) {sqrt(x + y)}))
# sdsNonlinear = rowMeans(outer(fitELKnonlinear$sigmasNoNugget^2, fitELKnonlinear$clusterVars, function(x, y) {sqrt(x + y)}))
# zlim = range(c(finalResults$standError, sds, sdsNonlinear))
# quilt.plot(cbind(datPred$Lon, datPred$Lat), sdsNonlinear, nx=500, ny=300, 
#            xlim=xrange, ylim=yrange, main="ELK SD (nonlinear covariate model)", 
#            xlab="Longitude", ylab="Latitude", zlim=zlim)
# dev.off()
# 
# pdf(paste0("Figures/applicationHeaton/ELKnonlinearResids", thisFileNameRoot, ".pdf"), width=5, height=5)
# ylim = range(c(finalResults$yHat-thisDataObject$TrueTempMissing, fitELK$preds-datPred$TrueTemp, fitELKnonlinear$preds-datPred$TrueTemp), na.rm=TRUE)
# plot(datPred$elev, fitELKnonlinear$preds-datPred$TrueTemp, 
#      xlab="Elevation (m)", ylab="ELK residuals", 
#      main="ELK residuals Vs. elevation (nonlinear covariate)", 
#      pch=19, cex=.1, col="blue", ylim=ylim)
# dev.off()
# 
# pdf(paste0("Figures/applicationHeaton/ELKnonlinearCovar", thisFileNameRoot, ".pdf"), width=5, height=5)
# xs = fitELKnonlinear$rwSummary$ID
# ys = xs*fitELKnonlinear$fixedEffectSummary$mean[4] + fitELKnonlinear$rwSummary$mean
# plot(xs, ys, 
#      xlab="NDVI", ylab="Effect (degrees F)", 
#      main="Estimated effect of NDVI on temperature", 
#      pch=19, cex=1, col="blue")
# dev.off()
# 
# # calculate distance to nearest observation
# # dists = rdist(cbind(datPred$Lon, datPred$Lat), thisDataObject$x)
# # nndists = apply(dists, 1, min)
# 
# pdf(paste0("Figures/applicationHeaton/LK-ELKbasicPredsComparison", thisFileNameRoot, ".pdf"), width=5, height=5)
# lims = range(c(finalResults$yHat, fitELK$preds, fitELKnonlinear$preds))
# plotWithColor(finalResults$yHat, fitELK$preds, nndists, colScale=yellowBlueCols, 
#               xlab="LatticeKrig", ylab="ELK", xlim=lims, ylim=lims, 
#               main="Prediction comparison (linear covariate models)", 
#               pch=19, cex=.2, ordering="increasing")
# dev.off()
# 
# pdf(paste0("Figures/applicationHeaton/LK-ELKbasicSDsComparison", thisFileNameRoot, ".pdf"), width=5, height=5)
# plotWithColor(finalResults$standError, sds, nndists, colScale=yellowBlueCols, 
#               xlab="LatticeKrig", ylab="ELK", 
#               main="SD comparison (linear covariate models)", 
#               pch=19, cex=.2, ordering="increasing")
# abline(0, 1, lty=2)
# dev.off()
# 
# # * Provide Tables ----
# tabLK = 
#   fitLK$MLE

# try fitting ELK with rw2d covariate interaction ----

# * fit model ----
comp.timeELKfitFinalInt = system.time(fitELKfinalInt <- fitLKINLAStandard2(thisDataObject$x, thisDataObject$y, 
                                                                   predCoords=cbind(datPred$Lon, datPred$Lat), 
                                                                   xObs=cbind(1, dat$Lon, dat$Lat, dat$elev, dat$ndvi, dat$elev*dat$ndvi), 
                                                                   xPred=cbind(1, datPred$Lon, datPred$Lat, datPred$elev, datPred$ndvi, datPred$elev*datPred$ndvi), 
                                                                   nonlinearCovariateInds=c(), 
                                                                   nonlinearCovariateInteractionInds=c(4, 5), 
                                                                   nu=nu, seed=234, nLayer=nlevel, NC=NC,
                                                                   nBuffer=nBuffer, priorPar=priorPar, normalize=TRUE, 
                                                                   intStrategy="eb", strategy="gaussian", fastNormalize=TRUE, 
                                                                   predictionType=c("mean", "median"), significanceCI=0.8, 
                                                                   printVerboseTimings=FALSE, nPostSamples=1000, family="normal",
                                                                   clusterEffect=TRUE, latInfo=latInfo, 
                                                                   initialEffectiveRange=3, 
                                                                   verbose=TRUE, separateRanges=separateRanges, 
                                                                   loadPrecomputationResults=TRUE, 
                                                                   precomputationFileNameRoot=precomputationFileNameRoot, 
                                                                   diagonal=thisDiagonal))

comp.timeELKfitFinalIntAll = comp.timeELKfitFinalInt + comp.timeELKprecomputation

# * get scoring rules ----
scoresELKfinalInt = getScores(datPred$TrueTemp, est=fitELKfinalInt$preds, var=fitELKfinalInt$sigmas^2, lower=NULL, upper=NULL, estMat=NULL, significance=.95, 
                      distances=NULL, breaks=30, doFuzzyReject=FALSE, getAverage=FALSE)

# calculate distance to nearest observation
if(!exists(nndists)) {
  if(nrow(thisDataObject$x) > 50000) {
    mean.neighbor = nrow(thisDataObject$x)/(5*2.5) *pi * .8
    dists = fields.rdist.near(cbind(datPred$Lon, datPred$Lat), thisDataObject$x, 1, mean.neighbor=mean.neighbor)
  } else {
    dists = rdist(cbind(datPred$Lon, datPred$Lat), thisDataObject$x)
  }
  nndists = apply(dists, 1, min)
}
scoresELK = cbind(NNDist=nndists, scoresELKfinalInt)

# * save results ----
fitELKfinalInt$mod$.args = NULL
fitELKfinalInt$mod$all.hyper = NULL
save(fitELKfinalInt, comp.timeELKfitFinalIntAll, comp.timeELKfitFinalInt, comp.timeELKprecomputation,
     file=paste0("savedOutput/heaton/resultsELKfinalInt", thisFileNameRoot, ".rda"))

# * plot results ----
# replot LK results with same scales
xrange = range(all.sat.temps$Lon)
yrange = range(all.sat.temps$Lat)

# LK plots
pdf(paste0("Figures/applicationHeaton/LKbasicPreds", thisFileNameRoot, ".pdf"), width=5, height=5)
zlim = range(c(finalResults$yHat, fitELK$preds, fitELKfinalInt$preds))
quilt.plot(cbind(thisDataObject$xMissing[,1], thisDataObject$xMissing[,2]), finalResults$yHat, nx=500, ny=300, 
           xlim=xrange, ylim=yrange, main="LK predictions (linear covariate model)", 
           xlab="Longitude", ylab="Latitude", zlim=zlim)
dev.off()

pdf(paste0("Figures/applicationHeaton/LKbasicSDs", thisFileNameRoot, ".pdf"), width=5, height=5)
sds = rowMeans(outer(fitELK$sigmasNoNugget^2, fitELK$clusterVars, function(x, y) {sqrt(x + y)}))
sdsFinalInt = rowMeans(outer(fitELKfinalInt$sigmasNoNugget^2, fitELKfinalInt$clusterVars, function(x, y) {sqrt(x + y)}))
zlim = range(c(finalResults$standError, sds, sdsFinalInt))
quilt.plot(cbind(thisDataObject$xMissing[,1], thisDataObject$xMissing[,2]), finalResults$standError, nx=500, ny=300, 
           xlim=xrange, ylim=yrange, main="LK SD (linear covariate model)", 
           xlab="Longitude", ylab="Latitude", zlim=zlim)
dev.off()

pdf(paste0("Figures/applicationHeaton/LKbasicResids", thisFileNameRoot, ".pdf"), width=5, height=5)
ylim = range(c(finalResults$yHat-thisDataObject$TrueTempMissing, 
               fitELK$preds-datPred$TrueTemp, 
               fitELKfinalInt$preds-datPred$TrueTemp), na.rm=TRUE)
plot(datPred$elev, finalResults$yHat-thisDataObject$TrueTempMissing, 
     xlab="Elevation (m)", ylab="LK residuals", 
     main="LK residuals Vs. elevation", 
     pch=19, cex=.1, col="blue", ylim=ylim)
dev.off()

# ELK plots
pdf(paste0("Figures/applicationHeaton/ELKbasicPredsNDVIxELEV", thisFileNameRoot, ".pdf"), width=5, height=5)
zlim = range(c(finalResults$yHat, fitELK$preds, fitELKfinalInt$preds))
quilt.plot(cbind(datPred$Lon, datPred$Lat), fitELK$preds, nx=500, ny=300, 
           xlim=xrange, ylim=yrange, main="ELK predictions (linear covariate model)", 
           xlab="Longitude", ylab="Latitude", zlim=zlim)
dev.off()

pdf(paste0("Figures/applicationHeaton/ELKbasicSDsNDVIxELEV", thisFileNameRoot, ".pdf"), width=5, height=5)
sds = rowMeans(outer(fitELK$sigmasNoNugget^2, fitELK$clusterVars, function(x, y) {sqrt(x + y)}))
sdsFinalInt = rowMeans(outer(fitELKfinalInt$sigmasNoNugget^2, fitELKfinalInt$clusterVars, function(x, y) {sqrt(x + y)}))
zlim = range(c(finalResults$standError, sds, sdsFinalInt))
quilt.plot(cbind(datPred$Lon, datPred$Lat), sds, nx=500, ny=300, 
           xlim=xrange, ylim=yrange, main="ELK SD (linear covariate model)", 
           xlab="Longitude", ylab="Latitude", zlim=zlim)
dev.off()

pdf(paste0("Figures/applicationHeaton/ELKbasicResidsNDVIxELEV", thisFileNameRoot, ".pdf"), width=5, height=5)
ylim = range(c(finalResults$yHat-thisDataObject$TrueTempMissing, 
               fitELK$preds-datPred$TrueTemp, 
               fitELKfinalInt$preds-datPred$TrueTemp), na.rm=TRUE)
plot(datPred$elev, fitELK$preds-datPred$TrueTemp, 
     xlab="Elevation (m)", ylab="LK residuals", 
     main="ELK residuals Vs. elevation", 
     pch=19, cex=.1, col="blue", ylim=ylim)
dev.off()

pdf(paste0("Figures/applicationHeaton/LK-ELKbasicPredsComparison", thisFileNameRoot, ".pdf"), width=5, height=5)
predLim = range(c(finalResults$yHat, fitELK$preds, fitELKfinalInt$preds))
plotWithColor(finalResults$yHat, fitELK$preds, nndists, colScale=yellowBlueCols, 
              xlab="LatticeKrig", ylab="ELK", xlim=predLim, ylim=predLim, 
              main="Prediction comparison", 
              pch=19, cex=.2, ordering="increasing")
dev.off()

pdf(paste0("Figures/applicationHeaton/LK-ELKfinalIntPredsComparison", thisFileNameRoot, ".pdf"), width=5, height=5)
predLim = range(c(finalResults$yHat, fitELK$preds, fitELKfinalInt$preds))
plotWithColor(finalResults$yHat, fitELKfinalInt$preds, nndists, colScale=yellowBlueCols, 
              xlab="LatticeKrig", ylab="ELK with 2D RW", xlim=predLim, ylim=predLim, 
              main="Prediction comparison", 
              pch=19, cex=.2, ordering="increasing")
dev.off()

pdf(paste0("Figures/applicationHeaton/ELKbasic-ELKfinalIntPredsComparison", thisFileNameRoot, ".pdf"), width=5, height=5)
predLim = range(c(finalResults$yHat, fitELK$preds, fitELKfinalInt$preds))
plotWithColor(fitELK$preds, fitELKfinalInt$preds, nndists, colScale=yellowBlueCols, 
              xlab="ELK (no 2D RW)", ylab="ELK (with 2D RW)", xlim=predLim, ylim=predLim, 
              main="Prediction comparison", 
              pch=19, cex=.2, ordering="increasing")
dev.off()

pdf(paste0("Figures/applicationHeaton/LK-ELKbasicSDsComparison", thisFileNameRoot, ".pdf"), width=5, height=5)
sdLim = range(c(finalResults$standError, sds, sdsFinalInt))
plotWithColor(finalResults$standError, sds, nndists, colScale=yellowBlueCols, 
              xlab="LatticeKrig", ylab="ELK (no 2D RW)", xlim=sdLim, ylim=sdLim, 
              main="SD comparison", 
              pch=19, cex=.2, ordering="increasing")
abline(0, 1, lty=2)
dev.off()

pdf(paste0("Figures/applicationHeaton/LK-ELKfinalIntSDsComparison", thisFileNameRoot, ".pdf"), width=5, height=5)
sdLim = range(c(finalResults$standError, sds, sdsFinalInt))
plotWithColor(finalResults$standError, sdsFinalInt, nndists, colScale=yellowBlueCols, 
              xlab="LatticeKrig", ylab="ELK (with 2D RW)", xlim=sdLim, ylim=sdLim, 
              main="SD comparison", 
              pch=19, cex=.2, ordering="increasing")
abline(0, 1, lty=2)
dev.off()

pdf(paste0("Figures/applicationHeaton/ELKbasic-ELKfinalIntSDsComparison", thisFileNameRoot, ".pdf"), width=5, height=5)
sdLim = range(c(finalResults$standError, sds, sdsFinalInt))
plotWithColor(sds, sdsFinalInt, nndists, colScale=yellowBlueCols, 
              xlab="ELK (no 2D RW)", ylab="ELK (with 2D RW)", xlim=sdLim, ylim=sdLim, 
              main="SD comparison", 
              pch=19, cex=.2, ordering="increasing")
abline(0, 1, lty=2)
dev.off()

# ELK nonlinear plots
xrange = range(all.sat.temps$Lon)
yrange = range(all.sat.temps$Lat)
pdf(paste0("Figures/applicationHeaton/ELKfinalIntPreds", thisFileNameRoot, ".pdf"), width=5, height=5)
zlim = range(c(finalResults$yHat, fitELK$preds, fitELKfinalInt$preds))
quilt.plot(cbind(datPred$Lon, datPred$Lat), fitELKfinalInt$preds, nx=500, ny=300, 
           xlim=xrange, ylim=yrange, main="ELK predictions (2d covariate interaction model)", 
           xlab="Longitude", ylab="Latitude", zlim=zlim)
dev.off()

pdf(paste0("Figures/applicationHeaton/ELKfinalIntSDs", thisFileNameRoot, ".pdf"), width=5, height=5)
sds = rowMeans(outer(fitELK$sigmasNoNugget^2, fitELK$clusterVars, function(x, y) {sqrt(x + y)}))
sdsNonlinearInt = rowMeans(outer(fitELKfinalInt$sigmasNoNugget^2, fitELKfinalInt$clusterVars, function(x, y) {sqrt(x + y)}))
zlim = range(c(finalResults$standError, sds, sdsNonlinearInt))
quilt.plot(cbind(datPred$Lon, datPred$Lat), sdsNonlinearInt, nx=500, ny=300, 
           xlim=xrange, ylim=yrange, main="ELK SD (2d covariate interaction model)", 
           xlab="Longitude", ylab="Latitude", zlim=zlim)
dev.off()

pdf(paste0("Figures/applicationHeaton/ELKfinalIntResids", thisFileNameRoot, ".pdf"), width=5, height=5)
ylim = range(c(finalResults$yHat-thisDataObject$TrueTempMissing, fitELK$preds-datPred$TrueTemp, 
               fitELKfinalInt$preds-datPred$TrueTemp), na.rm=TRUE)
plot(datPred$elev, fitELKfinalInt$preds-datPred$TrueTemp, 
     xlab="Elevation (m)", ylab="ELK residuals", 
     main="ELK residuals vs. elevation (final model)", 
     pch=19, cex=.1, col="blue", ylim=ylim)
dev.off()

# pdf(paste0("Figures/applicationHeaton/ELKfinalIntElev", thisFileNameRoot, ".pdf"), width=5, height=5)
# xs = fitELKfinalInt$rwSummary$nonlinearEffect1$ID
# ys = xs*fitELKfinalInt$fixedEffectSummary$mean[4] + fitELKfinalInt$rwSummary$nonlinearEffect1$mean
# plot(xs, ys, 
#      xlab="Elevation (m)", ylab="Effect (degrees F)", 
#      main="Main effect of elevation on temperature", 
#      pch=19, cex=1, col="blue")
# dev.off()
# 
# pdf(paste0("Figures/applicationHeaton/ELKfinalIntNDVI", thisFileNameRoot, ".pdf"), width=5, height=5)
# xs = fitELKfinalInt$rwSummary$nonlinearEffect2$ID
# ys = xs*fitELKfinalInt$fixedEffectSummary$mean[5] + fitELKfinalInt$rwSummary$nonlinearEffect2$mean
# plot(xs, ys, 
#      xlab="NDVI", ylab="Effect (degrees F)", 
#      main="Main effect of NDVI on temperature", 
#      pch=19, cex=1, col="blue")
# dev.off()

pdf(paste0("Figures/applicationHeaton/ELKfinalIntRW2D", thisFileNameRoot, ".pdf"), width=5, height=5)
knotCoords = fitELKfinalInt$rw2dKnotCoords[,1:2]
knotVals = fitELKfinalInt$rw2dSummary$mean
quilt.plot(knotCoords, knotVals, nx=30, ny=30, 
           xlab="Elevation (m)", ylab="NDVI", 
           main="RW2D Mean", col=makeGreenBlueDivergingColors(64, range(knotVals), 0, TRUE))
dev.off()

##### Final models ----

# * setup ----
# dataObject<- makeData(sat.temps$Lon,sat.temps$Lat, sat.temps$Temp)

# get elevation
# GMTED 2010 elevation data downloaded from https://topotools.cr.usgs.gov/gmted_viewer/viewer.htm
# temp = raster("gmted2010_30arcsec_US_lon_-90_-120_lat_30_50.tif")
# elevObs = extract(temp, SpatialPoints(dataObject$x, proj4string=CRS("+proj=longlat")),method="bilinear")
# elevPred = extract(temp, SpatialPoints(dataObject$xMissing, proj4string=CRS("+proj=longlat")),method="bilinear")
# 
# NC<- 40
# nlevel<- 3
# nlevel<- 2
# a.wght<-  10.25
# nu<- .1
# kappa = sqrt(a.wght - 4)
# 
# nSample = 1000
# sampleI = sample(1:nrow(dataObject$x), nSample, replace=FALSE)
# 
# # setup LKrig object
# LKinfoFinal<- LKrigSetup( dataObject$x[sampleI,],
#                           NC = NC,
#                           nlevel = nlevel,
#                           a.wght = a.wght,
#                           nu = nu)
# 
# # setup ELK object
# latInfo = LKinfo2ELKBasis(LKinfoFinal)
# Ms = getMs(latInfo=latInfo)
# Ms
# sapply(latInfo, function(x) {x$latWidth})

# * LK ----

# * * Fit model ----
require(spam64)
comp.timeLKfit <- system.time(fitFinal<- LatticeKrig(dataObject$x[sampleI,], dataObject$y[sampleI], 
                                                     LKinfo=LKinfoFinal, 
                                                     verbose=TRUE))

# * * Get SEs ----
set.seed(234)

M <- 100  # Number of conditional draws (used for standard errors) [Final model: 100]
comp.timeLKse = 
  system.time(outputSim<- 
                LKrig.sim.conditional(fitFinal,
                                      x.grid = dataObject$xMissing, 
                                      M = M))

# * * Get CIs ----
comp.timeLKci <- system.time(standardError<- sqrt( 
  apply( outputSim$g.draw, 1, "var") +
    fitFinal$sigma.MLE^2))


yHat<- outputSim$ghat
CI80Lower<- yHat + qnorm(.1) * standardError
CI80Upper<- yHat + qnorm(.9) * standardError

comp.timeLK = comp.timeLKfit + comp.timeLKse + comp.timeLKci

# * * Save results ----

finalResults<- list(x=dataObject$xMissing,
                    yHat=yHat, 
                    standError=standardError)

save( finalResults, fitFinal, comp.timeLK, comp.timeLKfit, comp.timeLKse, comp.timeLKci, file=paste0("savedOutput/heaton/finalResultsLK", thisFileNameRoot, ".rda") )

# * ELK ----

# * * set priors ----
priorPar = getPCPrior(max(c(latInfo[[1]]$xRangeDat, latInfo[[1]]$yRangeDat))/5, .01, 1, nLayer=nlevel, separateRanges=FALSE, latticeInfo=latInfo, useUrbanPrior=FALSE) # 37.06811/5

# * * fit model ----
sampleI = sample(1:nrow(dataObject$x), 1000)
thisDataObject = dataObject
thisDataObject$x = thisDataObject$x[sampleI,]
thisDataObject$y = thisDataObject$y[sampleI]
dat = data.frame(Lon=thisDataObject$x[,1], Lat=thisDataObject$x[,2], MaskedTemp=thisDataObject$y, 
                 elev=extract(elevRaster, SpatialPoints(cbind(thisDataObject$x[,1], thisDataObject$x[,2]), proj4string=CRS("+proj=longlat")),method="bilinear"), 
                 distToWater=extract(distToWaterRaster, SpatialPoints(cbind(thisDataObject$x[,1], thisDataObject$x[,2]), proj4string=CRS("+proj=longlat")),method="bilinear"), 
                 ndvi=extract(ndviRaster, SpatialPoints(cbind(thisDataObject$x[,1], thisDataObject$x[,2]), proj4string=CRS("+proj=longlat")),method="bilinear"))
datPred = data.frame(Lon=thisDataObject$xMissing[,1], Lat=thisDataObject$xMissing[,2], MaskedTemp=thisDataObject$yMissing, 
                     elev=extract(elevRaster, SpatialPoints(cbind(thisDataObject$xMissing[,1], thisDataObject$xMissing[,2]), proj4string=CRS("+proj=longlat")),method="bilinear"), 
                     distToWater=extract(distToWaterRaster, SpatialPoints(cbind(thisDataObject$xMissing[,1], thisDataObject$xMissing[,2]), proj4string=CRS("+proj=longlat")),method="bilinear"), 
                     ndvi=extract(ndviRaster, SpatialPoints(cbind(thisDataObject$xMissing[,1], thisDataObject$xMissing[,2]), proj4string=CRS("+proj=longlat")),method="bilinear"))
dat$Covar = dat$elev
datPred$Covar = datPred$elev
comp.timeELKall = system.time(fitFinal <- modHeaton(dat, 
                                                    predCoords=cbind(datPred$Lon, datPred$Lat), 
                                                    predCovar=datPred$elev, 
                                                    nu=nu, seed=234, nLayer=nlevel, NC=NC,
                                                    nBuffer=nBuffer, priorPar=priorPar, normalize=TRUE, 
                                                    intStrategy="eb", strategy="gaussian", fastNormalize=TRUE, 
                                                    predictionType=c("mean", "median"), significanceCI=0.8, 
                                                    printVerboseTimings=FALSE, nPostSamples=1000, family="normal",
                                                    clusterEffect=TRUE, latInfo=latInfo, 
                                                    initialEffectiveRange=3, 
                                                    verbose=TRUE, separateRanges=FALSE, 
                                                    loadPrecomputationResults=TRUE, 
                                                    precomputationFileNameRoot=precomputationFileNameRoot, 
                                                    diagonal=c(0.0)))


# try fitting ELK with multiple nonlinear effects
comp.timeELKfit = system.time(fitELKfinal2 <- fitLKINLAStandard2(thisDataObject$x, thisDataObject$y, 
                                                           predCoords=cbind(datPred$Lon, datPred$Lat), 
                                                           xObs=cbind(1, dat$Lon, dat$Lat, dat$elev, dat$ndvi, dat$elev*dat$ndvi), 
                                                           xPred=cbind(1, datPred$Lon, datPred$Lat, datPred$elev, datPred$ndvi, datPred$elev*datPred$ndvi), 
                                                           nonlinearCovariateInds=c(4, 5), 
                                                           nu=nu, seed=234, nLayer=nlevel, NC=NC,
                                                           nBuffer=nBuffer, priorPar=priorPar, normalize=TRUE, 
                                                           intStrategy="eb", strategy="gaussian", fastNormalize=TRUE, 
                                                           predictionType=c("mean", "median"), significanceCI=0.8, 
                                                           printVerboseTimings=FALSE, nPostSamples=1000, family="normal",
                                                           clusterEffect=TRUE, latInfo=latInfo, 
                                                           initialEffectiveRange=3, 
                                                           verbose=TRUE, separateRanges=FALSE, 
                                                           loadPrecomputationResults=TRUE, 
                                                           precomputationFileNameRoot=precomputationFileNameRoot, 
                                                           diagonal=c(0.0)))

# ELK nonlinear plots
pdf(paste0("Figures/applicationHeaton/ELKfinalPreds", thisFileNameRoot, ".pdf"), width=5, height=5)
zlim = range(c(finalResults$yHat, fitELK$preds, fitELKfinal$preds))
quilt.plot(cbind(datPred$Lon, datPred$Lat), fitELKfinal$preds, nx=500, ny=300, 
           xlim=xrange, ylim=yrange, main="ELK predictions (final covariate model)", 
           xlab="Longitude", ylab="Latitude", zlim=zlim)
dev.off()

pdf(paste0("Figures/applicationHeaton/ELKfinalSDs", thisFileNameRoot, ".pdf"), width=5, height=5)
sds = rowMeans(outer(fitELK$sigmasNoNugget^2, fitELK$clusterVars, function(x, y) {sqrt(x + y)}))
sdsFinal = rowMeans(outer(fitELKfinal$sigmasNoNugget^2, fitELKfinal$clusterVars, function(x, y) {sqrt(x + y)}))
zlim = range(c(finalResults$standError, sds, sdsFinal))
quilt.plot(cbind(datPred$Lon, datPred$Lat), sdsFinal, nx=500, ny=300, 
           xlim=xrange, ylim=yrange, main="ELK SD (final covariate model)", 
           xlab="Longitude", ylab="Latitude", zlim=zlim)
dev.off()

pdf(paste0("Figures/applicationHeaton/ELKfinalResids", thisFileNameRoot, ".pdf"), width=5, height=5)
ylim = range(c(finalResults$yHat-thisDataObject$TrueTempMissing, fitELK$preds-datPred$TrueTemp, 
               fitELKfinal$preds-datPred$TrueTemp), na.rm=TRUE)
plot(datPred$elev, fitELKfinal$preds-datPred$TrueTemp, 
     xlab="Elevation (m)", ylab="ELK residuals", 
     main="ELK residuals vs. elevation (final model)", 
     pch=19, cex=.1, col="blue", ylim=ylim)
dev.off()

pdf(paste0("Figures/applicationHeaton/ELKfinalElev", thisFileNameRoot, ".pdf"), width=5, height=5)
xs = fitELKfinal$rwSummary$nonlinearEffect1$ID
ys = xs*fitELKfinal$fixedEffectSummary$mean[4] + fitELKfinal$rwSummary$nonlinearEffect1$mean
plot(xs, ys, 
     xlab="Elevation (m)", ylab="Effect (degrees F)", 
     main="Estimated effect of elevation on temperature", 
     pch=19, cex=1, col="blue")
dev.off()

pdf(paste0("Figures/applicationHeaton/ELKfinalNDVI", thisFileNameRoot, ".pdf"), width=5, height=5)
xs = fitELKfinal$rwSummary$nonlinearEffect2$ID
ys = xs*fitELKfinal$fixedEffectSummary$mean[5] + fitELKfinal$rwSummary$nonlinearEffect2$mean
plot(xs, ys, 
     xlab="NDVI", ylab="Effect (degrees F)", 
     main="Estimated effect of NDVI on temperature", 
     pch=19, cex=1, col="blue")
dev.off()

# ELK nonlinear plots (2!)
pdf(paste0("Figures/applicationHeaton/ELKfinal2Preds", thisFileNameRoot, ".pdf"), width=5, height=5)
zlim = range(c(finalResults$yHat, fitELK$preds, fitELKfinal2$preds))
quilt.plot(cbind(datPred$Lon, datPred$Lat), fitELKfinal2$preds, nx=500, ny=300, 
           xlim=xrange, ylim=yrange, main="ELK predictions (final covariate model)", 
           xlab="Longitude", ylab="Latitude", zlim=zlim)
dev.off()

pdf(paste0("Figures/applicationHeaton/ELKfinal2SDs", thisFileNameRoot, ".pdf"), width=5, height=5)
sds = rowMeans(outer(fitELK$sigmasNoNugget^2, fitELK$clusterVars, function(x, y) {sqrt(x + y)}))
sdsNonlinear = rowMeans(outer(fitELKnonlinear$sigmasNoNugget^2, fitELKnonlinear$clusterVars, function(x, y) {sqrt(x + y)}))
sdsFinal = rowMeans(outer(fitELKfinal2$sigmasNoNugget^2, fitELKfinal2$clusterVars, function(x, y) {sqrt(x + y)}))
zlim = range(c(finalResults$standError, sds, sdsNonlinear, sdsFinal))
quilt.plot(cbind(datPred$Lon, datPred$Lat), sdsFinal, nx=500, ny=300, 
           xlim=xrange, ylim=yrange, main="ELK SD (final covariate model)", 
           xlab="Longitude", ylab="Latitude", zlim=zlim)
dev.off()

pdf(paste0("Figures/applicationHeaton/ELKfinal2Resids", thisFileNameRoot, ".pdf"), width=5, height=5)
ylim = range(c(finalResults$yHat-thisDataObject$TrueTempMissing, fitELK$preds-datPred$TrueTemp, 
               fitELKnonlinear$preds-datPred$TrueTemp, fitELKfinal2$preds-datPred$TrueTemp), na.rm=TRUE)
plot(datPred$elev, fitELKnonlinear$preds-datPred$TrueTemp, 
     xlab="Elevation (m)", ylab="ELK residuals", 
     main="ELK residuals vs. elevation (final model)", 
     pch=19, cex=.1, col="blue", ylim=ylim)
dev.off()

pdf(paste0("Figures/applicationHeaton/ELKfinal2Elev", thisFileNameRoot, ".pdf"), width=5, height=5)
xs = fitELKfinal2$rwSummary$nonlinearEffect1$ID
ys = xs*fitELKfinal2$fixedEffectSummary$mean[4] + fitELKfinal2$rwSummary$nonlinearEffect1$mean
plot(xs, ys, 
     xlab="Elevation (m)", ylab="Effect (degrees F)", 
     main="Estimated effect of elevation on temperature", 
     pch=19, cex=1, col="blue")
dev.off()

pdf(paste0("Figures/applicationHeaton/ELKfinal2NDVI", thisFileNameRoot, ".pdf"), width=5, height=5)
xs = fitELKfinal2$rwSummary$nonlinearEffect2$ID
ys = xs*fitELKfinal2$fixedEffectSummary$mean[5] + fitELKfinal2$rwSummary$nonlinearEffect2$mean
plot(xs, ys, 
     xlab="NDVI", ylab="Effect (degrees F)", 
     main="Estimated effect of NDVI on temperature", 
     pch=19, cex=1, col="blue")
dev.off()




