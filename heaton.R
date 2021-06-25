blueGreenCols = rev(makeGreenBlueSequentialColors(64))
yellowBlueCols = makeBlueGreenYellowSequentialColors(64)
infernoCols = inferno(64)
magmaCols = magma(64)
plasmaCols = plasma(64)

# analysis of the heaton et al. dataset
# https://link.springer.com/article/10.1007/s13253-018-00348-w
# https://github.com/finnlindgren/heatoncomparison
out = load("heaton/AllSatelliteTemps.RData")
names(all.sat.temps)
out = load("heaton/SatelliteTemps.RData")
names(sat.temps)

# project onto a reasonable coordinate system
# see https://gis.stackexchange.com/questions/141580/which-projection-is-best-for-mapping-the-contiguous-united-states 
# for coordinate systems for contiguous USA. Choose: EPSG 102004
# from lon/lat coords to easting/northing
lonLatCoords = SpatialPoints(cbind(all.sat.temps$Lon, all.sat.temps$Lat), proj4string=CRS("+proj=longlat"))
coordsUTM = spTransform(lonLatCoords, CRS("+proj=lcc +lat_1=33 +lat_2=45 +lat_0=39 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"))
eastNorth = attr(coordsUTM, "coords")
all.sat.temps$east = eastNorth[,1]
all.sat.temps$north = eastNorth[,2]

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

# temp = raster("elevation.nc")
# elev = extract(temp, SpatialPoints(cbind(all.sat.temps$Lon, all.sat.temps$Lat), proj4string=CRS("+proj=longlat")),method="bilinear")
# GMTED 2010 elevation data downloaded from https://topotools.cr.usgs.gov/gmted_viewer/viewer.htm
temp = raster("gmted2010_30arcsec_US_lon_-90_-120_lat_30_50.tif")
elev = extract(temp, SpatialPoints(cbind(all.sat.temps$Lon, all.sat.temps$Lat), proj4string=CRS("+proj=longlat")),method="bilinear")
# temp2 = raster("gmted2010_30arcsec_US_lon_-90_-120_lat_30_50_breaklineEmphasis.tif")
# breakline = extract(temp2, SpatialPoints(cbind(all.sat.temps$Lon, all.sat.temps$Lat), proj4string=CRS("+proj=longlat")),method="bilinear")

png("Figures/applicationHeaton/elevation.png", width=700, height=700)
quilt.plot(all.sat.temps$Lon, all.sat.temps$Lat, elev, 
           col=viridis(64), nx=500, ny=300, 
           xlab="Longitude", ylab="Latitude")
dev.off()

# png("Figures/applicationHeaton/breakline.png", width=700, height=700)
# quilt.plot(all.sat.temps$Lon, all.sat.temps$Lat, breakline, 
#            col=viridis(64), nx=500, ny=300, 
#            xlab="Longitude", ylab="Latitude")
# dev.off()

out = lm(TrueTemp ~ Lon + Lat + elev, data=all.sat.temps)
summary(out)

temp = raster("MCD13.A2014.unaccum.nc4")
ndvi = extract(temp, SpatialPoints(cbind(all.sat.temps$Lon, all.sat.temps$Lat), proj4string=CRS("+proj=longlat")),method="bilinear")

png("Figures/applicationHeaton/ndvi1.png", width=700, height=700)
quilt.plot(all.sat.temps$Lon, all.sat.temps$Lat, ndvi, 
           col=viridis(64), nx=500, ny=300, 
           xlab="Longitude", ylab="Latitude")
dev.off()

daysInMonth = c(31, 29, 31, 30, 31, 30, 31)
sum(daysInMonth) # 213, beginning of August
sum(daysInMonth) + 29 # 242, last day of August
temp = raster("MOD13Q1.006__250m_16_days_NDVI_doy2016225_aid0001.tif")
ndvi2 = extract(temp, SpatialPoints(cbind(all.sat.temps$Lon, all.sat.temps$Lat), proj4string=CRS("+proj=longlat")),method="bilinear")

png("Figures/applicationHeaton/ndvi2.png", width=700, height=700)
quilt.plot(all.sat.temps$Lon, all.sat.temps$Lat, ndvi2, 
           col=viridis(64), nx=500, ny=300, 
           xlab="Longitude", ylab="Latitude")
dev.off()

temp = raster("MOD13Q1.006__250m_16_days_EVI_doy2016225_aid0001.tif")
dvi = extract(temp, SpatialPoints(cbind(all.sat.temps$Lon, all.sat.temps$Lat), proj4string=CRS("+proj=longlat")),method="bilinear")

png("Figures/applicationHeaton/dvi.png", width=700, height=700)
quilt.plot(all.sat.temps$Lon, all.sat.temps$Lat, dvi, 
           col=viridis(64), nx=500, ny=300, 
           xlab="Longitude", ylab="Latitude")
dev.off()

library(ncdf4)
# nc_data <- nc_open("globolakes-static_distance_to_water_Map-300m-P5Y-2005-ESACCI_WB-fv1.0.nc")
nc_data <- nc_open("distanceFromWater.nc")
# names(nc_data$var)
distanceToWater = ncvar_get(nc_data, "distance_to_water") # 
lonGrid <- nc_data$dim$lon$vals
latGrid <- nc_data$dim$lat$vals
lon = rep(lonGrid, length(latGrid))
lat = rep(latGrid, each=length(lonGrid))
# temp = rasterFromXYZ(data.frame(x=lon, y=lat, z=c(distanceToWater)))
# dataExtent = extent(c(min(lonGrid), max(lonGrid), min(latGrid), max(latGrid)))
# temp <- raster(nrows=length(lonGrid), ncols=length(latGrid),
#                xmn=min(lonGrid), xmx=max(lonGrid),
#                ymn=min(latGrid), ymx=max(latGrid),
#                vals=c(distanceToWater), # latitude is backwards, so we must reverse matrix to compensate
#                ext = dataExtent, crs=CRS("+proj=longlat"))
temp <- raster(dataExtent, nrows=length(latGrid), ncols=length(lonGrid))
temp2 <- rasterize(cbind(lon, lat), temp, distanceToWater, fun=mean)
plot(temp2)

distToWater = extract(temp2, SpatialPoints(cbind(all.sat.temps$Lon, all.sat.temps$Lat), proj4string=CRS("+proj=longlat")),method="bilinear")



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



out = lm(TrueTemp ~ Lon + Lat + elev + ndvi + distToWater, data=all.sat.temps)
summary(out)
out = lm(TrueTemp ~ Lon + Lat + ndvi + elev, data=all.sat.temps)
summary(out)

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

preds = predict(out, fullPredFrame)
resids = fullPredFrame$TrueTemp - preds
png("Figures/applicationHeaton/linearResidsFull.png", width=700, height=700)
quilt.plot(fullPredFrame$Lon, fullPredFrame$Lat, resids, 
           col=viridis(64), nx=500, ny=300, 
           xlab="Longitude", ylab="Latitude")
dev.off()

##### try fitting LK ----
NC<- 40
nlevel<- 4
a.wght<-  10.25
kappa = sqrt(a.wght - 4)
nu<- .1

# create locations and responses (y)
dataObject<- makeData(sat.temps$Lon,sat.temps$Lat, sat.temps$Temp)

# setup LKrig object
LKinfoFinal<- LKrigSetup( dataObject$x,
                          NC = NC,
                          nlevel = nlevel,
                          a.wght = a.wght,
                          nu = nu)
comp.timeLKfit <- system.time(fitFinal<- LatticeKrig( dataObject$x, dataObject$y, LKinfo= LKinfoFinal))

set.seed(234)

M <- 100  # Number of conditional draws (used for standard errors) [Final model: 100]
comp.timeLKse = 
  system.time(outputSim<- 
                LKrig.sim.conditional(fitFinal,
                                      x.grid = dataObject$xMissing, 
                                      M = M))
comp.timeLK = comp.timeLKfit + comp.timeLKse

##### try fitting ELK ----
xRange = range(sat.temps$Lon)
yRange = range(sat.temps$Lat)

latInfo = LKinfo2ELKBasis(LKinfoFinal)
Ms = getMs(latInfo=latInfo)
Ms
sapply(latInfo, function(x) {x$latWidth})





