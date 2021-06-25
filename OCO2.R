# getting OCO2 data from Katzfuss MRA big data paper:
# https://arxiv.org/pdf/1805.03309.pdf
# get the dataset:
# https://disc.gsfc.nasa.gov/datasets/OCO2_L2_Lite_SIF_10r/summary?keywords=oco2
# download the dataset by first following these instructions:
# https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+cURL+And+Wget
# use the "-i" flag of wget to download list of links in the given file:
# wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --keep-session-cookies --content-disposition -i subset_OCO2_L2_Lite_SIF_10r_20210616_205947.txt
# data summarization
# ncview file_name
# ncdump -h file_name

library(ncdf4)
out = system("ls oco2/*.nc4", intern=TRUE)
oco2 = data.frame()
for(i in 1:length(out)) {
  thisFile = out[i]
  nc_data <- nc_open(thisFile)
  # names(nc_data$var)
  lon <- ncvar_get(nc_data, "Geolocation/longitude")
  lat <- ncvar_get(nc_data, "Geolocation/latitude")
  t <- ncvar_get(nc_data, "Delta_Time") # seconds since 1 January 1990
  # t <- ncvar_get(nc_data, "Geolocation/time_tai93") # seconds since 1 January 1993
  alt <- ncvar_get(nc_data, "Geolocation/altitude") # in meters according to ncdump -h
  sifDaily740 <- ncvar_get(nc_data, "Daily_SIF_740nm") # Daily Corrected Solar induced chlorophyll fluorescence at 740 nm: Daily_Average_SIF_740 = /Science/sif_740nm * /Science/daily_correction_factor. In W/m^2/sr/µm. 
  sifDaily757 <- ncvar_get(nc_data, "Daily_SIF_757nm")
  sifDaily771 <- ncvar_get(nc_data, "Daily_SIF_771nm")
  sif740 <- ncvar_get(nc_data, "SIF_740nm") # Solar induced chlorophyll fluorescence at retrieved wavelength: SIF_740 = 0.75 * (/Science/sif_757nm + 1.5*/Science/sif_771nm)
  sif740Sigma <- ncvar_get(nc_data, "SIF_Uncertainty_740nm") # Daily Corrected Solar induced chlorophyll fluorescence at 740 nm: Daily_Average_SIF_740 = /Science/sif_740nm * /Science/daily_correction_factor. In W/m^2/sr/µm. 
  # sifScienceDaily757 <- ncvar_get(nc_data, "Science/SIF_757nm")
  # sifScienceDaily771 <- ncvar_get(nc_data, "Science/SIF_771nm") # Offset-Adjusted Solar Induced Chlorophyll Fluorescence at 771nm
  landFraction <- ncvar_get(nc_data, "Science/sounding_land_fraction") # 0-100
  qualityFlag <- ncvar_get(nc_data, "Quality_Flag") # SIF Lite Quality Flag: 0 = best (passes quality control + cloud fraction = 0.0); 1 = good (passes quality control); 2 = bad (failed quality control); -1 = not investigated
  fileName = thisFile
  thisFrame = data.frame(lon, lat, t, alt, landFraction, sifDaily740, sifDaily757, sifDaily771, sif740, sif740Sigma, qualityFlag, fileName)
  oco2 = rbind(oco2, thisFrame)
  nc_close(nc_data) 
}

# use data that is only within the contiguous US
usmap = US()
usmap = map("usa", plot=FALSE)
temp = map("usa", plot=FALSE, fill=TRUE)
# usmap = polygon(usmap$x, usmap$y)
test = maps:::map.poly("usa")
nas = is.na(test$x)
usPolys = list()
startI = c(1, which(nas)+1)
endI = c(which(nas)-1, length(test$x))
for(i in 1:length(startI)) {
  theseIs = startI[i]:endI[i]
  thisPoly = Polygon(cbind(test$x[theseIs], test$y[theseIs]))
  usPolys = c(usPolys, list(thisPoly))
}
inPolys = do.call("cbind", lapply(usPolys, function(p) {in.poly(cbind(oco2$lon, oco2$lat), p@coords)}))
inUSA = apply(inPolys, 1, any)
# inUSA = in.poly(cbind(oco2$lon, oco2$lat), cbind(usmap$x, usmap$y))
# plot(cbind(usmap$x, usmap$y), col=rainbow(length(usmap$x)))
# plot(cbind(test$x, test$y), col=rainbow(length(usmap$x)))

# test2 = cbind(test$x[!is.na(test$x)], test$y[!is.na(test$y)])
# inUSA = maps:::in.one.polygon(test, cbind(oco2$lon, oco2$lat))
oco2good = oco2[inUSA,]
oco2 = oco2good
oco2 = oco2[oco2$qualityFlag != 2,] # remove bad quality data

dev.off()
plot.new()
quilt.plot(oco2$lon, oco2$lat, oco2$sifDaily740, col=viridis(64), 
           nx=400, ny=900, xlab="Longitude", ylab="Latitude")
US(add=TRUE)
quilt.plot(oco2$lon, oco2$lat, oco2$sifDaily757, col=viridis(64), 
           nx=400, ny=900, xlab="Longitude", ylab="Latitude")
US(add=TRUE)
quilt.plot(oco2$lon, oco2$lat, oco2$sifDaily771, col=viridis(64), 
           nx=400, ny=900, xlab="Longitude", ylab="Latitude")
US(add=TRUE)
plotWithColor(oco2$lon, oco2$lat, oco2$sifDaily740, colScale=viridis(64), 
              pch=19, cex=.01, xlab="Longitude", ylab="Latitude")
US(add=TRUE)

quilt.plot(oco2$lon, oco2$lat, oco2$t, col=viridis(64), 
           nx=150, ny=500, xlab="Longitude", ylab="Latitude")
US(add=TRUE)

quilt.plot(oco2$lon, oco2$lat, oco2$alt, col=viridis(64), 
           nx=150, ny=500, xlab="Longitude", ylab="Latitude")
US(add=TRUE)

png("Figures/applicationOCO2/SIF.png", width=500, height=500)
quilt.plot(oco2$lon, oco2$lat, oco2$sif740, #col=viridis(64), 
           nx=150, ny=500, xlab="Longitude", ylab="Latitude")
US(add=TRUE)
dev.off()

png("Figures/applicationOCO2/altitude.png", width=500, height=500)
quilt.plot(oco2$lon, oco2$lat, oco2$alt, col=viridis(64), 
           nx=150, ny=500, xlab="Longitude", ylab="Latitude")
US(add=TRUE)
dev.off()

png("Figures/applicationOCO2/landFraction.png", width=500, height=500)
quilt.plot(oco2$lon, oco2$lat, oco2$landFraction, col=viridis(64), 
           nx=150, ny=500, xlab="Longitude", ylab="Latitude")
US(add=TRUE)
dev.off()

png("Figures/applicationOCO2/SIF_SD.png", width=500, height=500)
quilt.plot(oco2$lon, oco2$lat, oco2$sif740Sigma, col=viridis(64), 
           nx=150, ny=500, xlab="Longitude", ylab="Latitude")
US(add=TRUE)
dev.off()

# test to make sure we're in 2018 and in a 31 day month (should be in August)
secondsPerDay = 60*60*24
leapYears = seq(1992, 2020, by=4)
years = 1990:2020
daysPerYear = rep(365, length(years))
daysPerYear[years %in% leapYears] = 366
yearStartInSeconds = c(0, cumsum(daysPerYear*secondsPerDay))[-32]
yearEndInSeconds = cumsum(daysPerYear*secondsPerDay)
yearStartInSeconds[years == 2018]
yearEndInSeconds[years == 2018]
mean((oco2$t >= yearStartInSeconds[years == 2018]) & (oco2$t <= yearEndInSeconds[years == 2018]))
diff(range(oco2$t))
secondsPerDay*31
(secondsPerDay*31-diff(range(oco2$t)))/60


# project onto a reasonable coordinate system
# see https://gis.stackexchange.com/questions/141580/which-projection-is-best-for-mapping-the-contiguous-united-states 
# for coordinate systems for contiguous USA. Choose: EPSG 102004
# from lon/lat coords to easting/northing
lonLatCoords = SpatialPoints(cbind(oco2$lon, oco2$lat), proj4string=CRS("+proj=longlat"))
coordsUTM = spTransform(lonLatCoords, CRS("+proj=lcc +lat_1=33 +lat_2=45 +lat_0=39 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"))
eastNorth = attr(coordsUTM, "coords")
oco2$east = eastNorth[,1]
oco2$north = eastNorth[,2]

# get elevation
temp = raster("elevation.nc")
elev = extract(temp, SpatialPoints(cbind(oco2$lon, oco2$lat)),method="bilinear")
plot(oco2$alt, elev, pch=".")
abline(0, 1, col="blue")
hist(oco2$alt - elev, breaks=70)
hist(abs(oco2$alt - elev), breaks=70)
mean(abs(oco2$alt - elev))
mean(abs(oco2$alt - elev)<50)
oco2$alt = elev

# save resulting dataset
save(oco2, file="savedOutput/OCO2/oco2.RData")
load("savedOutput/OCO2/oco2.RData")

# generate prediction points in lon/lat coords:
nas = is.na(usmap$x)
usmapProj = spTransform(SpatialPoints(cbind(usmap$x, usmap$y)[!nas,], proj4string=CRS("+proj=longlat")), CRS("+proj=lcc +lat_1=33 +lat_2=45 +lat_0=39 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"))
usmapProjCoords = attr(usmapProj, "coords")
eastRange = range(usmapProjCoords[,1])
northRange = range(usmapProjCoords[,2])
lonRange = usmap$range[1:2]
latRange = usmap$range[3:4]
res = 10/60 # 20 arcminutes
lonGrid = seq(lonRange[1], lonRange[2]+res, by=res)
latGrid = seq(latRange[1], latRange[2]+res, by=res)
predPoints = expand.grid(list(lon=lonGrid, lat=latGrid))
inPolys = do.call("cbind", lapply(usPolys, function(p) {in.poly(cbind(predPoints$lon, predPoints$lat), p@coords)}))
inUS = apply(inPolys, 1, any)
predPoints = predPoints[inUS,]
dim(predPoints)
predLonLat = predPoints

# project to easting/northing
predEastNorth = spTransform(SpatialPoints(predLonLat, proj4string=CRS("+proj=longlat")), CRS("+proj=lcc +lat_1=33 +lat_2=45 +lat_0=39 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"))
predEastNorth = attr(predEastNorth, "coords")

# get elevation at prediction points
predAlt = extract(temp, SpatialPoints(predLonLat),method="bilinear")

# plot elevation at prediction points
png("Figures/applicationOCO2/elevation.png", width=500, height=500)
quilt.plot(predEastNorth, predAlt, col=viridis(64), 
           nx=length(lonGrid)-70, ny=length(latGrid)-28, 
           xlab="Easting (km)", ylab="Northing (km)")
# polygon(usmapProjCoords[,1], usmapProjCoords[,2])
dev.off()

# plot SIF versus elevation
plot(oco2$alt, oco2$sif740, pch=".")

cuts = cut(oco2$alt, breaks=seq(-100, 3700, by=100), include.lowest=TRUE)
xs = seq(-50, 3650, by=100)
breakInd = as.integer(cuts)
xVals = xs[breakInd]
meanVals = aggregate(oco2$sifDaily740, by=list(Elevation=xVals), FUN=mean)$x

pdf("Figures/applicationOCO2/elevVsSIF_mean.pdf", width=5, height=5)
plot(oco2$alt, oco2$sifDaily740, pch=".", 
     xlab="Elevation (m)", ylab="Solar induced chlorophyll fluorescence")
points(xs, meanVals, col="red", pch=19, cex=.5)
dev.off()
