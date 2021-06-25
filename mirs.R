# script for MIRS TPW (total precipitable water in mm) application

library(ncdf4)
out = system("ls MIRS/allDat/*.nc", intern=TRUE)
mirs = data.frame()
for(i in 1:length(out)) {
  if(i %% 10 == 1) {
    print(paste0("reading file ", i, "/", length(out)))
  }
  
  thisFile = out[i]
  nc_data <- nc_open(thisFile)
  # names(nc_data$var)
  # names(nc_data$dim)
  # sapply(nc_data$dim, function(x){x$len})
  lon <- ncvar_get(nc_data, "Longitude") # (FOV x Scanline)
  lat <- ncvar_get(nc_data, "Latitude") # (FOV x Scanline)
  t <- ncvar_get(nc_data, "ScanTime_UTC") # Number of seconds since 00:00:00 UTC (Scanline)
  doy <- ncvar_get(nc_data, "ScanTime_doy") # julian day 1-366 in days (Scanline)
  
  tpw <- ncvar_get(nc_data, "TPW") # in mm (FOV x Scanline)
  
  # see https://www.star.nesdis.noaa.gov/mirs/documents/documentation/doc_v11r3/MIRS_Users_Manual.pdf
  # for detailed interpretation of Qc flag. In general, we only care about the first Qc dimension
  qualityFlag <- ncvar_get(nc_data, "Qc")[1,,] # Quality Flag: 0 = good; 1 = ok with some problems; 2 = bad (Qc_dim x FOV x Scanline)
  fileName = thisFile
  thisFrame = data.frame(c(lon), c(lat), rep(t, each=nrow(lat)), rep(doy, each=nrow(lat)), c(tpw), c(qualityFlag), fileName)
  mirs = rbind(mirs, thisFrame)
  nc_close(nc_data) 
}
names(mirs) = c("lon", "lat", "t", "doy", "tpw", "quality", "fileName")
dim(mirs)

mirsFull = mirs
mirs = mirsFull

# remove bad data
mirs = mirs[mirs$quality <= 1,]
dim(mirs)

# take only data near US
goodLon = (mirs$lon >= -130) & (mirs$lon <= -70)
goodLat = (mirs$lat >= 20) & (mirs$lat <= 50)
mirs = mirs[goodLon & goodLat,]
dim(mirs)

# take only data from February 1st:
mirs = mirs[mirs$doy == 32,]
dim(mirs)

mirsBest = mirs[mirs$quality == 0,]

# use data that is only within the contiguous US
# usmap = US()
# usmap = map("usa", plot=FALSE)
# temp = map("usa", plot=FALSE, fill=TRUE)
# # usmap = polygon(usmap$x, usmap$y)
# test = maps:::map.poly("usa")
# nas = is.na(test$x)
# usPolys = list()
# startI = c(1, which(nas)+1)
# endI = c(which(nas)-1, length(test$x))
# for(i in 1:length(startI)) {
#   theseIs = startI[i]:endI[i]
#   thisPoly = Polygon(cbind(test$x[theseIs], test$y[theseIs]))
#   usPolys = c(usPolys, list(thisPoly))
# }
# inPolys = do.call("cbind", lapply(usPolys, function(p) {in.poly(cbind(mirs$lon, mirs$lat), p@coords)}))
# inUSA = apply(inPolys, 1, any)
# inUSA = in.poly(cbind(oco2$lon, oco2$lat), cbind(usmap$x, usmap$y))
# plot(cbind(usmap$x, usmap$y), col=rainbow(length(usmap$x)))
# plot(cbind(test$x, test$y), col=rainbow(length(usmap$x)))

# test2 = cbind(test$x[!is.na(test$x)], test$y[!is.na(test$y)])
# inUSA = maps:::in.one.polygon(test, cbind(oco2$lon, oco2$lat))
# oco2good = oco2[inUSA,]
# oco2 = oco2good
# oco2 = oco2[oco2$qualityFlag != 2,] # remove bad quality data

png("Figures/applicationMIRS/tpw.png", width=500, height=500)
quilt.plot(mirs$lon, mirs$lat, mirs$tpw, col=viridis(64), 
           nx=200, ny=300, xlab="Longitude", ylab="Latitude")
world(add=TRUE)
dev.off()

png("Figures/applicationMIRS/tpwGoodQuality.png", width=500, height=500)
quilt.plot(mirsBest$lon, mirsBest$lat, mirsBest$tpw, col=viridis(64), 
           nx=200, ny=300, xlab="Longitude", ylab="Latitude")
world(add=TRUE)
dev.off()

quilt.plot(mirs$lon, mirs$lat, mirs$t/(24*60*60), col=viridis(64), 
           nx=200, ny=300, xlab="Longitude", ylab="Latitude")
US(add=TRUE)
