library(colorspace)
library(mgcv)

# makes plots as well as parameter estimate tables for the example applications in the manuscript
makeAllPlots = function(dataType=c("ed", "mort"), resultFilenames, modelClasses, modelVariations, 
                        areaLevels=c("Region", "County", "Pixel", "Cluster"), 
                        meanRange, meanTicks, meanTickLabels, widthRange, widthTicks, widthTickLabels, 
                        relWidthRange, relWidthTicks, relWidthTickLabels, 
                        plotNameRoot="Education", resultNameRoot="Ed", meanCols=makeRedBlueDivergingColors(64), 
                        widthCols=makeBlueYellowSequentialColors(64), popCols=makeBlueSequentialColors(64), 
                        ncols=29, urbCols=makeGreenBlueSequentialColors(ncols), 
                        plotUrbanMap=FALSE, kenyaLatRange=c(-4.6, 5), kenyaLonRange=c(33.5, 42.0), 
                        makeModelPredictions=TRUE, makeCovariograms=TRUE, makePairPlots=TRUE, 
                        loadResults=FALSE, saveResults=!loadResults, singleColorBar=TRUE, doModelClassPlots=FALSE, 
                        col=NULL, lty=NULL, pch=NULL) {
  plotNameRootLower = tolower(plotNameRoot)
  resultNameRootLower = tolower(resultNameRoot)
  
  dataType = match.arg(dataType)
  if(dataType == "mort") {
    out = load("../U5MR/kenyaData.RData")
    dat = mort
  }
  else {
    out = load("../U5MR/kenyaDataEd.RData")
    dat = ed
  }
  
  # print lattice resolutions
  print("SPDE mesh resolution: ")
  test <- getSPDEMeshKenya()
  # plot(test)
  print(meanSegmentLength(getSPDEMeshKenya(), 50))
  print("three layer coarse resolution:")
  test = makeLatGridsKenya(3,NC=14)
  print(test[[1]]$latWidth)
  print("three layer fine resolution:")
  print(test[[3]]$latWidth)
  print("two layer coarse resolution:")
  test = makeLatGridsKenya(2,NC=c(30, 214))
  print(test[[1]]$latWidth)
  print("two layer fine resolution:")
  print(test[[2]]$latWidth)
  
  if(makeModelPredictions) {
    plotModelPredictions(dat, resultFilenames, modelClasses, modelVariations, areaLevels, 
                         meanRange, meanTicks, meanTickLabels, widthRange, widthTicks, widthTickLabels, 
                         relWidthRange, relWidthTicks, relWidthTickLabels, 
                         varName, plotNameRoot, resultNameRoot, meanCols, widthCols, 
                         kenyaLatRange, kenyaLonRange, singleColorBar, doModelClassPlots)
  }
  
  if(makeCovariograms) {
    plotCovariograms(dat, resultFilenames, modelClasses, modelVariations, 
                     varName, plotNameRoot, resultNameRoot, 
                     cgramList=NULL, loadResults=loadResults, saveResults=saveResults, doModelClassPlots=doModelClassPlots, 
                     col=col, lty=lty, pch=pch)
  }
  
  if(makePairPlots) {
    makeFinalPercentResidualPlot(dat, resultFilenames, modelClasses, modelVariations, 
                                 areaLevels, meanRange, meanTicks, meanTickLabels, 
                                 plotNameRoot=plotNameRoot, resultNameRoot, 
                                 meanCols, ncols, urbCols, kenyaLatRange, kenyaLonRange, doModelClassPlots)
    
    makeFinalPairPlot(dat, resultFilenames, modelClasses, modelVariations, 
                      areaLevels, meanRange, meanTicks, meanTickLabels, 
                      plotNameRoot=plotNameRoot, resultNameRoot, 
                      meanCols, ncols, urbCols, kenyaLatRange, kenyaLonRange, doModelClassPlots)
    makePairPlots(dat, resultFilenames, modelClasses, modelVariations, 
                  areaLevels, meanRange, meanTicks, meanTickLabels, 
                  plotNameRoot=plotNameRoot, resultNameRoot, 
                  meanCols, ncols, urbCols, kenyaLatRange, kenyaLonRange, doModelClassPlots)
  }
  
  invisible(NULL)
}

plotDataVisualizations = function(dataType=c("ed", "mort"), dat=NULL, 
                                  varName="SEP", plotNameRoot="Education", resultNameRoot="Ed", meanCols=makeRedBlueDivergingColors(64), 
                                  sdCols=makeBlueYellowSequentialColors(64), popCols=makeBlueSequentialColors(64), 
                                  ncols=29, relativeCols=makeRedGreenDivergingColors(ncols), urbCols=makeGreenBlueSequentialColors(ncols), 
                                  plotUrbanMap=FALSE, kenyaLatRange=c(-4.6, 5), kenyaLonRange=c(33.5, 42.0)) {
  plotNameRootLower = tolower(plotNameRoot)
  resultNameRootLower = tolower(resultNameRoot)
  
  dataType = match.arg(dataType)
  if(dataType == "mort") {
    out = load("../U5MR/kenyaData.RData")
    dat = mort
  }
  else {
    out = load("../U5MR/kenyaDataEd.RData")
    dat = ed
  }
  
  out = load("../U5MR/adminMapData.RData")
  kenyaMap = adm0
  countyMap = adm1
  
  print("generating data visualizations...")
  
  # plot the actual data
  png(file=paste0("Figures/", resultNameRoot, "/clustersUrban", plotNameRoot, ".png"), width=800, height=800)
  par(oma=c( 0,0,0,0), mar=c(5.5, 6.1, 3.5, 6))
  urban = dat$urban
  plot(dat$lon[!urban], dat$lat[!urban], pch=19, col="green", main=paste0("Urban vs. rural clusters"), xlim=kenyaLonRange, 
       ylim=kenyaLatRange, xlab="Longitude", ylab="Latitude", cex=.2, asp=1, cex.axis=2, cex.lab=2, tck=-.03, mgp=c(4, 2, 0), cex.main=2)
  points(dat$lon[urban], dat$lat[urban], pch=19, col="blue", cex=.2)
  # world(add=TRUE)
  plotMapDat(mapDat=adm1, lwd=.5)
  dev.off()
  
  # plot a map of urbanicity if requested (this can take ~10 minutes)
  if(plotUrbanMap) {
    # inside this if statement since it takes around ten minutes to run
    makeUrbanMap(kmres=1, savePlot=TRUE, lonLim=kenyaLonRange, latLim=kenyaLatRange, main="")
  }
  browser()
  png(file=paste0("Figures/", resultNameRoot, "/empirical", plotNameRoot, ".png"), width=800, height=800)
  par(oma=c( 0,0,0,0), mar=c(5.5, 6.1, 3.5, 6))
  # plot(cbind(dat$lon, dat$lat), type="n", ylim=kenyaLatRange, xlim=kenyaLonRange, 
  #      xlab="Longitude", ylab="Latitude", main=paste0("Empirical ", varName), asp=1)
  plot(cbind(dat$lon, dat$lat), type="n", ylim=kenyaLatRange, xlim=kenyaLonRange, 
       xlab="Longitude", ylab="Latitude", main="", asp=1, cex.axis=2, cex.lab=2, tck=-.03, mgp=c(4, 2, 0))
  quilt.plot(dat$lon, dat$lat, dat$y / dat$n, nx=100, ny=100, col=meanCols, add=TRUE, axis.args=list(cex.axis=2, tck=-.7, hadj=-.1))
  # world(add=TRUE)
  plotMapDat(mapDat=adm1, lwd=.5)
  dev.off()
  
  shrunkLogit = function(x, shrink=.08) {
    x = (x - .5) * (1 - 2 * shrink) + .5
    logit(x)
  }
  shrunkExpit = function(x, shrink=.08) {
    x = expit(x)
    x = (x - 0.5) / (1 - 2 * shrink) + .5
    x
  }
  
  png(file=paste0("Figures/", resultNameRoot, "/empirical", plotNameRoot, "ModLogit.png"), width=500, height=500)
  par(oma=c( 0,0,0,2), mar=c(5.1, 4.1, 4.1, 6))
  # ticks = pretty(seq(0, max(dat$y / dat$n), l=10), n=10)
  # ticks = logit(ticks[-c(1, 11)])
  ticks = seq(0, 1, l=11)
  labels = as.character(ticks)
  labels[seq(2, 10, by=2)] = ""
  ticks = shrunkLogit(ticks)
  
  varRange = shrunkExpit(range(ticks))
  # par( oma=c( 0,0,0,5)) # save some room for the legend
  plot(cbind(dat$lon, dat$lat), type="n", main="", ylim=kenyaLatRange, 
       xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", asp=1)
  quilt.plot(cbind(dat$lon, dat$lat), shrunkLogit(dat$y / dat$n), col=meanCols, 
             nx=100, ny=100, add.legend=FALSE, add=TRUE, zlim=shrunkLogit(varRange))
  plotMapDat(mapDat=adm1, lwd=.5)
  # world(add=TRUE)
  # par( oma=c(0,0,0,2))
  image.plot(zlim=shrunkLogit(varRange), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
             col=meanCols, add = TRUE, axis.args=list(at=ticks, labels=labels))
  dev.off()
}

makeUrbanMap = function(popGrid=NULL, kmres=1, savePlot=FALSE, fileName=ifelse(whiteRural, "Figures/UrbanMapWhiteRural.png", "Figures/urbanMap.png"), 
                        nx=850, ny=1050, width=800, height=800, kenyaLatRange=c(-4.6, 5), kenyaLonRange=c(33.5, 42.0), 
                        lonLim=kenyaLonRange, latLim=kenyaLatRange, whiteRural=TRUE, main="") {
  # get prediction locations from population grid
  if(is.null(popGrid)) {
    if(kmres == 5)
      load("../U5MR/popGrid.RData")
    else
      popGrid = makeInterpPopGrid(kmres)
  }
  
  out = load("../U5MR/adminMapData.RData")
  countyMap = adm1
  
  # determine which points in Kenya are urban
  threshes = setThresholds()
  popThreshes = sapply(1:nrow(popGrid), function(i) {threshes$threshes[threshes$counties == popGrid$admin1[i]]})
  urban = popGrid$popOrig > popThreshes
  
  if(savePlot) {
    png(file=fileName, width=width, height=height)
    par(oma=c( 0,0,0,0), mar=c(5.5, 6.1, 3.5, 6))
  }
  plot(popGrid$lon, popGrid$lat, xlab="Longitude", ylab="Latitude", main=main, xlim=lonLim, ylim=latLim, asp=1, type="n", 
       cex.axis=2, cex.lab=2, tck=-.03, mgp=c(4, 2, 0))
  # quilt.plot(popGrid$lon, popGrid$lat, urban, col=c("green", "blue"), nx=850, ny=1050, add.legend = FALSE, 
  #            xlab="Longitude", ylab="Latitude", main=TeX("Urbanicity"), xlim=lonLim, ylim=latLim, asp=1)
  if(whiteRural)
    quilt.plot(popGrid$lon, popGrid$lat, urban, col=c(rgb(0, 0, 0, 0), "blue"), nx=850, ny=1050, add.legend = FALSE, add=TRUE)
  else {
    quilt.plot(popGrid$lon, popGrid$lat, urban, col=c("green", "blue"), nx=850, ny=1050, add.legend = FALSE, add=TRUE)
    legend("bottomleft", c("urban", "rural"), col=c("blue", "green"), pch=19)
  }
  # world(add=TRUE)
  plotMapDat(mapDat=adm1, lwd=.5)
  if(savePlot) {
    dev.off()
  }
  browser()
}

# this function plots central estimates and credible interval widths over a map of Kenya 
# 4 plots are made by default: one for each area level
plotModelPredictions = function(dat, resultFilenames, modelClasses, modelVariations, 
                                areaLevels=c("Region", "County", "Pixel", "Cluster"), 
                                meanRange, meanTicks, meanTickLabels, widthRange, widthTicks, widthTickLabels, 
                                relWidthRange, relWidthTicks, relWidthTickLabels, 
                                varName="education", plotNameRoot="Education", resultNameRoot="Ed", 
                                meanCols=makeRedBlueDivergingColors(64), 
                                widthCols=makeBlueYellowSequentialColors(64), 
                                kenyaLatRange=c(-4.6, 5), kenyaLonRange=c(33.5, 42.0), 
                                singleColorBar=TRUE, doModelClassPlots=FALSE) {
  plotNameRootLower = tolower(plotNameRoot)
  resultNameRootLower = tolower(resultNameRoot)
  numberModels = length(resultFilenames)
  uniqueModelClasses = unique(modelClasses)
  
  ##### if there are multiple unique model classes, call this function on each one individually 
  ##### as well as on the combination
  if(length(uniqueModelClasses) != 1 && doModelClassPlots) {
    for(i in 1:length(uniqueModelClasses)) {
      thisModelClass = uniqueModelClasses[i]
      thisI = modelClasses == thisModelClass
      
      if(sum(thisI) == 1)
        next
      
      thisResultFilenames = resultFilenames[thisI]
      thisModelClasses = modelClasses[thisI]
      thisModelVariations = modelVariations[thisI]
      plotModelPredictions(dat, thisResultFilenames, thisModelClasses, thisModelVariations, areaLevels, 
                           meanRange, meanTicks, meanTickLabels, widthRange, widthTicks, widthTickLabels, 
                           relWidthRange, relWidthTicks, relWidthTickLabels, 
                           varName, plotNameRoot, resultNameRoot, meanCols, widthCols, 
                           kenyaLatRange, kenyaLonRange, singleColorBar)
    }
  }
  
  # load the fine grid used to approximate continuous prediction
  print("Loading prediction grid and shapefiles...")
  load("../U5MR/popGrid.RData")
  
  # load shape files for plotting
  require(maptools)
  # regionMap = readShapePoly("../U5MR/mapData/kenya_region_shapefile/kenya_region_shapefile.shp", delete_null_obj=TRUE, force_ring=TRUE, repair=TRUE)
  out = load("regionMap.RData")
  out = load("../U5MR/adminMapData.RData")
  kenyaMap = adm0
  countyMap = adm1
  
  ##### central estimates and credible interval widths
  
  if(length(uniqueModelClasses) == 1)
    extraPlotNameRoot = uniqueModelClasses[1]
  else
    extraPlotNameRoot = ""
  
  for(i in 1:length(areaLevels)) {
    thisArea = areaLevels[i]
    
    # get map for this areal aggregation level
    if(thisArea == "Cluster") {
      next
    } else if(thisArea == "Region") {
      thisMap = regionMap
    } else {
      thisMap = countyMap
    }
    
    print("Plotting central estimates and credible interval widths...")
    width = 400 * numberModels
    
    png(file=paste0("Figures/", resultNameRoot, "/preds", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=1000)
    
    if(thisArea %in% c("Region", "County"))
      par(mfrow=c(2,numberModels), oma=c( 0,6,0,7), mar=c(5.1, 5.1, 4.1, 2.5))
    else
      par(mfrow=c(2,numberModels), oma=c( 0,6,0,7), mar=c(5.1, 5.1, 4.1, 2.5))
    
    # first load in the predictions
    predictionList = list()
    widthList = list()
    for(j in 1:numberModels) {
      modelName = bquote(.(modelClasses[j])[.(modelVariations[j])])
      
      # load this model and plot results in the given column
      out = load(resultFilenames[j])
      theseResults = results$aggregatedResults$predictions[[paste0(tolower(thisArea), "Predictions")]]
      predictionList = c(predictionList, list(theseResults$preds))
      widthList = c(widthList, list(theseResults$Q90 - theseResults$Q10))
    }
    
    # plot the predictions
    for(j in 1:numberModels) {
      modelName = bquote(.(modelClasses[j])[.(modelVariations[j])])
      thisAddColorBar = (j == numberModels) | !singleColorBar
      
      # load this model and plot results in the given column
      # out = load(resultFilenames[j])
      # theseResults = results$aggregatedResults$predictions[[paste0(tolower(thisArea), "Predictions")]]
      
      if(thisArea == "Region") {
        plotMapDat(plotVar=predictionList[[j]], new = TRUE, 
                   main="", scaleFun=logit, scaleFunInverse=expit, 
                   cols=meanCols, zlim=logit(meanRange), ticks=meanTicks, tickLabels=meanTickLabels, 
                   xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = thisAddColorBar, 
                   legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot= c(.97,1,.1,.9)), legend.width=3, 
                   plotArgs=list(cex.main=3, cex.axis=2, cex.lab=2), legend.mar=0)
      } else if(thisArea == "County") {
        plotMapDat(plotVar=predictionList[[j]], new = TRUE, 
                   main="", scaleFun=logit, scaleFunInverse=expit, 
                   cols=meanCols, zlim=logit(meanRange), ticks=meanTicks, tickLabels=meanTickLabels, 
                   xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = thisAddColorBar, 
                   legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot= c(.97,1,.1,.9)), legend.width=3, 
                   plotArgs=list(cex.main=3, cex.axis=2, cex.lab=2), legend.mar=0, lwd=.5, border=rgb(.4,.4,.4))
        plotMapDat(mapDat=regionMap, lwd=2.5)
      } else if(thisArea == "Pixel"){
        plot(cbind(popGrid$lon, popGrid$lat), type="n", main="", ylim=kenyaLatRange, 
             xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", asp=1, cex.main=3, cex.axis=2, cex.lab=2)
        quilt.plot(cbind(popGrid$lon, popGrid$lat), logit(predictionList[[j]]), 
                   nx=150, ny=150, add.legend=FALSE, add=TRUE, col=meanCols, zlim=range(logit(meanRange)))
        plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
        plotMapDat(mapDat=regionMap, lwd=2.5)
        points(dat$lon, dat$lat, pch=".")
        
        if(thisAddColorBar) {
          image.plot(zlim=range(logit(meanRange)), nlevel=length(meanCols), legend.only=TRUE, horizontal=FALSE,
                     col=meanCols, add = TRUE, axis.args=list(at=logit(meanTicks), labels=meanTickLabels, cex.axis=2, tck=-.7, hadj=-.1), 
                     legend.mar = 0, legend.cex=2, legend.width=3, smallplot= c(.97,1,.1,.9))
          
          # image.plot(legend.only = TRUE, zlim=c(0,1), nlevel=ncols, legend.mar=7, col=urbCols, add=TRUE, 
          #            legend.lab = "Urbanicity", legend.line=3.0, legend.width=1, legend.shrink=.9, 
          #            legend.cex=2, axis.args=list(cex.axis=1, tck=-1, hadj=-.1))
        }
      }
      
      mtext(side = 3, as.expression(modelName), line = 1, cex=2)
      if(j == 1)
        mtext(side = 2, "Estimates", line = 8, cex=2)
    }
    
    # plot the credible interval widths
    for(j in 1:numberModels) {
      modelName = bquote(.(modelClasses[j])[.(modelVariations[j])])
      thisAddColorBar = (j == numberModels) | !singleColorBar
      
      # load this model and plot results in the given column
      # out = load(resultFilenames[j])
      # theseResults = results$aggregatedResults$predictions[[paste0(tolower(thisArea), "Predictions")]]
      
      if(thisArea == "Region") {
        plotMapDat(plotVar=widthList[[j]], new = TRUE, 
                   main="", scaleFun=log, scaleFunInverse=exp, 
                   cols=widthCols, zlim=log(widthRange), ticks=widthTicks, tickLabels=widthTickLabels, 
                   xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = thisAddColorBar, 
                   legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot= c(.97,1,.1,.9)), legend.width=3, 
                   plotArgs=list(cex.main=3, cex.axis=2, cex.lab=2), legend.mar=0)
      } else if(thisArea == "County") {
        plotMapDat(plotVar=widthList[[j]], new = TRUE, 
                   main="", scaleFun=log, scaleFunInverse=exp, 
                   cols=widthCols, zlim=log(widthRange), ticks=widthTicks, tickLabels=widthTickLabels, 
                   xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = thisAddColorBar, 
                   legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot= c(.97,1,.1,.9)), legend.width=3, 
                   plotArgs=list(cex.main=3, cex.axis=2, cex.lab=2), legend.mar=0, lwd=.5, border=rgb(.4,.4,.4))
        plotMapDat(mapDat=regionMap, lwd=2.5)
      } else if(thisArea == "Pixel") {
        plot(cbind(popGrid$lon, popGrid$lat), type="n", main="", ylim=kenyaLatRange, 
             xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", asp=1, cex.main=3, cex.axis=2, cex.lab=2)
        quilt.plot(cbind(popGrid$lon, popGrid$lat), log(widthList[[j]]), 
                   nx=150, ny=150, add.legend=FALSE, add=TRUE, col=widthCols, zlim=range(log(widthRange)))
        plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
        plotMapDat(mapDat=regionMap, lwd=2.5)
        points(dat$lon, dat$lat, pch=".")
        
        if(thisAddColorBar) {
          image.plot(zlim=range(log(widthRange)), nlevel=length(widthCols), legend.only=TRUE, horizontal=FALSE,
                     col=widthCols, add = TRUE, axis.args=list(at=log(widthTicks), labels=widthTickLabels, cex.axis=2, tck=-.7, hadj=-.1), 
                     legend.mar = 0, legend.cex=2, legend.width=3, smallplot= c(.97,1,.1,.9))
          
          # image.plot(legend.only = TRUE, zlim=c(0,1), nlevel=ncols, legend.mar=7, col=urbCols, add=TRUE, 
          #            legend.lab = "Urbanicity", legend.line=3.0, legend.width=1, legend.shrink=.9, 
          #            legend.cex=2, axis.args=list(cex.axis=1, tck=-1, hadj=-.1))
        }
      }
      
      if(j == 1)
        mtext(side = 2, "80% CI Width", line = 8, cex=2)
    }
    dev.off()
    
    ##### Do the same thing but with relative credible interval widths
    png(file=paste0("Figures/", resultNameRoot, "/predsRelativeWidths", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=1000)
    
    if(thisArea %in% c("Region", "County"))
      par(mfrow=c(2,numberModels), oma=c( 0,6,0,7), mar=c(5.1, 5.1, 4.1, 2.5))
    else
      par(mfrow=c(2,numberModels), oma=c( 0,6,0,7), mar=c(5.1, 5.1, 4.1, 2.5))
    
    # first load in the predictions
    predictionList = list()
    widthList = list()
    for(j in 1:numberModels) {
      modelName = bquote(.(modelClasses[j])[.(modelVariations[j])])
      
      # load this model and plot results in the given column
      out = load(resultFilenames[j])
      theseResults = results$aggregatedResults$predictions[[paste0(tolower(thisArea), "Predictions")]]
      predictionList = c(predictionList, list(theseResults$preds))
      widthList = c(widthList, list((theseResults$Q90 - theseResults$Q10)/theseResults$preds))
    }
    
    # plot the predictions
    for(j in 1:numberModels) {
      modelName = bquote(.(modelClasses[j])[.(modelVariations[j])])
      thisAddColorBar = (j == numberModels) | !singleColorBar
      
      # load this model and plot results in the given column
      # out = load(resultFilenames[j])
      # theseResults = results$aggregatedResults$predictions[[paste0(tolower(thisArea), "Predictions")]]
      
      if(thisArea == "Region") {
        plotMapDat(plotVar=predictionList[[j]], new = TRUE, 
                   main="", scaleFun=logit, scaleFunInverse=expit, 
                   cols=meanCols, zlim=logit(meanRange), ticks=meanTicks, tickLabels=meanTickLabels, 
                   xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = thisAddColorBar, 
                   legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot= c(.97,1,.1,.9)), legend.width=3, 
                   plotArgs=list(cex.main=3, cex.axis=2, cex.lab=2), legend.mar=0)
      } else if(thisArea == "County") {
        plotMapDat(plotVar=predictionList[[j]], new = TRUE, 
                   main="", scaleFun=logit, scaleFunInverse=expit, 
                   cols=meanCols, zlim=logit(meanRange), ticks=meanTicks, tickLabels=meanTickLabels, 
                   xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = thisAddColorBar, 
                   legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot= c(.97,1,.1,.9)), legend.width=3, 
                   plotArgs=list(cex.main=3, cex.axis=2, cex.lab=2), legend.mar=0, lwd=.5, border=rgb(.4,.4,.4))
        plotMapDat(mapDat=regionMap, lwd=2.5)
      } else if(thisArea == "Pixel"){
        plot(cbind(popGrid$lon, popGrid$lat), type="n", main="", ylim=kenyaLatRange, 
             xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", asp=1, cex.main=3, cex.axis=2, cex.lab=2)
        quilt.plot(cbind(popGrid$lon, popGrid$lat), logit(predictionList[[j]]), 
                   nx=150, ny=150, add.legend=FALSE, add=TRUE, col=meanCols, zlim=range(logit(meanRange)))
        plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
        plotMapDat(mapDat=regionMap, lwd=2.5)
        points(dat$lon, dat$lat, pch=".")
        
        if(thisAddColorBar) {
          image.plot(zlim=range(logit(meanRange)), nlevel=length(meanCols), legend.only=TRUE, horizontal=FALSE,
                     col=meanCols, add = TRUE, axis.args=list(at=logit(meanTicks), labels=meanTickLabels, cex.axis=2, tck=-.7, hadj=-.1), 
                     legend.mar = 0, legend.cex=2, legend.width=3, smallplot= c(.97,1,.1,.9))
          
          # image.plot(legend.only = TRUE, zlim=c(0,1), nlevel=ncols, legend.mar=7, col=urbCols, add=TRUE, 
          #            legend.lab = "Urbanicity", legend.line=3.0, legend.width=1, legend.shrink=.9, 
          #            legend.cex=2, axis.args=list(cex.axis=1, tck=-1, hadj=-.1))
        }
      }
      
      mtext(side = 3, as.expression(modelName), line = 1, cex=2)
      if(j == 1)
        mtext(side = 2, "Estimates", line = 8, cex=2)
    }
    
    # plot the relative credible interval widths
    for(j in 1:numberModels) {
      modelName = bquote(.(modelClasses[j])[.(modelVariations[j])])
      thisAddColorBar = (j == numberModels) | !singleColorBar
      
      # load this model and plot results in the given column
      # out = load(resultFilenames[j])
      # theseResults = results$aggregatedResults$predictions[[paste0(tolower(thisArea), "Predictions")]]
      
      if(thisArea == "Region") {
        plotMapDat(plotVar=widthList[[j]], new = TRUE, 
                   main="", scaleFun=log, scaleFunInverse=exp, 
                   cols=widthCols, zlim=log(relWidthRange), ticks=relWidthTicks, tickLabels=relWidthTickLabels, 
                   xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = thisAddColorBar, 
                   legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot= c(.97,1,.1,.9)), legend.width=3, 
                   plotArgs=list(cex.main=3, cex.axis=2, cex.lab=2), legend.mar=0)
      } else if(thisArea == "County") {
        plotMapDat(plotVar=widthList[[j]], new = TRUE, 
                   main="", scaleFun=log, scaleFunInverse=exp, 
                   cols=widthCols, zlim=log(relWidthRange), ticks=relWidthTicks, tickLabels=relWidthTickLabels, 
                   xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = thisAddColorBar, 
                   legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot= c(.97,1,.1,.9)), legend.width=3, 
                   plotArgs=list(cex.main=3, cex.axis=2, cex.lab=2), legend.mar=0, lwd=.5, border=rgb(.4,.4,.4))
        plotMapDat(mapDat=regionMap, lwd=2.5)
      } else if(thisArea == "Pixel") {
        plot(cbind(popGrid$lon, popGrid$lat), type="n", main="", ylim=kenyaLatRange, 
             xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", asp=1, cex.main=3, cex.axis=2, cex.lab=2)
        quilt.plot(cbind(popGrid$lon, popGrid$lat), log(widthList[[j]]), 
                   nx=150, ny=150, add.legend=FALSE, add=TRUE, col=widthCols, zlim=range(log(relWidthRange)))
        plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
        plotMapDat(mapDat=regionMap, lwd=2.5)
        points(dat$lon, dat$lat, pch=".")
        
        if(thisAddColorBar) {
          image.plot(zlim=range(log(relWidthRange)), nlevel=length(widthCols), legend.only=TRUE, horizontal=FALSE,
                     col=widthCols, add = TRUE, axis.args=list(at=log(relWidthTicks), labels=relWidthTickLabels, cex.axis=2, tck=-.7, hadj=-.1), 
                     legend.mar = 0, legend.cex=2, legend.width=3, smallplot= c(.97,1,.1,.9))
          
          # image.plot(legend.only = TRUE, zlim=c(0,1), nlevel=ncols, legend.mar=7, col=urbCols, add=TRUE, 
          #            legend.lab = "Urbanicity", legend.line=3.0, legend.width=1, legend.shrink=.9, 
          #            legend.cex=2, axis.args=list(cex.axis=1, tck=-1, hadj=-.1))
        }
      }
      
      if(j == 1)
        mtext(side = 2, "Relative 80% CI Width", line = 8, cex=2)
    }
    dev.off()
  }
}

# this function plots central estimates and credible interval widths over a map of Kenya 
# 4 plots are made by default: one for each area level
makePairPlots = function(dat, resultFilenames, modelClasses, modelVariations, 
                         areaLevels=c("Region", "County", "Pixel"), 
                         meanRange, meanTicks, meanTickLabels, 
                         plotNameRoot="Education", resultNameRoot="Ed", meanCols=makeRedBlueDivergingColors(64), 
                         ncols=29, urbCols=makeGreenBlueSequentialColors(ncols), 
                         kenyaLatRange=c(-4.6, 5), kenyaLonRange=c(33.5, 42.0), doModelClassPlots=FALSE, 
                         areaMat=NULL, radiiMat=NULL) {
  plotNameRootLower = tolower(plotNameRoot)
  resultNameRootLower = tolower(resultNameRoot)
  numberModels = length(resultFilenames)
  uniqueModelClasses = unique(modelClasses)
  
  ##### if there are multiple unique model classes, call this function on each one individually 
  ##### as well as on the combination
  if(length(uniqueModelClasses) != 1 && doModelClassPlots) {
    for(i in 1:length(uniqueModelClasses)) {
      thisModelClass = uniqueModelClasses[i]
      thisI = modelClasses == thisModelClass
      
      if(sum(thisI) == 1)
        next
      
      thisResultFilenames = resultFilenames[thisI]
      thisModelClasses = modelClasses[thisI]
      thisModelVariations = modelVariations[thisI]
      makePairPlots(dat, thisResultFilenames, thisModelClasses, thisModelVariations, 
                    areaLevels, meanRange, meanTicks, meanTickLabels, 
                    plotNameRoot=paste0(plotNameRoot, thisModelClass), 
                    resultNameRoot, meanCols, ncols, urbCols, kenyaLatRange, kenyaLonRange)
    }
  }
  
  ##### central estimates and credible interval widths
  
  for(i in 1:length(areaLevels)) {
    thisArea = areaLevels[i]
    
    if(length(uniqueModelClasses) == 1)
      extraPlotNameRoot = uniqueModelClasses[1]
    else
      extraPlotNameRoot = ""
    
    if(thisArea %in% c("Region", "County")) {
      
    }
    else {
      
    }
    
    # collect predictions and proportion urban per area/point
    predsList = list()
    widthsList = list()
    # modelNames = c()
    modelNames = list()
    for(j in 1:numberModels) {
      # modelNames = c(modelNames, paste(modelClasses[j], modelVariations[j]))
      modelNames = c(modelNames, bquote(.(modelClasses[j])[.(modelVariations[j])]))
      
      # load this model and put results in the given column
      out = load(resultFilenames[j])
      theseResults = results$aggregatedResults$predictions[[paste0(tolower(thisArea), "Predictions")]]
      
      # get proportion urban (will be the same for each model)
      colI = cut(as.numeric(theseResults$urban), breaks=seq(0 - .0001, 1 + .0001, l=ncols+1), labels=FALSE)
      theseCols = urbCols[colI]
      
      predsList = c(predsList, list(theseResults$preds))
      widthsList = c(widthsList, list(theseResults$Q90 - theseResults$Q10))
    }
    
    valMat = do.call("cbind", predsList)
    zlim = range(valMat)
    
    if(thisArea %in% c("Region", "County")) {
      
      # valMat = rbind(1:5, valMat)
      my_line <- function(x,y,...){
        if(diff(range(x)) >= .04)
          xlim = zlim
        # else
        #   xlim = zlim2
        if(diff(range(y)) >= .04)
          ylim = zlim
        # else
        #   ylim = zlim2
        # if(diff(range(c(x, y))) > 0.04)
        #   par(usr = c(zlim, zlim))
        # else
        #   par(usr = c(zlim2, zlim2))
        # par(usr = c(xlim, ylim))
        # points(x,y,..., col="blue")
        abline(a = 0,b = 1,...)
        points(x,y,..., col=theseCols)
      }
      
      width = 400 * numberModels
      png(file=paste0("Figures/", resultNameRoot, "/pairPlot", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=width)
      
      # pairs(valMat, 
      #       modelNames, 
      #       pch=19, cex=.3, lower.panel=my_line, upper.panel = my_line, 
      #       main=paste0(thisArea, " estimate comparisons"))
      # lims = c(list(zlim), list(zlim), list(zlim2), list(zlim2), list(zlim2))
      lims = rep(list(zlim), numberModels)
      myPairs(valMat, 
              as.expression(unlist(modelNames)), 
              pch=19, cex=2, lower.panel=my_line, upper.panel = my_line, 
              main=paste0(thisArea, " estimate comparisons"), 
              lims=lims, oma=c(3,3,6,10), cex.main=2)
      image.plot(legend.only = TRUE, zlim=c(0,1), nlevel=ncols, legend.mar=7, col=urbCols, add=TRUE, 
                 legend.lab = "Urbanicity", legend.line=3.0, legend.width=1, legend.shrink=.9, 
                 legend.cex=2, axis.args=list(cex.axis=1, tck=-1, hadj=-.1))
      dev.off()
      
      valMat = do.call("cbind", widthsList)
      zlim = range(valMat)
      width = 400 * numberModels
      png(file=paste0("Figures/", resultNameRoot, "/pairPlotWidths", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=width)
      
      # pairs(valMat, 
      #       modelNames, 
      #       pch=19, cex=.3, lower.panel=my_line, upper.panel = my_line, 
      #       main=paste0(thisArea, " estimate comparisons"))
      # lims = c(list(zlim), list(zlim), list(zlim2), list(zlim2), list(zlim2))
      lims = rep(list(zlim), numberModels)
      myPairs(valMat, 
              as.expression(unlist(modelNames)), 
              pch=19, cex=2, lower.panel=my_line, upper.panel = my_line, 
              main=paste0(thisArea, " posterior 80% CI width comparisons"), 
              lims=lims, oma=c(3,3,6,10), cex.main=2)
      image.plot(legend.only = TRUE, zlim=c(0,1), nlevel=ncols, legend.mar=7, col=urbCols, add=TRUE, 
                 legend.lab = "Urbanicity", legend.line=3.0, legend.width=1, legend.shrink=.9, 
                 legend.cex=2, axis.args=list(cex.axis=1, tck=-1, hadj=-.1))
      dev.off()
    } else {
      my_line <- function(x,y,...){
        if(diff(range(x)) >= .04)
          xlim = zlim
        # else
        #   xlim = zlim2
        if(diff(range(y)) >= .04)
          ylim = zlim
        # else
        #   ylim = zlim2
        # if(diff(range(c(x, y))) > 0.04)
        #   par(usr = c(zlim, zlim))
        # else
        #   par(usr = c(zlim2, zlim2))
        # par(usr = c(xlim, ylim))
        # points(x,y,..., col="blue")
        abline(a = 0,b = 1,...)
        points(x,y,..., col=theseCols)
      }
      
      width = 2 * numberModels
      pdf(file=paste0("Figures/", resultNameRoot, "/pairPlot", plotNameRoot, extraPlotNameRoot, thisArea, ".pdf"), width=width, height=width)
      
      lims = rep(list(zlim), numberModels)
      myPairs(valMat, 
              as.expression(unlist(modelNames)), 
              pch=19, cex=.4, lower.panel=my_line, upper.panel = my_line, 
              main=paste0(thisArea, " estimate comparisons"), 
              lims=lims, oma=c(3,3,6,7))
      # legend("topleft", c("Urban", "Rural"), col=c(urbCols[ncols], urbCols[1]), pch=19)
      dev.off()
      
      valMat = do.call("cbind", widthsList)
      zlim = range(valMat)
      pdf(file=paste0("Figures/", resultNameRoot, "/pairPlotWidths", plotNameRoot, extraPlotNameRoot, thisArea, ".pdf"), width=width, height=width)
      
      lims = rep(list(zlim), numberModels)
      myPairs(valMat, 
              as.expression(unlist(modelNames)), 
              pch=19, cex=.4, lower.panel=my_line, upper.panel = my_line, 
              main=paste0(thisArea, " posterior 80% CI width comparisons"), 
              lims=lims, oma=c(3,3,6,7))
      legend("topleft", c("Urban", "Rural"), col=c(urbCols[ncols], urbCols[1]), pch=19)
      dev.off()
    }
  }
}

# this function plots central estimates and credible interval widths over a map of Kenya 
# 4 plots are made by default: one for each area level
makeFinalPairPlot = function(dat, resultFilenames, modelClasses, modelVariations, 
                         areaLevels=c("Region", "County", "Pixel"), 
                         meanRange, meanTicks, meanTickLabels, 
                         plotNameRoot="Education", resultNameRoot="Ed", meanCols=makeRedBlueDivergingColors(64), 
                         ncols=29, urbCols=makeGreenBlueSequentialColors(ncols), 
                         kenyaLatRange=c(-4.6, 5), kenyaLonRange=c(33.5, 42.0), doModelClassPlots=FALSE, 
                         areaMat=NULL, radiiMat=NULL, areaColors=makePurpleYellowSequentialColors(ncols)) {
  plotNameRootLower = tolower(plotNameRoot)
  resultNameRootLower = tolower(resultNameRoot)
  numberModels = length(resultFilenames)
  uniqueModelClasses = unique(modelClasses)
  
  ##### central estimates and credible interval widths
  
  for(i in 1:length(areaLevels)) {
    thisArea = areaLevels[i]
    
    if(length(uniqueModelClasses) == 1)
      extraPlotNameRoot = uniqueModelClasses[1]
    else
      extraPlotNameRoot = ""
    
    # collect predictions and proportion urban per area/point
    predsList = list()
    widthsList = list()
    # modelNames = c()
    modelNames = list()
    for(j in 1:numberModels) {
      # modelNames = c(modelNames, paste(modelClasses[j], modelVariations[j]))
      modelNames = c(modelNames, bquote(.(modelClasses[j])[.(modelVariations[j])]))
      
      # load this model and put results in the given column
      out = load(resultFilenames[j])
      theseResults = results$aggregatedResults$predictions[[paste0(tolower(thisArea), "Predictions")]]
      
      # get proportion urban (will be the same for each model)
      colI = cut(as.numeric(theseResults$urban), breaks=seq(0 - .0001, 1 + .0001, l=ncols+1), labels=FALSE)
      theseCols = urbCols[colI]
      
      predsList = c(predsList, list(theseResults$preds))
      widthsList = c(widthsList, list(theseResults$Q90 - theseResults$Q10))
    }
    
    valMat = do.call("cbind", predsList)
    zlim = range(valMat)
    
    if(thisArea %in% c("Region", "County")) {
      # determine colorings by area and radius
      theseAreas = getArea(thisArea)
      theseRadii = getRadius(thisArea)
      colI = cut(theseAreas, breaks=seq(min(theseAreas) - .0001, max(theseAreas) + .0001, l=length(areaColors)+1), labels=FALSE)
      theseAreaColors = areaColors[colI]
      colI = cut(theseRadii, breaks=seq(min(theseRadii) - .0001, max(theseRadii) + .0001, l=length(areaColors)+1), labels=FALSE)
      theseRadiusColors = areaColors[colI]
      theseUrbanColors = theseCols
      
      # valMat = rbind(1:5, valMat)
      if(thisArea == "Region")
        cex = 2.5
      else
        cex = 1
      my_line <- function(x,y,...){
        if(diff(range(x)) >= .04)
          xlim = zlim
        # else
        #   xlim = zlim2
        if(diff(range(y)) >= .04)
          ylim = zlim
        
        points(x,y,..., col=theseCols,cex=cex)
        abline(a = 0,b = 1,...)
      }
      
      width = 400 * (numberModels - 1)
      height = 400
      png(file=paste0("Figures/", resultNameRoot, "/finalPairPlot", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,5), mar=c(5.1, 4.1, 4.1, 6))
      
      for(i in 1:(numberModels-1)) {
        x = valMat[,i]
        y = valMat[,numberModels]
        ylab = ""
        
        plot(x, y, main="", type="n", xlab="", ylab="", 
             ylim=zlim, xlim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1)
        if(i == 1)
          mtext(side = 2, as.expression(modelNames[[numberModels]]), line = 4, cex=2)
        my_line(x, y, lwd=1, pch=19)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      image.plot(legend.only = TRUE, zlim=c(0,1), nlevel=ncols, legend.mar=4, col=urbCols, add=TRUE, 
                 legend.lab = "Urbanicity", legend.line=5.0, legend.width=3, legend.shrink=.9, 
                 legend.cex=1.5, axis.args=list(cex.axis=1.5, tck=-1, hadj=-.1))
      dev.off()
      
      png(file=paste0("Figures/", resultNameRoot, "/finalPairPlotAreaCol", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,5), mar=c(5.1, 4.1, 4.1, 6))
      theseCols = theseAreaColors
      
      for(i in 1:(numberModels-1)) {
        x = valMat[,i]
        y = valMat[,numberModels]
        ylab = ""
        
        plot(x, y, main="", type="n", xlab="", ylab="", 
             ylim=zlim, xlim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1)
        if(i == 1)
          mtext(side = 2, as.expression(modelNames[[numberModels]]), line = 4, cex=2)
        my_line(x, y, lwd=1, pch=19)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      image.plot(legend.only = TRUE, zlim=range(theseAreas), nlevel=ncols, legend.mar=4, col=areaColors, add=TRUE, 
                 legend.lab = "Area (sq. km)", legend.line=5.0, legend.width=3, legend.shrink=.9, 
                 legend.cex=1.5, axis.args=list(cex.axis=1.5, tck=-1, hadj=-.1))
      dev.off()
      
      png(file=paste0("Figures/", resultNameRoot, "/finalPairPlotRadiusCol", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,5), mar=c(5.1, 4.1, 4.1, 6))
      theseCols = theseRadiusColors
      
      for(i in 1:(numberModels-1)) {
        x = valMat[,i]
        y = valMat[,numberModels]
        ylab = ""
        
        plot(x, y, main="", type="n", xlab="", ylab="", 
             ylim=zlim, xlim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1)
        if(i == 1)
          mtext(side = 2, as.expression(modelNames[[numberModels]]), line = 4, cex=2)
        my_line(x, y, lwd=1, pch=19)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      image.plot(legend.only = TRUE, zlim=range(theseRadii), nlevel=ncols, legend.mar=4, col=areaColors, add=TRUE, 
                 legend.lab = "Radius (km)", legend.line=5.0, legend.width=3, legend.shrink=.9, 
                 legend.cex=1.5, axis.args=list(cex.axis=1.5, tck=-1, hadj=-.1))
      dev.off()
      
      theseCols = theseUrbanColors
      valMat = do.call("cbind", lapply(1:length(widthsList), function(x) {widthsList[[x]] / predsList[[x]]}))
      zlim = range(valMat)
      width = 400 * (numberModels - 1)
      height = 400
      png(file=paste0("Figures/", resultNameRoot, "/finalPairPlotRelWidths", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,5), mar=c(5.1, 4.1, 4.1, 6))
      
      for(i in 1:(numberModels-1)) {
        x = valMat[,i]
        y = valMat[,numberModels]
        ylab = ""
        
        plot(x, y, main="", type="n", xlab="", ylab="", 
             ylim=zlim, xlim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1)
        if(i == 1)
          mtext(side = 2, as.expression(modelNames[[numberModels]]), line = 4, cex=2)
        my_line(x, y, lwd=1, pch=19)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      image.plot(legend.only = TRUE, zlim=c(0,1), nlevel=ncols, legend.mar=4, col=urbCols, add=TRUE, 
                 legend.lab = "Urbanicity", legend.line=5.0, legend.width=3, legend.shrink=.9, 
                 legend.cex=1.5, axis.args=list(cex.axis=1.5, tck=-1, hadj=-.1))
      dev.off()
      
      theseCols = theseAreaColors
      png(file=paste0("Figures/", resultNameRoot, "/finalPairPlotRelWidthsAreaCol", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,5), mar=c(5.1, 4.1, 4.1, 6))
      
      for(i in 1:(numberModels-1)) {
        x = valMat[,i]
        y = valMat[,numberModels]
        ylab = ""
        
        plot(x, y, main="", type="n", xlab="", ylab="", 
             ylim=zlim, xlim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1)
        if(i == 1)
          mtext(side = 2, as.expression(modelNames[[numberModels]]), line = 4, cex=2)
        my_line(x, y, lwd=1, pch=19)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      image.plot(legend.only = TRUE, zlim=range(theseAreas), nlevel=ncols, legend.mar=4, col=areaColors, add=TRUE, 
                 legend.lab = "Area (sq. km)", legend.line=5.0, legend.width=3, legend.shrink=.9, 
                 legend.cex=1.5, axis.args=list(cex.axis=1.5, tck=-1, hadj=-.1))
      dev.off()
      
      theseCols = theseRadiusColors
      png(file=paste0("Figures/", resultNameRoot, "/finalPairPlotRelWidthsRadiusCol", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,5), mar=c(5.1, 4.1, 4.1, 6))
      
      for(i in 1:(numberModels-1)) {
        x = valMat[,i]
        y = valMat[,numberModels]
        ylab = ""
        
        plot(x, y, main="", type="n", xlab="", ylab="", 
             ylim=zlim, xlim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1)
        if(i == 1)
          mtext(side = 2, as.expression(modelNames[[numberModels]]), line = 4, cex=2)
        my_line(x, y, lwd=1, pch=19)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      image.plot(legend.only = TRUE, zlim=range(theseRadii), nlevel=ncols, legend.mar=4, col=areaColors, add=TRUE, 
                 legend.lab = "Radius (km)", legend.line=5.0, legend.width=3, legend.shrink=.9, 
                 legend.cex=1.5, axis.args=list(cex.axis=1.5, tck=-1, hadj=-.1))
      dev.off()
      
      valMat = do.call("cbind", widthsList)
      zlim = range(valMat)
      width = 400 * (numberModels - 1)
      height = 400
      theseCols = theseUrbanColors
      png(file=paste0("Figures/", resultNameRoot, "/finalPairPlotWidths", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,5), mar=c(5.1, 4.1, 4.1, 6))
      
      for(i in 1:(numberModels-1)) {
        x = valMat[,i]
        y = valMat[,numberModels]
        ylab = ""
        
        plot(x, y, main="", type="n", xlab="", ylab="", 
             ylim=zlim, xlim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1)
        if(i == 1)
          mtext(side = 2, as.expression(modelNames[[numberModels]]), line = 4, cex=2)
        my_line(x, y, lwd=1, pch=19)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      image.plot(legend.only = TRUE, zlim=c(0,1), nlevel=ncols, legend.mar=4, col=urbCols, add=TRUE, 
                 legend.lab = "Urbanicity", legend.line=5.0, legend.width=3, legend.shrink=.9, 
                 legend.cex=1.5, axis.args=list(cex.axis=1.5, tck=-1, hadj=-.1))
      dev.off()
      
      theseCols = theseAreaColors
      png(file=paste0("Figures/", resultNameRoot, "/finalPairPlotWidthsAreaCol", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,5), mar=c(5.1, 4.1, 4.1, 6))
      
      for(i in 1:(numberModels-1)) {
        x = valMat[,i]
        y = valMat[,numberModels]
        ylab = ""
        
        plot(x, y, main="", type="n", xlab="", ylab="", 
             ylim=zlim, xlim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1)
        if(i == 1)
          mtext(side = 2, as.expression(modelNames[[numberModels]]), line = 4, cex=2)
        my_line(x, y, lwd=1, pch=19)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      image.plot(legend.only = TRUE, zlim=range(theseAreas), nlevel=ncols, legend.mar=4, col=areaColors, add=TRUE, 
                 legend.lab = "Area (sq. km)", legend.line=5.0, legend.width=3, legend.shrink=.9, 
                 legend.cex=1.5, axis.args=list(cex.axis=1.5, tck=-1, hadj=-.1))
      dev.off()
      
      theseCols = theseRadiusColors
      png(file=paste0("Figures/", resultNameRoot, "/finalPairPlotWidthsRadiusCol", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,5), mar=c(5.1, 4.1, 4.1, 6))
      
      for(i in 1:(numberModels-1)) {
        x = valMat[,i]
        y = valMat[,numberModels]
        ylab = ""
        
        plot(x, y, main="", type="n", xlab="", ylab="", 
             ylim=zlim, xlim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1)
        if(i == 1)
          mtext(side = 2, as.expression(modelNames[[numberModels]]), line = 4, cex=2)
        my_line(x, y, lwd=1, pch=19)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      image.plot(legend.only = TRUE, zlim=range(theseRadii), nlevel=ncols, legend.mar=4, col=areaColors, add=TRUE, 
                 legend.lab = "Radius (km)", legend.line=5.0, legend.width=3, legend.shrink=.9, 
                 legend.cex=1.5, axis.args=list(cex.axis=1.5, tck=-1, hadj=-.1))
      dev.off()
    } else {
      thesePointTypes = sapply(theseResults$urban, function(x) {ifelse(x, 15, 17)})
      
      my_line <- function(x,y,...){
        if(diff(range(x)) >= .04)
          xlim = zlim
        
        if(diff(range(y)) >= .04)
          ylim = zlim
        
        # abline(a = 0,b = 1, col=rgb(.4, .4, .4),...)
        points(x,y,..., col=theseCols, pch=thesePointTypes)
        abline(a = 0,b = 1, ...)
      }
      
      width = 4 * (numberModels - 1)
      height = 4
      pdf(file=paste0("Figures/", resultNameRoot, "/finalPairPlot", plotNameRoot, extraPlotNameRoot, thisArea, ".pdf"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,0), mar=c(5.1, 4.1, 4.1, 4.1))
      
      zlim = range(valMat)
      for(i in 1:(numberModels-1)) {
        x = valMat[,i]
        y = valMat[,numberModels]
        ylab = ""
        
        plot(x, y, main="", type="n", xlab="", ylab="", 
             ylim=zlim, xlim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1)
        if(i == 1)
          mtext(side = 2, as.expression(modelNames[[numberModels]]), line = 4, cex=2)
        my_line(x, y, lwd=1, cex=.8)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      dev.off()
      
      valMat = do.call("cbind", lapply(1:length(widthsList), function(x) {widthsList[[x]] / predsList[[x]]}))
      pdf(file=paste0("Figures/", resultNameRoot, "/finalPairPlotRelWidths", plotNameRoot, extraPlotNameRoot, thisArea, ".pdf"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,0), mar=c(5.1, 4.1, 4.1, 4.1))
      
      zlim = range(valMat)
      for(i in 1:(numberModels-1)) {
        x = valMat[,i]
        y = valMat[,numberModels]
        ylab = ""
        
        plot(x, y, main="", type="n", xlab="", ylab="", 
             ylim=zlim, xlim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1)
        if(i == 1)
          mtext(side = 2, as.expression(modelNames[[numberModels]]), line = 4, cex=2)
        my_line(x, y, lwd=1, cex=.8)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      dev.off()
      
      valMat = do.call("cbind", widthsList)
      zlim = range(valMat)
      pdf(file=paste0("Figures/", resultNameRoot, "/finalPairPlotWidths", plotNameRoot, extraPlotNameRoot, thisArea, ".pdf"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,0), mar=c(5.1, 4.1, 4.1, 4.1))
      
      for(i in 1:(numberModels-1)) {
        x = valMat[,i]
        y = valMat[,numberModels]
        ylab = ""
        
        plot(x, y, main="", type="n", xlab="", ylab="", 
             ylim=zlim, xlim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1)
        if(i == 1)
          mtext(side = 2, as.expression(modelNames[[numberModels]]), line = 4, cex=2)
        my_line(x, y, lwd=1, cex=.8)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      dev.off()
    }
  }
}

# this function plots central estimates and credible interval widths over a map of Kenya 
# 4 plots are made by default: one for each area level
makeFinalPercentResidualPlot = function(dat, resultFilenames, modelClasses, modelVariations, 
                                        areaLevels=c("Region", "County", "Pixel"), 
                                        meanRange, meanTicks, meanTickLabels, 
                                        plotNameRoot="Education", resultNameRoot="Ed", meanCols=makeRedBlueDivergingColors(64), 
                                        ncols=29, urbCols=makeGreenBlueSequentialColors(ncols), 
                                        kenyaLatRange=c(-4.6, 5), kenyaLonRange=c(33.5, 42.0), doModelClassPlots=FALSE, 
                                        areaMat=NULL, radiiMat=NULL, areaColors=makePurpleYellowSequentialColors(ncols)) {
  plotNameRootLower = tolower(plotNameRoot)
  resultNameRootLower = tolower(resultNameRoot)
  numberModels = length(resultFilenames)
  uniqueModelClasses = unique(modelClasses)
  
  ##### central estimates and credible interval widths
  
  for(i in 1:length(areaLevels)) {
    thisArea = areaLevels[i]
    
    if(length(uniqueModelClasses) == 1)
      extraPlotNameRoot = uniqueModelClasses[1]
    else
      extraPlotNameRoot = ""
    
    # collect predictions and proportion urban per area/point
    predsList = list()
    widthsList = list()
    # modelNames = c()
    modelNames = list()
    for(j in 1:numberModels) {
      # modelNames = c(modelNames, paste(modelClasses[j], modelVariations[j]))
      modelNames = c(modelNames, bquote(.(modelClasses[j])[.(modelVariations[j])]))
      
      # load this model and put results in the given column
      out = load(resultFilenames[j])
      theseResults = results$aggregatedResults$predictions[[paste0(tolower(thisArea), "Predictions")]]
      
      # get proportion urban (will be the same for each model)
      colI = cut(as.numeric(theseResults$urban), breaks=seq(0 - .0001, 1 + .0001, l=ncols+1), labels=FALSE)
      theseCols = urbCols[colI]
      
      predsList = c(predsList, list(theseResults$preds))
      widthsList = c(widthsList, list(theseResults$Q90 - theseResults$Q10))
    }
    
    valMat = do.call("cbind", predsList)
    zlim = range(sweep(valMat[,-ncol(valMat)], 1, valMat[,ncol(valMat)], function(x,y) {100*(x-y)/y}))
    
    if(thisArea %in% c("Region", "County")) {
      # determine colorings by area, radius, and sparsity
      theseAreas = getArea(thisArea)
      theseRadii = getRadius(thisArea)
      theseSparsities = getAreaPerObservation(thisArea)
      
      colI = cut(theseAreas, breaks=seq(min(theseAreas) - .0001, max(theseAreas) + .0001, l=length(areaColors)+1), labels=FALSE)
      theseAreaColors = areaColors[colI]
      colI = cut(theseRadii, breaks=seq(min(theseRadii) - .0001, max(theseRadii) + .0001, l=length(areaColors)+1), labels=FALSE)
      theseRadiusColors = areaColors[colI]
      colI = cut(theseSparsities, breaks=seq(min(theseSparsities) - .0001, max(theseSparsities) + .0001, l=length(areaColors)+1), labels=FALSE)
      theseSparsityColors = areaColors[colI]
      theseUrbanColors = theseCols
      sparsityTicks = getAreaPerObservationTicks(thisArea)
      sparsityTickLabels = as.character(sparsityTicks)
      sparsityTicks = log(sparsityTicks)
      
      # valMat = rbind(1:5, valMat)
      if(thisArea == "Region")
        cex = 2.5
      else
        cex = 1
      my_line <- function(x,y,...){
        if(diff(range(x)) >= .04)
          xlim = zlim
        # else
        #   xlim = zlim2
        if(diff(range(y)) >= .04)
          ylim = zlim
        
        points(x,y,..., col=theseCols,cex=cex)
        abline(h=0, lty=2,...)
      }
      
      width = 400 * (numberModels - 1)
      height = 400
      png(file=paste0("Figures/", resultNameRoot, "/finalPercentResidualPlot", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,5), mar=c(5.1, 4.1, 4.1, 6))
      
      for(i in 1:(numberModels-1)) {
        y = valMat[,i]
        x = valMat[,numberModels]
        y = 100 * (x-y)/y
        ylab = ""
        
        plot(x, y, main="", type="n", xlab=bquote(.(modelNames[[numberModels]]) ~ " Estimates"), ylab="", 
             ylim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1, log="x")
        if(i == 1)
          mtext(side = 2, bquote("Pct. Diff. from " ~ .(modelNames[[numberModels]])), line = 4, cex=2)
        my_line(x, y, lwd=1, pch=19)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      image.plot(legend.only = TRUE, zlim=c(0,1), nlevel=ncols, legend.mar=4, col=urbCols, add=TRUE, 
                 legend.lab = "Urbanicity", legend.line=5.0, legend.width=3, legend.shrink=.9, 
                 legend.cex=1.5, axis.args=list(cex.axis=1.5, tck=-1, hadj=-.1))
      dev.off()
      
      png(file=paste0("Figures/", resultNameRoot, "/finalPercentResidualPlotAreaCol", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,5), mar=c(5.1, 4.1, 4.1, 6))
      theseCols = theseAreaColors
      
      for(i in 1:(numberModels-1)) {
        y = valMat[,i]
        x = valMat[,numberModels]
        y = 100 * (x-y)/y
        ylab = ""
        
        plot(x, y, main="", type="n", xlab=bquote(.(modelNames[[numberModels]]) ~ " Estimates"), ylab="", 
             ylim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1, log="x")
        if(i == 1)
          mtext(side = 2, bquote("Pct. Diff. from " ~ .(modelNames[[numberModels]])), line = 4, cex=2)
        my_line(x, y, lwd=1, pch=19)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      image.plot(legend.only = TRUE, zlim=range(theseAreas), nlevel=ncols, legend.mar=4, col=areaColors, add=TRUE, 
                 legend.lab = "Area (sq. km)", legend.line=5.0, legend.width=3, legend.shrink=.9, 
                 legend.cex=1.5, axis.args=list(cex.axis=1.5, tck=-1, hadj=-.1))
      dev.off()
      
      png(file=paste0("Figures/", resultNameRoot, "/finalPercentResidualPlotRadiusCol", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,5), mar=c(5.1, 4.1, 4.1, 6))
      theseCols = theseRadiusColors
      
      for(i in 1:(numberModels-1)) {
        y = valMat[,i]
        x = valMat[,numberModels]
        y = 100 * (x-y)/y
        ylab = ""
        
        plot(x, y, main="", type="n", xlab=bquote(.(modelNames[[numberModels]]) ~ " Estimates"), ylab="", 
             ylim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1, log="x")
        if(i == 1)
          mtext(side = 2, bquote("Pct. Diff. from " ~ .(modelNames[[numberModels]])), line = 4, cex=2)
        my_line(x, y, lwd=1, pch=19)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      image.plot(legend.only = TRUE, zlim=range(theseRadii), nlevel=ncols, legend.mar=4, col=areaColors, add=TRUE, 
                 legend.lab = "Radius (km)", legend.line=5.0, legend.width=3, legend.shrink=.9, 
                 legend.cex=1.5, axis.args=list(cex.axis=1.5, tck=-1, hadj=-.1))
      dev.off()
      
      png(file=paste0("Figures/", resultNameRoot, "/finalPercentResidualPlotSparsityCol", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,5), mar=c(5.1, 4.1, 4.1, 6))
      theseCols = theseSparsityColors
      
      for(i in 1:(numberModels-1)) {
        y = valMat[,i]
        x = valMat[,numberModels]
        y = 100 * (x-y)/y
        ylab = ""
        
        plot(x, y, main="", type="n", xlab=bquote(.(modelNames[[numberModels]]) ~ " Estimates"), ylab="", 
             ylim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1, log="x")
        if(i == 1)
          mtext(side = 2, bquote("Pct. Diff. from " ~ .(modelNames[[numberModels]])), line = 4, cex=2)
        my_line(x, y, lwd=1, pch=19)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      image.plot(legend.only = TRUE, zlim=range(theseSparsities), nlevel=ncols, legend.mar=4, col=areaColors, add=TRUE, 
                 legend.lab = expression("Sparsity (km"^2 ~" Per Cluster)"), legend.line=5.0, legend.width=3, legend.shrink=.9, 
                 legend.cex=1.5, axis.args=list(at=sparsityTicks, labels=sparsityTickLabels, cex.axis=1.5, tck=-1, hadj=-.1))
      dev.off()
      
      theseCols = theseUrbanColors
      valMat = do.call("cbind", lapply(1:length(widthsList), function(x) {widthsList[[x]] / predsList[[x]]}))
      zlim = range(sweep(valMat[,-ncol(valMat)], 1, valMat[,ncol(valMat)], function(x,y) {100*(x-y)/y}))
      width = 400 * (numberModels - 1)
      height = 400
      png(file=paste0("Figures/", resultNameRoot, "/finalPercentResidualPlotRelWidths", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,5), mar=c(5.1, 4.1, 4.1, 6))
      
      for(i in 1:(numberModels-1)) {
        x = valMat[,i]
        y = valMat[,numberModels]
        y = 100 * (x-y)/y
        ylab = ""
        
        plot(x, y, main="", type="n", xlab=bquote(.(modelNames[[numberModels]]) ~ " Rel. Width"), ylab="", 
             ylim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1, log="x")
        if(i == 1)
          mtext(side = 2, bquote("Pct. Diff. from " ~ .(modelNames[[numberModels]])), line = 4, cex=2)
        my_line(x, y, lwd=1, pch=19)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      image.plot(legend.only = TRUE, zlim=c(0,1), nlevel=ncols, legend.mar=4, col=urbCols, add=TRUE, 
                 legend.lab = "Urbanicity", legend.line=5.0, legend.width=3, legend.shrink=.9, 
                 legend.cex=1.5, axis.args=list(cex.axis=1.5, tck=-1, hadj=-.1))
      dev.off()
      
      theseCols = theseAreaColors
      png(file=paste0("Figures/", resultNameRoot, "/finalPercentResidualPlotRelWidthsAreaCol", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,5), mar=c(5.1, 4.1, 4.1, 6))
      
      for(i in 1:(numberModels-1)) {
        x = valMat[,i]
        y = valMat[,numberModels]
        y = 100 * (x-y)/y
        ylab = ""
        
        plot(x, y, main="", type="n", xlab=bquote(.(modelNames[[numberModels]]) ~ " Rel. Width"), ylab="", 
             ylim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1, log="x")
        if(i == 1)
          mtext(side = 2, bquote("Pct. Diff. from " ~ .(modelNames[[numberModels]])), line = 4, cex=2)
        my_line(x, y, lwd=1, pch=19)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      image.plot(legend.only = TRUE, zlim=range(theseAreas), nlevel=ncols, legend.mar=4, col=areaColors, add=TRUE, 
                 legend.lab = "Area (sq. km)", legend.line=5.0, legend.width=3, legend.shrink=.9, 
                 legend.cex=1.5, axis.args=list(cex.axis=1.5, tck=-1, hadj=-.1))
      dev.off()
      
      theseCols = theseRadiusColors
      png(file=paste0("Figures/", resultNameRoot, "/finalPercentResidualPlotRelWidthsRadiusCol", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,5), mar=c(5.1, 4.1, 4.1, 6))
      
      for(i in 1:(numberModels-1)) {
        x = valMat[,i]
        y = valMat[,numberModels]
        y = 100 * (x-y)/y
        ylab = ""
        
        plot(x, y, main="", type="n", xlab=bquote(.(modelNames[[numberModels]]) ~ " Rel. Width"), ylab="", 
             ylim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1, log="x")
        if(i == 1)
          mtext(side = 2, bquote("Pct. Diff. from " ~ .(modelNames[[numberModels]])), line = 4, cex=2)
        my_line(x, y, lwd=1, pch=19)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      image.plot(legend.only = TRUE, zlim=range(theseRadii), nlevel=ncols, legend.mar=4, col=areaColors, add=TRUE, 
                 legend.lab = "Radius (km)", legend.line=5.0, legend.width=3, legend.shrink=.9, 
                 legend.cex=1.5, axis.args=list(cex.axis=1.5, tck=-1, hadj=-.1))
      dev.off()
      
      theseCols = theseSparsityColors
      png(file=paste0("Figures/", resultNameRoot, "/finalPercentResidualPlotRelWidthsSparsityCol", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,5), mar=c(5.1, 4.1, 4.1, 6))
      
      for(i in 1:(numberModels-1)) {
        x = valMat[,i]
        y = valMat[,numberModels]
        y = 100 * (x-y)/y
        ylab = ""
        
        plot(x, y, main="", type="n", xlab=bquote(.(modelNames[[numberModels]]) ~ " Rel. Width"), ylab="", 
             ylim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1, log="x")
        if(i == 1)
          mtext(side = 2, bquote("Pct. Diff. from " ~ .(modelNames[[numberModels]])), line = 4, cex=2)
        my_line(x, y, lwd=1, pch=19)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      image.plot(legend.only = TRUE, zlim=range(theseSparsities), nlevel=ncols, legend.mar=4, col=areaColors, add=TRUE, 
                 legend.lab = expression("Sparsity (km"^2 ~" Per Cluster)"), legend.line=5.0, legend.width=3, legend.shrink=.9, 
                 legend.cex=1.5, axis.args=list(at=sparsityTicks, labels=sparsityTickLabels, cex.axis=1.5, tck=-1, hadj=-.1))
      dev.off()
      
      valMat = do.call("cbind", widthsList)
      zlim = range(sweep(valMat[,-ncol(valMat)], 1, valMat[,ncol(valMat)], function(x,y) {100*(x-y)/y}))
      width = 400 * (numberModels - 1)
      height = 400
      theseCols = theseUrbanColors
      png(file=paste0("Figures/", resultNameRoot, "/finalPercentResidualPlotWidths", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,5), mar=c(5.1, 4.1, 4.1, 6))
      
      for(i in 1:(numberModels-1)) {
        x = valMat[,i]
        y = valMat[,numberModels]
        y = 100 * (x-y)/y
        ylab = ""
        
        plot(x, y, main="", type="n", xlab=bquote(.(modelNames[[numberModels]]) ~ " Width"), ylab="", 
             ylim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1, log="x")
        if(i == 1)
          mtext(side = 2, bquote("Pct. Diff. from " ~ .(modelNames[[numberModels]])), line = 4, cex=2)
        my_line(x, y, lwd=1, pch=19)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      image.plot(legend.only = TRUE, zlim=c(0,1), nlevel=ncols, legend.mar=4, col=urbCols, add=TRUE, 
                 legend.lab = "Urbanicity", legend.line=5.0, legend.width=3, legend.shrink=.9, 
                 legend.cex=1.5, axis.args=list(cex.axis=1.5, tck=-1, hadj=-.1))
      dev.off()
      
      theseCols = theseAreaColors
      png(file=paste0("Figures/", resultNameRoot, "/finalPercentResidualPlotWidthsAreaCol", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,5), mar=c(5.1, 4.1, 4.1, 6))
      
      for(i in 1:(numberModels-1)) {
        x = valMat[,i]
        y = valMat[,numberModels]
        y = 100 * (x-y)/y
        ylab = ""
        
        plot(x, y, main="", type="n", xlab=bquote(.(modelNames[[numberModels]]) ~ " Width"), ylab="", 
             ylim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1, log="x")
        if(i == 1)
          mtext(side = 2, bquote("Pct. Diff. from " ~ .(modelNames[[numberModels]])), line = 4, cex=2)
        my_line(x, y, lwd=1, pch=19)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      image.plot(legend.only = TRUE, zlim=range(theseAreas), nlevel=ncols, legend.mar=4, col=areaColors, add=TRUE, 
                 legend.lab = "Area (sq. km)", legend.line=5.0, legend.width=3, legend.shrink=.9, 
                 legend.cex=1.5, axis.args=list(cex.axis=1.5, tck=-1, hadj=-.1))
      dev.off()
      
      theseCols = theseRadiusColors
      png(file=paste0("Figures/", resultNameRoot, "/finalPercentResidualPlotWidthsRadiusCol", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,5), mar=c(5.1, 4.1, 4.1, 6))
      
      for(i in 1:(numberModels-1)) {
        x = valMat[,i]
        y = valMat[,numberModels]
        y = 100 * (x-y)/y
        ylab = ""
        
        plot(x, y, main="", type="n", xlab=bquote(.(modelNames[[numberModels]]) ~ " Width"), ylab="", 
             ylim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1, log="x")
        if(i == 1)
          mtext(side = 2, bquote("Pct. Diff. from " ~ .(modelNames[[numberModels]])), line = 4, cex=2)
        my_line(x, y, lwd=1, pch=19)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      image.plot(legend.only = TRUE, zlim=range(theseRadii), nlevel=ncols, legend.mar=4, col=areaColors, add=TRUE, 
                 legend.lab = "Radius (km)", legend.line=5.0, legend.width=3, legend.shrink=.9, 
                 legend.cex=1.5, axis.args=list(cex.axis=1.5, tck=-1, hadj=-.1))
      dev.off()
      
      theseCols = theseSparsityColors
      png(file=paste0("Figures/", resultNameRoot, "/finalPercentResidualPlotWidthsSparsityCol", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,5), mar=c(5.1, 4.1, 4.1, 6))
      
      for(i in 1:(numberModels-1)) {
        x = valMat[,i]
        y = valMat[,numberModels]
        y = 100 * (x-y)/y
        ylab = ""
        
        plot(x, y, main="", type="n", xlab=bquote(.(modelNames[[numberModels]]) ~ " Width"), ylab="", 
             ylim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1, log="x")
        if(i == 1)
          mtext(side = 2, bquote("Pct. Diff. from " ~ .(modelNames[[numberModels]])), line = 4, cex=2)
        my_line(x, y, lwd=1, pch=19)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      image.plot(legend.only = TRUE, zlim=range(theseSparsities), nlevel=ncols, legend.mar=4, col=areaColors, add=TRUE, 
                 legend.lab = expression("Sparsity (km"^2 ~" Per Cluster)"), legend.line=5.0, legend.width=3, legend.shrink=.9, 
                 legend.cex=1.5, axis.args=list(at=sparsityTicks, labels=sparsityTickLabels, cex.axis=1.5, tck=-1, hadj=-.1))
      dev.off()
    } else {
      urbanPointTypes = sapply(theseResults$urban, function(x) {ifelse(x, 15, 17)})
      thesePointTypes = urbanPointTypes
      theseUrbanColors = theseCols
      
      # determine colorings by distance to observation
      if(thisArea == "Pixel") {
        theseDistances = getPredictionDistance(doLog=FALSE)
        
        colI = cut(theseDistances, breaks=seq(min(theseDistances) - .0001, max(theseDistances) + .0001, l=length(areaColors)+1), labels=FALSE)
        theseDistanceColors = areaColors[colI]
        theseUrbanColors = theseCols
        distanceTicks = getPredictionDistanceTicks()
        distanceTickLabels = as.character(distanceTicks)
        distanceTicks = distanceTicks
        distancePointTypes = 19
        cex = 0.1
      } else if(thisArea == "Cluster") {
        cex = 0.8
      }
      
      my_line <- function(x,y,...){
        if(diff(range(x)) >= .04)
          xlim = zlim
        
        if(diff(range(y)) >= .04)
          ylim = zlim
        
        # abline(a = 0,b = 1, col=rgb(.4, .4, .4),...)
        points(x,y,..., col=theseCols, pch=thesePointTypes)
        abline(h=0, lty=2, ...)
      }
      
      width = 4 * (numberModels - 1)
      height = 4
      pdf(file=paste0("Figures/", resultNameRoot, "/finalPercentResidualPlot", plotNameRoot, extraPlotNameRoot, thisArea, ".pdf"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,0), mar=c(5.1, 4.1, 4.1, 4.1))
      
      zlim = range(sweep(valMat[,-ncol(valMat)], 1, valMat[,ncol(valMat)], function(x,y) {100*(x-y)/y}))
      for(i in 1:(numberModels-1)) {
        x = valMat[,i]
        y = valMat[,numberModels]
        y = 100 * (x-y)/y
        ylab = ""
        
        plot(x, y, main="", type="n", xlab=bquote(.(modelNames[[numberModels]]) ~ " Estimates"), ylab="", 
             ylim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1, log="x")
        if(i == 1)
          mtext(side = 2, bquote("Pct. Diff. from " ~ .(modelNames[[numberModels]])), line = 4, cex=2)
        my_line(x, y, lwd=1, cex=cex)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      dev.off()
      
      if(thisArea == "Pixel") {
        thesePointTypes = distancePointTypes
        theseCols = theseDistanceColors
        pdf(file=paste0("Figures/", resultNameRoot, "/finalPercentResidualPlotDistCol", plotNameRoot, extraPlotNameRoot, thisArea, ".pdf"), width=width, height=height)
        par(mfrow=c(1, 3), oma=c(0,4,0,5), mar=c(5.1, 4.1, 4.1, 6))
        
        zlim = range(sweep(valMat[,-ncol(valMat)], 1, valMat[,ncol(valMat)], function(x,y) {100*(x-y)/y}))
        for(i in 1:(numberModels-1)) {
          x = valMat[,i]
          y = valMat[,numberModels]
          y = 100 * (x-y)/y
          ylab = ""
          
          plot(x, y, main="", type="n", xlab=bquote(.(modelNames[[numberModels]]) ~ " Estimates"), ylab="", 
               ylim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1, log="x")
          if(i == 1)
            mtext(side = 2, bquote("Pct. Diff. from " ~ .(modelNames[[numberModels]])), line = 4, cex=2)
          my_line(x, y, lwd=1, cex=cex)
          mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
        }
        image.plot(legend.only = TRUE, zlim=range(theseDistances), nlevel=ncols, legend.mar=4, col=areaColors, add=TRUE, 
                   legend.lab = "Distance to Observation (km)", legend.line=5.0, legend.width=3, legend.shrink=.9, 
                   legend.cex=1.5, axis.args=list(cex.axis=1.5, tck=-1, hadj=-.1))
        dev.off()
      }
      
      thesePointTypes = urbanPointTypes
      theseCols = theseUrbanColors
      valMat = do.call("cbind", lapply(1:length(widthsList), function(x) {widthsList[[x]] / predsList[[x]]}))
      pdf(file=paste0("Figures/", resultNameRoot, "/finalPercentResidualPlotRelWidths", plotNameRoot, extraPlotNameRoot, thisArea, ".pdf"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,0), mar=c(5.1, 4.1, 4.1, 4.1))
      
      zlim = range(sweep(valMat[,-ncol(valMat)], 1, valMat[,ncol(valMat)], function(x,y) {100*(x-y)/y}))
      for(i in 1:(numberModels-1)) {
        x = valMat[,i]
        y = valMat[,numberModels]
        y = 100 * (x-y)/y
        ylab = ""
        
        plot(x, y, main="", type="n", xlab=bquote(.(modelNames[[numberModels]]) ~ " Rel. Width"), ylab="", 
             ylim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1, log="x")
        if(i == 1)
          mtext(side = 2, bquote("Pct. Diff. from " ~ .(modelNames[[numberModels]])), line = 4, cex=2)
        my_line(x, y, lwd=1, cex=cex)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      dev.off()
      
      if(thisArea == "Pixel") {
        thesePointTypes = distancePointTypes
        theseCols = theseDistanceColors
        valMat = do.call("cbind", lapply(1:length(widthsList), function(x) {widthsList[[x]] / predsList[[x]]}))
        pdf(file=paste0("Figures/", resultNameRoot, "/finalPercentResidualPlotRelWidthsDistCol", plotNameRoot, extraPlotNameRoot, thisArea, ".pdf"), width=width, height=height)
        par(mfrow=c(1, 3), oma=c(0,4,0,5), mar=c(5.1, 4.1, 4.1, 6))
        
        zlim = range(sweep(valMat[,-ncol(valMat)], 1, valMat[,ncol(valMat)], function(x,y) {100*(x-y)/y}))
        for(i in 1:(numberModels-1)) {
          x = valMat[,i]
          y = valMat[,numberModels]
          y = 100 * (x-y)/y
          ylab = ""
          
          plot(x, y, main="", type="n", xlab=bquote(.(modelNames[[numberModels]]) ~ " Rel. Width"), ylab="", 
               ylim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1, log="x")
          if(i == 1)
            mtext(side = 2, bquote("Pct. Diff. from " ~ .(modelNames[[numberModels]])), line = 4, cex=2)
          my_line(x, y, lwd=1, cex=cex)
          mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
        }
        image.plot(legend.only = TRUE, zlim=range(theseDistances), nlevel=ncols, legend.mar=4, col=areaColors, add=TRUE, 
                   legend.lab = "Distance to Observation (km)", legend.line=5.0, legend.width=3, legend.shrink=.9, 
                   legend.cex=1.5, axis.args=list(cex.axis=1.5, tck=-1, hadj=-.1))
        dev.off()
      }
      
      thesePointTypes = urbanPointTypes
      theseCols = theseUrbanColors
      valMat = do.call("cbind", widthsList)
      zlim = range(sweep(valMat[,-ncol(valMat)], 1, valMat[,ncol(valMat)], function(x,y) {100*(x-y)/y}))
      pdf(file=paste0("Figures/", resultNameRoot, "/finalPercentResidualPlotWidths", plotNameRoot, extraPlotNameRoot, thisArea, ".pdf"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,0), mar=c(5.1, 4.1, 4.1, 4.1))
      
      for(i in 1:(numberModels-1)) {
        x = valMat[,i]
        y = valMat[,numberModels]
        y = 100 * (x-y)/y
        ylab = ""
        
        plot(x, y, main="", type="n", xlab=bquote(.(modelNames[[numberModels]]) ~ " Width"), ylab="", 
             ylim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1, log="x")
        if(i == 1)
          mtext(side = 2, bquote("Pct. Diff. from " ~ .(modelNames[[numberModels]])), line = 4, cex=2)
        my_line(x, y, lwd=1, cex=cex)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      dev.off()
      
      if(thisArea == "Pixel") {
        thesePointTypes = distancePointTypes
        theseCols = theseDistanceColors
        pdf(file=paste0("Figures/", resultNameRoot, "/finalPercentResidualPlotWidthsDistCol", plotNameRoot, extraPlotNameRoot, thisArea, ".pdf"), width=width, height=height)
        par(mfrow=c(1, 3), oma=c(0,4,0,5), mar=c(5.1, 4.1, 4.1, 6))
        
        for(i in 1:(numberModels-1)) {
          x = valMat[,i]
          y = valMat[,numberModels]
          y = 100 * (x-y)/y
          ylab = ""
          
          plot(x, y, main="", type="n", xlab=bquote(.(modelNames[[numberModels]]) ~ " Width"), ylab="", 
               ylim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1, log="x")
          if(i == 1)
            mtext(side = 2, bquote("Pct. Diff. from " ~ .(modelNames[[numberModels]])), line = 4, cex=2)
          my_line(x, y, lwd=1, cex=cex)
          mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
        }
        image.plot(legend.only = TRUE, zlim=range(theseDistances), nlevel=ncols, legend.mar=4, col=areaColors, add=TRUE, 
                   legend.lab = "Distance to Observation (km)", legend.line=5.0, legend.width=3, legend.shrink=.9, 
                   legend.cex=1.5, axis.args=list(cex.axis=1.5, tck=-1, hadj=-.1))
        dev.off()
      }
    }
  }
}

plotCovariograms = function(dat, resultFilenames, modelClasses, modelVariations, 
                            varName="education", plotNameRoot="Education", resultNameRoot="Ed", 
                            cgramList=NULL, loadResults=FALSE, saveResults=!loadResults, col=NULL, 
                            doModelClassPlots=FALSE, lty=1, pch=19, hNuggetShift=-3) {
  plotNameRootLower = tolower(plotNameRoot)
  resultNameRootLower = tolower(resultNameRoot)
  numberModels = length(resultFilenames)
  uniqueModelClasses = unique(modelClasses)
  
  if(is.null(col))
    col = rainbow(length(modelClasses))
  
  binomialFamily = grepl("LgtN", resultFilenames[1])
  
  # first get the covariograms if necessary
  if(is.null(cgramList)) {
    cgramList = list()
    for(j in 1:numberModels) {
      modelName = bquote(.(modelClasses[j])[.(modelVariations[j])])
      
      # load this model and get the covariogram if necessary
      print(paste0("Loading ", paste(modelClasses[j], modelVariations[j])))
      out = load(resultFilenames[j])
      if(!loadResults || !("cgram" %in% names(results))) {
        hyperDraws = results$fit$hyperMat
        browser()
        
        # determine if this model has a cluster effect
        thisModelClass = modelClasses[j]
        clusterEffect = grepl("includeClusterTRUE", resultFilenames[j])
        
        # hyperparameters will be drawn differently depending on the type of model
        if(thisModelClass == "SPDE") {
          # get hyperparameter draws
          effectiveRangeI = grepl("Range for field", rownames(hyperDraws))
          sdI = grepl("Stdev for field", rownames(hyperDraws))
          effectiveRangeVals = hyperDraws[effectiveRangeI,]
          varVals = hyperDraws[sdI,]^2
          if(clusterEffect)
            nuggetVarVals = 1/hyperDraws[3,]
          else
            nuggetVarVals = rep(0, ncol(hyperDraws))
          
          # get range of the data and the SPDE basis function mesh for which to compute the covariograms
          out = load(paste0("dataPointsKenya.RData"))
          xRangeDat = dataPointsKenya$xRange
          yRangeDat = dataPointsKenya$yRange
          mesh = results$fit$mesh
          
          # compute the covariance function for the different hyperparameter samples
          cgram = covarianceDistributionSPDE(effectiveRangeVals, varVals, nuggetVarVals, mesh, xRangeDat=xRangeDat, yRangeDat=yRangeDat)
        } else if(thisModelClass == "ELK") {
          # get lattice information object, determine whether it's a separateRange model
          latInfo = results$fit$latInfo
          nLayer = length(latInfo)
          separateRanges = grepl("separateRangesTRUE", resultFilenames[j])
          if(separateRanges)
            alphaI = (1 + nLayer+1 + 1):(1 + nLayer+1 + nLayer-1)
          else
            alphaI = 4:(3+nLayer-1)
          zSamples = matrix(hyperDraws[alphaI,], ncol=length(alphaI))
          xSamples = t(matrix(apply(zSamples, 1, multivariateExpit), ncol=length(alphaI)))
          xSamples = rbind(xSamples, 1-colSums(xSamples))
          
          # get hyperparameter draws
          nuggetVarVals = rep(0, ncol(hyperDraws))
          if(separateRanges) {
            kappaVals = sweep(sqrt(8)/exp(hyperDraws[2:(nLayer+1),]), 1, sapply(latInfo, function(x) {x$latWidth}), "*")
            rhoVals = exp(hyperDraws[nLayer+2,])
          } else {
            kappaVals = sqrt(8)/exp(hyperDraws[2,]) * latInfo[[1]]$latWidth
            rhoVals = exp(hyperDraws[3,])
          }
          alphaMat = xSamples
          
          # compute the covariance function for many different hyperparameter samples
          out = covarianceDistributionLKINLA(latInfo, kappaVals, rhoVals, nuggetVarVals, alphaMat)
        } else {
          stop(paste0("Unrecognized model class: ", thisModelClass))
        }
        
        # save results if necessary
        if(saveResults) {
          results$cgram = cgram
          save(results, file=resultFilenames[j])
        }
      } else {
        cgram = results$cgram
      }
      
      # append to our list of covariograms
      cgramList = c(cgramList, list(cgram))
    }
  }
  
  ##### if there are multiple unique model classes, call this function on each one individually 
  ##### as well as on the combination
  if(length(uniqueModelClasses) != 1 && doModelClassPlots) {
    for(i in 1:length(uniqueModelClasses)) {
      thisModelClass = uniqueModelClasses[i]
      thisI = modelClasses == thisModelClass
      if(sum(thisI) == 1)
        next
      
      thisResultFilenames = resultFilenames[thisI]
      thisModelClasses = modelClasses[thisI]
      thisModelVariations = modelVariations[thisI]
      thiscgramList = cgramList[thisI]
      thisColors = col[thisI]
      thispch = pch[thisI]
      thislty = lty[thisI]
      plotCovariograms(dat, thisResultFilenames, thisModelClasses, thisModelVariations, 
                       varName, plotNameRoot, resultNameRoot, thiscgramList, 
                       loadResults=TRUE, saveResults=FALSE, col=thisColors, lty=thislty, 
                       pch=thispch, hNuggetShift=hNuggetShift)
    }
  }
  
  # reorder the models if necessary so that we have 2 blocks of 4 (note that these permutations are 
  # the inverse permutations of themselves)
  # if(numberModels == 8) {
  #   reordering = c(t(rbind(c(1, 2, 5, 6), 
  #                          c(3, 4, 7, 8))))
  #   resultFilenames = resultFilenames[reordering]
  #   modelClasses = modelClasses[reordering]
  #   modelVariations = modelVariations[reordering]
  #   cgramList = cgramList[reordering]
  # } else if(numberModels == 4) {
  #   reordering = c(t(rbind(c(1, 3), 
  #                          c(2, 4))))
  #   resultFilenames = resultFilenames[reordering]
  #   modelClasses = modelClasses[reordering]
  #   modelVariations = modelVariations[reordering]
  #   cgramList = cgramList[reordering]
  # }
  # if(numberModels == 6) {
  #   reordering = c(t(rbind(c(1, 3, 4), 
  #                          c(2, 5, 6)))) # urban effects on bottom, first column is SPDE
  #   resultFilenames = resultFilenames[reordering]
  #   modelClasses = modelClasses[reordering]
  #   modelVariations = modelVariations[reordering]
  #   cgramList = cgramList[reordering]
  # } else {
  #   reordering = 1:numberModels # no need to reorder
  #   resultFilenames = resultFilenames[reordering]
  #   modelClasses = modelClasses[reordering]
  #   modelVariations = modelVariations[reordering]
  #   cgramList = cgramList[reordering]
  # }
  
  # get range of covariogram values
  yRange = c()
  yRangeNoCIs = c()
  for(j in 1:numberModels) {
    modelName = bquote(.(modelClasses[j])[.(modelVariations[j])])
    
    thiscgram = cgramList[[j]]
    yRange = range(c(yRange, thiscgram$cov, thiscgram$lowerCov, thiscgram$upperCov))
    yRangeNoCIs = range(c(yRangeNoCIs, thiscgram$cov))
  }
  
  # modify plot file names if necessary
  if(length(uniqueModelClasses) == 1)
    extraPlotNameRoot = uniqueModelClasses[1]
  else
    extraPlotNameRoot = ""
  
  allModelNames = list()
  for(j in 1:numberModels)
    allModelNames = c(allModelNames, bquote(.(modelClasses[j])[.(modelVariations[j])]))
  
  ##### Plot everything
  print("Plotting covariograms...")
  # col = rainbow(numberModels)
  # width = 4 * ceiling(numberModels/2)
  
  pdf(file=paste0("Figures/", resultNameRoot, "/covariogramsAll", plotNameRoot, extraPlotNameRoot, ".pdf"), width=5, height=5)
  
  # plot the covariograms together
  # allModelNames = paste(modelClasses, modelVariations)
  for(j in 1:numberModels) {
    modelName = bquote(.(modelClasses[j])[.(modelVariations[j])])
    thiscgram = cgramList[[j]]
    d = thiscgram$d
    sortI = sort(d, index.return=TRUE)$ix
    d = d[sortI]
    
    if(modelClasses[j] == "SPDE") {
      covMean = thiscgram$cov[sortI]
      upperCov=thiscgram$upperCov[1,sortI] # second row is the 95% CI, while the first is the 80% CI
      lowerCov=thiscgram$lowerCov[1,sortI]
    } else {
      covMean = thiscgram$cov[sortI]
      upperCov=thiscgram$upperCov[sortI]
      lowerCov=thiscgram$lowerCov[sortI]
    }
    
    if(binomialFamily) {
      d0 = d == 0
      if(j == 1) {
        plot(hNuggetShift, mean(covMean[d0]), pch=pch[j], cex=.4, main="", xlab="Distance (km)", ylab="Covariance", 
             ylim=yRange, xlim=range(d), col=col[j])
        lines(d[!d0], covMean[!d0], col=col[j], lwd=2)
      } else {
        points(hNuggetShift, mean(covMean[d0]), pch=pch[j], cex=.4, col=col[j])
        lines(d[!d0], covMean[!d0], col=col[j], lwd=2)
      }
      
      points(hNuggetShift, mean(lowerCov[d0]), pch=pch[j], cex=.2, col=col[j], lty=2)
      lines(d[!d0], lowerCov[!d0], col=col[j], lty=2, lwd=2)
      points(hNuggetShift, mean(upperCov[d0]), pch=pch[j], cex=.2, col=col[j], lty=2)
      lines(d[!d0], upperCov[!d0], col=col[j], lty=2, lwd=2)
    } else {
      if(j == 1) {
        plot(d, covMean, type="l", main="", xlab="Distance (km)", ylab="Covariance", 
             ylim=yRange, col=col[j], lwd=2)
      } else {
        lines(d, covMean, col=col[j], lwd=2)
      }
      
      lines(d, lowerCov, lty=2, col=col[j], lwd=2)
      lines(d, upperCov, lty=2, col=col[j], lwd=2)
    }
    
    # lines(d, mixtureCovFun(d), col="green")
    # legend("topright", c("Truth", "Estimate", "80% CI"), lty=c(1, 1, 2), col=c("green", "black", "black"))
    
  }
  legend("topright", as.expression(allModelNames), lty=1, col=col, cex=ifelse(numberModels >= 5, .7, 1))
  dev.off()
  
  pdf(file=paste0("Figures/", resultNameRoot, "/covariogramsAllNoCIs", plotNameRoot, extraPlotNameRoot, ".pdf"), width=5, height=5)
  
  # plot the covariograms together
  for(j in 1:numberModels) {
    modelName = bquote(.(modelClasses[j])[.(modelVariations[j])])
    thiscgram = cgramList[[j]]
    d = thiscgram$d
    sortI = sort(d, index.return=TRUE)$ix
    d = d[sortI]
    
    if(modelClasses[j] == "SPDE") {
      covMean = thiscgram$cov[sortI]
      # upperCov=thiscgram$upperCov[1,sortI] # second row is the 95% CI, while the first is the 80% CI
      # lowerCov=thiscgram$lowerCov[1,sortI]
    } else {
      covMean = thiscgram$cov[sortI]
      # upperCov=thiscgram$upperCov[sortI]
      # lowerCov=thiscgram$lowerCov[sortI]
    }
    
    if(binomialFamily) {
      d0 = d == 0
      if(j == 1) {
        plot(hNuggetShift, mean(covMean[d0]), pch=pch[j], cex=.8, main="", xlab="Distance (km)", ylab="Covariance", 
             ylim=yRangeNoCIs, xlim=range(d), col=col[j], lty=lty[j])
        lines(d[!d0], covMean[!d0], col=col[j], lty=lty[j], lwd=2)
      } else {
        points(hNuggetShift, mean(covMean[d0]), pch=pch[j], cex=.8, col=col[j], lty=lty[j])
        lines(d[!d0], covMean[!d0], col=col[j], lty=lty[j], lwd=2)
        }
    } else {
      if(j == 1) {
        plot(d, covMean, type="l", main="", xlab="Distance (km)", ylab="Covariance", 
             ylim=yRangeNoCIs, col=col[j], lty=lty[j], lwd=2)
      } else {
        lines(d, covMean, col=col[j], lty=lty[j], lwd=2)
      }
    }
    
    # lines(d, lowerCov, lty=2, col=col[j])
    # lines(d, upperCov, lty=2, col=col[j])
    # lines(d, mixtureCovFun(d), col="green")
    # legend("topright", c("Truth", "Estimate", "80% CI"), lty=c(1, 1, 2), col=c("green", "black", "black"))
    
  }
  legend("topright", as.expression(allModelNames), lty=lty, col=col, cex=ifelse(numberModels >= 5, .5, 1))
  dev.off()
  
  print("Plotting correlograms...")
  
  pdf(file=paste0("Figures/", resultNameRoot, "/correlogramsAll", plotNameRoot, extraPlotNameRoot, ".pdf"), width=5, height=5)
  
  # plot the correlograms together
  for(j in 1:numberModels) {
    modelName = bquote(.(modelClasses[j])[.(modelVariations[j])])
    thiscgram = cgramList[[j]]
    d = thiscgram$d
    sortI = sort(d, index.return=TRUE)$ix
    d = d[sortI]
    
    if(modelClasses[j] == "SPDE") {
      corMean = thiscgram$cor[sortI]
      upperCor=thiscgram$upperCor[1,sortI] # second row is the 95% CI, while the first is the 80% CI
      lowerCor=thiscgram$lowerCor[1,sortI]
    } else {
      corMean = thiscgram$cor[sortI]
      upperCor=thiscgram$upperCor[sortI]
      lowerCor=thiscgram$lowerCor[sortI]
    }
    
    if(binomialFamily) {
      d0 = d == 0
      if(j == 1) {
        plot(hNuggetShift, mean(corMean[d0]), pch=pch[j], cex=.4, main="", xlab="Distance (km)", ylab="Correlation", 
             ylim=c(0,1), xlim=range(d), col=col[j])
        lines(d[!d0], corMean[!d0], col=col[j], lwd=2)
      } else {
        points(hNuggetShift, mean(corMean[d0]), pch=pch[j], cex=.4, col=col[j])
        lines(d[!d0], corMean[!d0], col=col[j], lwd=2)
      }
      
      points(hNuggetShift, mean(lowerCor[d0]), pch=pch[j], cex=.2, col=col[j], lty=2)
      lines(d[!d0], lowerCor[!d0], col=col[j], lty=2, lwd=2)
      points(hNuggetShift, mean(upperCor[d0]), pch=pch[j], cex=.2, col=col[j], lty=2)
      lines(d[!d0], upperCor[!d0], col=col[j], lty=2, lwd=2)
    } else {
      if(j == 1) {
        plot(d, corMean, type="l", main="", xlab="Distance (km)", ylab="Correlation", 
             ylim=c(0,1), col=col[j], lwd=2)
      } else {
        lines(d, corMean, col=col[j], lwd=2)
      }
      
      lines(d, lowerCor, lty=2, col=col[j], lwd=2)
      lines(d, upperCor, lty=2, col=col[j], lwd=2)
    }
    # lines(d, mixtureCovFun(d), col="green")
    # legend("topright", c("Truth", "Estimate", "80% CI"), lty=c(1, 1, 2), col=c("green", "black", "black"))
    
  }
  legend("topright", as.expression(allModelNames), lty=1, col=col, cex=ifelse(numberModels >= 5, .5, 1))
  dev.off()
  
  pdf(file=paste0("Figures/", resultNameRoot, "/correlogramsAllNoCIs", plotNameRoot, extraPlotNameRoot, ".pdf"), width=5, height=5)
  
  # plot the correlograms together
  for(j in 1:numberModels) {
    modelName = bquote(.(modelClasses[j])[.(modelVariations[j])])
    thiscgram = cgramList[[j]]
    d = thiscgram$d
    sortI = sort(d, index.return=TRUE)$ix
    d = d[sortI]
    
    if(modelClasses[j] == "SPDE") {
      corMean = thiscgram$cor[sortI]
      # upperCor=thiscgram$upperCor[1,sortI] # second row is the 95% CI, while the first is the 80% CI
      # lowerCor=thiscgram$lowerCor[1,sortI]
    } else {
      corMean = thiscgram$cor[sortI]
      # upperCor=thiscgram$upperCor[sortI]
      # lowerCor=thiscgram$lowerCor[sortI]
    }
    
    if(binomialFamily) {
      d0 = d == 0
      if(j == 1) {
        plot(hNuggetShift, mean(corMean[d0]), pch=pch[j], cex=.8, main="", xlab="Distance (km)", ylab="Correlation", 
             ylim=c(0,1), xlim=range(d), col=col[j], lty=lty[j])
        lines(d[!d0], corMean[!d0], col=col[j], lty=lty[j], lwd=2)
      } else {
        points(hNuggetShift, mean(corMean[d0]), pch=pch[j], cex=.8, col=col[j], lty=lty[j])
        lines(d[!d0], corMean[!d0], col=col[j], lty=lty[j], lwd=2)
      }
    } else {
      if(j == 1) {
        plot(d, corMean, type="l", main="", xlab="Distance (km)", ylab="Correlation", 
             ylim=c(0,1), col=col[j], lty=lty[j], lwd=2)
      } else {
        lines(d, corMean, col=col[j], lty=lty[j], lwd=2)
      }
    }
    
    # lines(d, lowerCor, lty=2, col=col[j])
    # lines(d, upperCor, lty=2, col=col[j])
    # lines(d, mixtureCovFun(d), col="green")
    # legend("topright", c("Truth", "Estimate", "80% CI"), lty=c(1, 1, 2), col=c("green", "black", "black"))
    
  }
  legend("topright", as.expression(allModelNames), lty=lty, col=col, cex=ifelse(numberModels >= 5, .5, 1))
  dev.off()
  
}

printModelPredictionTables = function(dataType=c("ed", "mort"), resultFilenames=NULL, modelClasses=c("SPDE", "ELK"), modelVariations=c("U", "U"), 
                                      nDigitsPredictions=3, nDigitsParameters=3, byRow=FALSE, areaLevels=c("Region", "County")) {
  
  dataType = match.arg(dataType)
  if(dataType == "mort") {
    out = load("../U5MR/kenyaData.RData")
    dat = mort
    resultNameRoot = "Mort"
  }
  else {
    out = load("../U5MR/kenyaDataEd.RData")
    dat = ed
    resultNameRoot = "Ed"
  }
  resultNameRootLower = tolower(resultNameRoot)
  
  if(is.null(resultFilenames)) {
    resultFilenames = c(
      paste0("savedOutput/", resultNameRoot, "/resultsSPDE", resultNameRootLower, "_urbanEffectTRUE_LgtN.RData"), 
      paste0("savedOutput/", resultNameRoot, "/resultsLKINLA", resultNameRootLower, "_separateRangesTRUE_urbanEffectTRUE_LgtN_noUrbanPrior.RData")
    )
  }
  
  ##### load estimates
  # construct a list of all estimates. First list index is model class, second index is model variation, 
  # third is aggregation level and parameters
  print("getting SPDE estimates...")
  out = load(resultFilenames[1])
  spdeResultsRegion = results$aggregatedResults$predictions$regionPredictions
  spdeResultsCounty = results$aggregatedResults$predictions$countyPredictions
  spdeParameterTable = results$aggregatedResults$parameterSummary
  
  ##### ELK estimates
  print("getting ELK estimates...")
  out = load(resultFilenames[2])
  lkinlaResultsRegion = results$aggregatedResults$predictions$regionPredictions
  lkinlaResultsCounty = results$aggregatedResults$predictions$countyPredictions
  lkinlaParameterTable = results$aggregatedResults$parameterSummary
  
  ##### now we construct the latex table for each areal level
  for(i in 1:length(areaLevels)) {
    thisArea = areaLevels[i]
    
    # get the model results for this level of areal aggregation
    if(thisArea == "Region") {
      lkinlaResults = lkinlaResultsRegion
      spdeResults = spdeResultsRegion
    } else if(thisArea == "County") {
      lkinlaResults = lkinlaResultsCounty
      spdeResults = spdeResultsCounty
    } else
      stop(paste0("Unrecognized area name: ", thisArea))
    
    # construct result table
    tab = cbind(Estimates=spdeResults$preds, Q10=spdeResults$Q10, Q90=spdeResults$Q90)
    tab = cbind(tab, Estimates=lkinlaResults$preds, Q10=lkinlaResults$Q10, Q90=lkinlaResults$Q90)
    tab = format(tab, digits=nDigitsPredictions) # round to the nearest 1 child per 1000
    tab = cbind(as.character(spdeResults$areaName), tab)
    colnames(tab) = c(thisArea, rep(c("Est", "Q10", "Q90"), 2))
    rownames(tab) = NULL
    
    # print out relevant summary statistics about the predictions
    print(paste0("SPDE U range of ", thisArea, " predictions: ", diff(range(as.numeric(tab[,2])))))
    print(paste0("SPDE U median ", thisArea, " 80% CI width: ", median(as.numeric(tab[,4])-as.numeric(tab[,3]))))
    print(paste0("ELK U range of ", thisArea, " predictions: ", diff(range(as.numeric(tab[,5])))))
    print(paste0("ELK U median 80% ", thisArea, " CI width: ", median(as.numeric(tab[,7])-as.numeric(tab[,6]))))
    
    # print out the reformatted table
    require(kableExtra)
    print(paste0("Predictions at ", tolower(thisArea), " level:"))
    fullTab = tab %>%
      kable("latex", escape = F, booktabs = TRUE, format.args=list(drop0trailing=FALSE, scientific=FALSE), 
            align=c("l", rep("r", ncol(tab) - 1)), longtable=thisArea == "County", caption = thisArea) %>% 
      kable_styling(latex_options =c("repeat_header", "scale_down"))
    numberColumns = c(" "=1, "SPDE U"=3, "ELK U"=3)
    print(add_header_above(fullTab, numberColumns, italic=TRUE, bold=TRUE, escape=FALSE, line=TRUE))
  }
  
  
  ##### now construct the parameter tables
  ##### now print the parameter estimates:
  print("printing parameter estimates...")
  
  # SPDE
  
  parameters = round(spdeParameterTable, digits=nDigitsParameters)
  
  # modify the row names do not include the word "Summary"
  allNames = rownames(parameters)
  # rownames(parameters) = unlist(sapply(allNames, strsplit, split="Summary"))
  browser()
  byRow = FALSE
  if(byRow == FALSE)
    nonRangePar = format(parameters[-nrow(parameters),], digits=nDigitsParameters, scientific=FALSE)
  else
    nonRangePar = t(apply(parameters[-nrow(parameters),], 1, format, digits=2, scientific=FALSE))
  formattedParameters = rbind(nonRangePar, format(parameters[nrow(parameters),], digits=0))
  rownames(formattedParameters) = c("Intercept", "Urban", "Total Var", "Spatial Var", "Cluster Var", "Total SD", "Spatial SD", "Cluster SD", "Range")
  
  # rename the columns
  colnames(formattedParameters) = c("Est", "SD", "Q10", "Q50", "Q90")
  
  # print(paste0("Parameter summary table for SPDE ", typeText, " model:"))
  # print(xtable(parameters, digit=3))
  parTable = formattedParameters
  parTable = cbind(rownames(parTable), parTable)
  rownames(parTable) = NULL
  
  # ELK
  
  parameters = round(lkinlaParameterTable, digits=nDigitsParameters)
  
  # modify the row names do not include the word "Summary"
  allNames = rownames(parameters)
  # rownames(parameters) = unlist(sapply(allNames, strsplit, split="Summary"))
  
  if(byRow == FALSE)
    nonRangePar = format(parameters[-(9:10),], digits=nDigitsParameters, scientific=FALSE)
  else
    nonRangePar = t(apply(parameters[-(9:10),], 1, format, digits=2, scientific=FALSE))
  formattedParameters = rbind(nonRangePar[1:8,], format(parameters[9:10,], digits=0), nonRangePar[9:nrow(nonRangePar),])
  rownames(formattedParameters) = c("Intercept", "Urban", "Total Var", "Spatial Var", "Cluster Var", "Total SD", "Spatial SD", "Cluster SD", "Range1", "Range2", 
                                    "alpha1", "alpha2")
  
  # rename the columns
  colnames(formattedParameters) = c("Est", "SD", "Q10", "Q50", "Q90")
  
  # print(paste0("Parameter summary table for SPDE ", typeText, " model:"))
  # print(xtable(parameters, digit=3))
  
  formattedParameters = cbind(rownames(formattedParameters), formattedParameters)
  rownames(formattedParameters) = NULL
  colnames(parTable)[1] = "Parameter"
  colnames(formattedParameters)[1] = "Parameter"
  parTable = rbind(parTable, formattedParameters)
  
  fullTab = kable(parTable, "latex", booktabs = T, escape=FALSE, format.args=list(drop0trailing=FALSE, scientific=FALSE), 
                  align=c("l", rep("r", ncol(parTable) - 1)), longtable=FALSE, caption = "Parameter") %>% 
    kable_styling(latex_options="repeat_header")
  
  fullTab = fullTab %>% 
    pack_rows("SPDE U", 1, 9, escape=FALSE, bold=TRUE, italic=TRUE) %>% 
    pack_rows("ELK U", 10, nrow(parTable), escape=FALSE, bold=TRUE, italic=TRUE)
  print(fullTab)
}

# this is mostly a test function for plotting the predictions of a single model
# this function plots central estimates and credible interval widths over a map of Kenya 
# 4 plots are made by default: one for each area level
plotSingleModelPredictions = function(dat=NULL, results, modelName="", targetPop=c("women", "children"), 
                                      areaLevels=c("Region", "County", "Pixel", "Cluster"), 
                                      plotNameRoot="", inExampleFolder=FALSE, 
                                      meanRange=NULL, meanTicks=NULL, meanTickLabels=NULL, widthRange=NULL, widthTicks=NULL, widthTickLabels=NULL, 
                                      meanCols=makeRedBlueDivergingColors(64), 
                                      widthCols=makeBlueYellowSequentialColors(64), 
                                      kenyaLatRange=c(-4.6, 5), kenyaLonRange=c(33.5, 42.0)) {
  targetPop = match.arg(targetPop)
  if(targetPop == "women") {
    varName="education"
    plotNameRoot=paste0("Education", plotNameRoot)
    resultNameRoot="Ed"
    if(is.null(dat)) {
      out = load("../U5MR/kenyaDataEd.RData")
      dat = mort
    }
  } else if(targetPop == "children") {
    varName="mort"
    plotNameRoot=paste0("Mort", plotNameRoot)
    resultNameRoot="Mort"
    if(is.null(dat)) {
      out = load("../U5MR/kenyaData.RData")
      dat = mort
    }
  }
  if(!inExampleFolder)
    resultNameRoot = ""
  plotNameRootLower = tolower(plotNameRoot)
  resultNameRootLower = tolower(resultNameRoot)
  
  # load the fine grid used to approximate continuous prediction (the actual population density values don't matter)
  print("Loading prediction grid and shapefiles...")
  load("../U5MR/popGrid.RData")
  
  # load shape files for plotting
  require(maptools)
  # regionMap = readShapePoly("../U5MR/mapData/kenya_region_shapefile/kenya_region_shapefile.shp", delete_null_obj=TRUE, force_ring=TRUE, repair=TRUE)
  out = load("regionMap.RData")
  out = load("../U5MR/adminMapData.RData")
  kenyaMap = adm0
  countyMap = adm1
  
  ##### central estimates and credible interval widths
  
  for(i in 1:length(areaLevels)) {
    thisArea = areaLevels[i]
    
    # get map for this areal aggregation level
    if(thisArea == "Cluster") {
      next
    } else if(thisArea == "Pixel") {
      thisMap = countyMap
    } else if(thisArea == "Region") {
      thisMap = regionMap
    } else if(thisArea == "County") {
      thisMap = countyMap
    }
    
    print("Plotting central estimates and credible interval widths...")
    numberModels=1
    width = 400 * numberModels
    png(file=paste0("Figures/", resultNameRoot, "/", plotNameRoot, thisArea, ".png"), width=width, height=1000)
    
    if(thisArea %in% c("Region", "County"))
      par(mfrow=c(2,numberModels))
    else
      par(mfrow=c(2,numberModels), oma=c( 0,0,0,2), mar=c(5.1, 4.1, 4.1, 6))
    
    # first load in the predictions
    predictionList = list()
    widthList = list()
    
    # load this model and plot results in the given column
    theseResults = results$aggregatedResults$predictions[[paste0(tolower(thisArea), "Predictions")]]
    predictionList = c(predictionList, list(theseResults$preds))
    widthList = c(widthList, list(theseResults$Q90 - theseResults$Q10))
    
    # plot the predictions
    
    # load this model and plot results in the given column
    # out = load(resultFilenames[j])
    # theseResults = results$aggregatedResults$predictions[[paste0(tolower(thisArea), "Predictions")]]
    
    if(is.null(meanRange))
      thisMeanRange = logit(range(predictionList[[1]]))
    else
      thisMeanRange = logit(meanRange)
    if(is.null(meanTicks))
      thisMeanTicks = NULL
    else
      thisMeanTicks = logit(meanTicks)
    
    if(thisArea %in% c("Region", "County")) {
      
      plotMapDat(plotVar=predictionList[[1]], new = TRUE, 
                 main=bquote(.(modelName) ~ " estimates"), scaleFun=logit, scaleFunInverse=expit, 
                 cols=meanCols, zlim=thisMeanRange, ticks=meanTicks, tickLabels=meanTickLabels, 
                 xlim=kenyaLonRange, ylim=kenyaLatRange)
    } else if(thisArea == "Pixel"){
      plot(cbind(popGrid$lon, popGrid$lat), type="n", main=bquote(.(modelName) ~ " estimates"), ylim=kenyaLatRange, 
           xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", asp=1)
      quilt.plot(cbind(popGrid$lon, popGrid$lat), logit(predictionList[[1]]), 
                 nx=150, ny=150, add.legend=FALSE, add=TRUE, col=meanCols, zlim=thisMeanRange)
      plotMapDat(mapDat=thisMap, lwd=.5)
      points(dat$lon, dat$lat, pch=".")
      if(is.null(thisMeanTicks) || is.null(meanTickLabels)) {
        thisMeanTicks = logit(pretty(expit(thisMeanRange), n=5))
        meanTickLabels = as.character(expit(thisMeanTicks))
        image.plot(zlim=thisMeanRange, nlevel=length(meanCols), legend.only=TRUE, horizontal=FALSE,
                   col=meanCols, add = TRUE, axis.args=list(at=thisMeanTicks, labels=meanTickLabels), legend.mar = 0)
      } else {
        image.plot(zlim=thisMeanRange, nlevel=length(meanCols), legend.only=TRUE, horizontal=FALSE,
                   col=meanCols, add = TRUE, axis.args=list(at=thisMeanTicks, labels=meanTickLabels), legend.mar = 0)
      }
      # image.plot(zlim=range(logit(meanRange)), nlevel=length(meanCols), legend.only=TRUE, horizontal=FALSE,
      #            col=meanCols, add = TRUE, axis.args=list(at=logit(meanTicks), labels=meanTickLabels), legend.mar = 0)
    }
    
    # plot the credible interval widths
    
    # plot results in the given column
    # theseResults = results$aggregatedResults$predictions[[paste0(tolower(thisArea), "Predictions")]]
    
    if(is.null(widthRange))
      thisWidthRange = log(range(widthList[[1]]))
    else
      thisWidthRange = log(widthRange)
    if(is.null(widthTicks))
      thisWidthTicks = NULL
    else
      thisWidthTicks = log(widthTicks)
    
    if(thisArea %in% c("Region", "County")) {
      plotMapDat(plotVar=widthList[[1]], new = TRUE, 
                 main=bquote(.(modelName) ~ " 80% CI width"), scaleFun=log, scaleFunInverse=exp, 
                 cols=widthCols, zlim=thisWidthRange, ticks=widthTicks, tickLabels=widthTickLabels, 
                 xlim=kenyaLonRange, ylim=kenyaLatRange)
    } else if(thisArea == "Pixel") {
      
      plot(cbind(popGrid$lon, popGrid$lat), type="n", main=bquote(.(modelName) ~ " 80% CI width"), ylim=kenyaLatRange, 
           xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", asp=1)
      quilt.plot(cbind(popGrid$lon, popGrid$lat), log(widthList[[1]]), 
                 nx=150, ny=150, add.legend=FALSE, add=TRUE, col=widthCols, zlim=thisWidthRange)
      plotMapDat(mapDat=thisMap, lwd=.5)
      points(dat$lon, dat$lat, pch=".")
      if(is.null(thisWidthTicks) || is.null(widthTickLabels)) {
        thisWidthTicks = log(pretty(exp(thisWidthRange), n=5))
        widthTickLabels = as.character(exp(thisWidthTicks))
        image.plot(zlim=thisWidthRange, nlevel=length(widthCols), legend.only=TRUE, horizontal=FALSE,
                   col=widthCols, add = TRUE, axis.args=list(at=thisWidthTicks, labels=widthTickLabels), legend.mar = 0)
      } else {
        image.plot(zlim=thisWidthRange, nlevel=length(widthCols), legend.only=TRUE, horizontal=FALSE,
                   col=widthCols, add = TRUE, axis.args=list(at=thisWidthTicks, labels=widthTickLabels), legend.mar = 0)
      }
    }
    dev.off()
  }
}













##### END OF MAIN PLOTTING FUNCTIONS #####

##### The rest of the functions in this script are utility functions for plotting

makeRedBlueSequentialColors = function(n, ggplot=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  if(!ggplot)
    sequential_hcl(n, h1=10, h2=-115, c1=100, c2=100, l1=44, l2=59, p1=0, p2=2.3)
  else
    scale_colour_continuous_sequential(h1=10, h2=-115, c1=100, c2=100, l1=44, l2=59, p1=0, p2=2.3, n_interp=n)
}

makePurpleRedSequentialColors = function(n, ggplot=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  if(!ggplot)
    sequential_hcl(n, h1=-87, h2=10, c1=65, c2=79, l1=13, l2=60, p1=3, p2=1.6)
  else
    scale_colour_continuous_sequential(h1=-87, h2=10, c1=65, c2=79, l1=13, l2=60, p1=3, p2=1.6, n_interp=n)
}

makePurpleRedDivergingColors = function(n, valRange=NULL, center=NULL, rev=FALSE, ggplot=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  # if(!ggplot)
  #   sequential_hcl(n, h1=-87, h2=10, c1=65, c2=79, l1=13, l2=60, p1=3, p2=1.6)
  # else
  #   scale_colour_continuous_sequential(h1=-87, h2=10, c1=65, c2=79, l1=13, l2=60, p1=3, p2=1.6, n_interp=n)
  
  if(is.null(valRange) && is.null(center)) {
    # diverging_hcl(n, h1=10, h2=-115, c1=90, l1=40, l2=100, p1=0.9, p2=0.6)
    if(!ggplot) {
      divergingx_hcl(n, h1=-87, h3=10, c1=65, c3=79, l1=13, l3=60, p1=3, p2=1.6, h2=0, l2=100, c2=100, p3=3, p4=1.6, rev=rev)
      # diverging_hcl(n, h1=-87, h2=10, c1=65, c2=79, l1=13, l2=60, p1=3, p2=1.6, rev=rev)
    }
    else {
      stop("ggplot doesn't work with makePurpleRedDivergingColors")
      scale_colour_continuous_diverging(h1=-87, h2=10, c1=65, c2=79, l1=13, l2=60, p1=3, p2=1.6, rev=rev, n_interp=n)
    }
  }
  else {
    # in this case we want white to be at the center of valRange if center is NULL
    if(!ggplot) {
      propUp = (valRange[2] - center) / diff(valRange)
      propDown = 1 - propUp
      totalColors = ceiling(2 * max(propUp, propDown) * n)
      tempColors = makePurpleRedDivergingColors(totalColors, rev=rev)
      totalMissingColors = totalColors - n
      
      if(propUp >= propDown)
        tempColors[-(1:totalMissingColors)]
      else
        tempColors[1:n]
    } else {
      stop("ggplot doesn't work with makePurpleRedDivergingColors")
      if(is.null(center))
        center = min(valRange) + abs(diff(valRange))/2
      scale_colour_continuous_diverging(h1=-87, h2=10, c1=65, c2=79, l1=13, l2=60, p1=3, p2=1.6, rev=rev, n_interp=n, mid=center)
    }
  }
}

makeGreenBlueSequentialColors = function(n, ggplot=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  if(!ggplot)
    sequential_hcl(n, h1=128, h2=250, c1=117, cmax=74, c2=107, l1=71, l2=55, p1=2, p2=2)
  else
    scale_colour_continuous_sequential(h1=128, h2=250, c1=117, cmax=74, c2=107, l1=71, l2=55, p1=2, p2=2, n_interp=n)
}

makeGreenBlueDivergingColors = function(n, valRange=NULL, center=NULL, rev=FALSE, ggplot=FALSE, p1=1) {
  # library("colorspace")
  # pal <-choose_palette()
  # if(!ggplot)
  #   sequential_hcl(n, h1=128, h2=250, c1=117, cmax=74, c2=107, l1=71, l2=55, p1=2, p2=2)
  # else
  #   scale_colour_continuous_sequential(h1=128, h2=250, c1=117, cmax=74, c2=107, l1=71, l2=55, p1=2, p2=2, n_interp=n)
  
  if(is.null(valRange) && is.null(center)) {
    if(!ggplot)
      diverging_hcl(n, h1=128, h2=250, c1=100, l1=71, l2=95, p1=p1, rev=rev)
    else
      scale_colour_continuous_diverging(h1=128, h2=250, c1=100, l1=71, l2=95, p1=p1, rev=rev, n_interp=n)
  }
  else {
    # in this case we want white to be at the center of valRange if center is NULL
    if(!ggplot) {
      propUp = (valRange[2] - center) / diff(valRange)
      propDown = 1 - propUp
      totalColors = ceiling(2 * max(propUp, propDown) * n)
      tempColors = makeGreenBlueDivergingColors(totalColors, rev=rev, p1=p1)
      totalMissingColors = totalColors - n
      
      if(propUp >= propDown)
        tempColors[-(1:totalMissingColors)]
      else
        tempColors[1:n]
    } else {
      if(is.null(center))
        center = min(valRange) + abs(diff(valRange))/2
      scale_colour_continuous_diverging(h1=128, h2=250, c1=100, l1=71, l2=95, p1=p1, rev=rev, n_interp=n, mid=center)
    }
  }
}

makePurpleYellowSequentialColors = function(n, rev=FALSE, ggplot=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  if(!ggplot)
    sequential_hcl(n, h1=-100, h2=100, c1=60, cmax=74, c2=100, l1=15, l2=95, p1=2, p2=0.9, rev=rev)
  else
    scale_colour_continuous_sequential(h1=-100, h2=100, c1=60, cmax=74, c2=100, l1=15, l2=95, p1=2, p2=0.9, rev=rev, n_interp=n)
}

makeRedBlueDivergingColors = function(n, valRange=NULL, center=NULL, rev=FALSE, ggplot=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  if(is.null(valRange) && is.null(center)) {
    # diverging_hcl(n, h1=10, h2=-115, c1=90, l1=40, l2=100, p1=0.9, p2=0.6)
    if(!ggplot)
      diverging_hcl(n, h1=10, h2=-115, c1=90, l1=40, l2=100, p1=0.9, rev=rev)
    else
      scale_colour_continuous_diverging(h1=10, h2=-115, c1=90, l1=40, l2=100, p1=0.9, rev=rev, n_interp=n)
  }
  else {
    # in this case we want white to be at the center of valRange if center is NULL
    if(!ggplot) {
      propUp = (valRange[2] - center) / diff(valRange)
      propDown = 1 - propUp
      totalColors = ceiling(2 * max(propUp, propDown) * n)
      tempColors = makeRedBlueDivergingColors(totalColors, rev=rev)
      totalMissingColors = totalColors - n
      
      if(propUp >= propDown)
        tempColors[-(1:totalMissingColors)]
      else
        tempColors[1:n]
    } else {
      if(is.null(center))
        center = min(valRange) + abs(diff(valRange))/2
      scale_colour_continuous_diverging(h1=10, h2=-115, c1=90, l1=40, l2=100, p1=0.9, rev=rev, n_interp=n, mid=center)
    }
  }
}

makeRedGrayBlueDivergingColors = function(n, valRange=NULL, center=NULL, rev=FALSE, ggplot=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  if(is.null(valRange) && is.null(center)) {
    if(!ggplot)
      diverging_hcl(n, h1=10, h2=-115, c1=90, l1=40, l2=90, p1=0.9, rev=rev)
    if(!ggplot)
      scale_colour_continuous_diverging(n.interp, h1=10, h2=-115, c1=90, l1=40, l2=90, p1=0.9, rev=rev)
    # diverging_hcl(n, h1=10, h2=-115, c1=90, l1=40, l2=100, p1=0.9, p2=0.6)
  }
  else {
    # in this case we want white to be at the center of valRange if center is NULL
    if(!ggplot) {
      propUp = (valRange[2] - center) / diff(valRange)
      propDown = 1 - propUp
      totalColors = ceiling(2 * max(propUp, propDown) * n)
      tempColors = makeRedGrayBlueDivergingColors(totalColors, rev=rev)
      totalMissingColors = totalColors - n
      
      if(propUp >= propDown)
        tempColors[-(1:totalMissingColors)]
      else
        tempColors[1:n]
    } else {
      if(is.null(center))
        center = min(valRange) + abs(diff(valRange))/2
      scale_colour_continuous_diverging(n.interp, h1=10, h2=-115, c1=90, l1=40, l2=90, p1=0.9, rev=rev, mid=center)
    }
  }
}

makeBlueSequentialColors = function(n, ggplot=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  # sequential_hcl(n, h1=260, c1=80, l1=30, l2=90, p1=1.5, rev=TRUE)
  if(!ggplot)
    sequential_hcl(n, h1=245, c1=50, cmax=75, l1=20, l2=98, p1=0.8, rev=TRUE)
  else
    scale_colour_continuous_sequential(h1=245, c1=50, cmax=75, l1=20, l2=98, p1=0.8, rev=TRUE, n_interp=n)
}

makeGreenSequentialColors = function(n, ggplot=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  # sequential_hcl(n, h1=260, c1=80, l1=30, l2=90, p1=1.5, rev=TRUE)
  if(!ggplot)
    sequential_hcl(n, h1=128, c1=100, l1=72, l2=95, p1=1.0, rev=TRUE)
  else
    scale_colour_continuous_sequential(h1=128, c1=100, l1=72, l2=95, p1=1.0, rev=TRUE, n_interp=n)
}

makeYellowSequentialColors = function(n, ggplot=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  # sequential_hcl(n, h1=260, c1=80, l1=30, l2=90, p1=1.5, rev=TRUE)
  if(!ggplot)
    sequential_hcl(n, h1=86, c1=100, l1=70, l2=95, p1=1.2, rev=TRUE)
  else
    scale_colour_continuous_sequential(h1=86, c1=100, l1=70, l2=95, p1=1.2, rev=TRUE, n_interp=n)
}

makeBlueGreenYellowSequentialColors = function(n, ggplot=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  if(!ggplot)
    sequential_hcl(n, h1=300, h2=75, c1=40, c2=95, l1=15, l2=90, p1=1.0, p2=1.1)
  else
    scale_colour_continuous_sequential(h1=300, h2=75, c1=40, c2=95, l1=15, l2=90, p1=1.0, p2=1.1, n_interp=n)
}

makeRedYellowBlueColors = function(n, ggplot=FALSE) {
  if(!ggplot)
    divergingx_hcl(n, palette="RdYlBu")
  else
    scale_colour_continuous_sequential(palette="RdYlBu", n_interp=n)
}

makeRedYellowBlueColorsMod = function(n, ggplot=FALSE) {
  # Original parameters:
  # h1 h2  h3  c1 c2 c3 l1 l2 l3  p1
  # 10 85 260 105 45 70 35 98 35 1.5
  if(!ggplot) {
    divergingx_hcl(n,
                   h1=10, h2=85, h3=260,
                   c1=105, c2=45, c3=70,
                   l1=35, l2=98, l3=35, 
                   p1=1.5, p2=1.2, p3=.6, p4=1.2, 
                   cmax1=150, cmax2=10)
  }
  else
    scale_colour_continuous_sequential(palette="RdYlBu", n_interp=n)
}

makeRedYellowBlueColorsMod2 = function(n, ggplot=FALSE) {
  # Original parameters:
  # h1=10; h2=85; h3=260
  # c1=105; c2=45; c3=70
  # l1=35; l2=98; l3=35
  # p1=1.5; p2=1.2; p3=.6; p4=1.2
  # cmax1=150; cmax2=10
  
  h1=15; h2=79; h3=250
  c1=100; c2=100; c3=40
  l1=60; l2=70; l3=33
  p1=1.2; p3=.5; p4=1
  # 2: p1=.5; p2=1
  
  if(!ggplot) {
    # divergingx_hcl(n, 
    #                h1=10, h2=85, h3=260, 
    #                c1=105, c2=45, c3=70, 
    #                l1=80, l2=55, l3=30)
    if((n %% 2) == 0) {
      lower = sequential_hcl(h1=h1, h2=h2, c1=c1, c2=c2, 
                             l1=l1, l2=l2, p1=p1, n)[seq(1, n, by=2)]
      upper = sequential_hcl(h1=h2, h2=h3, c1=c2, c2=c3, 
                             l1=l2, l2=l3, p1=p3, p2=p4, n)[seq(2, n, by=2)]
    } else {
      lower = sequential_hcl(h1=h1, h2=h2, c1=c1, c2=c2, 
                             l1=l1, l2=l2, p1=p1, ceiling(n/2))
      upper = sequential_hcl(h1=h2, h2=h3, c1=c2, c2=c3, 
                             l1=l2, l2=l3, p1=p3, p2=p4, ceiling(n/2))[-1]
    }
    c(lower, upper)
  }
  else
    stop("ggplot version not supported")
}

makeRedYellowSequentialColors = function(n, ggplot=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  if(!ggplot)
    sequential_hcl(n, h1=15, h2=79, c1=100, c2=72, l1=40, l2=90, p1=1.2)
  else
    scale_colour_continuous_sequential(h1=15, h2=79, c1=100, c2=52, l1=55, l2=95, p1=1.2, n_interp=n)
}

makeRedGreenDivergingColors = function(n, ggplot=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  if(!ggplot)
    sequential_hcl(n, h1=265, h2=101, c1=100, l1=50, l2=92, p1=0.6, p2=1.5)
  else
    scale_colour_continuous_sequential(h1=265, h2=101, c1=100, l1=50, l2=92, p1=0.6, p2=1.5, n_interp=n)
}

# given continuous color scale and range, chooses colors based on a set of values
getColorsFromScale = function(vals, valRange=NULL, cols, scaleFun=function(x) {x}, 
                              forceValuesInRange=FALSE) {
  
  if(is.null(valRange)) {
    nas = !is.finite(scaleFun(vals))
    valRange = range(vals[!nas])
  }
  
  if(forceValuesInRange) {
    vals[vals < valRange[1]] = valRange[1]
    vals[vals > valRange[2]] = valRange[2]
  }
  
  valRange = scaleFun(valRange)
  vals = scaleFun(vals)
  vals = vals - valRange[1]
  vals = vals/(valRange[2] - valRange[1])
  col = cols[round(vals*(length(cols)-1))+1]
  
  col
}

plotWithColor = function(x, y, z, zlim=NULL, colScale=tim.colors(), xlab=NULL, ylab=NULL, 
                         legend.mar=7, new=TRUE, scaleFun = function(x) {x}, scaleFunInverse = function(x) {x}, 
                         n.ticks=5, min.n=5, ticks=NULL, tickLabels=NULL, legend.width=1.2, addColorBar=TRUE, 
                         legendArgs=list(), leaveRoomForLegend=TRUE, forceColorsInRange=FALSE, orderI=NULL, 
                         ordering=c("none", "increasing", "decreasing"), colorName = c("col", "bg"), ...) {
  ordering = match.arg(ordering)
  colorName = match.arg(colorName)
  
  # remove NA points
  nas = is.na(x) | is.na(y) | is.na(z)
  if(any(nas)) {
    warning("Removing NAs")
    x = x[!nas]
    y = y[!nas]
    z = z[!nas]
  }
  
  # do setup for ploting data if necessary
  if(is.null(zlim)) {
    nas = !is.finite(scaleFun(z))
    zlim = range(z[!nas])
  }
  
  # order the plotting of the points
  if(is.null(orderI)) {
    if(ordering == "increasing") {
      orderI = sort(z, index.return=TRUE)$ix
    } else if(ordering == "decreasing") {
      orderI = sort(z, decreasing=TRUE, index.return=TRUE)$ix
    } else {
      orderI = 1:length(z)
    }
  }
  x = x[orderI]
  y = y[orderI]
  z = z[orderI]
  
  # if(forceColorsInRange) {
  #   z[z > zlim[2]] = zlim[2]
  #   z[z < zlim[1]] = zlim[1]
  # }
  
  # get colors of points
  cols = getColorsFromScale(z, zlim, colScale, scaleFun, forceColorsInRange)
  
  # generate new plot if necessary
  # browser()
  if(new) {
    # set graphical parameters so the legend won't overlap with plot
    currPar = par()
    newPar = currPar
    newMar = newPar$mar
    newMar[4] = max(newMar[4], legend.mar)
    newPar$mar = newMar
    if(currPar$mar[4] != newMar[4])
      suppressWarnings({par(newPar)})
    
    if(is.null(xlab)) {
      xlab = "X"
    }
    if(is.null(ylab)) {
      ylab = "Y"
    }
    
    # par( oma=c( 0,0,0,6)) # leave room for the legend
    if(colorName == "col") {
      do.call("plot", c(list(x=x, y=y, col=cols, xlab=xlab, ylab=ylab), list(...)))
    } else {
      do.call("plot", c(list(x=x, y=y, bg=cols, xlab=xlab, ylab=ylab), list(...)))
    }
  } else {
    if(colorName == "col") {
      do.call("points", c(list(x=x, y=y, col=cols), list(...)))
    } else {
      do.call("points", c(list(x=x, y=y, bg=cols), list(...)))
    }
  }
  
  if(addColorBar) {
    # add legend
    # par( oma=c(0,0,0,2))
    if(is.null(ticks))
      ticks = scaleFun(pretty(zlim, n=n.ticks, min.n=min.n))
    else
      ticks = scaleFun(ticks)
    if(is.null(tickLabels))
      tickLabels = scaleFunInverse(ticks)
    
    # par( oma=c( 0,0,0,3))
    
    # set list of arguments to image.plot
    legendArgs$zlim=scaleFun(zlim)
    legendArgs$nlevel=length(colScale)
    legendArgs$legend.only=TRUE
    legendArgs$horizontal=FALSE
    legendArgs$col=colScale
    legendArgs$add = TRUE
    if(is.null(legendArgs$axis.args))
      legendArgs$axis.args=list(at=ticks, labels=tickLabels)
    else {
      legendArgs$axis.args$at=ticks
      legendArgs$axis.args$labels=tickLabels
    }
    legendArgs$legend.mar=legend.mar
    legendArgs$legend.width=legend.width
    
    do.call("image.plot", legendArgs)
    
    # image.plot(zlim=zlim, nlevel=length(cols), legend.only=TRUE, horizontal=FALSE, 
    #            col=cols, add = TRUE)
  }
  invisible(NULL)
}

myPairs = function(x, labels, panel = points, ..., horInd = 1:nc, verInd = 1:nc, 
                   lower.panel = panel, upper.panel = panel, diag.panel = NULL, 
                   text.panel = textPanel, label.pos = 0.5 + has.diag/3, line.main = 3, 
                   cex.labels = NULL, font.labels = 1, row1attop = TRUE, gap = 1, 
                   log = "", lims=NULL) 
{
  if (doText <- missing(text.panel) || is.function(text.panel)) 
    textPanel <- function(x = 0.5, y = 0.5, txt, cex, font) text(x, 
                                                                 y, txt, cex = cex, font = font)
  localAxis <- function(side, x, y, i, j, xpd, bg, col = NULL, main, 
                        oma, ...) {
    if(!is.null(lims)) {
      x = lims[[i]]
      y = lims[[j]]
    }
    xpd <- NA
    if (side%%2L == 1L && xl[j]) 
      xpd <- FALSE
    if (side%%2L == 0L && yl[i]) 
      xpd <- FALSE
    if (side%%2L == 1L) 
      Axis(x, side = side, xpd = xpd, ...)
    else Axis(y, side = side, xpd = xpd, ...)
  }
  localPlot <- function(..., main, oma, font.main, cex.main) plot(...)
  localLowerPanel <- function(..., main, oma, font.main, cex.main) lower.panel(...)
  localUpperPanel <- function(..., main, oma, font.main, cex.main) upper.panel(...)
  localDiagPanel <- function(..., main, oma, font.main, cex.main) diag.panel(...)
  dots <- list(...)
  nmdots <- names(dots)
  if (!is.matrix(x)) {
    x <- as.data.frame(x)
    for (i in seq_along(names(x))) {
      if (is.factor(x[[i]]) || is.logical(x[[i]])) 
        x[[i]] <- as.numeric(x[[i]])
      if (!is.numeric(unclass(x[[i]]))) 
        stop("non-numeric argument to 'pairs'")
    }
  }
  else if (!is.numeric(x)) 
    stop("non-numeric argument to 'pairs'")
  panel <- match.fun(panel)
  if ((has.lower <- !is.null(lower.panel)) && !missing(lower.panel)) 
    lower.panel <- match.fun(lower.panel)
  if ((has.upper <- !is.null(upper.panel)) && !missing(upper.panel)) 
    upper.panel <- match.fun(upper.panel)
  if ((has.diag <- !is.null(diag.panel)) && !missing(diag.panel)) 
    diag.panel <- match.fun(diag.panel)
  if (row1attop) {
    tmp <- lower.panel
    lower.panel <- upper.panel
    upper.panel <- tmp
    tmp <- has.lower
    has.lower <- has.upper
    has.upper <- tmp
  }
  nc <- ncol(x)
  if (nc < 2L) 
    stop("only one column in the argument to 'pairs'")
  if (!all(horInd >= 1L && horInd <= nc)) 
    stop("invalid argument 'horInd'")
  if (!all(verInd >= 1L && verInd <= nc)) 
    stop("invalid argument 'verInd'")
  if (doText) {
    if (missing(labels)) {
      labels <- colnames(x)
      if (is.null(labels)) 
        labels <- paste("var", 1L:nc)
    }
    else if (is.null(labels)) 
      doText <- FALSE
  }
  oma <- if ("oma" %in% nmdots) 
    dots$oma
  main <- if ("main" %in% nmdots) 
    dots$main
  if (is.null(oma)) 
    oma <- c(4, 4, if (!is.null(main)) 6 else 4, 4)
  opar <- par(mfcol = c(length(horInd), length(verInd)), mar = rep.int(gap/2, 
                                                                       4), oma = oma)
  on.exit(par(opar))
  dev.hold()
  on.exit(dev.flush(), add = TRUE)
  xl <- yl <- logical(nc)
  if (is.numeric(log)) 
    xl[log] <- yl[log] <- TRUE
  else {
    xl[] <- grepl("x", log)
    yl[] <- grepl("y", log)
  }
  ni <- length(iSet <- if (row1attop) horInd else rev(horInd))
  nj <- length(jSet <- verInd)
  for (j in jSet) for (i in iSet) {
    l <- paste0(if (xl[j]) 
      "x"
      else "", if (yl[i]) 
        "y"
      else "")
    if(is.null(lims)) {
      localPlot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE, 
                type = "n", ..., log = l)
    } else {
      localPlot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE, 
                type = "n", xlim=lims[[j]], ylim=lims[[i]], ..., log = l)
    }
    if (i == j || (i < j && has.lower) || (i > j && has.upper)) {
      box()
      j.odd <- (match(j, jSet) + !row1attop)%%2L
      i.odd <- (match(i, iSet) + !row1attop)%%2L
      if (i == iSet[1L] && (!j.odd || !has.upper || !has.lower)) 
        localAxis(3L, x[, j], x[, i], j, i, ...)
      if (i == iSet[ni] && (j.odd || !has.upper || !has.lower)) 
        localAxis(1L, x[, j], x[, i], j, i, ...)
      if (j == jSet[1L] && (!i.odd || !has.upper || !has.lower)) 
        localAxis(2L, x[, j], x[, i], j, i, ...)
      if (j == jSet[nj] && (i.odd || !has.upper || !has.lower)) 
        localAxis(4L, x[, j], x[, i], j, i, ...)
      mfg <- par("mfg")
      if (i == j) {
        if (has.diag) 
          localDiagPanel(as.vector(x[, i]), ...)
        if (doText) {
          par(usr = c(0, 1, 0, 1))
          if (is.null(cex.labels)) {
            l.wid <- strwidth(labels, "user")
            cex.labels <- max(0.8, min(2, 0.9/max(l.wid)))
          }
          xlp <- if (xl[i]) 
            10^0.5
          else 0.5
          ylp <- if (yl[j]) 
            10^label.pos
          else label.pos
          text.panel(xlp, ylp, labels[i], cex = cex.labels, 
                     font = font.labels)
        }
      }
      else if (i < j) 
        localLowerPanel(as.vector(x[, j]), as.vector(x[, 
                                                       i]), ...)
      else localUpperPanel(as.vector(x[, j]), as.vector(x[, 
                                                          i]), ...)
      if (any(par("mfg") != mfg)) 
        stop("the 'panel' function made a new plot")
    }
    else par(new = FALSE)
  }
  if (!is.null(main)) {
    font.main <- if ("font.main" %in% nmdots) 
      dots$font.main
    else par("font.main")
    cex.main <- if ("cex.main" %in% nmdots) 
      dots$cex.main
    else par("cex.main")
    mtext(main, 3, line.main, outer = TRUE, at = 0.5, cex = cex.main, 
          font = font.main)
  }
  invisible(NULL)
}

# plotExampleGaussianProcess(extraPlotName="Nugget", phi=.05, sigma2=0.5^2)
# plotExampleGaussianProcess(sigma2=0)
plotExampleGaussianProcess = function(resGP=512, mu=0, marginalVariance=1^2, phi=.25, kappa=1, sigma2=.1^2, seed=1, extraPlotName="") {
  require(fields)
  require(RandomFields)
  require(spatstat)
  
  set.seed(seed)
  
  # genMaternGP generates a Matern covariance GP on the unit square
  # with the parameters given the in text
  genMaternGP = function(nsim=1, nx=resGP, ny=resGP, asList=TRUE, coords=NULL, method="circulant") {
    #mu = 4, sigma^2=1.5, phi=0.15, kappa=1, beta=2, tau^2 = 0
    
    # use RFmatern and RFsimulate
    obj = RMmatern(nu=kappa, var=marginalVariance, scale=phi)
    # obj = RMwhittle(nu=kappa, var=sigmasq, scale=phi)
    
    if(is.null(coords)) {
      coordsSet=TRUE
      xs = seq(-1, 1, length=nx)
      ys = seq(-1, 1, length=ny)
      coords = make.surface.grid(list(x=xs, y=ys))
    }
    else
      coordsSet = FALSE
    
    if(method == "instrinsic")
      sims = as.matrix(RFsimulate(RPintrinsic(obj), x=coords[,1], y=coords[,2], n=nsim)) + rnorm(nrow(coords)*nsim, sd=sqrt(sigma2)) + mu
    else if(method == "circulant")
      sims = as.matrix(RFsimulate(RPcirculant(obj), x=coords[,1], y=coords[,2], n=nsim)) + rnorm(nrow(coords)*nsim, sd=sqrt(sigma2)) + mu
    else if(method == "cutoff")
      sims = as.matrix(RFsimulate(RPcutoff(obj), x=coords[,1], y=coords[,2], n=nsim)) + rnorm(nrow(coords)*nsim, sd=sqrt(sigma2)) + mu
    
    list(coords=coords, sims=sims)
  }
    
  GP = genMaternGP()
  GPCoords = GP$coords
  # par(mfrow=c(1,1), family="serif")
  png(paste0("Figures/Illustrations/exampleGP", extraPlotName, ".png"), width=800, height=800)
  quilt.plot(GPCoords, GP$sims, nx=resGP, ny=resGP)
  # axis(1, at=seq(-1, 1, l=3))
  # axis(2, at=seq(-1, 1, l=3))
  dev.off()
}

plotExampleMaternCorrelation = function(effectiveScales = c(.1, .5, 1), sigma2=.1^2) {
  cols = rainbow(length(effectiveScales))
  
  pdf("Figures/Illustrations/matern.pdf", width=5, height=5)
  xs = seq(0, 1, l=200)
  plot(xs, (1/(1 + sqrt(sigma2))) * Matern(xs, effectiveScales[1] / sqrt(8), smoothness = 1), type="l", col=cols[1], 
       main="Matern correlation functions", xlab="Distance", ylab="Correlation", ylim=c(0,1))
  if(length(effectiveScales) >= 1) {
    for(i in 2:length(effectiveScales)) {
      lines(xs, (1/(1 + sqrt(sigma2))) * Matern(xs, effectiveScales[i] / sqrt(8), smoothness = 1), col=cols[i])
    }
    legend("topright", paste0("Effective range=", effectiveScales), lty=1, col=cols)
  }
  points(0, 1, pch=19, cex=.7)
  dev.off()
  
  # pdf("Figures/Illustrations/maternINLA.pdf", width=5, height=5)
  # xs = seq(0, 1, l=200)
  # plot(xs, (1 - sqrt(sigma2)) * inla.matern.cov(x=xs, kappa=1 / (effectiveScales[1] / sqrt(8)), nu = 1, corr=TRUE, d=2), type="l", col=cols[1], 
  #      main="Matern correlation functions", xlab="Distance", ylab="Correlation", ylim=c(0,1))
  # if(length(effectiveScales) >= 1) {
  #   for(i in 2:length(effectiveScales)) {
  #     lines(xs, (1 - sqrt(sigma2)) * inla.matern.cov(x=xs, kappa=1 / (effectiveScales[i] / sqrt(8)), nu = 1, corr=TRUE, d=2), col=cols[i])
  #   }
  #   legend("topright", paste0("Effective range=", effectiveScales), lty=1, col=cols)
  # }
  # points(0, 1, pch=19, cex=.7)
  # dev.off()
}

plotBasisKnots = function(NC=14, nBuffer=5, main="Basis Knot Locations") {
  # get lattice points, prediction points
  nLayer=3
  xRangeDat = c(-1, 1)
  yRangeDat = c(-1, 1)
  latInfo = makeLatGrids(xRangeDat, yRangeDat, NC, nBuffer, nLayer=nLayer)
  
  # get first, coarsest layer
  pdf("Figures/illustrations/basisKnots.pdf", width=5, height=5)
  pts = latInfo[[1]]$latCoords
  border = rbind(c(xRangeDat[1], yRangeDat[1]), 
                 c(xRangeDat[1], yRangeDat[2]), 
                 c(xRangeDat[2], yRangeDat[2]), 
                 c(xRangeDat[2], yRangeDat[1]))
  plot(pts[,1], pts[,2], type="n", xlab="x", ylab="y", main=main)
  polygon(border[,1], border[,2], col=rgb(.5, .5, .5, .5), border=rgb(0, 0, 0, 0))
  points(pts[,1], pts[,2], pch="+", cex=1)
  
  # plot second layer
  pts = latInfo[[2]]$latCoords
  points(pts[,1], pts[,2], pch="+", cex=.6, col="darkGreen")
  
  # plot third layer
  pts = latInfo[[3]]$latCoords
  points(pts[,1], pts[,2], pch="+", cex=.25, col="red")
  dev.off()
}

plotDiskDistanceDistribution = function() {
  r=1
  dDiskDist = function(d, r = 1) {
    out = 4 * d / (pi*r^2) * (acos(d / (2*r)) - d / (2*r) * sqrt(1 - (d / (2*r))^2))
    out[d<0] = 0
    out[d>2*r] = 0
    out
  }
  
  ds = seq(0, 2 * r, l=200)
  pdf("Figures/Illustrations/diskDistanceDistribution.pdf", width=5, height=5)
  plot(ds, dDiskDist(ds, r), type="l", main="Distances between points in disk (radius R)", ylab="Probability density", xlab="Distance", xaxt="n", yaxt="n")
  xtick<-seq(0, 2, by=1)
  axis(side=1, at=xtick, labels = FALSE)
  text(x=xtick,  par("usr")[3]-.02, 
       labels = c("0", "R", "2R"), pos = 1, xpd = TRUE)
  
  ytick<-seq(0, max(dDiskDist(ds, r)), l=5)
  axis(side=2, at=ytick, labels = FALSE)
  text(par("usr")[1]-.01, ytick, 
       labels = c("0", TeX("$\\frac{1}{5R}$"), TeX("$\\frac{2}{5R}$"), TeX("$\\frac{3}{5R}$"), TeX("$\\frac{4}{5R}$")), pos = 2, xpd = TRUE)
  dev.off()
  
  print(ds[which.max(dDiskDist(ds, r))])
}

generatePlottingSymbols = function() {
  urbCols = makeGreenBlueSequentialColors(29)
  pdf('Figures/Illustrations/greenTriangle.pdf', width=.2, height=.2)
  par(mar=c(0,0,0,0))
  plot.new()
  plot.window(xlim=c(-1,1),ylim=c(-1,1), xaxs="i", yaxs="i")
  points(0, -0.5, pch=17, cex=1, col=urbCols[1])
  dev.off()
  
  pdf('Figures/Illustrations/blueSquare.pdf', width=.2, height=.2)
  par(mar=c(0,0,0,0))
  plot.new()
  plot.window(xlim=c(-1,1),ylim=c(-1,1), xaxs="i", yaxs="i")
  points(0, -0.5, pch=15, cex=1, col=urbCols[29])
  dev.off()
  
  pdf('Figures/Illustrations/blueDisc.pdf', width=.2, height=.2)
  par(mar=c(0,0,0,0))
  plot.new()
  plot.window(xlim=c(-1,1),ylim=c(-1,1), xaxs="i", yaxs="i")
  points(0, -0.5, pch=19, cex=1, col="skyblue")
  dev.off()
  
  pdf('Figures/Illustrations/blueCircle', width=.2, height=.2)
  par(mar=c(0,0,0,0))
  plot.new()
  plot.window(xlim=c(-1,1),ylim=c(-1,1), xaxs="i", yaxs="i")
  points(0, -0.5, pch=1, cex=1, col="skyblue")
  dev.off()
  
  pdf('Figures/Illustrations/greenDisc.pdf', width=.2, height=.2)
  par(mar=c(0,0,0,0))
  plot.new()
  plot.window(xlim=c(-1,1),ylim=c(-1,1), xaxs="i", yaxs="i")
  points(0, -0.5, pch=19, cex=1, col="darkgreen")
  dev.off()
  
  pdf('Figures/Illustrations/greenCircle', width=.2, height=.2)
  par(mar=c(0,0,0,0))
  plot.new()
  plot.window(xlim=c(-1,1),ylim=c(-1,1), xaxs="i", yaxs="i")
  points(0, -0.5, pch=1, cex=1, col="darkgreen")
  dev.off()
}

plotPlusAtAngle = function(x, y=NULL, angle=0, size=1, ...) {
  if(is.null(y)) {
    y = x[,2]
    x = x[,1]
  }
  
  rotationMat = rbind(c(cos(angle), -sin(angle)), 
                      c(sin(angle), cos(angle)))
  shiftX = c(0, -1, 0, 1) * size
  shiftY = c(-1, 0, 1, 0) * size
  thesePoints = rbind(shiftX, shiftY)
  thesePoints = rotationMat %*% thesePoints
  thesePoints1 = t(cbind(thesePoints[,1], thesePoints[,3]))
  thesePoints2 = t(cbind(thesePoints[,2], thesePoints[,4]))
  
  for(i in 1:length(x)) {
    thisX = x[i]
    thisY = y[i]
    lines(thisX + thesePoints1[,1], thisY + thesePoints1[,2], ...)
    lines(thisX + thesePoints2[,1], thisY + thesePoints2[,2], ...)
  }
}





