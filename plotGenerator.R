library(colorspace)
library(mgcv)

# makes plots as well as parameter estimate tables for the example applications in the manuscript
makeAllPlots = function(dataType=c("ed", "mort"), resultFilenames, modelClasses, modelVariations, 
                        areaLevels=c("Region", "County", "Pixel", "Cluster"), 
                        meanRange, meanTicks, meanTickLabels, widthRange, widthTicks, widthTickLabels, 
                        plotNameRoot="Education", resultNameRoot="Ed", meanCols=makeRedBlueDivergingColors(64), 
                        widthCols=makeBlueYellowSequentialColors(64), popCols=makeBlueSequentialColors(64), 
                        ncols=29, urbCols=makeGreenBlueSequentialColors(ncols), 
                        plotUrbanMap=FALSE, kenyaLatRange=c(-4.6, 5), kenyaLonRange=c(33.5, 42.0), 
                        makeModelPredictions=TRUE, makeCovariograms=TRUE, makePairPlots=TRUE, 
                        loadResults=FALSE, saveResults=!loadResults) {
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
  
  if(makeModelPredictions) {
    plotModelPredictions(dat, resultFilenames, modelClasses, modelVariations, areaLevels, 
                         meanRange, meanTicks, meanTickLabels, widthRange, widthTicks, widthTickLabels, 
                         varName, plotNameRoot, resultNameRoot, meanCols, widthCols, 
                         kenyaLatRange, kenyaLonRange)
  }
  
  if(makeCovariograms) {
    plotCovariograms(dat, resultFilenames, modelClasses, modelVariations, 
                     varName, plotNameRoot, resultNameRoot, 
                     cgramList=NULL, loadResults=loadResults, saveResults=saveResults)
  }
  
  if(makePairPlots) {
    makePairPlots(dat, resultFilenames, modelClasses, modelVariations, 
                  areaLevels, meanRange, meanTicks, meanTickLabels, 
                  plotNameRoot=plotNameRoot, resultNameRoot, 
                  meanCols, ncols, urbCols, kenyaLatRange, kenyaLonRange)
  }
  
  invisible(NULL)
}

# this function plots central estimates and credible interval widths over a map of Kenya 
# 4 plots are made by default: one for each area level
plotModelPredictions = function(dat, resultFilenames, modelClasses, modelVariations, 
                                areaLevels=c("Region", "County", "Pixel", "Cluster"), 
                                meanRange, meanTicks, meanTickLabels, widthRange, widthTicks, widthTickLabels, 
                                varName="education", plotNameRoot="Education", resultNameRoot="Ed", 
                                meanCols=makeRedBlueDivergingColors(64), 
                                widthCols=makeBlueYellowSequentialColors(64), 
                                kenyaLatRange=c(-4.6, 5), kenyaLonRange=c(33.5, 42.0)) {
  plotNameRootLower = tolower(plotNameRoot)
  resultNameRootLower = tolower(resultNameRoot)
  numberModels = length(resultFilenames)
  uniqueModelClasses = unique(modelClasses)
  
  ##### if there are multiple unique model classes, call this function on each one individually 
  ##### as well as on the combination
  if(length(uniqueModelClasses) != 1) {
    for(i in 1:length(uniqueModelClasses)) {
      thisModelClass = uniqueModelClasses[i]
      thisI = modelClasses == thisModelClass
      
      thisResultFilenames = resultFilenames[thisI]
      thisModelClasses = modelClasses[thisI]
      thisModelVariations = modelVariations[thisI]
      plotModelPredictions(dat, thisResultFilenames, thisModelClasses, thisModelVariations, areaLevels, 
                           meanRange, meanTicks, meanTickLabels, widthRange, widthTicks, widthTickLabels, 
                           varName, plotNameRoot, resultNameRoot, meanCols, widthCols, 
                           kenyaLatRange, kenyaLonRange)
    }
  }
  
  # load the fine grid used to approximate continuous prediction
  print("Loading prediction grid and shapefiles...")
  load("../U5MR/popGrid.RData")
  
  # load shape files for plotting
  require(maptools)
  regionMap = readShapePoly("../U5MR/mapData/kenya_region_shapefile/kenya_region_shapefile.shp", delete_null_obj=TRUE, force_ring=TRUE, repair=TRUE)
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
    } else if(thisArea == "County") {
      thisMap = countyMap
    }
    
    print("Plotting central estimates and credible interval widths...")
    width = 400 * numberModels
    png(file=paste0("Figures/", resultNameRoot, "/preds", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=1000)
    
    if(thisArea %in% c("Region", "County"))
      par(mfrow=c(2,numberModels))
    else
      par(mfrow=c(2,numberModels), oma=c( 0,0,0,1.5), mar=c(5.1, 4.1, 4.1, 6))
    
    # first load in the predictions
    predictionList = list()
    widthList = list()
    for(j in 1:numberModels) {
      modelName = paste(modelClasses[j], modelVariations[j])
      
      # load this model and plot results in the given column
      out = load(resultFilenames[j])
      theseResults = results$aggregatedResults$predictions[[paste0(tolower(thisArea), "Predictions")]]
      predictionList = c(predictionList, list(theseResults$preds))
      widthList = c(widthList, list(theseResults$Q90 - theseResults$Q10))
    }
    
    # plot the predictions
    for(j in 1:numberModels) {
      modelName = paste(modelClasses[j], modelVariations[j])
      
      # load this model and plot results in the given column
      # out = load(resultFilenames[j])
      # theseResults = results$aggregatedResults$predictions[[paste0(tolower(thisArea), "Predictions")]]
      
      if(thisArea %in% c("Region", "County")) {
        plotMapDat(plotVar=predictionList[[j]], new = TRUE, 
                   main=paste0(modelName, " estimates"), scaleFun=logit, scaleFunInverse=expit, 
                   cols=meanCols, zlim=logit(meanRange), ticks=meanTicks, tickLabels=meanTickLabels, 
                   xlim=kenyaLonRange, ylim=kenyaLatRange)
      } else if(thisArea == "Pixel"){
        plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0(modelName, " estimates"), ylim=kenyaLatRange, 
             xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", asp=1)
        quilt.plot(cbind(popGrid$lon, popGrid$lat), logit(predictionList[[j]]), 
                   nx=150, ny=150, add.legend=FALSE, add=TRUE, col=meanCols, zlim=range(logit(meanRange)))
        plotMapDat(mapDat=thisMap, lwd=.5)
        points(dat$lon, dat$lat, pch=".")
        image.plot(zlim=range(logit(meanRange)), nlevel=length(meanCols), legend.only=TRUE, horizontal=FALSE,
                   col=meanCols, add = TRUE, axis.args=list(at=logit(meanTicks), labels=meanTickLabels), legend.mar = 0)
      }
    }
    
    # plot the credible interval widths
    for(j in 1:numberModels) {
      modelName = paste(modelClasses[j], modelVariations[j])
      
      # load this model and plot results in the given column
      # out = load(resultFilenames[j])
      # theseResults = results$aggregatedResults$predictions[[paste0(tolower(thisArea), "Predictions")]]
      
      if(thisArea %in% c("Region", "County")) {
        plotMapDat(plotVar=widthList[[j]], new = TRUE, 
                   main=paste0(modelName, " 80% CI width"), scaleFun=log, scaleFunInverse=exp, 
                   cols=widthCols, zlim=log(widthRange), ticks=widthTicks, tickLabels=widthTickLabels, 
                   xlim=kenyaLonRange, ylim=kenyaLatRange)
      } else if(thisArea == "Pixel") {
        
        plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0(modelName, " 80% CI width"), ylim=kenyaLatRange, 
             xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", asp=1)
        quilt.plot(cbind(popGrid$lon, popGrid$lat), log(widthList[[j]]), 
                   nx=150, ny=150, add.legend=FALSE, add=TRUE, col=widthCols, zlim=range(log(widthRange)))
        plotMapDat(mapDat=thisMap, lwd=.5)
        points(dat$lon, dat$lat, pch=".")
        image.plot(zlim=range(log(widthRange)), nlevel=length(widthCols), legend.only=TRUE, horizontal=FALSE,
                   col=widthCols, add = TRUE, axis.args=list(at=log(widthTicks), labels=widthTickLabels), legend.mar = 0)
      }
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
                         kenyaLatRange=c(-4.6, 5), kenyaLonRange=c(33.5, 42.0)) {
  plotNameRootLower = tolower(plotNameRoot)
  resultNameRootLower = tolower(resultNameRoot)
  numberModels = length(resultFilenames)
  uniqueModelClasses = unique(modelClasses)
  
  ##### if there are multiple unique model classes, call this function on each one individually 
  ##### as well as on the combination
  if(length(uniqueModelClasses) != 1) {
    for(i in 1:length(uniqueModelClasses)) {
      thisModelClass = uniqueModelClasses[i]
      thisI = modelClasses == thisModelClass
      
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
      width = 200 * numberModels
      png(file=paste0("Figures/", resultNameRoot, "/pairPlot", plotNameRoot, extraPlotNameRoot, thisArea, ".pdf"), width=width, height=width)
    }
    else {
      width = 2 * numberModels
      pdf(file=paste0("Figures/", resultNameRoot, "/pairPlot", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=width)
    }
    
    # collect predictions and proportion urban per area/point
    predsList = list()
    modelNames = c()
    for(j in 1:numberModels) {
      modelNames = c(modelNames, paste(modelClasses[j], modelVariations[j]))
      
      # load this model and put results in the given column
      out = load(resultFilenames[j])
      theseResults = results$aggregatedResults$predictions[[paste0(tolower(thisArea), "Predictions")]]
      
      # get proportion urban (will be the same for each model)
      colI = cut(as.numeric(theseResults$urban), breaks=seq(0 - .0001, 1 + .0001, l=ncols+1), labels=FALSE)
      theseCols = urbCols[colI]
      
      predsList = c(predsList, list(theseResults$preds))
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
      
      # pairs(valMat, 
      #       modelNames, 
      #       pch=19, cex=.3, lower.panel=my_line, upper.panel = my_line, 
      #       main=paste0(thisArea, " estimate comparisons"))
      # lims = c(list(zlim), list(zlim), list(zlim2), list(zlim2), list(zlim2))
      lims = rep(list(zlim), numberModels)
      myPairs(valMat, 
              modelNames, 
              pch=19, cex=.4, lower.panel=my_line, upper.panel = my_line, 
              main=paste0(thisArea, " estimate comparisons"), 
              lims=lims, oma=c(3,3,6,7))
      image.plot(legend.only = TRUE, zlim=c(0,1), nlevel=ncols, legend.mar=3.3, col=urbCols, add=TRUE, 
                 legend.lab = "Urbanicity", legend.line=1.2, legend.width=.5, legend.shrink=.8, 
                 legend.cex=.8, axis.args=list(cex.axis=.5, tck=-1, hadj=.8))
      dev.off()
    } else {
      lims = rep(list(zlim), numberModels)
      myPairs(valMat, 
              modelNames, 
              pch=19, cex=.4, lower.panel=my_line, upper.panel = my_line, 
              main=paste0(thisArea, " estimate comparisons"), 
              lims=lims, oma=c(3,3,6,7))
      legend("topleft", c("Urban", "Rural"), col=c(urbCols[ncols], urbCols[1]), pch=19)
    }
    dev.off()
  }
}

plotCovariograms = function(dat, resultFilenames, modelClasses, modelVariations, 
                            varName="education", plotNameRoot="Education", resultNameRoot="Ed", 
                            cgramList=NULL, loadResults=FALSE, saveResults=!loadResults) {
  plotNameRootLower = tolower(plotNameRoot)
  resultNameRootLower = tolower(resultNameRoot)
  numberModels = length(resultFilenames)
  uniqueModelClasses = unique(modelClasses)
  
  # first get the covariograms if necessary
  if(is.null(cgramList)) {
    cgramList = list()
    for(j in 1:numberModels) {
      modelName = paste(modelClasses[j], modelVariations[j])
      
      # load this model and get the covariogram if necessary
      print(paste0("Loading ", modelName))
      out = load(resultFilenames[j])
      if(!loadResults || !("cgram" %in% names(results))) {
        hyperDraws = results$fit$hyperMat
        
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
        } else if(thisModelClass == "LK-INLA") {
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
          nuggetVarVals = rep(0, nrow(hyperDraws))
          if(separateRanges) {
            kappaVals = t(sweep(2.3/exp(hyperDraws[2:(nLayer+1),]), 2, sapply(latInfo, function(x) {x$latWidth}), "*"))
            rhoVals = exp(hyperDraws[nLayer+2,])
          } else {
            kappaVals = 2.3/exp(hyperDraws[2,]) * latInfo[[1]]$latWidth
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
  if(length(uniqueModelClasses) != 1) {
    for(i in 1:length(uniqueModelClasses)) {
      thisModelClass = uniqueModelClasses[i]
      thisI = modelClasses == thisModelClass
      
      thisResultFilenames = resultFilenames[thisI]
      thisModelClasses = modelClasses[thisI]
      thisModelVariations = modelVariations[thisI]
      thiscgramList = cgramList[thisI]
      plotCovariograms(dat, thisResultFilenames, thisModelClasses, thisModelVariations, 
                       varName, plotNameRoot, resultNameRoot, thiscgramList, 
                       loadResults=TRUE, saveResults=FALSE)
    }
  }
  
  # reorder the models if necessary so that we have 2 blocks of 4 (note that these permutations are 
  # the inverse permutations of themselves)
  if(numberModels == 8) {
    reordering = c(t(rbind(c(1, 2, 5, 6), 
                           c(3, 4, 7, 8))))
    resultFilenames = resultFilenames[reordering]
    modelClasses = modelClasses[reordering]
    modelVariations = modelVariations[reordering]
    cgramList = cgramList[reordering]
  } else if(numberModels == 4) {
    reordering = c(t(rbind(c(1, 3), 
                           c(2, 4))))
    resultFilenames = resultFilenames[reordering]
    modelClasses = modelClasses[reordering]
    modelVariations = modelVariations[reordering]
    cgramList = cgramList[reordering]
  }
  
  # get range of covariogram values
  yRange = c()
  for(j in 1:numberModels) {
    modelName = paste(modelClasses[j], modelVariations[j])
    
    thiscgram = cgramList[[j]]
    yRange = range(c(yRange, thiscgram$covMean, thiscgram$lowerCov, thiscgram$upperCov))
  }
  
  # modify plot file names if necessary
  if(length(uniqueModelClasses) == 1)
    extraPlotNameRoot = uniqueModelClasses[1]
  else
    extraPlotNameRoot = ""
  
  ##### Plot everything
  print("Plotting covariograms...")
  cols = rainbow(numberModels)
  browser()
  width = 4 * ceiling(numberModels/2)
  pdf(file=paste0("Figures/", resultNameRoot, "/covariograms", plotNameRoot, extraPlotNameRoot, ".pdf"), width=width, height=8)
  par(mfrow=c(2,ceiling(numberModels/2)))
  
  # plot the covariograms separately
  for(j in 1:numberModels) {
    modelName = paste(modelClasses[j], modelVariations[j])
    
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
    
    plot(d, covMean, type="l", main=paste0("Posterior of ", modelName, " covariance function"), xlab="Distance", ylab="Covariance", 
         ylim=yRange)
    lines(d, lowerCov, lty=2)
    lines(d, upperCov, lty=2)
    # lines(d, mixtureCovFun(d), col="green")
    # legend("topright", c("Truth", "Estimate", "80% CI"), lty=c(1, 1, 2), col=c("green", "black", "black"))
    legend("topright", c("Estimate", "80% CI"), lty=c(1, 2), col="black")
  }
  dev.off()
  
  pdf(file=paste0("Figures/", resultNameRoot, "/covariogramsAll", plotNameRoot, extraPlotNameRoot, ".pdf"), width=5, height=5)
  
  # plot the covariograms together
  allModelNames = paste(modelClasses, modelVariations)
  for(j in 1:numberModels) {
    modelName = paste(modelClasses[j], modelVariations[j])
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
    
    if(j == 1) {
      plot(d, covMean, type="l", main=paste0("Covariance estimates and 80% CIs"), xlab="Distance (km)", ylab="Covariance", 
           ylim=yRange, col=cols[j])
    } else {
      lines(d, covMean, col=cols[j])
    }
    
    lines(d, lowerCov, lty=2, col=cols[j])
    lines(d, upperCov, lty=2, col=cols[j])
    # lines(d, mixtureCovFun(d), col="green")
    # legend("topright", c("Truth", "Estimate", "80% CI"), lty=c(1, 1, 2), col=c("green", "black", "black"))
    
  }
  legend("topright", allModelNames[reordering], lty=1, col=cols[reordering], cex=ifelse(numberModels >= 5, .5, 1))
  dev.off()
  
  print("Plotting correlograms...")
  width = 4 * ceiling(numberModels/2)
  pdf(file=paste0("Figures/", resultNameRoot, "/correlograms", plotNameRoot, extraPlotNameRoot, ".pdf"), width=width, height=8)
  par(mfrow=c(2,ceiling(numberModels/2)))
  
  # plot the correlograms
  for(j in 1:numberModels) {
    modelName = paste(modelClasses[j], modelVariations[j])
    
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
    
    plot(d, corMean, type="l", main=paste0("Posterior of ", modelName, " correlation function"), xlab="Distance", ylab="Correlation", 
         ylim=c(0,1))
    lines(d, lowerCor, lty=2)
    lines(d, upperCor, lty=2)
    # lines(d, mixtureCorFun(d), col="green")
    # legend("topright", c("Truth", "Estimate", "80% CI"), lty=c(1, 1, 2), col=c("green", "black", "black"))
    legend("topright", c("Estimate", "80% CI"), lty=c(1, 2), col="black")
  }
  dev.off()
  
  pdf(file=paste0("Figures/", resultNameRoot, "/correlogramsAll", plotNameRoot, extraPlotNameRoot, ".pdf"), width=5, height=5)
  
  # plot the correlograms together
  allModelNames = paste(modelClasses, modelVariations)
  for(j in 1:numberModels) {
    modelName = paste(modelClasses[j], modelVariations[j])
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
    
    if(j == 1) {
      plot(d, corMean, type="l", main=paste0("Correlation estimates and 80% CIs"), xlab="Distance (km)", ylab="Correlation", 
           ylim=yRange, col=cols[j])
    } else {
      lines(d, corMean, col=cols[j])
    }
    
    lines(d, lowerCor, lty=2, col=cols[j])
    lines(d, upperCor, lty=2, col=cols[j])
    # lines(d, mixtureCovFun(d), col="green")
    # legend("topright", c("Truth", "Estimate", "80% CI"), lty=c(1, 1, 2), col=c("green", "black", "black"))
    
  }
  legend("topright", allModelNames[reordering], lty=1, col=cols[reordering], cex=ifelse(numberModels >= 5, .5, 1))
  dev.off()
}

printModelPredictionTables = function(dataType=c("mort", "ed"), resultNameRoot="Ed", areaLevels=c("Region", "County"), 
                                      nDigitsPredictions=3, nDigitsParameters=3, byRow=FALSE) {
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
  
  ##### SPDE estimates
  print("getting SPDE estimates...")
  
  includeUrban = TRUE
  includeCluster = TRUE
  nameRoot = paste0("SPDE", resultNameRootLower, "_clusterEffect", includeCluster, 
                    "_urbanEffect", includeUrban)
  out = load(paste0("savedOutput/results", nameRoot, '.RData'))
  spdeResultsRegion = results$aggregatedResults$predictions$regionPredictions
  spdeResultsCounty = results$aggregatedResults$predictions$countyPredictions
  spdeParameterTable = results$aggregatedResults$parameterSummary
  
  ##### LK-INLA estimates
  print("getting LK-INLA estimates...")
  
  includeUrban = TRUE
  includeCluster = TRUE
  nameRoot = paste0("LKINLA", resultNameRootLower, "_clusterEffect", includeCluster, 
                    "_urbanEffect", includeUrban)
  out = load(paste0("savedOutput/results", nameRoot, '.RData'))
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
    print(paste0("SPDE UC range of predictions: ", diff(range(as.numeric(tab[,2])))))
    print(paste0("SPDE UC median 80% CI width: ", median(as.numeric(tab[,4])-as.numeric(tab[,3]))))
    print(paste0("LK-INLA UC range of predictions: ", diff(range(as.numeric(tab[,5])))))
    print(paste0("LK-INLA UC median 80% CI width: ", median(as.numeric(tab[,7])-as.numeric(tab[,6]))))
    
    # print out the reformatted table
    require(kableExtra)
    fullTab = tab %>%
      kable("latex", escape = F, booktabs = TRUE, format.args=list(drop0trailing=FALSE, scientific=FALSE), 
            align=c("l", rep("r", ncol(tab) - 1)), longtable=TRUE, caption = "Longtable") %>% 
      kable_styling(latex_options =c("repeat_header", "scale_down"))
    numberColumns = c(" "=1, "SPDE UC"=3, "LK-INLA UC"=3)
    print(add_header_above(fullTab, numberColumns, italic=TRUE, bold=TRUE, escape=FALSE, line=TRUE))
  }
  
  
  ##### now construct the parameter tables
  ##### now print the parameter estimates:
  print("printing parameter estimates...")
  browser()
  # SPDE
  
  parameters = round(spdeParameterTable, digits=nDigitsParameters)
  
  # modify the row names do not include the word "Summary"
  allNames = rownames(parameters)
  # rownames(parameters) = unlist(sapply(allNames, strsplit, split="Summary"))
  
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
  
  # LK-INLA
  
  parameters = round(lkinlaParameterTable, digits=nDigitsParameters)
  
  # modify the row names do not include the word "Summary"
  allNames = rownames(parameters)
  # rownames(parameters) = unlist(sapply(allNames, strsplit, split="Summary"))
  
  if(byRow == FALSE)
    nonRangePar = format(parameters[-9,], digits=nDigitsParameters, scientific=FALSE)
  else
    nonRangePar = t(apply(parameters[-9,], 1, format, digits=2, scientific=FALSE))
  formattedParameters = rbind(nonRangePar[1:8,], format(parameters[3,], digits=0), nonRangePar[9:nrow(nonRangePar),])
  rownames(formattedParameters) = c("Intercept", "Urban", "Total Var", "Spatial Var", "Cluster Var", "Total SD", "Spatial SD", "Cluster SD", "Range", 
                                    "alpha1", "alpha2", "alpha3")
  
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
                  align=c("l", rep("r", ncol(parTable) - 1)), longtable=TRUE, caption = "Longtable") %>% 
    kable_styling(latex_options="repeat_header")
  nrow(parTable)
  fullTab = fullTab %>% 
    pack_rows("SPDE UC", 1, 9, escape=FALSE, bold=TRUE, italic=TRUE) %>% 
    pack_rows("LK-INLA UC", 10, nrow(parTable), escape=FALSE, bold=TRUE, italic=TRUE)
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
  regionMap = readShapePoly("../U5MR/mapData/kenya_region_shapefile/kenya_region_shapefile.shp", delete_null_obj=TRUE, force_ring=TRUE, repair=TRUE)
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
                 main=paste0(modelName, " estimates"), scaleFun=logit, scaleFunInverse=expit, 
                 cols=meanCols, zlim=thisMeanRange, ticks=meanTicks, tickLabels=meanTickLabels, 
                 xlim=kenyaLonRange, ylim=kenyaLatRange)
    } else if(thisArea == "Pixel"){
      plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0(modelName, " estimates"), ylim=kenyaLatRange, 
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
                 main=paste0(modelName, " 80% CI width"), scaleFun=log, scaleFunInverse=exp, 
                 cols=widthCols, zlim=thisWidthRange, ticks=widthTicks, tickLabels=widthTickLabels, 
                 xlim=kenyaLonRange, ylim=kenyaLatRange)
    } else if(thisArea == "Pixel") {
      
      plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0(modelName, " 80% CI width"), ylim=kenyaLatRange, 
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

makeRedBlueSequentialColors = function(n) {
  # library("colorspace")
  # pal <-choose_palette()
  sequential_hcl(n, h1=10, h2=-115, c1=100, c2=100, l1=44, l2=59, p1=0, p2=2.3)
}

makeGreenBlueSequentialColors = function(n) {
  # library("colorspace")
  # pal <-choose_palette()
  sequential_hcl(n, h1=128, h2=250, c1=117, cmax=74, c2=107, l1=71, l2=55, p1=2, p2=2)
}

makeRedBlueDivergingColors = function(n, valRange=NULL, center=NULL, rev=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  if(is.null(valRange)  || is.null(center)) {
    diverging_hcl(n, h1=10, h2=-115, c1=90, l1=40, l2=100, p1=0.9, rev=rev)
    # diverging_hcl(n, h1=10, h2=-115, c1=90, l1=40, l2=100, p1=0.9, p2=0.6)
  }
  else {
    # in this case we want white to be at the center of valRange
    propUp = (valRange[2] - center) / diff(valRange)
    propDown = 1 - propUp
    totalColors = ceiling(2 * max(propUp, propDown) * n)
    tempColors = makeRedBlueDivergingColors(totalColors, rev=rev)
    totalMissingColors = totalColors - n
    
    if(propUp >= propDown)
      tempColors[-(1:totalMissingColors)]
    else
      tempColors[1:n]
  }
}

makeRedGrayBlueDivergingColors = function(n, valRange=NULL, center=NULL, rev=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  if(is.null(valRange)  || is.null(center)) {
    diverging_hcl(n, h1=10, h2=-115, c1=90, l1=40, l2=90, p1=0.9, rev=rev)
    # diverging_hcl(n, h1=10, h2=-115, c1=90, l1=40, l2=100, p1=0.9, p2=0.6)
  }
  else {
    # in this case we want white to be at the center of valRange
    propUp = (valRange[2] - center) / diff(valRange)
    propDown = 1 - propUp
    totalColors = ceiling(2 * max(propUp, propDown) * n)
    tempColors = makeRedGrayBlueDivergingColors(totalColors, rev=rev)
    totalMissingColors = totalColors - n
    
    if(propUp >= propDown)
      tempColors[-(1:totalMissingColors)]
    else
      tempColors[1:n]
  }
}

makeBlueSequentialColors = function(n) {
  # library("colorspace")
  # pal <-choose_palette()
  # sequential_hcl(n, h1=260, c1=80, l1=30, l2=90, p1=1.5, rev=TRUE)
  sequential_hcl(n, h1=245, c1=50, cmax=75, l1=20, l2=98, p1=0.8, rev=TRUE)
}

makeBlueYellowSequentialColors = function(n) {
  # library("colorspace")
  # pal <-choose_palette()
  sequential_hcl(n, h1=300, h2=75, c1=40, c2=95, l1=15, l2=90, p1=1.0, p2=1.1)
}

makeRedGreenDivergingColors = function(n) {
  # library("colorspace")
  # pal <-choose_palette()
  diverging_hcl(n, h1=265, h2=101, c1=100, l1=50, l2=92, p1=0.6, p2=1.5)
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









