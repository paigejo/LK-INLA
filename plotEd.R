# script for plotting predictions for secondary education completion in Kenya
# source("plotGenerator.R")
resultName = "Ed"
resultNameRootLower = tolower(resultName)

modelClasses = c(rep("SPDE", 2), rep("LK-INLA", 4))
modelVariations = c("u", "U", "ui", "uI", "Ui", "UI")
groupPlot = rep(c(FALSE, TRUE), 4)

##### before we make any plots, get the scale on which to put all of them
# first get the file names for the results to load in later
filenames = c()

argList = list(list(urbanEffect = FALSE), 
               list(urbanEffect = TRUE))
for(i in 1:length(argList)) {
  args = argList[[i]]
  urbanEffect = args$urbanEffect
  filenames = c(filenames, paste0("savedOutput/", resultName, "/resultsSPDE", resultNameRootLower, "_urbanEffect", urbanEffect, ".RData"))
}

argList = list(list(urbanEffect = FALSE, separateRanges = FALSE), 
               list(urbanEffect = FALSE, separateRanges = TRUE), 
               list(urbanEffect = TRUE, separateRanges = FALSE), 
               list(urbanEffect = TRUE, separateRanges = TRUE))
for(i in 1:length(argList)) {
  args = argList[[i]]
  separateRanges = args$separateRanges
  urbanEffect = args$urbanEffect
  filenames = c(filenames, paste0("savedOutput/", resultName, "/resultsLKINLA", resultNameRootLower, "_separateRanges", separateRanges, 
                    "_urbanEffect", urbanEffect, ".RData"))
}

# now load in the files, and calculate the range of the predictions
regionPredRange = c()
countyPredRange = c()
pixelPredRange = c()
clusterPredRange = c()
regionWidthRange = c()
countyWidthRange = c()
pixelWidthRange = c()
clusterWidthRange = c()
for(i in 1:length(filenames)) {
  load(filenames[i])
  
  allAreas = c("Region", "County", "Pixel", "Cluster")
  for(j in 1:length(allAreas)) {
    thisArea = allAreas[j]
    theseResults = results$aggregatedResults$predictions[[paste0(tolower(thisArea), "Predictions")]]
     
    if(thisArea == "Region") {
      regionPredRange = range(c(regionPredRange, theseResults$preds))
      regionWidthRange = range(c(regionWidthRange, theseResults$Q90 - theseResults$Q10))
    } else if(thisArea == "County") {
      countyPredRange = range(c(countyPredRange, theseResults$preds))
      countyWidthRange = range(c(countyWidthRange, theseResults$Q90 - theseResults$Q10))
    } else if(thisArea == "Pixel") {
      pixelPredRange = range(c(pixelPredRange, theseResults$preds))
      pixelWidthRange = range(c(pixelWidthRange, theseResults$Q90 - theseResults$Q10))
    } else if(thisArea == "Cluster") {
      clusterPredRange = range(c(clusterPredRange, theseResults$preds))
      clusterWidthRange = range(c(clusterWidthRange, theseResults$Q90 - theseResults$Q10))
    }
  }
}
fullPredRange = range(c(regionPredRange, countyPredRange, pixelPredRange, clusterPredRange))
fullWidthRange = range(c(regionWidthRange, countyWidthRange, pixelWidthRange, clusterWidthRange))
meanRange = fullPredRange
widthRange = fullWidthRange

# get correlograms/covariograms
if(FALSE) {
  recalculate = FALSE
  cgramList = list()
  for(j in 1:length(filenames)) {
    modelName = paste(modelClasses[j], modelVariations[j])
    
    # load this model and get the covariogram if necessary
    print(paste0("Loading ", modelName))
    out = load(filenames[j])
    if(!("cgram" %in% names(results)) || recalculate) {
      print("Calculating covariogram...")
      hyperDraws = results$fit$hyperMat
      
      # determine if this model has a cluster effect
      thisModelClass = modelClasses[j]
      
      # hyperparameters will be drawn differently depending on the type of model
      if(thisModelClass == "SPDE") {
        # get hyperparameter draws
        effectiveRangeVals = hyperDraws[2,]
        varVals = hyperDraws[3,]^2
        nuggetVarVals = rep(0, ncol(hyperDraws))
        
        # get range of the data and the SPDE basis function mesh for which to compute the covariograms
        out = load(paste0("dataPointsKenya.RData"))
        xRangeDat = dataPointsKenya$xRange
        yRangeDat = dataPointsKenya$yRange
        mesh = results$fit$mesh
        
        # compute the covariance function for the different hyperparameter samples
        cgram = covarianceDistributionSPDE(effectiveRangeVals, varVals, nuggetVarVals, mesh, xRangeDat=xRangeDat, yRangeDat=yRangeDat)
      } else if(thisModelClass == "LK-INLA") {
        # get lattice information object
        latInfo = results$fit$latInfo
        
        # get hyperparameter draws
        separateRanges = grepl("separateRangesTRUE", filenames[j])
        nLayer = length(latInfo)
        if(separateRanges)
          alphaI = (1 + nLayer+1 + 1):(1 + nLayer+1 + nLayer-1)
        else
          alphaI = 4:(3+nLayer-1)
        zSamples = matrix(hyperDraws[alphaI,], ncol=length(alphaI))
        xSamples = t(matrix(apply(zSamples, 1, multivariateExpit), ncol=length(alphaI)))
        xSamples = rbind(xSamples, 1-colSums(xSamples))
        
        
        nuggetVarVals = rep(0, ncol(hyperDraws))
        if(separateRanges) {
          kappaVals = sweep(2.3/exp(hyperDraws[2:(nLayer+1),]), 2, sapply(latInfo, function(x) {x$latWidth}), "*")
          rhoVals = exp(hyperDraws[nLayer+2,])
        } else {
          latticeWidth = latInfo[[1]]$latWidth
          kappaVals = 2.3/exp(hyperDraws[2,]) * latticeWidth
          rhoVals = exp(hyperDraws[3,])
        }
        alphaMat = xSamples
        
        # compute the covariance function for many different hyperparameter samples
        cgram = covarianceDistributionLKINLA(latInfo, kappaVals, rhoVals, nuggetVarVals, alphaMat)
      } else {
        stop(paste0("Unrecognized model class: ", thisModelClass))
      }
      
      # save results
      results$cgram = cgram
      print("Saving covariogram...")
      save(results, file=filenames[j])
    } else {
      print("Loading covariogram...")
      cgram = results$cgram
    }
    
    # append to our list of covariograms
    cgramList = c(cgramList, list(cgram))
  }
}

# set plot tick marks to be reasonable on logit and log scales
meanTicks = pretty(c(.01, meanRange[2]), n=10)
meanTicks = meanTicks[-c(1, 6, 8, 10, 12)]
meanTickLabels = as.character(meanTicks)
widthTicks = pretty(widthRange, n=10)
widthTickLabels = as.character(widthTicks)
widthTicks = widthTicks[-1]
widthTickLabels = widthTickLabels[-1]

# meanTickLabels[c(5, 7, 9, 11, 13)] = ""

# add in a few extra tick marks
meanTicks = c(.01, meanTicks)
meanTickLabels = c("0.01", meanTickLabels)

# rather than plot everything with the same color scales, we should use a different color scale for each areal level
# makeAllPlots(dataType="ed", filenames, modelClasses, modelVariations, 
#              c("Region", "County", "Pixel", "Cluster"), 
#              meanRange, meanTicks, meanTickLabels, widthRange, widthTicks, widthTickLabels, 
#              plotNameRoot="Education", resultNameRoot="Ed", meanCols=makeRedBlueDivergingColors(64), 
#              widthCols=makeBlueYellowSequentialColors(64), popCols=makeBlueSequentialColors(64), 
#              ncols=29, urbCols=makeGreenBlueSequentialColors(ncols), 
#              plotUrbanMap=FALSE, makeModelPredictions=TRUE)
ncols = 29
makeAllPlots(dataType="ed", filenames, modelClasses, modelVariations, 
             "Region", 
             regionPredRange, meanTicks, meanTickLabels, regionWidthRange, widthTicks, widthTickLabels, 
             plotNameRoot="Education", resultNameRoot="Ed", meanCols=makeRedBlueDivergingColors(64), 
             widthCols=makeBlueYellowSequentialColors(64), popCols=makeBlueSequentialColors(64), 
             ncols=29, urbCols=makeGreenBlueSequentialColors(ncols), loadResults=TRUE, saveResults=TRUE, 
             plotUrbanMap=FALSE, makeModelPredictions=FALSE, makeCovariograms=TRUE, makePairPlots=TRUE)

makeAllPlots(dataType="ed", filenames, modelClasses, modelVariations, 
             "County", 
             countyPredRange, meanTicks, meanTickLabels, countyWidthRange, widthTicks, widthTickLabels, 
             plotNameRoot="Education", resultNameRoot="Ed", meanCols=makeRedBlueDivergingColors(64), 
             widthCols=makeBlueYellowSequentialColors(64), popCols=makeBlueSequentialColors(64), 
             ncols=29, urbCols=makeGreenBlueSequentialColors(ncols), loadResults=TRUE, saveResults=FALSE, 
             plotUrbanMap=FALSE, makeModelPredictions=TRUE, makeCovariograms=FALSE, makePairPlots=TRUE)

makeAllPlots(dataType="ed", filenames, modelClasses, modelVariations, 
             "Pixel", 
             pixelPredRange, meanTicks, meanTickLabels, pixelWidthRange, widthTicks, widthTickLabels, 
             plotNameRoot="Education", resultNameRoot="Ed", meanCols=makeRedBlueDivergingColors(64), 
             widthCols=makeBlueYellowSequentialColors(64), popCols=makeBlueSequentialColors(64), 
             ncols=29, urbCols=makeGreenBlueSequentialColors(29), loadResults=TRUE, saveResults=FALSE, 
             plotUrbanMap=FALSE, makeModelPredictions=TRUE, makeCovariograms=FALSE, makePairPlots=TRUE)

makeAllPlots(dataType="ed", filenames, modelClasses, modelVariations, 
             "Cluster", 
             clusterPredRange, meanTicks, meanTickLabels, clusterWidthRange, widthTicks, widthTickLabels, 
             plotNameRoot="Education", resultNameRoot="Ed", meanCols=makeRedBlueDivergingColors(64), 
             widthCols=makeBlueYellowSequentialColors(64), popCols=makeBlueSequentialColors(64), 
             ncols=29, urbCols=makeGreenBlueSequentialColors(29), loadResults=TRUE, saveResults=FALSE, 
             plotUrbanMap=FALSE, makeModelPredictions=TRUE, makeCovariograms=FALSE, makePairPlots=TRUE)

# now make the same plots, but only using models including a cluster effect to better highlight 
# differences between models caused by including or leaving out an urban effect
makeAllPlots(dataType="ed", filenames[groupPlot], modelClasses[groupPlot], modelVariations[groupPlot], 
             "Region", 
             regionPredRange, meanTicks, meanTickLabels, regionWidthRange, widthTicks, widthTickLabels, 
             plotNameRoot="Education", resultNameRoot="Ed", meanCols=makeRedBlueDivergingColors(64), 
             widthCols=makeBlueYellowSequentialColors(64), popCols=makeBlueSequentialColors(64), 
             ncols=29, urbCols=makeGreenBlueSequentialColors(ncols), loadResults=TRUE, saveResults=FALSE, 
             plotUrbanMap=FALSE, makeModelPredictions=TRUE, makeCovariograms=TRUE, makePairPlots=TRUE)

makeAllPlots(dataType="ed", filenames[groupPlot], modelClasses[groupPlot], modelVariations[groupPlot], 
             "County", 
             countyPredRange, meanTicks, meanTickLabels, countyWidthRange, widthTicks, widthTickLabels, 
             plotNameRoot="Education", resultNameRoot="Ed", meanCols=makeRedBlueDivergingColors(64), 
             widthCols=makeBlueYellowSequentialColors(64), popCols=makeBlueSequentialColors(64), 
             ncols=29, urbCols=makeGreenBlueSequentialColors(ncols), loadResults=TRUE, saveResults=FALSE, 
             plotUrbanMap=FALSE, makeModelPredictions=TRUE, makeCovariograms=FALSE, makePairPlots=TRUE)

makeAllPlots(dataType="ed", filenames[groupPlot], modelClasses[groupPlot], modelVariations[groupPlot], 
             "Pixel", 
             pixelPredRange, meanTicks, meanTickLabels, pixelWidthRange, widthTicks, widthTickLabels, 
             plotNameRoot="Education", resultNameRoot="Ed", meanCols=makeRedBlueDivergingColors(64), 
             widthCols=makeBlueYellowSequentialColors(64), popCols=makeBlueSequentialColors(64), 
             ncols=29, urbCols=makeGreenBlueSequentialColors(ncols), loadResults=TRUE, saveResults=FALSE, 
             plotUrbanMap=FALSE, makeModelPredictions=TRUE, makeCovariograms=FALSE, makePairPlots=TRUE)

makeAllPlots(dataType="ed", filenames[groupPlot], modelClasses[groupPlot], modelVariations[groupPlot], 
             "Cluster", 
             clusterPredRange, meanTicks, meanTickLabels, clusterWidthRange, widthTicks, widthTickLabels, 
             plotNameRoot="Education", resultNameRoot="Ed", meanCols=makeRedBlueDivergingColors(64), 
             widthCols=makeBlueYellowSequentialColors(64), popCols=makeBlueSequentialColors(64), 
             ncols=29, urbCols=makeGreenBlueSequentialColors(ncols), loadResults=TRUE, saveResults=FALSE, 
             plotUrbanMap=FALSE, makeModelPredictions=TRUE, makeCovariograms=FALSE, makePairPlots=TRUE)

# printModelPredictionTables("ed", resultNameRoot="Ed", nDigitsPredictions=2)