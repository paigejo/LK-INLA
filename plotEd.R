# script for plotting predictions for secondary education completion in Kenya
# source("plotGenerator.R")
resultName = "Ed"
resultNameRootLower = tolower(resultName)

family = "binomial" # logit normal binomial model versus beta binomial
urbanPrior = TRUE # whether to place fine scale urban effective range prior on the fine scale layer

familyText=""
if(family == "binomial")
  familyText = "_LgtN"

# modelClasses = c(rep("SPDE", 2), rep("ELK", 4))
# modelVariations = c("u", "U", "ui", "uI", "Ui", "UI")
modelClasses = c(rep("SPDE", 2), rep("ELK", 2))
modelVariations = rep(c("u", "U"), 2)
# groupPlot = rep(c(FALSE, TRUE), 4)

##### before we make any plots, get the scale on which to put all of them
# first get the file names for the results to load in later
filenames = c()

argList = list(list(urbanEffect = FALSE), 
               list(urbanEffect = TRUE))
for(i in 1:length(argList)) {
  args = argList[[i]]
  urbanEffect = args$urbanEffect
  filenames = c(filenames, paste0("savedOutput/", resultName, "/resultsSPDE", resultNameRootLower, "_urbanEffect", urbanEffect, familyText, ".RData"))
}

# argList = list(list(urbanEffect = FALSE, separateRanges = FALSE), 
#                list(urbanEffect = FALSE, separateRanges = TRUE), 
#                list(urbanEffect = TRUE, separateRanges = FALSE), 
#                list(urbanEffect = TRUE, separateRanges = TRUE))
argList = list(list(urbanEffect = FALSE, separateRanges = TRUE), 
               list(urbanEffect = TRUE, separateRanges = TRUE))
for(i in 1:length(argList)) {
  args = argList[[i]]
  separateRanges = args$separateRanges
  urbanEffect = args$urbanEffect
  
  urbanPriorText = ""
  if(!urbanPrior && separateRanges)
    urbanPriorText = "_noUrbanPrior"
  
  filenames = c(filenames, paste0("savedOutput/", resultName, "/resultsLKINLA", resultNameRootLower, "_separateRanges", separateRanges, 
                    "_urbanEffect", urbanEffect, familyText, urbanPriorText, ".RData"))
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
  recalculate = TRUE
  maxRadius = 500
  cgramList = list()
  # for(j in 5:length(filenames)) {
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
        if(family == "betabinomial") {
          effectiveRangeVals = hyperDraws[2,]
          varVals = hyperDraws[3,]^2
          nuggetVarVals = rep(0, ncol(hyperDraws))
        } else {
          effectiveRangeVals = hyperDraws[1,]
          varVals = hyperDraws[2,]^2
          nuggetVarVals = 1/hyperDraws[3,]
        }
        
        # get range of the data and the SPDE basis function mesh for which to compute the covariograms
        out = load(paste0("dataPointsKenya.RData"))
        xRangeDat = dataPointsKenya$xRange
        yRangeDat = dataPointsKenya$yRange
        mesh = results$fit$mesh
        
        # compute the covariance function for the different hyperparameter samples
        cgram = covarianceDistributionSPDE(effectiveRangeVals, varVals, nuggetVarVals, mesh, xRangeDat=xRangeDat, yRangeDat=yRangeDat, maxRadius=maxRadius)
      } else if(thisModelClass == "ELK") {
        # get lattice information object
        latInfo = results$fit$latInfo
        
        # get hyperparameter draws
        separateRanges = grepl("separateRangesTRUE", filenames[j])
        nLayer = length(latInfo)
        
        if(family == "betabinomial") {
          if(separateRanges)
            alphaI = (1 + nLayer+1 + 1):(1 + nLayer+1 + nLayer-1)
          else
            alphaI = 4:(3+nLayer-1)
          
          nuggetVarVals = rep(0, ncol(hyperDraws))
          if(separateRanges) {
            kappaVals = sweep(2.3/exp(hyperDraws[2:(nLayer+1),]), 2, sapply(latInfo, function(x) {x$latWidth}), "*")
            rhoVals = exp(hyperDraws[nLayer+2,])
          } else {
            latticeWidth = latInfo[[1]]$latWidth
            kappaVals = 2.3/exp(hyperDraws[2,]) * latticeWidth
            rhoVals = exp(hyperDraws[3,])
          }
        } else {
          if(separateRanges)
            alphaI = (nLayer+1 + 1):(nLayer+1 + nLayer-1)
          else
            alphaI = 3:(2+nLayer-1)
          
          nuggetVarVals = 1/hyperDraws[nrow(hyperDraws),]
          if(separateRanges) {
            kappaVals = sweep(2.3/exp(hyperDraws[1:nLayer,]), 2, sapply(latInfo, function(x) {x$latWidth}), "*")
            rhoVals = exp(hyperDraws[nLayer+1,])
          } else {
            latticeWidth = latInfo[[1]]$latWidth
            kappaVals = 2.3/exp(hyperDraws[1,]) * latticeWidth
            rhoVals = exp(hyperDraws[2,])
          }
        }
        
        zSamples = matrix(hyperDraws[alphaI,], ncol=length(alphaI))
        xSamples = t(matrix(apply(zSamples, 1, multivariateExpit), ncol=length(alphaI)))
        xSamples = rbind(xSamples, 1-colSums(xSamples))
        alphaMat = xSamples
        
        # compute the covariance function for many different hyperparameter samples
        cgram = covarianceDistributionLKINLA(latInfo, kappaVals, rhoVals, nuggetVarVals, alphaMat, maxRadius=maxRadius)
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

regionWidthTicks = pretty(regionWidthRange, n=6)
countyWidthTicks = pretty(countyWidthRange, n=6)
regionWidthTickLabels = as.character(regionWidthTicks)
countyWidthTickLabels = as.character(countyWidthTicks)

# meanTickLabels[c(5, 7, 9, 11, 13)] = ""

# add in a few extra tick marks
meanTicks = c(.01, meanTicks)
meanTickLabels = c("0.01", meanTickLabels)

plotNameRoot = paste0("Education", familyText, urbanPriorText)

# first plot data visualizations
plotDataVisualizations(plotUrbanMap=FALSE)

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
             regionPredRange, meanTicks, meanTickLabels, regionWidthRange, regionWidthTicks, regionWidthTickLabels, 
             plotNameRoot=plotNameRoot, resultNameRoot="Ed", meanCols=makeRedBlueDivergingColors(64), 
             widthCols=makeBlueYellowSequentialColors(64), popCols=makeBlueSequentialColors(64), 
             ncols=29, urbCols=makeGreenBlueSequentialColors(ncols), loadResults=TRUE, saveResults=TRUE, 
             plotUrbanMap=FALSE, makeModelPredictions=TRUE, makeCovariograms=TRUE, makePairPlots=TRUE)

makeAllPlots(dataType="ed", filenames, modelClasses, modelVariations, 
             "County", 
             countyPredRange, meanTicks, meanTickLabels, countyWidthRange, countyWidthTicks, countyWidthTickLabels, 
             plotNameRoot=plotNameRoot, resultNameRoot="Ed", meanCols=makeRedBlueDivergingColors(64), 
             widthCols=makeBlueYellowSequentialColors(64), popCols=makeBlueSequentialColors(64), 
             ncols=29, urbCols=makeGreenBlueSequentialColors(ncols), loadResults=TRUE, saveResults=FALSE, 
             plotUrbanMap=FALSE, makeModelPredictions=TRUE, makeCovariograms=FALSE, makePairPlots=TRUE)

makeAllPlots(dataType="ed", filenames, modelClasses, modelVariations, 
             "Pixel", 
             pixelPredRange, meanTicks, meanTickLabels, pixelWidthRange, widthTicks, widthTickLabels, 
             plotNameRoot=plotNameRoot, resultNameRoot="Ed", meanCols=makeRedBlueDivergingColors(64), 
             widthCols=makeBlueYellowSequentialColors(64), popCols=makeBlueSequentialColors(64), 
             ncols=29, urbCols=makeGreenBlueSequentialColors(29), loadResults=TRUE, saveResults=FALSE, 
             plotUrbanMap=FALSE, makeModelPredictions=TRUE, makeCovariograms=FALSE, makePairPlots=TRUE)

makeAllPlots(dataType="ed", filenames, modelClasses, modelVariations, 
             "Cluster", 
             clusterPredRange, meanTicks, meanTickLabels, clusterWidthRange, widthTicks, widthTickLabels, 
             plotNameRoot=plotNameRoot, resultNameRoot="Ed", meanCols=makeRedBlueDivergingColors(64), 
             widthCols=makeBlueYellowSequentialColors(64), popCols=makeBlueSequentialColors(64), 
             ncols=29, urbCols=makeGreenBlueSequentialColors(29), loadResults=TRUE, saveResults=FALSE, 
             plotUrbanMap=FALSE, makeModelPredictions=TRUE, makeCovariograms=FALSE, makePairPlots=TRUE)

# show reduction in oversmoothing
group = c(1, 3, 2, 4)
groupPlotName = paste0(plotNameRoot, "oversmoothing")
lty = c(2, 2, 1, 1)
col = c(do.call("rgb", as.list(rep(.6, 3))), "black", do.call("rgb", as.list(rep(.6, 3))), "black")
makeAllPlots(dataType="ed", filenames[group], modelClasses[group], modelVariations[group], 
             "Region", 
             regionPredRange, meanTicks, meanTickLabels, regionWidthRange, regionWidthTicks, regionWidthTickLabels, 
             plotNameRoot=groupPlotName, resultNameRoot="Ed", meanCols=makeRedBlueDivergingColors(64), 
             widthCols=makeBlueYellowSequentialColors(64), popCols=makeBlueSequentialColors(64), 
             ncols=29, urbCols=makeGreenBlueSequentialColors(ncols), loadResults=TRUE, saveResults=TRUE, 
             plotUrbanMap=FALSE, makeModelPredictions=TRUE, makeCovariograms=TRUE, makePairPlots=TRUE, 
             doModelClassPlots=FALSE, col=col, lty=lty)

makeAllPlots(dataType="ed", filenames[group], modelClasses[group], modelVariations[group], 
             "County", 
             countyPredRange, meanTicks, meanTickLabels, countyWidthRange, countyWidthTicks, countyWidthTickLabels, 
             plotNameRoot=groupPlotName, resultNameRoot="Ed", meanCols=makeRedBlueDivergingColors(64), 
             widthCols=makeBlueYellowSequentialColors(64), popCols=makeBlueSequentialColors(64), 
             ncols=29, urbCols=makeGreenBlueSequentialColors(ncols), loadResults=TRUE, saveResults=FALSE, 
             plotUrbanMap=FALSE, makeModelPredictions=TRUE, makeCovariograms=FALSE, makePairPlots=TRUE)

makeAllPlots(dataType="ed", filenames[group], modelClasses[group], modelVariations[group], 
             "Pixel", 
             pixelPredRange, meanTicks, meanTickLabels, pixelWidthRange, widthTicks, widthTickLabels, 
             plotNameRoot=groupPlotName, resultNameRoot="Ed", meanCols=makeRedBlueDivergingColors(64), 
             widthCols=makeBlueYellowSequentialColors(64), popCols=makeBlueSequentialColors(64), 
             ncols=29, urbCols=makeGreenBlueSequentialColors(29), loadResults=TRUE, saveResults=FALSE, 
             plotUrbanMap=FALSE, makeModelPredictions=TRUE, makeCovariograms=FALSE, makePairPlots=TRUE)

makeAllPlots(dataType="ed", filenames[group], modelClasses[group], modelVariations[group], 
             "Cluster", 
             clusterPredRange, meanTicks, meanTickLabels, clusterWidthRange, widthTicks, widthTickLabels, 
             plotNameRoot=groupPlotName, resultNameRoot="Ed", meanCols=makeRedBlueDivergingColors(64), 
             widthCols=makeBlueYellowSequentialColors(64), popCols=makeBlueSequentialColors(64), 
             ncols=29, urbCols=makeGreenBlueSequentialColors(29), loadResults=TRUE, saveResults=FALSE, 
             plotUrbanMap=FALSE, makeModelPredictions=FALSE, makeCovariograms=FALSE, makePairPlots=TRUE)

# printModelPredictionTables("ed", resultNameRoot="Ed", nDigitsPredictions=2)