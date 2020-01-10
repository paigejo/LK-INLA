# script for plotting predictions for neonatal mortality rates in Kenya
# source("plotGenerator.R")
resultNameRootLower = "mort"

modelClasses = c(rep("SPDE", 4), rep("LK-INLA", 4))
modelVariations = rep(c("uc", "uC", "Uc", "UC"), 2)
groupPlot = rep(c(FALSE, TRUE), 4)

##### before we make any plots, get the scale on which to put all of them
# first get the file names for the results to load in later
filenames = c()
argList = list(list(urbanEffect = FALSE, clusterEffect = FALSE), 
               list(urbanEffect = FALSE, clusterEffect = TRUE), 
               list(urbanEffect = TRUE, clusterEffect = FALSE), 
               list(urbanEffect = TRUE, clusterEffect = TRUE))

for(i in 1:length(argList)) {
  args = argList[[i]]
  clusterEffect = args$clusterEffect
  urbanEffect = args$urbanEffect
  filenames = c(filenames, paste0("savedOutput/resultsSPDE", resultNameRootLower, "_clusterEffect", clusterEffect, 
                                  "_urbanEffect", urbanEffect, ".RData"))
}
argList = list(list(urbanEffect = FALSE, clusterEffect = FALSE), 
               list(urbanEffect = FALSE, clusterEffect = TRUE), 
               list(urbanEffect = TRUE, clusterEffect = FALSE), 
               list(urbanEffect = TRUE, clusterEffect = TRUE))

for(i in 1:length(argList)) {
  args = argList[[i]]
  clusterEffect = args$clusterEffect
  urbanEffect = args$urbanEffect
  filenames = c(filenames, paste0("savedOutput/resultsLKINLA", resultNameRootLower, "_clusterEffect", clusterEffect, 
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

# set plot tick marks to be reasonable on logit and log scales
meanTicks = pretty(c(.01, meanRange[2]), n=10)
meanTicks = meanTicks[-c(5, 7, 9, 11, 13)]
meanTickLabels = as.character(meanTicks)
widthTicks = pretty(widthRange, n=10)
widthTickLabels = as.character(widthTicks)
# meanTickLabels[c(5, 7, 9, 11, 13)] = ""

# add in a few extra tick marks
# meanTicks = c(.01, meanTicks)
# meanTickLabels = c("0.01", meanTickLabels)

# rather than plot everything with the same color scales, we should use a different color scale for each areal level
# makeAllPlots(dataType="mort", filenames, modelClasses, modelVariations, 
#              c("Region", "County", "Pixel", "Cluster"), 
#              meanRange, meanTicks, meanTickLabels, widthRange, widthTicks, widthTickLabels, 
#              plotNameRoot="Mort", resultNameRoot="Mort", meanCols=rev(makeRedBlueDivergingColors(64)), 
#              widthCols=makeBlueYellowSequentialColors(64), popCols=makeBlueSequentialColors(64), 
#              ncols=29, urbCols=makeGreenBlueSequentialColors(ncols), 
#              plotUrbanMap=FALSE, makeModelPredictions=TRUE)

makeAllPlots(dataType="mort", filenames, modelClasses, modelVariations, 
             c("Region"), 
             regionPredRange, meanTicks, meanTickLabels, regionWidthRange, widthTicks, widthTickLabels, 
             plotNameRoot="Mort", resultNameRoot="Mort", meanCols=rev(makeRedBlueDivergingColors(64)), 
             widthCols=makeBlueYellowSequentialColors(64), popCols=makeBlueSequentialColors(64), 
             ncols=29, urbCols=makeGreenBlueSequentialColors(ncols), 
             plotUrbanMap=FALSE, makeModelPredictions=TRUE)

makeAllPlots(dataType="mort", filenames, modelClasses, modelVariations, 
             c("County"), 
             countyPredRange, meanTicks, meanTickLabels, countyWidthRange, widthTicks, widthTickLabels, 
             plotNameRoot="Mort", resultNameRoot="Mort", meanCols=rev(makeRedBlueDivergingColors(64)), 
             widthCols=makeBlueYellowSequentialColors(64), popCols=makeBlueSequentialColors(64), 
             ncols=29, urbCols=makeGreenBlueSequentialColors(ncols), 
             plotUrbanMap=FALSE, makeModelPredictions=TRUE)

makeAllPlots(dataType="mort", filenames, modelClasses, modelVariations, 
             c("Pixel"), 
             pixelPredRange, meanTicks, meanTickLabels, pixelWidthRange, widthTicks, widthTickLabels, 
             plotNameRoot="Mort", resultNameRoot="Mort", meanCols=rev(makeRedBlueDivergingColors(64)), 
             widthCols=makeBlueYellowSequentialColors(64), popCols=makeBlueSequentialColors(64), 
             ncols=29, urbCols=makeGreenBlueSequentialColors(ncols), 
             plotUrbanMap=FALSE, makeModelPredictions=TRUE)

makeAllPlots(dataType="mort", filenames, modelClasses, modelVariations, 
             c("Cluster"), 
             clusterPredRange, meanTicks, meanTickLabels, clusterWidthRange, widthTicks, widthTickLabels, 
             plotNameRoot="Mort", resultNameRoot="Mort", meanCols=rev(makeRedBlueDivergingColors(64)), 
             widthCols=makeBlueYellowSequentialColors(64), popCols=makeBlueSequentialColors(64), 
             ncols=29, urbCols=makeGreenBlueSequentialColors(ncols), 
             plotUrbanMap=FALSE, makeModelPredictions=TRUE)

# now make the same plots, but only using models including a cluster effect to better highlight 
# differences between models caused by including or leaving out an urban effect
makeAllPlots(dataType="mort", filenames[groupPlot], modelClasses[groupPlot], modelVariations[groupPlot], 
             c("Region"), 
             regionPredRange, meanTicks, meanTickLabels, regionWidthRange, widthTicks, widthTickLabels, 
             plotNameRoot="MortClustModels", resultNameRoot="Mort", meanCols=rev(makeRedBlueDivergingColors(64)), 
             widthCols=makeBlueYellowSequentialColors(64), popCols=makeBlueSequentialColors(64), 
             ncols=29, urbCols=makeGreenBlueSequentialColors(ncols), 
             plotUrbanMap=FALSE, makeModelPredictions=TRUE)

makeAllPlots(dataType="mort", filenames[groupPlot], modelClasses[groupPlot], modelVariations[groupPlot], 
             c("County"), 
             countyPredRange, meanTicks, meanTickLabels, countyWidthRange, widthTicks, widthTickLabels, 
             plotNameRoot="MortClustModels", resultNameRoot="Mort", meanCols=rev(makeRedBlueDivergingColors(64)), 
             widthCols=makeBlueYellowSequentialColors(64), popCols=makeBlueSequentialColors(64), 
             ncols=29, urbCols=makeGreenBlueSequentialColors(ncols), 
             plotUrbanMap=FALSE, makeModelPredictions=TRUE)

makeAllPlots(dataType="mort", filenames[groupPlot], modelClasses[groupPlot], modelVariations[groupPlot], 
             c("Pixel"), 
             pixelPredRange, meanTicks, meanTickLabels, pixelWidthRange, widthTicks, widthTickLabels, 
             plotNameRoot="MortClustModels", resultNameRoot="Mort", meanCols=rev(makeRedBlueDivergingColors(64)), 
             widthCols=makeBlueYellowSequentialColors(64), popCols=makeBlueSequentialColors(64), 
             ncols=29, urbCols=makeGreenBlueSequentialColors(ncols), 
             plotUrbanMap=FALSE, makeModelPredictions=TRUE)

makeAllPlots(dataType="mort", filenames[groupPlot], modelClasses[groupPlot], modelVariations[groupPlot], 
             c("Cluster"), 
             clusterPredRange, meanTicks, meanTickLabels, clusterWidthRange, widthTicks, widthTickLabels, 
             plotNameRoot="MortClustModels", resultNameRoot="Mort", meanCols=rev(makeRedBlueDivergingColors(64)), 
             widthCols=makeBlueYellowSequentialColors(64), popCols=makeBlueSequentialColors(64), 
             ncols=29, urbCols=makeGreenBlueSequentialColors(ncols), 
             plotUrbanMap=FALSE, makeModelPredictions=TRUE)


# printModelPredictionTables("mort", resultNameRoot="Mort", nDigitsPredictions=2)
