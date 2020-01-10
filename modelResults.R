# this script combines the results for any given model when fit to the data sets

# this function fits a given model to the (100) data sets given. Note that fitModelFun must return 
# the following variables (all of them being vectors):
## preds
## sigmas
## lower
## upper
## interceptSummary
## rangeSummary
## sdSummary
## varSummary
# otherVariableNames: if the user wants any other variables to be summarized they should pass the 
#                     variable names in the vector `otherVariableNames'
# otherArgs: a list of other arguments to the input fitModelFun aside from the standard. Examples 
#            include "betas" for the GPpreds function
fitModelToDataSets = function(fitModelFun, dataSets, randomSeeds=NULL, otherVariableNames=NULL, otherArgs=NULL, 
                              maxDataSets=NULL, parClust=cl, includeIntercept=TRUE, seed=123) {
  
  # remove data sets past maxDataSets
  if(!is.null(maxDataSets)) {
    dataSets$xTrain = dataSets$xTrain[,1:maxDataSets]
    dataSets$yTrain = dataSets$yTrain[,1:maxDataSets]
    dataSets$zTrain = dataSets$zTrain[,1:maxDataSets]
    dataSets$xTest = dataSets$xTest[,1:maxDataSets]
    dataSets$yTest = dataSets$yTest[,1:maxDataSets]
    dataSets$zTest = dataSets$zTest[,1:maxDataSets]
  }
  nsim = ncol(dataSets$xTrain)
  
  # set random number seed for generating random number seeds...
  if(!is.null(seed))
    set.seed(seed)
  
  # generate random seeds for each data set
  if(is.null(randomSeeds))
    randomSeeds = sample(1:2000000, ncol(dataSets$xTrain), replace=FALSE)
  else if(length(randomSeeds) != length(dataSets)) {
    warning(paste0("length(randomSeeds) (", length(randomSeeds), ") is not equal to length(dataSets), (", ncol(dataSets$xTrain), "). Regenerating randomSeeds"))
    randomSeeds = sample(1:2000000, ncol(dataSets$xTrain), replace=FALSE)
  }
  
  # this function combines the results from the model fit to each data set, and can be called in parallel
  combineResults = function(...) {
    print("Combining results...")
    
    # not sure why the [[1]] is necessary in the non-parallel case
    if(!doParallel)
      results = list(...)[[1]]
    else
      results = list(...)
    
    # predictive distribution summary statistics
    preds = do.call("cbind", lapply(results, function(x) {c(x$preds)}))
    sigmas = do.call("cbind", lapply(results, function(x) {c(x$sigmas)}))
    lower = do.call("cbind", lapply(results, function(x) {c(x$lower)}))
    upper = do.call("cbind", lapply(results, function(x) {c(x$upper)}))
    
    # parameter estimate summary statistics
    interceptSummary = do.call("rbind", lapply(results, function(x) {c(x$interceptSummary)}))
    interceptSummary = colMeans(matrix(unlist(interceptSummary), ncol=5))
    rangeSummary = do.call("rbind", lapply(results, function(x) {c(x$rangeSummary)}))
    if(!is.null(rangeSummary))
      rangeSummary = colMeans(matrix(unlist(rangeSummary), ncol=5))
    sdSummary = do.call("rbind", lapply(results, function(x) {c(x$sdSummary)}))
    if(!is.null(sdSummary))
      sdSummary = colMeans(matrix(unlist(sdSummary), ncol=5))
    varSummary = do.call("rbind", lapply(results, function(x) {c(x$varSummary)}))
    if(!is.null(varSummary))
      varSummary = colMeans(matrix(unlist(varSummary), ncol=5))
    # if(includeClustEffect) {
    #   nuggetVarSummary = do.call("rbind", lapply(results, function(x) {x$nuggetVarSummary}))
    #   nuggetVarSummary = colMeans(nuggetVarSummary)
    #   nuggetSDSummary = do.call("rbind", lapply(results, function(x) {x$nuggetSDSummary}))
    #   nuggetSDSummary = colMeans(nuggetSDSummary)
    # } else {
    #   nuggetVarSummary = NULL
    #   nuggetSDSummary = NULL
    # }
    
    # summarize compute time
    computeTime = sapply(results, function(x) {c(x$computeTime)})
    computeSummary = data.frame(matrix(c(mean(computeTime), sd(computeTime), quantile(probs=c(0.1, 0.5, 0.9), computeTime)), nrow=1))
    names(computeSummary) = names(interceptSummary)
    
    # summarize other variables if necessary, such as latticeKrig layer weights and other hyperparameters
    otherVariableSummaries = list()
    if(!is.null(otherVariableNames)) {
      for(i in 1:length(otherVariableNames)) {
        thisVariableName = otherVariableNames[i]
        thisVariableSummary = do.call("rbind", lapply(results, function(x) {c(x[[thisVariableName]])}))
        thisVariableSummary = colMeans(matrix(unlist(thisVariableSummary), ncol=5))
        otherVariableSummaries[[thisVariableName]] = thisVariableSummary
      }
    }
    
    list(preds=preds, sigmas=sigmas, lower=lower, upper=upper, 
         interceptSummary=interceptSummary, rangeSummary=rangeSummary, 
         varSummary=varSummary, sdSummary=sdSummary, computeSummary=computeSummary, 
         otherVariableSummaries=otherVariableSummaries)
  }
  
  # this function is a wrapper for the input model fitting function that could, if the user 
  # desires, be run in parallel
  fitModelWrapper = function(i, doSink=FALSE) {
    thisSeed = randomSeeds[i]
    set.seed(thisSeed)
    
    # get this dataset and associated
    xTrain = dataSets$xTrain[,i]
    yTrain = dataSets$yTrain[,i]
    zTrain = dataSets$zTrain[,i]
    xTest = dataSets$xTest[,i]
    yTest = dataSets$yTest[,i]
    # zTest = dataSets$zTest[,i]
    
    if(includeIntercept) {
      xObs = matrix(rep(1, length(zTrain)), ncol=1)
      xPred = matrix(rep(1, length(xTest)), ncol=1)
    }
    else {
      xObs = NULL
      xPred = NULL
    }
    
    standardArgList = list(obsCoords=cbind(xTrain, yTrain), obsValues=zTrain, xObs=xObs, 
             predCoords=cbind(xTest, yTest), xPred=xPred, significanceCI=.8)
    fullArgList = c(standardArgList, otherArgs)
    computeTime = system.time(out <- do.call("fitModelFun", fullArgList))[3]
    c(out, list(computeTime=computeTime))
  }
  
  # compute Bias & MSE & mean(Var) & 80\% coverage for each simulation
  if(is.null(parClust)) {
    ## sequential version
    
    surveyResults = list()
    for(i in 1:nsim) {
      surveyResults = c(surveyResults, list(fitModelWrapper(i, FALSE)))
    }
    results = combineResults(surveyResults)
  } else {
    # parallel version
    
    results = foreach(i = 1:nsim, .combine=combineResults, .verbose=TRUE, .multicombine=TRUE, .export=ls()) %dopar% {
      fitModelWrapper(i, FALSE)
    }
  }
  
  list(results=results, randomSeeds=randomSeeds)
}

# this function fits a given model to the (100) data sets given. Note that fitModelFun must return 
# the following variables (all of them being vectors):
## preds
## sigmas
## lower
## upper
## interceptSummary
## rangeSummary
## sdSummary
## varSummary
# otherVariableNames: if the user wants any other variables to be summarized they should pass the 
#                     variable names in the vector `otherVariableNames'
# otherArgs: a list of other arguments to the input fitModelFun aside from the standard. Examples 
#            include "betas" for the GPpreds function
fitModelToKenyaData = function(fitModelFun, dat, randomSeeds=NULL, otherVariableNames=NULL, otherArgs=NULL, 
                               parClust=cl, urbanEffect=TRUE, clusterEffect=TRUE, seed=123) {
  includeIntercept=TRUE
  
  # set random number seed for generating random number seeds...
  if(!is.null(seed))
    set.seed(seed)
  
  
  list(results=results, randomSeeds=randomSeeds)
}

# given model output, aggregates predictions to the requested levels
# NOTE: for validation, all "obs____" variables must be modified to include the left out cluster 
#       predictions in the correct order, 
aggregateModelResultsKenya = function(results, clusterLevel=TRUE, pixelLevel=TRUE, countyLevel=TRUE, 
                                      regionLevel=TRUE, targetPop=c("women", "children")) {
  targetPop = match.arg(targetPop)
  
  # get information on all the observations
  if(targetPop == "children") {
    out = load("../U5MR/kenyaData.RData")
    dat = mort
  }
  else {
    out = load("../U5MR/kenyaDataEd.RData")
    dat = ed
  }
  # obsValues = dat$y/dat$n
  # obsCoords = cbind(dat$east, dat$north)
  # obsNs = dat$n
  # xObs = matrix(rep(1, length(obsValues)), ncol=1)
  obsUrban = dat$urban
  
  # Get population density grid, adjusted for the target population
  kmres = results$kmres
  if(kmres == 5) {
    if(targetPop == "children") {
      load("../U5MR/popGridAdjusted.RData")
    }
    else {
      load("../U5MR/popGridAdjustedWomen.RData")
    }
  }
  else {
    popGrid = makeInterpPopGrid(kmres, adjustPopSurface=TRUE, targetPop)
  }
  predPts = cbind(popGrid$east, popGrid$north)
  predsUrban = popGrid$urban
  predsCounty = popGrid$admin1
  predsRegion = countyToRegion(predsCounty)
  
  # Cluster level predictions (no aggregation required)
  if(clusterLevel) {
    clusterPredictions = data.frame(areaName=1:length(results$obsPreds), 
                                    urban=results$obsUrban, 
                                    preds=results$obsPreds, 
                                    SDs=results$obsSDs, 
                                    Q10=results$obsLower, 
                                    Q50=results$obsMedian, 
                                    Q90=results$obsUpper)
  } else {
    clusterPredictions = NULL
  }
  
  # Pixel level predictions (no aggregation required)
  if(pixelLevel) {
    pixelPredictions = data.frame(areaName=1:length(results$preds), 
                                  urban=predsUrban, 
                                  preds=results$preds, 
                                  SDs=results$sigmas, 
                                  Q10=results$lower, 
                                  Q50=results$median, 
                                  Q90=results$upper)
  } else {
    pixelPredictions = NULL
  }
  
  # From here onwards, we will need to aggregate predictions over the 
  # population density grid. Use the following function to get  numerical 
  # integration matrix for a given level of areal aggregation
  getIntegrationMatrix = function(areaNames) {
    densities = popGrid$popOrig
    uniqueNames = sort(unique(areaNames))
    
    integrationMatrix = sapply(1:length(uniqueNames), function(i) {
      areaI = areaNames == uniqueNames[i]
      theseDensities = densities
      theseDensities[!areaI] = 0
      theseDensities * (1/sum(theseDensities))
    })
    
    integrationMatrix
  }
  
  # Use the following function to perform the
  # aggregations
  getIntegratedPredictions = function(areaNames, urbanProportions) {
    # get numerical integration matrix
    A = getIntegrationMatrix(areaNames)
    
    # aggregate the prediction matrix
    newPredMat = t(A) %*% results$predMat
    
    # calculate relevant summary statistics
    data.frame(areaName=sort(unique(areaNames)), 
               urban=urbanProportions, 
               preds=rowMeans(newPredMat), 
               SDs=apply(newPredMat, 1, sd), 
               Q10=apply(newPredMat, 1, quantile, probs=.1), 
               Q50=apply(newPredMat, 1, quantile, probs=.5), 
               Q90=apply(newPredMat, 1, quantile, probs=.9))
  }
  
  # County level predictions
  if(countyLevel) {
    load("../U5MR/poppc.RData")
    urbanProportions = poppc$popUrb / poppc$popTotal
    sortI = sort(poppc$County, index.return=TRUE)$ix
    urbanProportions = urbanProportions[sortI]
    
    countyPredictions = getIntegratedPredictions(predsCounty, urbanProportions)
  } else {
    countyPredictions = NULL
  }
  
  # Region level predictions
  if(regionLevel) {
    load("../U5MR/poppr.RData")
    urbanProportions = poppr$popUrb / poppr$popTotal
    sortI = sort(as.character(poppr$Region), index.return=TRUE)$ix
    urbanProportions = urbanProportions[sortI]
    
    regionPredictions = getIntegratedPredictions(predsRegion, urbanProportions)
  } else {
    regionPredictions = NULL
  }
  
  # combine fixed effects and hyperparameter summary tables into one single table
  names(results$fixedEffectSummary) = colnames(results$parameterSummaryTable)
  parameterSummary = rbind(results$fixedEffectSummary, 
                           results$parameterSummaryTable)
  
  # return results
  list(predictions = list(clusterPredictions=clusterPredictions, 
                          pixelPredictions=pixelPredictions, 
                          countyPredictions=countyPredictions, 
                          regionPredictions=regionPredictions), 
       parameterSummary = parameterSummary)
}











