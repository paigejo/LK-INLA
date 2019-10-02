# this script combines the results for any given model when fit to the simulated data sets

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
                              maxDataSets=NULL, parClust=cl) {
  
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
  
  # generate random seeds for each data set
  if(is.null(randomSeeds))
    randomSeeds = sample(1:2000000, length(dataSets), replace=FALSE)
  else if(length(randomSeeds) != length(dataSets)) {
    warning(paste0("length(randomSeeds) (", length(randomSeeds), ") is not equal to length(dataSets), (", length(dataSets), "). Regenerating randomSeeds"))
    randomSeeds = sample(1:2000000, length(dataSets), replace=FALSE)
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
    preds = do.call("rbind", lapply(results, function(x) {x$preds}))
    sigmas = do.call("rbind", lapply(results, function(x) {x$sigmas}))
    lower = do.call("rbind", lapply(results, function(x) {x$lower}))
    upper = do.call("rbind", lapply(results, function(x) {x$upper}))
    
    # parameter estimate summary statistics
    interceptSummary = do.call("rbind", lapply(results, function(x) {x$interceptSummary}))
    interceptSummary = colMeans(interceptSummary)
    rangeSummary = do.call("rbind", lapply(results, function(x) {x$rangeSummary}))
    rangeSummary = colMeans(rangeSummary)
    sdSummary = do.call("rbind", lapply(results, function(x) {x$sdSummary}))
    sdSummary = colMeans(sdSummary)
    varSummary = do.call("rbind", lapply(results, function(x) {x$varSummary}))
    varSummary = colMeans(varSummary)
    # if(includeClustEffect) {
    #   nuggetVarSummary = do.call("rbind", lapply(results, function(x) {x$nuggetVarSummary}))
    #   nuggetVarSummary = colMeans(nuggetVarSummary)
    #   nuggetSDSummary = do.call("rbind", lapply(results, function(x) {x$nuggetSDSummary}))
    #   nuggetSDSummary = colMeans(nuggetSDSummary)
    # } else {
    #   nuggetVarSummary = NULL
    #   nuggetSDSummary = NULL
    # }
    
    # summarize other variables if necessary, such as latticeKrig layer weights and other hyperparameters
    otherVariableSummaries = list()
    if(!is.null(otherVariableNames)) {
      for(i in 1:length(otherVariableNames)) {
        thisVariableName = otherVariableNames[i]
        thisVariableSummary = do.call("rbind", lapply(results, function(x) {x[[thisVariableName]]}))
        thisVariableSummary = colMeans(thisVariableSummary)
        otherVariableSummaries[[thisVariableName]] = thisVariableSummary
      }
    }
    
    list(preds=preds, sigmas=sigmas, lower=lower, upper=upper, 
         interceptSummary=interceptSummary, rangeSummary=rangeSummary, 
         varSummary=varSummary, sdSummary=sdSummary, 
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
    
    standardArgList = list(obsCoords=cbind(xTrain, yTrain), obsValues=zTrain, xObs=matrix(rep(1, length(obsValues)), ncol=1), 
             predCoords=cbind(xTest, yTest), xPred = matrix(rep(1, nrow(predCoords)), ncol=1), significanceCI=.8)
    fullArgList = c(standardArgList, otherArgs)
    do.call("fitModelFun", fullArgList)
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