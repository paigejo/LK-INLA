library(kableExtra)

# validate the smoothing models by leaving out data from one county at a time
validateExample = function(dat=NULL, targetPop=c("women", "children"), leaveOutRegion=FALSE, 
                           startI=ifelse(urbanPrior, 5, 1), loadPreviousFit=FALSE, verbose=TRUE, endI=Inf, loadPreviousResults=FALSE, 
                           urbanPrior=TRUE, family=c("betabinomial", "binomial")) {
  targetPop = match.arg(targetPop)
  family = match.arg(family)
  
  # load in relevant data for the given example
  if(targetPop == "women") {
    resultNameRoot="Ed"
    if(is.null(dat)) {
      out = load("../U5MR/kenyaDataEd.RData")
      dat = ed
    }
    load("../U5MR/popGridAdjustedWomen.RData")
  } else if(targetPop == "children") {
    resultNameRoot="Mort"
    if(is.null(dat)) {
      out = load("../U5MR/kenyaData.RData")
      dat = mort
    }
    load("../U5MR/popGridAdjusted.RData")
  }
  resultNameRootLower = tolower(resultNameRoot)
  dataType = resultNameRootLower
  
  familyText=""
  if(family == "binomial")
    familyText = "_LgtN"
  clusterEffect = family == "binomial"
  
  # get region names
  regions = sort(unique(countyToRegion(as.character(dat$admin1))))
  
  ##### Do validation
  resultsListSPDE = list()
  resultsListLKINLA = list()
  
  ##### run SPDE
  argList = list(list(dat = dat, urbanEffect = FALSE), 
                 list(dat = dat, urbanEffect = TRUE))
  otherArguments = list(dataType=dataType, verbose=verbose, loadPreviousFit=loadPreviousFit, family=family, 
                        loadPreviousResults=loadPreviousResults, stratifiedValidation=!leaveOutRegion)
  
  modelNames = c()
  for(i in 1:length(argList)) {
    args = argList[[i]]
    separateRanges = args$separateRanges
    urbanEffect = args$urbanEffect
    fileName = paste0("savedOutput/validation/resultsSPDE", resultNameRootLower, "_urbanEffect", urbanEffect, familyText, "_LORegion", leaveOutRegion)
    if(urbanEffect)
      urbanText = "U"
    else
      urbanText = "u"
    modelNames = c(modelNames, paste0("spde", urbanText))
    
    if(i > endI)
      return(invisible(NULL))
    
    if(startI <= i) {
      print(paste0("Fitting SPDE model with urbanEffect=", urbanEffect, "..."))
      results = do.call("validateSPDEKenyaDat", c(args, otherArguments))
      
      results = list(fit=results)
      # save(results, file=paste0(fileName, ".RData")) # too much space
      
      results$fit$fullModelFit = NULL
      save(results, file=paste0(fileName, "compact.RData"))
    } else {
      print(paste0("Loading SPDE model results with urbanEffect=", urbanEffect, "..."))
      if(file.exists(paste0(fileName, "compact.RData")))
        load(paste0(fileName, "compact.RData"))
      else {
        load(paste0(fileName, ".RData"))
        results$fit$fullModelFit = NULL
        save(results, file=paste0(fileName, "compact.RData"))
      }
    }
    
    resultsListSPDE = c(resultsListSPDE, list(results))
  }
  names(resultsListSPDE) = modelNames
  
  
  ##### run LK-INLA
  argList = list(list(dat = dat, separateRanges = FALSE, urbanEffect = FALSE), 
                 list(dat = dat, separateRanges = FALSE, urbanEffect = TRUE), 
                 list(dat = dat, separateRanges = TRUE, urbanEffect = FALSE), 
                 list(dat = dat, separateRanges = TRUE, urbanEffect = TRUE))
  otherArguments = list(dataType=dataType, loadPreviousFit=loadPreviousFit, family=family, 
                        loadPreviousResults=loadPreviousResults, stratifiedValidation=!leaveOutRegion)
  
  for(i in 1:length(argList)) {
    args = argList[[i]]
    separateRanges = args$separateRanges
    urbanEffect = args$urbanEffect
    
    urbanPriorText = ""
    if(!urbanPrior && separateRanges)
      urbanPriorText = "_noUrbanPrior"
    
    fileName = paste0("savedOutput/validation/resultsLKINLA", resultNameRootLower, "_urbanEffect", urbanEffect, 
                      "_separateRanges", separateRanges, familyText, urbanPriorText, "_LORegion", leaveOutRegion)
    if(separateRanges)
      separateText = "S"
    else
      separateText = "s"
    if(urbanEffect)
      urbanText = "U"
    else
      urbanText = "u"
    modelNames = c(modelNames, paste0("lkinla", urbanText, separateText))
    
    if(i+2 > endI)
      return(invisible(NULL))
    
    if(startI <= i + 2) {
      print(paste0("Fitting LK-INLA model with separateRanges=", separateRanges, " and urbanEffect=", urbanEffect, "..."))
      results = do.call("validateLKINLAKenyaDat", c(args, otherArguments))
      
      results = list(fit=results)
      # save(results, file=paste0(fileName, ".RData")) # too much space
      
      results$fit$fullModelFit = NULL
      save(results, file=paste0(fileName, "compact.RData"))
    } else {
      print(paste0("Loading LK-INLA model results with separateRanges=", separateRanges, " and urbanEffect=", urbanEffect, "..."))
      if(file.exists(paste0(fileName, "compact.RData")))
        load(paste0(fileName, "compact.RData"))
      else {
        load(paste0(fileName, ".RData"))
        results$fit$fullModelFit = NULL
        save(results, file=paste0(fileName, "compact.RData"))
      }
      
      if(names(results$fit)[34] == "") {
        # leftover error from previous run that we can fix without rerunning the whole thing
        names(results$fit)[33:34] = c("singleScores", "singleScoresBinomial")
        save(results, file=paste0(fileName, "compact.RData"))
      }
    }
    
    resultsListLKINLA = c(resultsListLKINLA, list(results))
  }
  names(resultsListLKINLA) = modelNames[3:6]
  
  ##### concatenate all scoring rules into a single table (one for in sample, and one for out of sample)
  allModelNames = c(names(resultsListSPDE), names(resultsListLKINLA))
  
  ## concatenate all in sample results
  scoresInSample = rbind(t(sapply(resultsListSPDE, function(x) {colMeans(x$fit$inSamplePooledScores)})), 
                         t(sapply(resultsListSPDE, function(x) {colMeans(x$fit$inSampleUrbanScores)})), 
                         t(sapply(resultsListSPDE, function(x) {colMeans(x$fit$inSampleRuralScores, na.rm=TRUE)})), 
                         t(sapply(resultsListLKINLA, function(x) {colMeans(x$fit$inSamplePooledScores)})), 
                         t(sapply(resultsListLKINLA, function(x) {colMeans(x$fit$inSampleUrbanScores)})), 
                         t(sapply(resultsListLKINLA, function(x) {colMeans(x$fit$inSampleRuralScores, na.rm=TRUE)})))
  
  # reorder them in order to group them by model rather than by urbanicity
  tempModelISPDE = seq(1, 1+4, by=2)
  tempModelILKINLA = seq(7, 7+4*2, by=4)
  scoresInSample = scoresInSample[c(tempModelISPDE, tempModelISPDE+1, 
                                    tempModelILKINLA, tempModelILKINLA+1, tempModelILKINLA+2, tempModelILKINLA+3),]
  scoresInSample = data.frame(scoresInSample)
  
  # add in column saying how the scores were aggregated (across all clusters, urban clusters, or rural clusters). 
  # Also add in the model name as the name of the row
  scoresInSample = cbind("Subset"=rep(c("Avg", "Urban", "Rural"), 6), scoresInSample)
  # rownames(scoresInSample) = rep(allModelNames, each=3)
  
  ## concatenate all in sample results (accounting for binomial variation)
  scoresInSampleBinomial = rbind(t(sapply(resultsListSPDE, function(x) {colMeans(x$fit$inSamplePooledScoresBinomial)})), 
                                 t(sapply(resultsListSPDE, function(x) {colMeans(x$fit$inSampleUrbanScoresBinomial)})), 
                                 t(sapply(resultsListSPDE, function(x) {colMeans(x$fit$inSampleRuralScoresBinomial, na.rm=TRUE)})), 
                                 t(sapply(resultsListLKINLA, function(x) {colMeans(x$fit$inSamplePooledScoresBinomial)})), 
                                 t(sapply(resultsListLKINLA, function(x) {colMeans(x$fit$inSampleUrbanScoresBinomial)})), 
                                 t(sapply(resultsListLKINLA, function(x) {colMeans(x$fit$inSampleRuralScoresBinomial, na.rm=TRUE)})))
  
  # reorder them in order to group them by model rather than by urbanicity
  tempModelISPDE = seq(1, 1+4, by=2)
  tempModelILKINLA = seq(7, 7+4*2, by=4)
  scoresInSampleBinomial = scoresInSampleBinomial[c(tempModelISPDE, tempModelISPDE+1, 
                                                    tempModelILKINLA, tempModelILKINLA+1, tempModelILKINLA+2, tempModelILKINLA+3),]
  scoresInSampleBinomial = data.frame(scoresInSampleBinomial)
  
  # add in column saying how the scores were aggregated (across all clusters, urban clusters, or rural clusters). 
  # Also add in the model name as the name of the row
  scoresInSampleBinomial = cbind("Subset"=rep(c("Avg", "Urban", "Rural"), 6), scoresInSampleBinomial)
  # rownames(scoresInSample) = rep(allModelNames, each=3)
  
  ## concatenate all leave out region results
  scoresLeaveOutRegion = rbind(t(sapply(resultsListSPDE, function(x) {colMeans(x$fit$pooledScoreTable[,-1])})), 
                               t(sapply(resultsListSPDE, function(x) {colMeans(x$fit$urbanScoreTable[,-1])})), 
                               t(sapply(resultsListSPDE, function(x) {colMeans(x$fit$ruralScoreTable[,-1], na.rm=TRUE)})), 
                               t(sapply(resultsListLKINLA, function(x) {colMeans(x$fit$pooledScoreTable[,-1])})), 
                               t(sapply(resultsListLKINLA, function(x) {colMeans(x$fit$urbanScoreTable[,-1])})), 
                               t(sapply(resultsListLKINLA, function(x) {colMeans(x$fit$ruralScoreTable[,-1], na.rm=TRUE)})))
  
  # reorder them in order to group them by model rather than by urbanicity
  tempModelISPDE = seq(1, 1+4, by=2)
  tempModelILKINLA = seq(7, 7+4*2, by=4)
  scoresLeaveOutRegion = scoresLeaveOutRegion[c(tempModelISPDE, tempModelISPDE+1, 
                                                tempModelILKINLA, tempModelILKINLA+1, tempModelILKINLA+2, tempModelILKINLA+3),]
  scoresLeaveOutRegion = data.frame(scoresLeaveOutRegion)
  
  # add in column saying how the scores were aggregated (across all clusters, urban clusters, or rural clusters). 
  # Also add in the model name as the name of the row
  scoresLeaveOutRegion = cbind("Subset"=rep(c("Avg", "Urban", "Rural"), 6), scoresLeaveOutRegion)
  # rownames(scoresInSample) = rep(allModelNames, each=3)
  
  ## concatenate all leave out region results (accounting for binomial variation)
  scoresLeaveOutRegionBinomial = rbind(t(sapply(resultsListSPDE, function(x) {colMeans(x$fit$pooledScoreTableBinomial[,-1])})), 
                               t(sapply(resultsListSPDE, function(x) {colMeans(x$fit$urbanScoreTableBinomial[,-1])})), 
                               t(sapply(resultsListSPDE, function(x) {colMeans(x$fit$ruralScoreTableBinomial[,-1], na.rm=TRUE)})), 
                               t(sapply(resultsListLKINLA, function(x) {colMeans(x$fit$pooledScoreTableBinomial[,-1])})), 
                               t(sapply(resultsListLKINLA, function(x) {colMeans(x$fit$urbanScoreTableBinomial[,-1])})), 
                               t(sapply(resultsListLKINLA, function(x) {colMeans(x$fit$ruralScoreTableBinomial[,-1], na.rm=TRUE)})))
  
  # reorder them in order to group them by model rather than by urbanicity
  tempModelISPDE = seq(1, 1+4, by=2)
  tempModelILKINLA = seq(7, 7+4*2, by=4)
  scoresLeaveOutRegionBinomial = scoresLeaveOutRegionBinomial[c(tempModelISPDE, tempModelISPDE+1, 
                                                tempModelILKINLA, tempModelILKINLA+1, tempModelILKINLA+2, tempModelILKINLA+3),]
  scoresLeaveOutRegionBinomial = data.frame(scoresLeaveOutRegionBinomial)
  
  # add in column saying how the scores were aggregated (across all clusters, urban clusters, or rural clusters). 
  # Also add in the model name as the name of the row
  scoresLeaveOutRegionBinomial = cbind("Subset"=rep(c("Avg", "Urban", "Rural"), 6), scoresLeaveOutRegionBinomial)
  # rownames(scoresInSample) = rep(allModelNames, each=3)
  
  ## get all leave one out results, remove them from the end sample results
  leaveOneOutI = sapply(c("CPO", "WAIC", "DIC"), function(x) {which(grepl(x, names(scoresInSample)))})
  scoresLeaveOneOut = scoresInSample[,c(1, leaveOneOutI)]
  scoresInSample = scoresInSample[,-leaveOneOutI]
  scoresLeaveOneOutBinomial = scoresInSampleBinomial[,c(1, leaveOneOutI)]
  scoresInSampleBinomial = scoresInSampleBinomial[,-leaveOneOutI]
  
  ##### concatenate binned and single scoring rule results
  # binnedScoringRulesuuAll=averageBinnedScores(binnedScoringRulesuuAll), binnedScoringRulesuUAll=averageBinnedScores(binnedScoringRulesuUAll), 
  # binnedScoringRulesUuAll=averageBinnedScores(binnedScoringRulesUuAll), binnedScoringRulesUUAll=averageBinnedScores(binnedScoringRulesUUAll), 
  # binnedScoringRulesAuAll=averageBinnedScores(binnedScoringRulesAuAll), binnedScoringRulesAUAll=averageBinnedScores(binnedScoringRulesAUAll), 
  # binnedScoringRulesuAAll=averageBinnedScores(binnedScoringRulesuAAll), binnedScoringRulesUAAll=averageBinnedScores(binnedScoringRulesUAAll), 
  # binnedScoringRulesAAAll=averageBinnedScores(binnedScoringRulesAAAll), 
  # binnedScoringRulesuuBinomialAll=averageBinnedScores(binnedScoringRulesuuBinomialAll), binnedScoringRulesuUBinomialAll=averageBinnedScores(binnedScoringRulesuUBinomialAll), 
  # binnedScoringRulesUuBinomialAll=averageBinnedScores(binnedScoringRulesUuBinomialAll), binnedScoringRulesUUBinomialAll=averageBinnedScores(binnedScoringRulesUUBinomialAll), 
  # binnedScoringRulesAuBinomialAll=averageBinnedScores(binnedScoringRulesAuBinomialAll), binnedScoringRulesAUBinomialAll=averageBinnedScores(binnedScoringRulesAUBinomialAll), 
  # binnedScoringRulesuABinomialAll=averageBinnedScores(binnedScoringRulesuABinomialAll), binnedScoringRulesUABinomialAll=averageBinnedScores(binnedScoringRulesUABinomialAll), 
  # binnedScoringRulesAABinomialAll=averageBinnedScores(binnedScoringRulesAABinomialAll), 
  # 
  # singleScores=singleScores, singleScoresBinomial=singleScoresBinomial
  binnedScoringRulesuuAll = c(lapply(resultsListSPDE, function(x) {x$fit$binnedScoringRulesuuAll}), lapply(resultsListLKINLA, function(x) {x$fit$binnedScoringRulesuuAll}))
  names(binnedScoringRulesuuAll) = modelNames
  binnedScoringRulesuUAll = c(lapply(resultsListSPDE, function(x) {x$fit$binnedScoringRulesuUAll}), lapply(resultsListLKINLA, function(x) {x$fit$binnedScoringRulesuUAll}))
  names(binnedScoringRulesuuAll) = modelNames
  binnedScoringRulesUuAll = c(lapply(resultsListSPDE, function(x) {x$fit$binnedScoringRulesUuAll}), lapply(resultsListLKINLA, function(x) {x$fit$binnedScoringRulesUuAll}))
  names(binnedScoringRulesuuAll) = modelNames
  binnedScoringRulesUUAll = c(lapply(resultsListSPDE, function(x) {x$fit$binnedScoringRulesUUAll}), lapply(resultsListLKINLA, function(x) {x$fit$binnedScoringRulesUUAll}))
  names(binnedScoringRulesuuAll) = modelNames
  binnedScoringRulesAuAll = c(lapply(resultsListSPDE, function(x) {x$fit$binnedScoringRulesAuAll}), lapply(resultsListLKINLA, function(x) {x$fit$binnedScoringRulesAuAll}))
  names(binnedScoringRulesuuAll) = modelNames
  binnedScoringRulesAUAll = c(lapply(resultsListSPDE, function(x) {x$fit$binnedScoringRulesAUAll}), lapply(resultsListLKINLA, function(x) {x$fit$binnedScoringRulesAUAll}))
  names(binnedScoringRulesuuAll) = modelNames
  binnedScoringRulesuAAll = c(lapply(resultsListSPDE, function(x) {x$fit$binnedScoringRulesuAAll}), lapply(resultsListLKINLA, function(x) {x$fit$binnedScoringRulesuAAll}))
  names(binnedScoringRulesuuAll) = modelNames
  binnedScoringRulesUAAll = c(lapply(resultsListSPDE, function(x) {x$fit$binnedScoringRulesUAAll}), lapply(resultsListLKINLA, function(x) {x$fit$binnedScoringRulesUAAll}))
  names(binnedScoringRulesuuAll) = modelNames
  binnedScoringRulesAAAll = c(lapply(resultsListSPDE, function(x) {x$fit$binnedScoringRulesAAAll}), lapply(resultsListLKINLA, function(x) {x$fit$binnedScoringRulesAAAll}))
  names(binnedScoringRulesuuAll) = modelNames
  
  binnedScoringRulesuuBinomialAll = c(lapply(resultsListSPDE, function(x) {x$fit$binnedScoringRulesuuBinomialAll}), lapply(resultsListLKINLA, function(x) {x$fit$binnedScoringRulesuuBinomialAll}))
  names(binnedScoringRulesuuAll) = modelNames
  binnedScoringRulesuUBinomialAll = c(lapply(resultsListSPDE, function(x) {x$fit$binnedScoringRulesuUBinomialAll}), lapply(resultsListLKINLA, function(x) {x$fit$binnedScoringRulesuUBinomialAll}))
  names(binnedScoringRulesuuAll) = modelNames
  binnedScoringRulesUuBinomialAll = c(lapply(resultsListSPDE, function(x) {x$fit$binnedScoringRulesUuBinomialAll}), lapply(resultsListLKINLA, function(x) {x$fit$binnedScoringRulesUuBinomialAll}))
  names(binnedScoringRulesuuAll) = modelNames
  binnedScoringRulesUUBinomialAll = c(lapply(resultsListSPDE, function(x) {x$fit$binnedScoringRulesUUBinomialAll}), lapply(resultsListLKINLA, function(x) {x$fit$binnedScoringRulesUUBinomialAll}))
  names(binnedScoringRulesuuAll) = modelNames
  binnedScoringRulesAuBinomialAll = c(lapply(resultsListSPDE, function(x) {x$fit$binnedScoringRulesAuBinomialAll}), lapply(resultsListLKINLA, function(x) {x$fit$binnedScoringRulesAuBinomialAll}))
  names(binnedScoringRulesuuAll) = modelNames
  binnedScoringRulesAUBinomialAll = c(lapply(resultsListSPDE, function(x) {x$fit$binnedScoringRulesAUBinomialAll}), lapply(resultsListLKINLA, function(x) {x$fit$binnedScoringRulesAUBinomialAll}))
  names(binnedScoringRulesuuAll) = modelNames
  binnedScoringRulesuABinomialAll = c(lapply(resultsListSPDE, function(x) {x$fit$binnedScoringRulesuABinomialAll}), lapply(resultsListLKINLA, function(x) {x$fit$binnedScoringRulesuABinomialAll}))
  names(binnedScoringRulesuuAll) = modelNames
  binnedScoringRulesUABinomialAll = c(lapply(resultsListSPDE, function(x) {x$fit$binnedScoringRulesUABinomialAll}), lapply(resultsListLKINLA, function(x) {x$fit$binnedScoringRulesUABinomialAll}))
  names(binnedScoringRulesuuAll) = modelNames
  binnedScoringRulesAABinomialAll = c(lapply(resultsListSPDE, function(x) {x$fit$binnedScoringRulesAABinomialAll}), lapply(resultsListLKINLA, function(x) {x$fit$binnedScoringRulesAABinomialAll}))
  names(binnedScoringRulesuuAll) = modelNames
  
  singleScoresAll = c(lapply(resultsListSPDE, function(x) {x$fit$singleScores}), lapply(resultsListLKINLA, function(x) {x$fit$singleScores}))
  names(singleScoresAll) = modelNames
  singleScoresBinomialAll = c(lapply(resultsListSPDE, function(x) {x$fit$singleScoresBinomial}), lapply(resultsListLKINLA, function(x) {x$fit$singleScoresBinomial}))
  names(singleScoresBinomialAll) = modelNames
  
  ##### Save all scoring rule tables
  if(!urbanPrior)
    urbanPriorText = "_noUrbanPrior"
  fileName = paste0("savedOutput/validation/validationResults", resultNameRoot, familyText, urbanPriorText, "_LORegion", leaveOutRegion, ".RData")
  allScores = list(scoresInSample=scoresInSample, scoresLeaveOneOut=scoresLeaveOneOut, scoresLeaveOutRegion=scoresLeaveOutRegion, 
                   scoresInSampleBinomial=scoresInSampleBinomial, scoresLeaveOneOutBinomial=scoresLeaveOneOutBinomial, scoresLeaveOutRegionBinomial=scoresLeaveOutRegionBinomial, 
                   
                   binnedScoringRulesuuAll=binnedScoringRulesuuAll, 
                   binnedScoringRulesuUAll=binnedScoringRulesuUAll, 
                   binnedScoringRulesUuAll=binnedScoringRulesUuAll, 
                   binnedScoringRulesUUAll=binnedScoringRulesUUAll, 
                   binnedScoringRulesAuAll=binnedScoringRulesAuAll, 
                   binnedScoringRulesAUAll=binnedScoringRulesAUAll, 
                   binnedScoringRulesuAAll=binnedScoringRulesuAAll, 
                   binnedScoringRulesUAAll=binnedScoringRulesUAAll, 
                   binnedScoringRulesAAAll=binnedScoringRulesAAAll, 
                   binnedScoringRulesuuBinomialAll=binnedScoringRulesuuBinomialAll, 
                   binnedScoringRulesuUBinomialAll=binnedScoringRulesuUBinomialAll, 
                   binnedScoringRulesUuBinomialAll=binnedScoringRulesUuBinomialAll, 
                   binnedScoringRulesUUBinomialAll=binnedScoringRulesUUBinomialAll, 
                   binnedScoringRulesAuBinomialAll=binnedScoringRulesAuBinomialAll, 
                   binnedScoringRulesAUBinomialAll=binnedScoringRulesAUBinomialAll, 
                   binnedScoringRulesuABinomialAll=binnedScoringRulesuABinomialAll, 
                   binnedScoringRulesUABinomialAll=binnedScoringRulesUABinomialAll, 
                   binnedScoringRulesAABinomialAll=binnedScoringRulesAABinomialAll, 
                   
                   singleScoresAll=singleScoresAll, 
                   singleScoresBinomialAll=singleScoresBinomialAll)
  save(allScores, file=fileName)
  
  allScores
}

# print out validation results and plot the PITs
printValidationResults = function(resultNameRoot="Ed") {
  require(stringr)
  require(dplyr)
  require(kableExtra)
  
  # first load the validation results
  out = load(paste0("resultsValidationAll", resultNameRoot, ".RData"))
  modelNames = rownames(validationResults$scoresInSample)
  
  # change the modelNames to the adjusted ones:
  badNames = c("IV", "III", "II", "I")
  goodNames = c("UC", "Uc", "uC", "uc")
  for(i in 1:length(badNames)) {
    badName = badNames[i]
    goodName = goodNames[i]
    whichI = grepl(badName, modelNames)
    modelNames[whichI] = gsub(badName, goodName, modelNames[whichI])
  }
  rownames(validationResults$scoresInSample) = modelNames
  rownames(validationResults$scoresLeaveOutCounty) = modelNames
  rownames(validationResults$scoresLeaveOutCluster) = modelNames[-1]
  
  modelCodes = gsub(" ", "", modelNames, fixed = TRUE)
  modelCodes[3] = "BYM2uClust"
  modelCodes[4] = "BYM2uClust'"
  modelCodes[5] = "BYM2Urbc"
  modelCodes[6] = "BYM2UrbClust"
  modelCodes[7] = "BYM2UrbClust'"
  modelCodes[9] = "SPDEuClust"
  modelCodes[10] = "SPDEUrbc"
  modelCodes[11] = "SPDEUrbClust"
  ## plot the PITs
  for(i in 1:length(modelNames)) {
    thisModelCode = modelCodes[i]
    pdf(paste0("figures/validationPITHistogram", resultNameRoot, "_", thisModelCode, ".pdf"), width=9, height=5)
    par(mfrow=c(1,2))
    hist(unlist(validationResults$pitInSample[[i]]), main=paste0(modelNames[i], " In Sample PITs"), 
         xlab="PIT", freq=FALSE, breaks=30, col="skyblue")
    
    hist(unlist(validationResults$pitLeaveOutCounty[[i]]), main=paste0(modelNames[i], " Leave Out County PITs"), 
         xlab="PIT", freq=FALSE, breaks=30, col="skyblue")
    dev.off()
    
    # only CPO is calculated when just leaving out individual clusters
    # if(i != 1) {
    #   hist(unlist(validationResults$pitLeaveOutCluster[[i]]), main=paste0("Histogram of ", modelNames[i], " Left Out Cluster PITs"), 
    #        xlab="PIT", freq=FALSE, breaks=30, col="skyblue")
    # }
  }
  
  ## modify the BYM2 model names 
  tab = validationResults$scoresInSample
  
  ## print out the scoring rule tables
  print("In sample results:")
  displayRow = "s"
  displayIC = rep("f", 2)
  displayMSE = rep("f", 3)
  displayVar = rep("f", 3)
  displayBias = rep("f", 3)
  displayCRPS = "f"
  display = c(displayRow, displayIC, displayMSE, displayVar, displayBias, displayCRPS)
  scalingPower = c(MSE=2, Var=3, Bias=3)
  if(resultNameRoot == "Mort") {
    scalingPower = c(MSE=4, Var=4, Bias=4)
  }
  colScaleIC = rep(1, 2)
  colScaleMSE = rep(10^scalingPower[1], 3)
  colScaleVar = rep(10^scalingPower[2], 3)
  colScaleBias = rep(10^scalingPower[3], 3)
  colScaleCRPS = 1
  colScale = c(colScaleIC, colScaleMSE, colScaleVar, colScaleBias, colScaleCRPS)
  colUnitsIC = rep("", 2)
  colUnitsMSE = rep(paste0(" ($\\times 10^{-", as.character(scalingPower[1]), "}$)"), 3)
  colUnitsVar = rep(paste0(" ($\\times 10^{-", as.character(scalingPower[2]), "}$)"), 3)
  colUnitsBias = rep(paste0(" ($\\times 10^{-", as.character(scalingPower[3]), "}$)"), 3)
  colUnitsCRPS = ""
  colUnits = c(colUnitsIC, colUnitsMSE, colUnitsVar, colUnitsBias, colUnitsCRPS)
  digitsIC = rep(0, 2)
  digitsMSE = rep(1, 3)
  digitsVar = rep(1, 3)
  digitsBias = rep(1, 3)
  digitsCRPS = 2
  colDigits = c(digitsIC, digitsMSE, digitsVar, digitsBias, digitsCRPS)
  
  tab = validationResults$scoresInSample
  for(i in 1:ncol(tab)) {
    tab[,i] = as.numeric(round(tab[,i] * colScale[i], digits=colDigits[i]))
    colnames(tab)[i] = paste0(colnames(tab)[i], colUnits[i])
  }
  colnames(tab) = gsub(".", " ", colnames(tab), fixed=TRUE)
  colnames(tab) = gsub("Urban", "Urb", colnames(tab), fixed=TRUE)
  colnames(tab) = gsub("Rural", "Rur", colnames(tab), fixed=TRUE)
  tab = tab[,-1] # filter out WAIC since it doesn't work for spatially correlated data
  colDigits=colDigits[-1]
  display=display[-1]
  colUnits = colUnits[-1]
  # print(xtable(tab, digits=c(1,colDigits), display=display), 
  #       include.colnames=TRUE,
  #       hline.after=0, 
  #       math.style.exponents=TRUE, 
  #       sanitize.text.function=function(x){x}, 
  #       scalebox=0.5)
  tab = t(tab)
  colnames(tab) = gsub("SPDE ", "", colnames(tab), fixed=TRUE)
  colnames(tab) = gsub("BYM2 ", "", colnames(tab), fixed=TRUE)
  colnames(tab) = gsub("Smoothed Direct", " ", colnames(tab), fixed=TRUE)
  tab = tab[c(2:11, 1), ]
  colUnits = colUnits[c(2:11, 1)]
  # temp = xtable2kable(xtable(tab), 
  #                     include.colnames=TRUE,
  #                     hline.after=0, 
  #                     math.style.exponents=TRUE, 
  #                     sanitize.text.function=function(x){x})
  # print(add_header_above(kable_styling(temp), c(" " = 1, "Smoothed Direct" = 1, "BYM2"=4, "SPDE"=4)))
  
  # remove the smoothed direct results
  tab = tab[,-1]
  
  # bold the best entries of each row, italicize worst entries of each row
  centers = c(rep(0, nrow(tab)))
  rowBest = apply(abs(sweep(tab, 1, centers, "-")), 1, min, na.rm=TRUE)
  rowWorst = apply(abs(sweep(tab, 1, centers, "-")), 1, max, na.rm=TRUE)
  boldFun = function(i) {
    vals = tab[i,]
    isNA = is.na(vals)
    out = rep(FALSE, length(vals))
    out[!isNA] = abs(tab[i,!isNA] - centers[i]) <= rowBest[i]
    out
  }
  italicFun = function(i) {
    vals = tab[i,]
    isNA = is.na(vals)
    out = rep(FALSE, length(vals))
    out[!isNA] = abs(tab[i,!isNA] - centers[i]) >= rowWorst[i]
    out
  }
  test = t(sapply(1:nrow(tab), function(i) {cell_spec(tab[i,], "latex", bold=boldFun(i), italic=italicFun(i), 
                                                      monospace=FALSE, underline=FALSE, strikeout=FALSE)}))
  
  # revert the column names to their true values
  colnames(test) = colnames(tab)
  scoreVariations = c("Avg", "Urban", "Rural")
  scoreVariations = c(rep(scoreVariations, 3), "Avg", "Avg")
  test = cbind(" "=scoreVariations, test)
  rownames(test)=NULL
  
  # group the rows by urbanicity
  fullTab = test %>%
    kable("latex", escape = F, booktabs = T, format.args=list(drop0trailing=FALSE, scientific=FALSE), 
          align=c("l", rep("r", ncol(test) - 1))) %>% kable_styling()
  uniqueScoreTypes = c("MSE", "Var", "Bias", "CRPS", "DIC")
  uniqueColUnits = colUnits[c(1, 4, 7, 10:11)]
  for(i in 1:length(uniqueScoreTypes)) {
    uniqueScoreTypes[i] = paste0(uniqueScoreTypes[i], uniqueColUnits[i])
  }
  scoreTypeGroups = c(rep(uniqueScoreTypes[1:3], each=3), uniqueScoreTypes[4:5])
  
  for(i in 1:length(uniqueScoreTypes)) {
    startR = match(TRUE, scoreTypeGroups == uniqueScoreTypes[i])
    endR = nrow(tab) - match(TRUE, rev(scoreTypeGroups == uniqueScoreTypes[i])) + 1
    fullTab = fullTab %>% pack_rows(uniqueScoreTypes[i], startR, endR, latex_gap_space = "0.3em", escape=FALSE)
  }
  print(add_header_above(fullTab, c(" " = 1, "BYM2"=4, "SPDE"=4), italic=TRUE, bold=TRUE, escape=FALSE))
  
  # print(xtable(tab, digits=digits, display=display), 
  #       include.colnames=TRUE,
  #       hline.after=0, 
  #       math.style.exponents=TRUE, 
  #       sanitize.text.function=function(x){x}, 
  #       scalebox=0.5)
  
  print("")
  print("Leave out county results:")
  
  scalingPower = c(MSE=2, Var=3, Bias=3)
  if(resultNameRoot == "Mort") {
    scalingPower = c(MSE=4, Var=4, Bias=4)
  }
  displayRow = "s"
  displayMSE = rep("f", 3)
  displayVar = rep("f", 3)
  displayBias = rep("f", 3)
  displayCPO = "f"
  displayCRPS = "f"
  display = c(displayRow, displayMSE, displayVar, displayBias, displayCPO, displayCRPS)
  colScaleMSE = rep(10^scalingPower[1], 3)
  colScaleVar = rep(10^scalingPower[2], 3)
  colScaleBias = rep(10^scalingPower[3], 3)
  colScaleCPO = 1
  colScaleCRPS = 1
  colScale = c(colScaleMSE, colScaleVar, colScaleBias, colScaleCPO, colScaleCRPS)
  colUnitsMSE = rep(paste0(" ($\\times 10^{-", as.character(scalingPower[1]), "}$)"), 3)
  colUnitsVar = rep(paste0(" ($\\times 10^{-", as.character(scalingPower[2]), "}$)"), 3)
  colUnitsBias = rep(paste0(" ($\\times 10^{-", as.character(scalingPower[3]), "}$)"), 3)
  colUnitsCPO = ""
  colUnitsCRPS = ""
  colUnits = c(colUnitsMSE, colUnitsVar, colUnitsBias, colUnitsCPO, colUnitsCRPS)
  digitsMSE = rep(1, 3)
  digitsVar = rep(1, 3)
  digitsBias = rep(1, 3)
  digitsCPO = 2
  digitsCRPS = 2
  if(resultNameRoot == "Mort") {
    digitsCPO = 3
    digitsCRPS = 3
  }
  colDigits = c(digitsMSE, digitsVar, digitsBias, digitsCPO, digitsCRPS)
  
  tab = data.frame(validationResults$scoresLeaveOutCounty)
  for(i in 1:ncol(tab)) {
    tab[,i] = as.numeric(round(unlist(tab[,i]) * colScale[i], digits=colDigits[i]))
    colnames(tab)[i] = paste0(colnames(tab)[i], colUnits[i])
  }
  colnames(tab) = gsub(".", " ", colnames(tab), fixed=TRUE)
  colnames(tab) = gsub("Urban", "Urb", colnames(tab), fixed=TRUE)
  colnames(tab) = gsub("Rural", "Rur", colnames(tab), fixed=TRUE)
  # print(xtable(tab, digits=c(1,colDigits), display=display), 
  #       include.colnames=TRUE,
  #       hline.after=0, 
  #       math.style.exponents=TRUE, 
  #       sanitize.text.function=function(x){x}, 
  #       scalebox=0.5)
  # print(xtable(tab, digits=c(1,colDigits), display=display), 
  #       include.colnames=TRUE,
  #       hline.after=0, 
  #       math.style.exponents=TRUE, 
  #       sanitize.text.function=function(x){x}, 
  #       scalebox=0.5)
  tab = t(tab)
  colnames(tab) = gsub("SPDE ", "", colnames(tab), fixed=TRUE)
  colnames(tab) = gsub("BYM2 ", "", colnames(tab), fixed=TRUE)
  colnames(tab) = gsub("Smoothed Direct", " ", colnames(tab), fixed=TRUE)
  # temp = xtable2kable(xtable(tab), 
  #                     include.colnames=TRUE,
  #                     hline.after=0, 
  #                     math.style.exponents=TRUE, 
  #                     sanitize.text.function=function(x){x})
  # print(add_header_above(kable_styling(temp), c(" " = 1, "Smoothed Direct" = 1, "BYM2"=4, "SPDE"=4)))
  
  # remove the smoothed direct results
  tab = tab[,-1]
  
  # bold the best entries of each row, italicize worst entries of each row
  centers = c(rep(0, nrow(tab)))
  rowWorst = apply(abs(sweep(tab, 1, centers, "-")), 1, max, na.rm=TRUE)
  rowBest = apply(abs(sweep(tab, 1, centers, "-")), 1, min, na.rm=TRUE)
  temp = rowWorst[10]
  rowWorst[10] = rowBest[10] # highest CPO is best, not the lowest
  rowBest[10] = temp
  boldFun = function(i) {
    vals = tab[i,]
    isNA = is.na(vals)
    out = rep(FALSE, length(vals))
    out[!isNA] = abs(tab[i,!isNA] - centers[i]) == rowBest[i]
    out
  }
  italicFun = function(i) {
    vals = tab[i,]
    isNA = is.na(vals)
    out = rep(FALSE, length(vals))
    out[!isNA] = abs(tab[i,!isNA] - centers[i]) == rowWorst[i]
    out
  }
  if(resultNameRoot == "Ed") {
    # bold the best model for each scoring rule, italicize the worst
    test = t(sapply(1:nrow(tab), function(i) {cell_spec(tab[i,], "latex", bold=boldFun(i), italic=italicFun(i), 
                                                        monospace=FALSE, underline=FALSE, strikeout=FALSE)}))
  }
  else {
    # don't bother bolding or italicizing, since all models perform essentially equally
    test = tab
  }
  
  # revert the column names to their true values
  colnames(test) = colnames(tab)
  scoreVariations = c("Avg", "Urban", "Rural")
  scoreVariations = c(rep(scoreVariations, 3), "Avg", "Avg")
  test = cbind(" "=scoreVariations, test)
  rownames(test)=NULL
  
  # group the rows by urbanicity
  fullTab = test %>%
    kable("latex", escape = F, booktabs = T, format.args=list(drop0trailing=FALSE, scientific=FALSE), 
          align=c("l", rep("r", ncol(test) - 1))) %>% kable_styling()
  uniqueScoreTypes = c("MSE", "Var", "Bias", "CPO", "CRPS")
  uniqueColUnits = colUnits[c(1, 4, 7, 10:11)]
  for(i in 1:length(uniqueScoreTypes)) {
    uniqueScoreTypes[i] = paste0(uniqueScoreTypes[i], uniqueColUnits[i])
  }
  scoreTypeGroups = c(rep(uniqueScoreTypes[1:3], each=3), uniqueScoreTypes[4:5])
  
  for(i in 1:length(uniqueScoreTypes)) {
    startR = match(TRUE, scoreTypeGroups == uniqueScoreTypes[i])
    endR = nrow(tab) - match(TRUE, rev(scoreTypeGroups == uniqueScoreTypes[i])) + 1
    fullTab = fullTab %>% pack_rows(uniqueScoreTypes[i], startR, endR, latex_gap_space = "0.3em", escape=FALSE)
  }
  print(add_header_above(fullTab, c(" " = 1, "BYM2"=4, "SPDE"=4), italic=TRUE, bold=TRUE, escape=FALSE))
  
  print("")
  print("Leave out cluster results:")
  colScale = colScaleCPO
  colDigits = digitsCPO
  colUnits = colUnitsCPO
  display = c(displayRow, displayCPO)
  tab = data.frame(validationResults$scoresLeaveOutCluster)
  for(i in 1:ncol(tab)) {
    tab[,i] = as.numeric(round(unlist(tab[,i]) * colScale[i], digits=colDigits[i]))
    colnames(tab)[i] = paste0(colnames(tab)[i], colUnits[i])
  }
  # print(xtable(tab, digits=c(1,colDigits), display=display), 
  #       include.colnames=TRUE,
  #       hline.after=0, 
  #       math.style.exponents=TRUE, 
  #       sanitize.text.function=function(x){x}, 
  #       scalebox=0.5)
  
  tab = t(tab)
  
  # bold the best entries of each row, italicize worst entries of each row
  centers = c(rep(0, nrow(tab)))
  rowWorst = apply(abs(sweep(tab, 1, centers, "-")), 1, min, na.rm=TRUE)
  rowBest = apply(abs(sweep(tab, 1, centers, "-")), 1, max, na.rm=TRUE)
  boldFun = function(i) {
    vals = tab[i,]
    isNA = is.na(vals)
    out = rep(FALSE, length(vals))
    out[!isNA] = abs(tab[i,!isNA] - centers[i]) == rowBest[i]
    out
  }
  italicFun = function(i) {
    vals = tab[i,]
    isNA = is.na(vals)
    out = rep(FALSE, length(vals))
    out[!isNA] = abs(tab[i,!isNA] - centers[i]) == rowWorst[i]
    out
  }
  test = t(sapply(1:nrow(tab), function(i) {cell_spec(tab[i,], "latex", bold=boldFun(i), italic=italicFun(i), 
                                                      monospace=FALSE, underline=FALSE, strikeout=FALSE)}))
  
  # revert the column names to their true values
  colnames(test) = word(colnames(tab), 2)
  rownames(test)= "CPO"
  
  # construct the kable version of the table
  fullTab = test %>%
    kable("latex", escape = F, booktabs = T, format.args=list(drop0trailing=FALSE, scientific=FALSE), 
          align=c("l", rep("r", ncol(test) - 1))) %>% kable_styling()
  
  # group the rows by urbanicity
  print(add_header_above(fullTab, c(" "=1, "BYM2"=4, "SPDE"=4), italic=TRUE, bold=TRUE, escape=FALSE))
}

plotValidationResults = function(dat=NULL, targetPop=c("women", "children"), leaveOutRegion=FALSE, 
                                 startI=0, loadPreviousFit=TRUE, verbose=TRUE, endI=Inf, loadPreviousResults=FALSE, 
                                 family = c("betabinomial", "binomial")) {
  targetPop = match.arg(targetPop)
  family = match.arg(family)
  
  familyText=""
  if(family == "binomial")
    familyText = "_LgtN"
  
  # load in relevant data for the given example
  if(targetPop == "women") {
    resultNameRoot="Ed"
    if(is.null(dat)) {
      out = load("../U5MR/kenyaDataEd.RData")
      dat = ed
    }
    load("../U5MR/popGridAdjustedWomen.RData")
  } else if(targetPop == "children") {
    resultNameRoot="Mort"
    if(is.null(dat)) {
      out = load("../U5MR/kenyaData.RData")
      dat = mort
    }
    load("../U5MR/popGridAdjusted.RData")
  }
  resultNameRootLower = tolower(resultNameRoot)
  dataType = resultNameRootLower
  
  # get region names
  regions = sort(unique(countyToRegion(as.character(dat$admin1))))
  
  # Load all scoring rule tables
  fileName = paste0("savedOutput/validation/validationResults", resultNameRoot, "_LORegion", leaveOutRegion, ".RData")
  # allScores = list(scoresInSample=scoresInSample, scoresLeaveOneOut=scoresLeaveOneOut, scoresLeaveOutRegion=scoresLeaveOutRegion, 
  #                  scoresInSampleBinomial=scoresInSampleBinomial, scoresLeaveOneOutBinomial=scoresLeaveOneOutBinomial, scoresLeaveOutRegionBinomial=scoresLeaveOutRegionBinomial, 
  #                  
  #                  binnedScoringRulesuuAll=binnedScoringRulesuuAll, 
  #                  binnedScoringRulesuUAll=binnedScoringRulesuUAll, 
  #                  binnedScoringRulesUuAll=binnedScoringRulesUuAll, 
  #                  binnedScoringRulesUUAll=binnedScoringRulesUUAll, 
  #                  binnedScoringRulesAuAll=binnedScoringRulesAuAll, 
  #                  binnedScoringRulesAUAll=binnedScoringRulesAUAll, 
  #                  binnedScoringRulesuAAll=binnedScoringRulesuAAll, 
  #                  binnedScoringRulesUAAll=binnedScoringRulesUAAll, 
  #                  binnedScoringRulesAAAll=binnedScoringRulesAAAll, 
  #                  binnedScoringRulesuuBinomialAll=binnedScoringRulesuuBinomialAll, 
  #                  binnedScoringRulesuUBinomialAll=binnedScoringRulesuUBinomialAll, 
  #                  binnedScoringRulesUuBinomialAll=binnedScoringRulesUuBinomialAll, 
  #                  binnedScoringRulesUUBinomialAll=binnedScoringRulesUUBinomialAll, 
  #                  binnedScoringRulesAuBinomialAll=binnedScoringRulesAuBinomialAll, 
  #                  binnedScoringRulesAUBinomialAll=binnedScoringRulesAUBinomialAll, 
  #                  binnedScoringRulesuABinomialAll=binnedScoringRulesuABinomialAll, 
  #                  binnedScoringRulesUABinomialAll=binnedScoringRulesUABinomialAll, 
  #                  binnedScoringRulesAABinomialAll=binnedScoringRulesAABinomialAll, 
  #                  
  #                  singleScoresAll=singleScoresAll, 
  #                  singleScoresBinomialAll=singleScoresBinomialAll)
  # save(allScores, file=fileName)
  out = load(fileName)
  # modelNames = names(allScores$singleScoresBinomialAll)
  modelNames = c(expression("SPDE"[u]), expression("SPDE"[U]), 
                 expression("LK-INLA"[ui]), expression("LK-INLA"[Ui]), 
                 expression("LK-INLA"[uI]), expression("LK-INLA"[UI]))
  
  ##### Start with single observation scoring rules
  # construct information for separating out prediction types
  sortI = sort(allScores$singleScoresBinomialAll$spdeu$dataI, index.return=TRUE)$ix
  
  ## pair plots
  # plot prediction standard deviations
  my_line <- function(x,y,...){
    # if(diff(range(x)) >= .04)
    xlim = zlim
    # else
    #   xlim = zlim2
    # if(diff(range(y)) >= .04)
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
    points(x[!dat$urban],y[!dat$urban],..., col="green", pch=19, cex=.1)
    points(x[dat$urban],y[dat$urban],..., col="blue", pch=19, cex=.1)
  }
  
  values = data.frame(allScores$singleScoresBinomialAll$spdeu$Width[sortI], allScores$singleScoresBinomialAll$lkinlaus$Width[sortI], allScores$singleScoresBinomialAll$lkinlauS$Width[sortI])
  zlim = range(c(as.matrix(values)))
  lims = rep(list(zlim), length(values))
  myPairs(values, 
          labels=modelNames, 
          lower.panel=my_line, upper.panel = my_line, 
          main=paste0("80% CI Width"), log="xy",
          lims=lims, oma=c(3,3,6,7))
  
  my_line <- function(x,y,...){
    # if(diff(range(x)) >= .04)
    xlim = zlim
    # else
    #   xlim = zlim2
    # if(diff(range(y)) >= .04)
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
    points(x[!dat$urban],y[!dat$urban],..., col="green", pch=19, cex=.1)
    points(x[dat$urban],y[dat$urban],..., col="blue", pch=19, cex=.1)
  }
  
  values = data.frame(allScores$singleScoresAll$spdeu$Width[sortI], allScores$singleScoresAll$lkinlaus$Width[sortI], allScores$singleScoresAll$lkinlauS$Width[sortI])
  zlim = range(c(as.matrix(values)))
  lims = rep(list(zlim), length(values))
  myPairs(values, 
          labels=modelNames, 
          lower.panel=my_line, upper.panel = my_line, 
          main=paste0("80% CI Width"), 
          lims=lims, oma=c(3,3,6,7))
  
  ## plot scores versus nearest neighbor distance
  plot(allScores$singleScoresBinomialAll$spdeu$NNDist[sortI][!dat$urban], 
       allScores$singleScoresBinomialAll$spdeu$RMSE[sortI][!dat$urban], 
       pch=19, cex=.1, col="green", 
       main="RMSE vs. NN Distance", xlab="Distance (km)", ylab="RMSE", ylim=zlim)
  points(allScores$singleScoresBinomialAll$spdeu$NNDist[sortI][dat$urban], 
         allScores$singleScoresBinomialAll$spdeu$RMSE[sortI][dat$urban], 
         pch=19, cex=.1, col="blue")
  
  plot(allScores$singleScoresBinomialAll$lkinlaus$NNDist[sortI][!dat$urban], 
         allScores$singleScoresBinomialAll$lkinlaus$RMSE[sortI][!dat$urban], 
         pch=19, cex=.1, col="green", 
       main="RMSE vs. NN Distance", xlab="Distance (km)", ylab="RMSE", ylim=zlim)
  points(allScores$singleScoresBinomialAll$lkinlaus$NNDist[sortI][dat$urban], 
         allScores$singleScoresBinomialAll$lkinlaus$RMSE[sortI][dat$urban], 
         pch=19, cex=.1, col="blue")
  
  ## plot scores spatially
  out = load("../U5MR/adminMapData.RData")
  kenyaMap = adm0
  countyMap = adm1
  latRange=c(-4.6, 5)
  lonRange=c(33.5, 42.0)
  coords = cbind(dat$lon, dat$lat)
  popCoords = cbind(popGrid$lon, popGrid$lat)
  popUrban = popGrid$urban
  
  values = data.frame(allScores$singleScoresBinomialAll$spdeu$RMSE[sortI], 
                      allScores$singleScoresBinomialAll$lkinlaus$RMSE[sortI], 
                      allScores$singleScoresBinomialAll$lkinlauS$RMSE[sortI])
  zlim = range(c(as.matrix(values)))
  modelNames = c(expression("SPDE"[u]), expression("LK-INLA"[ui]), expression("LK-INLA"[uI]))
  
  width = 1000
  png(file=paste0("Figures/", resultNameRoot, "/spatialRMSE_LORegion", leaveOutRegion, ".png"), width=width, height=1000)
  par(mfrow=c(2,2), oma=c( 0,0,0,2), mar=c(5.1, 4.1, 4.1, 6))
  for(i in 1:ncol(values)) {
    quilt.plot(coords, values[,i], xlim=lonRange, ylim=latRange, main=bquote(.(modelNames[[i]]) ~ " RMSE"), 
               nx=150, ny=150, col=makeRedGrayBlueDivergingColors(n=29, valRange=zlim, center=0), zlim=zlim)
    plotMapDat(mapDat=countyMap, lwd=.5)
  }
  dev.off()
  
  values = data.frame(allScores$singleScoresBinomialAll$spdeu$Bias[sortI], 
                      allScores$singleScoresBinomialAll$lkinlaus$Bias[sortI], 
                      allScores$singleScoresBinomialAll$lkinlauS$Bias[sortI])
  zlim = range(c(as.matrix(values)))
  modelNames = c(expression("SPDE"[u]), expression("LK-INLA"[ui]), expression("LK-INLA"[uI]))
  
  png(file=paste0("Figures/", resultNameRoot, "/spatialBias_LORegion", leaveOutRegion, ".png"), width=width, height=1000)
  par(mfrow=c(2,2), oma=c( 0,0,0,2), mar=c(5.1, 4.1, 4.1, 6))
  for(i in 1:ncol(values)) {
    quilt.plot(coords, values[,i], xlim=lonRange, ylim=latRange, main=bquote(.(modelNames[[i]]) ~ " Bias"), 
               nx=150, ny=150, col=makeRedGrayBlueDivergingColors(n=29, valRange=zlim, center=0), zlim=zlim)
    plotMapDat(mapDat=countyMap, lwd=.5)
  }
  dev.off()
  
  values = data.frame(allScores$singleScoresBinomialAll$spdeu$Width[sortI], 
                      allScores$singleScoresBinomialAll$lkinlaus$Width[sortI], 
                      allScores$singleScoresBinomialAll$lkinlauS$Width[sortI])
  zlim = range(c(as.matrix(values)))
  modelNames = c(expression("SPDE"[u]), expression("LK-INLA"[ui]), expression("LK-INLA"[uI]))
  
  png(file=paste0("Figures/", resultNameRoot, "/spatialWidthBinomial_LORegion", leaveOutRegion, ".png"), width=width, height=1000)
  par(mfrow=c(2,2), oma=c( 0,0,0,2), mar=c(5.1, 4.1, 4.1, 6))
  for(i in 1:ncol(values)) {
    quilt.plot(coords, values[,i], xlim=lonRange, ylim=latRange, main=bquote(.(modelNames[[i]]) ~ " Width"), 
               nx=120, ny=120, col=makeBlueYellowSequentialColors(n=64), zlim=zlim)
    plotMapDat(mapDat=countyMap, lwd=.5)
  }
  dev.off()
  
  values = data.frame(allScores$singleScoresAll$spdeu$Width[sortI], 
                      allScores$singleScoresAll$lkinlaus$Width[sortI], 
                      allScores$singleScoresAll$lkinlauS$Width[sortI])
  zlim = range(c(as.matrix(values)))
  modelNames = c(expression("SPDE"[u]), expression("LK-INLA"[ui]), expression("LK-INLA"[uI]))
  
  png(file=paste0("Figures/", resultNameRoot, "/spatialWidth_LORegion", leaveOutRegion, ".png"), width=width, height=1000)
  par(mfrow=c(2,2), oma=c( 0,0,0,2), mar=c(5.1, 4.1, 4.1, 6))
  quilt.plot(popCoords, popUrban, xlim=lonRange, ylim=latRange, main="Urbanicity", 
             nx=170, ny=170, col=c(rgb(0,0,0,0), "blue"), add.legend=FALSE)
  plotMapDat(mapDat=countyMap, lwd=.5)
  for(i in 1:ncol(values)) {
    quilt.plot(coords, values[,i], xlim=lonRange, ylim=latRange, main=bquote(.(modelNames[[i]]) ~ " Width"), 
               nx=120, ny=120, col=makeBlueYellowSequentialColors(n=64), zlim=zlim)
    plotMapDat(mapDat=countyMap, lwd=.5)
  }
  dev.off()
  
  png(file=paste0("Figures/", resultNameRoot, "/womenPerCluster_LORegion", leaveOutRegion, ".png"), width=500, height=500)
  quilt.plot(coords, dat$n, xlim=lonRange, ylim=latRange, main="Women Per Cluster", 
             nx=120, ny=120, col=makeBlueYellowSequentialColors(n=64))
  plotMapDat(mapDat=countyMap, lwd=.5)
  dev.off()
  
  png(file=paste0("Figures/", resultNameRoot, "/womenPerClusterLog_LORegion", leaveOutRegion, ".png"), width=500, height=500)
  quilt.plot(coords, dat$n, xlim=lonRange, ylim=latRange, main="Women Per Cluster", 
             nx=120, ny=120, col=makeBlueYellowSequentialColors(n=64), FUN=function(x){mean(log10(x))})
  plotMapDat(mapDat=countyMap, lwd=.5)
  dev.off()
  
  png(file=paste0("Figures/", resultNameRoot, "/womenPerCell_LORegion", leaveOutRegion, ".png"), width=500, height=500)
  quilt.plot(coords, dat$n, xlim=lonRange, ylim=latRange, main="Women Per Grid Cell", 
             nx=120, ny=120, col=makeBlueYellowSequentialColors(n=64), FUN = function(x) {log10(sum(x)+1)})
  plotMapDat(mapDat=countyMap, lwd=.5)
  dev.off()
  
  png(file=paste0("Figures/", resultNameRoot, "/womenPerCellLog_LORegion", leaveOutRegion, ".png"), width=500, height=500)
  quilt.plot(coords, dat$n, xlim=lonRange, ylim=latRange, main="Women Per Grid Cell", 
             nx=120, ny=120, col=makeBlueYellowSequentialColors(n=64), FUN = function(x) {log10(sum(x)+1)})
  plotMapDat(mapDat=countyMap, lwd=.5)
  dev.off()
  
  png(file=paste0("Figures/", resultNameRoot, "/Variance_LORegion", leaveOutRegion, ".png"), width=500, height=500)
  quilt.plot(coords, dat$y/dat$n, xlim=lonRange, ylim=latRange, main="Variance", 
             nx=50, ny=50, col=makeBlueYellowSequentialColors(n=64), FUN = function(x) {var(x)})
  plotMapDat(mapDat=countyMap, lwd=.5)
  dev.off()
  
  ##### Binned scoring rules versus distance
  ## RMSE
  # AA
  values = data.frame(allScores$binnedScoringRulesAABinomialAll$spdeu$RMSE, allScores$binnedScoringRulesAABinomialAll$lkinlaus$RMSE, allScores$binnedScoringRulesAABinomialAll$lkinlauS$RMSE)
  ns = allScores$binnedScoringRulesAABinomialAll$spdeu$nPerBin
  sizes = sqrt(ns) / 5
  zlim = range(c(as.matrix(values)))
  cols = rainbow(3)
  d = allScores$binnedScoringRulesAABinomialAll$spdeu$NNDist
  
  for(i in 1:ncol(values)) {
    if(i == 1) {
      plot(d, 
           values[,i], 
           cex = 1/sizes, col=cols[i], 
           main="RMSE vs. Distance From Observation", xlab="Distance (km)", ylab="RMSE", ylim=zlim)
    } else
      points(d, values[,i], col=cols[i], cex = 1/sizes)
    
    lines(d, values[,i], col=cols[i])
  }
  legend("topright", c(expression("SPDE"[u]), expression("LK-INLA"[ui]), expression("LK-INLA"[uI])), col=cols, lty=1, pch=1)
  
  # Au
  values = data.frame(allScores$binnedScoringRulesAuBinomialAll$spdeu$RMSE, allScores$binnedScoringRulesAuBinomialAll$lkinlaus$RMSE, allScores$binnedScoringRulesAuBinomialAll$lkinlauS$RMSE)
  ns = allScores$binnedScoringRulesAuBinomialAll$spdeu$nPerBin
  sizes = sqrt(ns) / 5
  zlim = range(c(as.matrix(values)))
  cols = rainbow(3)
  d = allScores$binnedScoringRulesAuBinomialAll$spdeu$NNDist
  
  for(i in 1:ncol(values)) {
    if(i == 1) {
      plot(d, 
           values[,i], 
           cex = 1/sizes, col=cols[i], 
           main="RMSE vs. Distance From Observation to Rural Pt.", xlab="Distance (km)", ylab="RMSE", ylim=zlim)
    } else
      points(d, values[,i], col=cols[i], cex = 1/sizes)
    
    lines(d, values[,i], col=cols[i])
  }
  legend("topright", c(expression("SPDE"[u]), expression("LK-INLA"[ui]), expression("LK-INLA"[uI])), col=cols, lty=1, pch=1)
  
  # Uu
  values = data.frame(allScores$binnedScoringRulesUuBinomialAll$spdeu$RMSE, allScores$binnedScoringRulesUuBinomialAll$lkinlaus$RMSE, allScores$binnedScoringRulesUuBinomialAll$lkinlauS$RMSE)
  ns = allScores$binnedScoringRulesUuBinomialAll$spdeu$nPerBin
  sizes = sqrt(ns) / 5
  zlim = range(c(as.matrix(values)))
  cols = rainbow(3)
  d = allScores$binnedScoringRulesUuBinomialAll$spdeu$NNDist
  
  for(i in 1:ncol(values)) {
    if(i == 1) {
      plot(d, 
           values[,i], 
           cex = 1/sizes, col=cols[i], 
           main="RMSE vs. Distance From Urban to Rural Pt.", xlab="Distance (km)", ylab="RMSE", ylim=zlim)
    } else
      points(d, values[,i], col=cols[i], cex = 1/sizes)
    
    lines(d, values[,i], col=cols[i])
  }
  legend("topright", c(expression("SPDE"[u]), expression("LK-INLA"[ui]), expression("LK-INLA"[uI])), col=cols, lty=1, pch=1)
  
  # uu
  values = data.frame(allScores$binnedScoringRulesuuBinomialAll$spdeu$RMSE, allScores$binnedScoringRulesuuBinomialAll$lkinlaus$RMSE, allScores$binnedScoringRulesuuBinomialAll$lkinlauS$RMSE)
  ns = allScores$binnedScoringRulesuuBinomialAll$spdeu$nPerBin
  sizes = sqrt(ns) / 5
  zlim = range(c(as.matrix(values)))
  cols = rainbow(3)
  d = allScores$binnedScoringRulesuuBinomialAll$spdeu$NNDist
  
  for(i in 1:ncol(values)) {
    if(i == 1) {
      plot(d, 
           values[,i], 
           cex = 1/sizes, col=cols[i], 
           main="RMSE vs. Distance From Rural to Rural", xlab="Distance (km)", ylab="RMSE", ylim=zlim)
    } else
      points(d, values[,i], col=cols[i], cex = 1/sizes)
    
    lines(d, values[,i], col=cols[i])
  }
  legend("topright", c(expression("SPDE"[u]), expression("LK-INLA"[ui]), expression("LK-INLA"[uI])), col=cols, lty=1, pch=1)
  
  ## CRPS
  # AA
  values = data.frame(allScores$binnedScoringRulesAABinomialAll$spdeu$CRPS, allScores$binnedScoringRulesAABinomialAll$lkinlaus$CRPS, allScores$binnedScoringRulesAABinomialAll$lkinlauS$CRPS)
  ns = allScores$binnedScoringRulesAABinomialAll$spdeu$nPerBin
  sizes = sqrt(ns) / 5
  zlim = range(c(as.matrix(values)))
  cols = rainbow(3)
  d = allScores$binnedScoringRulesAABinomialAll$spdeu$NNDist
  
  for(i in 1:ncol(values)) {
    if(i == 1) {
      plot(d, 
           values[,i], 
           cex = 1/sizes, col=cols[i], 
           main="CRPS vs. Distance From Observation", xlab="Distance (km)", ylab="CRPS", ylim=zlim)
    } else
      points(d, values[,i], col=cols[i], cex = 1/sizes)
    
    lines(d, values[,i], col=cols[i])
  }
  legend("topright", c(expression("SPDE"[u]), expression("LK-INLA"[ui]), expression("LK-INLA"[uI])), col=cols, lty=1, pch=1)
  
  # Au
  values = data.frame(allScores$binnedScoringRulesAuBinomialAll$spdeu$CRPS, allScores$binnedScoringRulesAuBinomialAll$lkinlaus$CRPS, allScores$binnedScoringRulesAuBinomialAll$lkinlauS$CRPS)
  ns = allScores$binnedScoringRulesAuBinomialAll$spdeu$nPerBin
  sizes = sqrt(ns) / 5
  zlim = range(c(as.matrix(values)))
  cols = rainbow(3)
  d = allScores$binnedScoringRulesAuBinomialAll$spdeu$NNDist
  
  for(i in 1:ncol(values)) {
    if(i == 1) {
      plot(d, 
           values[,i], 
           cex = 1/sizes, col=cols[i], 
           main="CRPS vs. Distance From Observation to Rural Pt.", xlab="Distance (km)", ylab="CRPS", ylim=zlim)
    } else
      points(d, values[,i], col=cols[i], cex = 1/sizes)
    
    lines(d, values[,i], col=cols[i])
  }
  legend("topright", c(expression("SPDE"[u]), expression("LK-INLA"[ui]), expression("LK-INLA"[uI])), col=cols, lty=1, pch=1)
  
  # Uu
  values = data.frame(allScores$binnedScoringRulesUuBinomialAll$spdeu$CRPS, allScores$binnedScoringRulesUuBinomialAll$lkinlaus$CRPS, allScores$binnedScoringRulesUuBinomialAll$lkinlauS$CRPS)
  ns = allScores$binnedScoringRulesUuBinomialAll$spdeu$nPerBin
  sizes = sqrt(ns) / 5
  zlim = range(c(as.matrix(values)))
  cols = rainbow(3)
  d = allScores$binnedScoringRulesUuBinomialAll$spdeu$NNDist
  
  for(i in 1:ncol(values)) {
    if(i == 1) {
      plot(d, 
           values[,i], 
           cex = 1/sizes, col=cols[i], 
           main="CRPS vs. Distance From Urban to Rural Pt.", xlab="Distance (km)", ylab="CRPS", ylim=zlim)
    } else
      points(d, values[,i], col=cols[i], cex = 1/sizes)
    
    lines(d, values[,i], col=cols[i])
  }
  legend("topright", c(expression("SPDE"[u]), expression("LK-INLA"[ui]), expression("LK-INLA"[uI])), col=cols, lty=1, pch=1)
  
  # uu
  values = data.frame(allScores$binnedScoringRulesuuBinomialAll$spdeu$CRPS, allScores$binnedScoringRulesuuBinomialAll$lkinlaus$CRPS, allScores$binnedScoringRulesuuBinomialAll$lkinlauS$CRPS)
  ns = allScores$binnedScoringRulesuuBinomialAll$spdeu$nPerBin
  sizes = sqrt(ns) / 5
  zlim = range(c(as.matrix(values)))
  cols = rainbow(3)
  d = allScores$binnedScoringRulesuuBinomialAll$spdeu$NNDist
  
  for(i in 1:ncol(values)) {
    if(i == 1) {
      plot(d, 
           values[,i], 
           cex = 1/sizes, col=cols[i], 
           main="CRPS vs. Distance From Rural to Rural", xlab="Distance (km)", ylab="CRPS", ylim=zlim)
    } else
      points(d, values[,i], col=cols[i], cex = 1/sizes)
    
    lines(d, values[,i], col=cols[i])
  }
  legend("topright", c(expression("SPDE"[u]), expression("LK-INLA"[ui]), expression("LK-INLA"[uI])), col=cols, lty=1, pch=1)
  
  ## Width
  # AA
  values = data.frame(allScores$binnedScoringRulesAABinomialAll$spdeu$Width, allScores$binnedScoringRulesAABinomialAll$lkinlaus$Width, allScores$binnedScoringRulesAABinomialAll$lkinlauS$Width)
  ns = allScores$binnedScoringRulesAABinomialAll$spdeu$nPerBin
  sizes = sqrt(ns) / 5
  zlim = range(c(as.matrix(values)))
  cols = rainbow(3)
  d = allScores$binnedScoringRulesAABinomialAll$spdeu$NNDist
  
  for(i in 1:ncol(values)) {
    if(i == 1) {
      plot(d, 
           values[,i], 
           cex = 1/sizes, col=cols[i], 
           main="Width vs. Distance From Observation", xlab="Distance (km)", ylab="Width", ylim=zlim)
    } else
      points(d, values[,i], col=cols[i], cex = 1/sizes)
    
    lines(d, values[,i], col=cols[i])
  }
  legend("topright", c(expression("SPDE"[u]), expression("LK-INLA"[ui]), expression("LK-INLA"[uI])), col=cols, lty=1, pch=1)
  
  # Au
  values = data.frame(allScores$binnedScoringRulesAuBinomialAll$spdeu$Width, allScores$binnedScoringRulesAuBinomialAll$lkinlaus$Width, allScores$binnedScoringRulesAuBinomialAll$lkinlauS$Width)
  ns = allScores$binnedScoringRulesAuBinomialAll$spdeu$nPerBin
  sizes = sqrt(ns) / 5
  zlim = range(c(as.matrix(values)))
  cols = rainbow(3)
  d = allScores$binnedScoringRulesAuBinomialAll$spdeu$NNDist
  
  for(i in 1:ncol(values)) {
    if(i == 1) {
      plot(d, 
           values[,i], 
           cex = 1/sizes, col=cols[i], 
           main="Width vs. Distance From Observation to Rural Pt.", xlab="Distance (km)", ylab="Width", ylim=zlim)
    } else
      points(d, values[,i], col=cols[i], cex = 1/sizes)
    
    lines(d, values[,i], col=cols[i])
  }
  legend("top", c(expression("SPDE"[u]), expression("LK-INLA"[ui]), expression("LK-INLA"[uI])), col=cols, lty=1, pch=1)
  
  # Uu
  values = data.frame(allScores$binnedScoringRulesUuBinomialAll$spdeu$Width, allScores$binnedScoringRulesUuBinomialAll$lkinlaus$Width, allScores$binnedScoringRulesUuBinomialAll$lkinlauS$Width)
  ns = allScores$binnedScoringRulesUuBinomialAll$spdeu$nPerBin
  sizes = sqrt(ns) / 5
  zlim = range(c(as.matrix(values)))
  cols = rainbow(3)
  d = allScores$binnedScoringRulesUuBinomialAll$spdeu$NNDist
  
  for(i in 1:ncol(values)) {
    if(i == 1) {
      plot(d, 
           values[,i], 
           cex = 1/sizes, col=cols[i], 
           main="Width vs. Distance From Urban to Rural Pt.", xlab="Distance (km)", ylab="Width", ylim=zlim)
    } else
      points(d, values[,i], col=cols[i], cex = 1/sizes)
    
    lines(d, values[,i], col=cols[i])
  }
  legend("topright", c(expression("SPDE"[u]), expression("LK-INLA"[ui]), expression("LK-INLA"[uI])), col=cols, lty=1, pch=1)
  
  # uu
  values = data.frame(allScores$binnedScoringRulesuuBinomialAll$spdeu$Width, allScores$binnedScoringRulesuuBinomialAll$lkinlaus$Width, allScores$binnedScoringRulesuuBinomialAll$lkinlauS$Width)
  ns = allScores$binnedScoringRulesuuBinomialAll$spdeu$nPerBin
  sizes = sqrt(ns) / 5
  zlim = range(c(as.matrix(values)))
  cols = rainbow(3)
  d = allScores$binnedScoringRulesuuBinomialAll$spdeu$NNDist
  
  for(i in 1:ncol(values)) {
    if(i == 1) {
      plot(d, 
           values[,i], 
           cex = 1/sizes, col=cols[i], 
           main="Width vs. Distance From Rural to Rural", xlab="Distance (km)", ylab="Width", ylim=zlim)
    } else
      points(d, values[,i], col=cols[i], cex = 1/sizes)
    
    lines(d, values[,i], col=cols[i])
  }
  legend("top", c(expression("SPDE"[u]), expression("LK-INLA"[ui]), expression("LK-INLA"[uI])), col=cols, lty=1, pch=1)
  
  ## Bias
  # AA
  values = data.frame(allScores$binnedScoringRulesAABinomialAll$spdeu$Bias, allScores$binnedScoringRulesAABinomialAll$lkinlaus$Bias, allScores$binnedScoringRulesAABinomialAll$lkinlauS$Bias)
  ns = allScores$binnedScoringRulesAABinomialAll$spdeu$nPerBin
  sizes = sqrt(ns) / 5
  zlim = range(c(as.matrix(values)))
  cols = rainbow(3)
  d = allScores$binnedScoringRulesAABinomialAll$spdeu$NNDist
  
  for(i in 1:ncol(values)) {
    if(i == 1) {
      plot(d, 
           values[,i], 
           cex = 1/sizes, col=cols[i], 
           main="Bias vs. Distance From Observation", xlab="Distance (km)", ylab="Bias", ylim=zlim)
    } else
      points(d, values[,i], col=cols[i], cex = 1/sizes)
    
    lines(d, values[,i], col=cols[i])
  }
  abline(h=0, lty=2)
  legend("bottomleft", c(expression("SPDE"[u]), expression("LK-INLA"[ui]), expression("LK-INLA"[uI])), col=cols, lty=1, pch=1)
  
  # Au
  values = data.frame(allScores$binnedScoringRulesAuBinomialAll$spdeu$Bias, allScores$binnedScoringRulesAuBinomialAll$lkinlaus$Bias, allScores$binnedScoringRulesAuBinomialAll$lkinlauS$Bias)
  ns = allScores$binnedScoringRulesAuBinomialAll$spdeu$nPerBin
  sizes = sqrt(ns) / 5
  zlim = range(c(as.matrix(values)))
  cols = rainbow(3)
  d = allScores$binnedScoringRulesAuBinomialAll$spdeu$NNDist
  
  for(i in 1:ncol(values)) {
    if(i == 1) {
      plot(d, 
           values[,i], 
           cex = 1/sizes, col=cols[i], 
           main="Bias vs. Distance From Observation to Rural Pt.", xlab="Distance (km)", ylab="Bias", ylim=zlim)
    } else
      points(d, values[,i], col=cols[i], cex = 1/sizes)
    
    lines(d, values[,i], col=cols[i])
  }
  abline(h=0, lty=2)
  legend("topleft", c(expression("SPDE"[u]), expression("LK-INLA"[ui]), expression("LK-INLA"[uI])), col=cols, lty=1, pch=1)
  
  # Uu
  values = data.frame(allScores$binnedScoringRulesUuBinomialAll$spdeu$Bias, allScores$binnedScoringRulesUuBinomialAll$lkinlaus$Bias, allScores$binnedScoringRulesUuBinomialAll$lkinlauS$Bias)
  ns = allScores$binnedScoringRulesUuBinomialAll$spdeu$nPerBin
  sizes = sqrt(ns) / 5
  zlim = range(c(as.matrix(values)))
  cols = rainbow(3)
  d = allScores$binnedScoringRulesUuBinomialAll$spdeu$NNDist
  
  for(i in 1:ncol(values)) {
    if(i == 1) {
      plot(d, 
           values[,i], 
           cex = 1/sizes, col=cols[i], 
           main="Bias vs. Distance From Urban to Rural Pt.", xlab="Distance (km)", ylab="Bias", ylim=zlim)
    } else
      points(d, values[,i], col=cols[i], cex = 1/sizes)
    
    lines(d, values[,i], col=cols[i])
  }
  abline(h=0, lty=2)
  legend("topleft", c(expression("SPDE"[u]), expression("LK-INLA"[ui]), expression("LK-INLA"[uI])), col=cols, lty=1, pch=1)
  
  # uu
  values = data.frame(allScores$binnedScoringRulesuuBinomialAll$spdeu$Bias, allScores$binnedScoringRulesuuBinomialAll$lkinlaus$Bias, allScores$binnedScoringRulesuuBinomialAll$lkinlauS$Bias)
  ns = allScores$binnedScoringRulesuuBinomialAll$spdeu$nPerBin
  sizes = sqrt(ns) / 5
  zlim = range(c(as.matrix(values)))
  cols = rainbow(3)
  d = allScores$binnedScoringRulesuuBinomialAll$spdeu$NNDist
  
  for(i in 1:ncol(values)) {
    if(i == 1) {
      plot(d, 
           values[,i], 
           cex = 1/sizes, col=cols[i], 
           main="Bias vs. Distance From Rural to Rural", xlab="Distance (km)", ylab="Bias", ylim=zlim)
    } else
      points(d, values[,i], col=cols[i], cex = 1/sizes)
    
    lines(d, values[,i], col=cols[i])
  }
  abline(h=0, lty=2)
  legend("bottomleft", c(expression("SPDE"[u]), expression("LK-INLA"[ui]), expression("LK-INLA"[uI])), col=cols, lty=1, pch=1)
  
  
  
  
  
  values = data.frame(allScores$binnedScoringRulesuuBinomialAll$spdeu$CRPS, allScores$binnedScoringRulesuuBinomialAll$lkinlaus$CRPS, allScores$binnedScoringRulesuuBinomialAll$lkinlauS$CRPS)
  ns = allScores$binnedScoringRulesuuBinomialAll$spdeu$nPerBin
  sizes = sqrt(ns) / 5
  zlim = range(c(as.matrix(values)))
  cols = rainbow(3)
  d = allScores$binnedScoringRulesuuBinomialAll$spdeu$NNDist
  
  for(i in 1:ncol(values)) {
    if(i == 1) {
      plot(d, 
           values[,i], 
           cex = 1/sizes, col=cols[i], 
           main="CRPS vs. Distance From Rural to Rural", xlab="Distance (km)", ylab="CRPS", ylim=zlim)
    } else
      points(d, values[,i], col=cols[i], cex = 1/sizes)
    
    lines(d, values[,i], col=cols[i])
  }
  legend("topright", c(expression("SPDE"[u]), expression("LK-INLA"[ui]), expression("LK-INLA"[uI])), col=cols, lty=1, pch=1)
  
  values = data.frame(allScores$binnedScoringRulesUuBinomialAll$spdeu$Bias, allScores$binnedScoringRulesUuBinomialAll$lkinlaus$Bias, allScores$binnedScoringRulesUuBinomialAll$lkinlauS$Bias)
  ns = allScores$binnedScoringRulesUuBinomialAll$spdeu$nPerBin
  sizes = sqrt(ns) / 5
  zlim = range(c(as.matrix(values)))
  cols = rainbow(3)
  d = allScores$binnedScoringRulesUuBinomialAll$spdeu$NNDist
  
  for(i in 1:ncol(values)) {
    if(i == 1) {
      plot(d, 
           values[,i], 
           cex = 1/sizes, col=cols[i], 
           main="Bias vs. Distance From Urban", xlab="Distance (km)", ylab="Bias", ylim=zlim)
    } else
      points(d, values[,i], col=cols[i], cex = 1/sizes)
    
    lines(d, values[,i], col=cols[i])
  }
  abline(h=0, lty= 2)
  legend("topleft", c(expression("SPDE"[u]), expression("LK-INLA"[ui]), expression("LK-INLA"[uI])), col=cols, lty=1, pch=1)
  
}

plotValidationSamples = function(dat=NULL, targetPop=c("women", "children")) {
  targetPop = match.arg(targetPop)
  
  # load in relevant data for the given example
  if(targetPop == "women") {
    resultNameRoot="Ed"
    if(is.null(dat)) {
      out = load("../U5MR/kenyaDataEd.RData")
      dat = ed
    }
  } else if(targetPop == "children") {
    resultNameRoot="Mort"
    if(is.null(dat)) {
      out = load("../U5MR/kenyaData.RData")
      dat = mort
    }
  }
  dataType = tolower(resultNameRoot)
  
  sampleTable = getValidationI(dat=dat, dataType=dataType)
  out = load("../U5MR/adminMapData.RData")
  kenyaMap = adm0
  countyMap = adm1
  latRange=c(-4.6, 5)
  lonRange=c(33.5, 42.0)
  coords = cbind(dat$lon, dat$lat)
  
  for(i in 1:ncol(sampleTable)) {
    leftOutI = sampleTable[,i]
    plot(coords[!leftOutI,], pch=19, cex=.4, main=paste0("Fold ", i), xlim=lonRange, ylim=latRange)
    points(coords[leftOutI,], pch=19, cex=.4, col="red")
  }
}











