library(kableExtra)

# validate the smoothing models by leaving out data from one county at a time
validateExample = function(dat=NULL, targetPop=c("women", "children"), 
                           startI=0, loadPreviousFit=TRUE, verbose=TRUE, endI=Inf, loadPreviousResults=FALSE) {
  targetPop = match.arg(targetPop)
  
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
  
  ##### Do validation
  resultsListSPDE = list()
  resultsListLKINLA = list()
  
  ##### run SPDE
  argList = list(list(dat = dat, urbanEffect = FALSE), 
                 list(dat = dat, urbanEffect = TRUE), 
                 list(dat = dat, urbanEffect = FALSE), 
                 list(dat = dat, urbanEffect = TRUE))
  otherArguments = list(dataType=dataType, verbose=verbose, loadPreviousFit=loadPreviousFit, family="betabinomial", 
                        loadPreviousResults=loadPreviousResults)
  
  modelNames = c()
  for(i in 1:length(argList)) {
    args = argList[[i]]
    separateRanges = args$separateRanges
    urbanEffect = args$urbanEffect
    fileName = paste0("savedOutput/validation/resultsSPDE", resultNameRootLower, "_urbanEffect", urbanEffect)
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
        load(paste0(fileName, ".RData"))
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
  otherArguments = list(dataType=dataType, loadPreviousFit=loadPreviousFit, family="betabinomial", 
                        loadPreviousResults=loadPreviousResults)
  
  for(i in 1:length(argList)) {
    args = argList[[i]]
    separateRanges = args$separateRanges
    urbanEffect = args$urbanEffect
    fileName = paste0("savedOutput/validation/resultsLKINLA", resultNameRootLower, "_urbanEffect", urbanEffect, 
                      "_separateRanges", separateRanges)
    if(separateRanges)
      separateText = "S"
    else
      separateText = "s"
    if(urbanEffect)
      urbanText = "U"
    else
      urbanText = "u"
    modelNames = c(modelNames, paste0("lkinla", urbanText, separateText))
    
    if(i+4 > endI)
      return(invisible(NULL))
    
    if(startI <= i + 4) {
      print(paste0("Fitting LK-INLA model with separateRanges=", separateRanges, " and urbanEffect=", urbanEffect, "..."))
      results = do.call("validateLKINLAKenyaDat", c(args, otherArguments))
      
      results = list(fit=results)
      # save(results, file=paste0(fileName, ".RData")) # too much space
      
      results$fit$fullModelFit = NULL
      save(results, file=paste0(fileName, "compact.RData"))
    } else {
      print(paste0("Loading LK-INLA model results with separateRanges=", separateRanges, " and urbanEffect=", urbanEffect, "..."))
      if(file.exists(paste0(fileName, "compact.RData")))
        load(paste0(fileName, ".RData"))
      else {
        load(paste0(fileName, ".RData"))
        results$fit$fullModelFit = NULL
        save(results, file=paste0(fileName, "compact.RData"))
      }
    }
    
    resultsListLKINLA = c(resultsListLKINLA, list(results))
  }
  names(resultsListLKINLA) = modelNames[5:8]
  
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
  tempModelISPDE = seq(1, 1+4*2, by=4)
  tempModelILKINLA = seq(13, 13+4*2, by=4)
  scoresInSample = scoresInSample[c(tempModelISPDE, tempModelISPDE+1, tempModelISPDE+2, tempModelISPDE+3, 
                                    tempModelILKINLA, tempModelILKINLA+1, tempModelILKINLA+2, tempModelILKINLA+3),]
  scoresInSample = data.frame(scoresInSample)
  
  # add in column saying how the scores were aggregated (across all clusters, urban clusters, or rural clusters). 
  # Also add in the model name as the name of the row
  scoresInSample = cbind("Subset"=rep(c("Avg", "Urban", "Rural"), 8), scoresInSample)
  # rownames(scoresInSample) = rep(allModelNames, each=3)
  
  ## concatenate all in sample results (accounting for binomial variation)
  scoresInSampleBinomial = rbind(t(sapply(resultsListSPDE, function(x) {colMeans(x$fit$inSamplePooledScoresBinomial)})), 
                                 t(sapply(resultsListSPDE, function(x) {colMeans(x$fit$inSampleUrbanScoresBinomial)})), 
                                 t(sapply(resultsListSPDE, function(x) {colMeans(x$fit$inSampleRuralScoresBinomial, na.rm=TRUE)})), 
                                 t(sapply(resultsListLKINLA, function(x) {colMeans(x$fit$inSamplePooledScoresBinomial)})), 
                                 t(sapply(resultsListLKINLA, function(x) {colMeans(x$fit$inSampleUrbanScoresBinomial)})), 
                                 t(sapply(resultsListLKINLA, function(x) {colMeans(x$fit$inSampleRuralScoresBinomial, na.rm=TRUE)})))
  
  # reorder them in order to group them by model rather than by urbanicity
  tempModelISPDE = seq(1, 1+4*2, by=4)
  tempModelILKINLA = seq(13, 13+4*2, by=4)
  scoresInSampleBinomial = scoresInSampleBinomial[c(tempModelISPDE, tempModelISPDE+1, tempModelISPDE+2, tempModelISPDE+3, 
                                                    tempModelILKINLA, tempModelILKINLA+1, tempModelILKINLA+2, tempModelILKINLA+3),]
  scoresInSampleBinomial = data.frame(scoresInSampleBinomial)
  
  # add in column saying how the scores were aggregated (across all clusters, urban clusters, or rural clusters). 
  # Also add in the model name as the name of the row
  scoresInSampleBinomial = cbind("Subset"=rep(c("Avg", "Urban", "Rural"), 8), scoresInSampleBinomial)
  # rownames(scoresInSample) = rep(allModelNames, each=3)
  
  ## concatenate all leave out region results
  scoresLeaveOutRegion = rbind(t(sapply(resultsListSPDE, function(x) {colMeans(x$fit$pooledScoreTable[,-1])})), 
                               t(sapply(resultsListSPDE, function(x) {colMeans(x$fit$urbanScoreTable[,-1])})), 
                               t(sapply(resultsListSPDE, function(x) {colMeans(x$fit$ruralScoreTable[,-1], na.rm=TRUE)})), 
                               t(sapply(resultsListLKINLA, function(x) {colMeans(x$fit$pooledScoreTable[,-1])})), 
                               t(sapply(resultsListLKINLA, function(x) {colMeans(x$fit$urbanScoreTable[,-1])})), 
                               t(sapply(resultsListLKINLA, function(x) {colMeans(x$fit$ruralScoreTable[,-1], na.rm=TRUE)})))
  
  # reorder them in order to group them by model rather than by urbanicity
  tempModelISPDE = seq(1, 1+4*2, by=4)
  tempModelILKINLA = seq(13, 13+4*2, by=4)
  scoresLeaveOutRegion = scoresLeaveOutRegion[c(tempModelISPDE, tempModelISPDE+1, tempModelISPDE+2, tempModelISPDE+3, 
                                                tempModelILKINLA, tempModelILKINLA+1, tempModelILKINLA+2, tempModelILKINLA+3),]
  scoresLeaveOutRegion = data.frame(scoresLeaveOutRegion)
  
  # add in column saying how the scores were aggregated (across all clusters, urban clusters, or rural clusters). 
  # Also add in the model name as the name of the row
  scoresLeaveOutRegion = cbind("Subset"=rep(c("Avg", "Urban", "Rural"), 8), scoresLeaveOutRegion)
  # rownames(scoresInSample) = rep(allModelNames, each=3)
  
  ## concatenate all leave out region results (accounting for binomial variation)
  scoresLeaveOutRegionBinomial = rbind(t(sapply(resultsListSPDE, function(x) {colMeans(x$fit$pooledScoreTableBinomial[,-1])})), 
                               t(sapply(resultsListSPDE, function(x) {colMeans(x$fit$urbanScoreTableBinomial[,-1])})), 
                               t(sapply(resultsListSPDE, function(x) {colMeans(x$fit$ruralScoreTableBinomial[,-1], na.rm=TRUE)})), 
                               t(sapply(resultsListLKINLA, function(x) {colMeans(x$fit$pooledScoreTableBinomial[,-1])})), 
                               t(sapply(resultsListLKINLA, function(x) {colMeans(x$fit$urbanScoreTableBinomial[,-1])})), 
                               t(sapply(resultsListLKINLA, function(x) {colMeans(x$fit$ruralScoreTableBinomial[,-1], na.rm=TRUE)})))
  
  # reorder them in order to group them by model rather than by urbanicity
  tempModelISPDE = seq(1, 1+4*2, by=4)
  tempModelILKINLA = seq(13, 13+4*2, by=4)
  scoresLeaveOutRegionBinomial = scoresLeaveOutRegionBinomial[c(tempModelISPDE, tempModelISPDE+1, tempModelISPDE+2, tempModelISPDE+3, 
                                                tempModelILKINLA, tempModelILKINLA+1, tempModelILKINLA+2, tempModelILKINLA+3),]
  scoresLeaveOutRegionBinomial = data.frame(scoresLeaveOutRegionBinomial)
  
  # add in column saying how the scores were aggregated (across all clusters, urban clusters, or rural clusters). 
  # Also add in the model name as the name of the row
  scoresLeaveOutRegionBinomial = cbind("Subset"=rep(c("Avg", "Urban", "Rural"), 8), scoresLeaveOutRegionBinomial)
  # rownames(scoresInSample) = rep(allModelNames, each=3)
  
  ## get all leave one out results, remove them from the end sample results
  browser() # remove CPO, WAIC, DIC columns from in sample scores, add to LOO cluster scores
  leaveOneOutI = sapply(c("CPO", "WAIC", "DIC"), function(x) {which(grepl(x, names(scoresInSample)))})
  scoresLeaveOneOut = scoresInSample[,c(1, leaveOneOutI)]
  scoresInSample = scoresInSample[,-leaveOneOutI]
  scoresLeaveOneOutBinomial = scoresInSampleBinomial[,c(1, leaveOneOutI)]
  scoresInSampleBinomial = scoresInSampleBinomial[,-leaveOneOutI]
  
  ##### Save all scoring rule tables
  fileName = paste0("savedOutput/validation/validationResults", resultNameRoot, ".RData")
  allScores = list(scoresInSample=scoresInSample, scoresLeaveOneOut=scoresLeaveOneOut, scoresLeaveOutRegion=scoresLeaveOutRegion, 
                   scoresInSampleBinomial=scoresInSampleBinomial, scoresLeaveOneOutBinomial=scoresLeaveOneOutBinomial, scoresLeaveOutRegionBinomial=scoresLeaveOutRegionBinomial)
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





