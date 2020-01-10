# function for analyzing example datasets

# first name elements of ed to be the same as the corresponding elements of the simulated datasets

generateExampleResults = function(targetPop=c("women", "children"), verbose=TRUE, startI=1) {
  targetPop = match.arg(targetPop)
  if(targetPop == "women") {
    load("../U5MR/kenyaDataEd.RData")
    dat=ed
    dataType="ed"
    resultNameRoot="Ed"
  } else {
    load("../U5MR/kenyaData.RData")
    dat=mort
    dataType="mort"
    resultNameRoot="Mort"
  }
  resultNameRootLower = tolower(resultNameRoot)
  
  ##### run SPDE 
  argList = list(list(dat = dat, clusterEffect = FALSE, urbanEffect = FALSE), 
                 list(dat = dat, clusterEffect = FALSE, urbanEffect = TRUE), 
                 list(dat = dat, clusterEffect = TRUE, urbanEffect = FALSE), 
                 list(dat = dat, clusterEffect = TRUE, urbanEffect = TRUE))
  otherArguments = list(dataType=dataType, verbose=verbose, dataType=dataType)
  
  for(i in 1:length(argList)) {
    if(startI <= i) {
      args = argList[[i]]
      clusterEffect = args$clusterEffect
      urbanEffect = args$urbanEffect
      fileName = paste0("savedOutput/resultsSPDE", resultNameRootLower, "_clusterEffect", clusterEffect, 
                        "_urbanEffect", urbanEffect, ".RData")
      
      
      print(paste0("Fitting SPDE model with clusterEffect=", clusterEffect, " and urbanEffect=", urbanEffect, "..."))
      spdeResults = do.call("fitSPDEKenyaDat", c(args, otherArguments))
      
      print(paste0("Aggregating SPDE model with clusterEffect=", clusterEffect, " and urbanEffect=", urbanEffect, "..."))
      aggregatedSPDEresults = aggregateModelResultsKenya(spdeResults, clusterLevel=TRUE, pixelLevel=TRUE, 
                                                         countyLevel=TRUE, regionLevel=TRUE, targetPop=targetPop)
      results = list(fit=spdeResults, aggregatedResults=aggregatedSPDEresults)
      save(results, file=fileName)
    }
  }
  
  ##### run LK-INLA
  argList = list(list(dat = dat, clusterEffect = FALSE, urbanEffect = FALSE), 
                 list(dat = dat, clusterEffect = FALSE, urbanEffect = TRUE), 
                 list(dat = dat, clusterEffect = TRUE, urbanEffect = FALSE), 
                 list(dat = dat, clusterEffect = TRUE, urbanEffect = TRUE))
  otherArguments = list(dataType=dataType, verbose=verbose, dataType=dataType)
  
  for(i in 1:length(argList)) {
    if(startI <= i + 4) {
      args = argList[[i]]
      clusterEffect = args$clusterEffect
      urbanEffect = args$urbanEffect
      fileName = paste0("savedOutput/resultsLKINLA", resultNameRootLower, "_clusterEffect", clusterEffect, 
                        "_urbanEffect", urbanEffect, ".RData")
      
      print(paste0("Fitting LK-INLA model with clusterEffect=", clusterEffect, " and urbanEffect=", urbanEffect, "..."))
      lkinlaResults = do.call("fitLKINLAKenyaDat", c(args, otherArguments))
      
      print(paste0("Aggregating LK-INLA model with clusterEffect=", clusterEffect, " and urbanEffect=", urbanEffect, "..."))
      aggregatedLKINLAresults = aggregateModelResultsKenya(lkinlaResults, clusterLevel=TRUE, pixelLevel=TRUE, 
                                                         countyLevel=TRUE, regionLevel=TRUE, targetPop=targetPop)
      results = list(fit=spdeResults, aggregatedResults=aggregatedLKINLAresults)
      save(results, file=fileName)
    }
  }
  invisible(NULL)
}