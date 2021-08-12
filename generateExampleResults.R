# function for analyzing example datasets

# first name elements of ed to be the same as the corresponding elements of the simulated datasets

generateExampleResults = function(targetPop=c("women", "children"), verbose=TRUE, startI=7, endI=Inf, family=c("binomial", "betabinomial"), 
                                  urbanPrior=FALSE, nBuffer=15, skipThreeLayer=FALSE) {
  family = match.arg(family)
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
  
  familyText=""
  if(family == "binomial")
    familyText = "_LgtN"
  clusterEffect = family == "binomial"
  
  ##### run SPDE 
  argList = list(list(dat = dat, urbanEffect = FALSE), 
                 list(dat = dat, urbanEffect = TRUE))
  otherArguments = list(dataType=dataType, verbose=verbose, family=family, clusterEffect=clusterEffect)
  
  for(i in 1:length(argList)) {
    if(startI <= i && i <= endI) {
      args = argList[[i]]
      urbanEffect = args$urbanEffect
      fileName = paste0("savedOutput/", resultNameRoot, "/resultsSPDE", resultNameRootLower, 
                        "_urbanEffect", urbanEffect, familyText, ".RData")
      
      print(paste0("Fitting SPDE model with urbanEffect=", urbanEffect, "..."))
      spdeResults = do.call("fitSPDEKenyaDat", c(args, otherArguments))
      
      print(paste0("Aggregating SPDE model with urbanEffect=", urbanEffect, "..."))
      aggregatedSPDEresults = aggregateModelResultsKenya(spdeResults, clusterLevel=TRUE, pixelLevel=TRUE, 
                                                         countyLevel=TRUE, regionLevel=TRUE, targetPop=targetPop)
      
      spdeResults$mlik = spdeResults$mod$mlik
      spdeResults$mod = NULL # make sure saved file isn't too large
      results = list(fit=spdeResults, aggregatedResults=aggregatedSPDEresults)
      save(results, file=fileName)
    }
  }
  
  ##### run LK-INLA
  dirichletConcentration1=1.5
  dirichletConcentration2=3
  dirichletConcentration3=5
  argList = list(list(dat = dat, separateRanges = FALSE, urbanEffect = FALSE), 
                 list(dat = dat, separateRanges = FALSE, urbanEffect = TRUE), 
                 list(dat = dat, separateRanges = TRUE, urbanEffect = FALSE), 
                 list(dat = dat, separateRanges = TRUE, urbanEffect = TRUE), 
                 list(dat = dat, separateRanges = FALSE, urbanEffect = FALSE, dirichletConcentration=dirichletConcentration2), 
                 list(dat = dat, separateRanges = FALSE, urbanEffect = TRUE, dirichletConcentration=dirichletConcentration2), 
                 list(dat = dat, separateRanges = TRUE, urbanEffect = FALSE, dirichletConcentration=dirichletConcentration2), 
                 list(dat = dat, separateRanges = TRUE, urbanEffect = TRUE, dirichletConcentration=dirichletConcentration2), 
                 list(dat = dat, separateRanges = FALSE, urbanEffect = FALSE, dirichletConcentration=dirichletConcentration3), 
                 list(dat = dat, separateRanges = FALSE, urbanEffect = TRUE, dirichletConcentration=dirichletConcentration3), 
                 list(dat = dat, separateRanges = TRUE, urbanEffect = FALSE, dirichletConcentration=dirichletConcentration3), 
                 list(dat = dat, separateRanges = TRUE, urbanEffect = TRUE, dirichletConcentration=dirichletConcentration3))
  otherArguments = list(dataType=dataType, verbose=verbose, family=family, clusterEffect=clusterEffect, 
                        useUrbanPrior=urbanPrior, nBuffer=nBuffer)
  
  for(i in 1:length(argList)) {
    if(startI <= i + 2 && i + 2 <= endI) {
      args = argList[[i]]
      separateRanges = args$separateRanges
      urbanEffect = args$urbanEffect
      dirichletConcentration = args$dirichletConcentration
      
      if(skipThreeLayer && !separateRanges) {
        next
      }
      
      urbanPriorText = ""
      if(!urbanPrior && separateRanges)
        urbanPriorText = "_noUrbanPrior"
      
      concentrationText = paste0("_conc", dirichletConcentration)
      
      fileName = paste0("savedOutput/", resultNameRoot, "/resultsLKINLA", resultNameRootLower, "_separateRanges", separateRanges, 
                        "_urbanEffect", urbanEffect, familyText, urbanPriorText, concentrationText, ".RData")
      
      print(paste0("Fitting LK-INLA model with separateRanges=", separateRanges, ", urbanEffect=", urbanEffect, ", and concentration parameter ", dirichletConcentration, "..."))
      lkinlaResults = do.call("fitLKINLAKenyaDat", c(args, otherArguments))
      
      print(paste0("Aggregating LK-INLA model with separateRanges=", separateRanges, " and urbanEffect=", urbanEffect, "..."))
      aggregatedLKINLAresults = aggregateModelResultsKenya(lkinlaResults, clusterLevel=TRUE, pixelLevel=TRUE, 
                                                         countyLevel=TRUE, regionLevel=TRUE, targetPop=targetPop)
      
      lkinlaResults$mlik = lkinlaResults$mod$mlik
      lkinlaResults$mod = NULL # make sure saved file isn't too large
      results = list(fit=lkinlaResults, aggregatedResults=aggregatedLKINLAresults)
      save(results, file=fileName)
    }
  }
  invisible(NULL)
}