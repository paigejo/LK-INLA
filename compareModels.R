# script for comparing models

# Outline of comparison:
# n simulations have been generated from a Gaussian process with a certain covariance structure
# Each simulation has been fitted with the following models:
#   Lattice Krig
#       3, 4 (?), 5 (?) layers
#   Lattice Krig INLA
#       3, 4 (?), 5 (?) layers
#   SPDE
#   GP (true covariance assumed)
# For each simulation calculate relevant scoring rules
# 


# This model produces an xtable summarizing the scoring rules aggregated at the given level of interest. 
# Same as runCompareModels, except uses the updated scoring rule function getScores, so it produces 
# scoring rules based on the predictive distribution with and without binomial variation.
# test: whether or not to use the test data set based results or the true
#       simulation study results. The test data set results are based on fewer
#       observations, requiring fewer computations and should be used to test the code
# tausq: the spatial nugget in the simulated data, can only be the default or 0
# resultType: the level of aggregation over which scoring rules should be computed
# sampling: whether or not to use the simple random sample or urban over sampled results
# recomputeTruth: by default the should be set to TRUE, unless the user calls this function 
#                 twice in a row with the same set of inputs
# modelsI: which model results should be included. Models come in the order:
#          c("GP", ""), 
#          so the default value of this argument is to include results for the naive, 
#          direct, and mercer models. 1:11 is all models
# produceFigures: whether or not to produce figures based on the results
# printIEvery: how often to print progress
# maxDataSets: if not null, set this to a small integer to only score this many datasets 
#              for testing purposes
# nsim: the number of points with which to approximate the distribution of probability on the logit scale
# saveResults: whether or not to save a disk image of the results
# loadResults: if TRUE loads the saved disk image rather than recomputing results
# xtable.args: the arguments passed to xtable for printing results
# tableFormat: if "1", binomial scores are considered the same model as non-binomial scores. 
#             i "2", binomial scores are put on extra rows of the printed table
runCompareModels2 = function(test=FALSE, tausq=.1^2, margVar=.15^2, gamma=-1, 
                             beta0=-1.75, resultType=c("county", "pixel", "EA"), 
                             sampling=c("SRS", "oversamp"), recomputeTruth=TRUE, modelsI=1:21, 
                             produceFigures=FALSE, big=FALSE, printIEvery=50, 
                             maxDataSets=NULL, nsim=10, saveResults=TRUE, loadResults=TRUE, 
                             xtable.args=list(digits=c(0, 2, 2, 2, 2, 1, 2), display=rep("f", 7), auto=TRUE), 
                             tableFormat=c("2", "1"), colScale=c(10^4, 10^5, 10^5, 10^3, 100, 100), 
                             colUnits=c(" ($\\times 10^{-4}$)", " ($\\times 10^{-5}$)", " ($\\times 10^{-5}$)", 
                                        " ($\\times 10^{-3}$)", " ($\\times 10^{-2}$)", " ($\\times 10^{-2}$)"), 
                             colDigits=c(2, 2, 2, 2, 1, 2), counties=sort(unique(poppc$admin1)), 
                             loadTempProgress=FALSE, includeBVarResults=FALSE, continuousSPDEonly=TRUE, 
                             strictPriors=FALSE, doFancyTables=FALSE, printScoreTable=TRUE, printParTable=TRUE) {
  
  # match the arguments with their correct values
  resultType = match.arg(resultType)
  sampling = match.arg(sampling)
  tableFormat = match.arg(tableFormat)
  
  # test to make sure only the naive and direct models are included if big is set to TRUE
  if(!identical(modelsI, 1:2) && big)
    stop("if big is set to TRUE, only naive and direct results can be included")
  
  # load the true superpopulation
  testText = ifelse(test, "Test", "")
  bigText = ifelse(big, "Big", "")
  strictPriorText = ifelse(strictPriors, "strictPrior", "")
  if(!test)
    load(paste0("simDataMultiBeta-1.75margVar", round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 
                "HHoldVar0urbanOverSamplefrac0", bigText, ".RData"))
  else
    load(paste0("simDataMultiBeta-1.75margVar", round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 
                "HHoldVar0urbanOverSamplefrac0Test", bigText, ".RData"))
  eaDat = SRSDat$eaDat
  
  if(sampling == "SRS") {
    clustDat = SRSDat
  } else {
    clustDat = overSampDat
  }
  maxDataSets = ifelse(is.null(maxDataSets), length(clustDat$clustDat), maxDataSets)
  
  # allModels = c("naive", "direct", "mercer", "bym", "bymMod", "bymNoUrb", "bymNoUrbMod", "bymNoClust", "bymNoUrbClust", "spde", "spdeNoUrb")
  # allNames = c("Naive", "Direct ", "Smoothed Direct", "BYM (no urban/cluster)", "BYM (no urban)", "BYM (no cluster)", "BYM", "SPDE (no urban)", "SPDE")
  # allNamesBinomial = c("Naive Binom.", "Direct Binom.", "Mercer et al. Binom.", "BYM Binom. (no urb/clust)", "BYM Binom. (no urb)", "BYM Binom. (no clust)", "BYM Binom.", "SPDE Binom. (no urb)", "SPDE Binom.")
  # BYM models are in order of complexity: no urban/cluster, no urban, no cluster, full
  allNames = c("Naive", "Direct", "Smoothed Direct", "BYM2 ucA", "BYM2 uCA", "BYM2 uCA'", "BYM2 UcA", "BYM2 UCA", "BYM2 UCA'", 
               "BYM2 uca", "BYM2 uCa", "BYM2 uCa'", "BYM2 Uca", "BYM2 UCa", "BYM2 UCa'", 
               "SPDE uc", "SPDE uC", "SPDE Uc", "SPDE UC", "SPDE uC'", "SPDE UC'")
  allNamesBinomial = paste0(allNames, " Bin.")
  allModels = allNames
  models = allModels[modelsI]
  
  # this string carries all the information about the run
  runId = paste0("Beta-1.75margVar", round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 
                 "HHoldVar0urbanOverSamplefrac0", strictPriorText, testText, bigText, sampling, 
                 "models", do.call("paste0", as.list(modelsI)), "nsim", nsim, "MaxDataSetI", maxDataSets)
  
  # compute the results if necessary
  if(!loadResults && !loadTempProgress) {
    
    ## generate the true values at the county level for the 
    ## given settings if these are new settings
    counties = sort(unique(eaDat$admin1))
    if(recomputeTruth){
      truth = getTruth(resultType, eaDat)
      if(resultType == "county") {
        save(truth, file="truthbycounty.RData")
      } else if(resultType == "pixel") {
        save(truth, file="truthbyPixel.RData")
      } else if(resultType == "EA") {
        save(truth, file="truthbyEA.RData")
      }
    } else {
      # load the truth
      load(paste0("truthby", resultType, ".RData"))
    }
    
    ## Now pick which models to include in the results
    
    # load data
    tauText = ifelse(tausq == 0, "0", "0.01")
    testText = ifelse(test, "Test", "")
    if("Naive" %in% models || "Direct" %in% models) {
      out = load(paste0("resultsDirectNaiveBeta-1.75margVar", round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 
                        "HHoldVar0urbanOverSamplefrac0", testText, bigText, ".RData"))
      if(sampling == "SRS") {
        directEst = directEstSRS
        naive = naiveSRS
      }
      else {
        directEst = directEstoverSamp
        naive = naiveoverSamp
      }
    }
    if("Smoothed Direct" %in% models) {
      if(!test)
        load(paste0("resultsMercerBeta-1.75margVar", round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 
                    "HHoldVar0urbanOverSamplefrac0", strictPriorText, ".RData"))
      else
        load(paste0("resultsMercerBeta-1.75margVar", round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 
                    "HHoldVar0urbanOverSamplefrac0", strictPriorText, "Test.RData"))
      if(sampling == "SRS") {
        mercer = mercerSRS
        mercerPar = mercerSRSPar
      }
      else {
        mercer = merceroverSamp
        mercerPar = merceroverSampPar
      }
    }
    if("BYM2 ucA" %in% models) {
      includeUrbanRural = FALSE
      includeCluster = FALSE
      aggregateByPopulation = FALSE
      load(paste0('bym2Beta-1.75margVar', round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 'UrbRur',
                  includeUrbanRural, 'Cluster', includeCluster, "aggByPop", aggregateByPopulation, "maxDataSets", 100, strictPriorText, testText, '.RData'))
      # out = load(paste0("kenyaSpatialDesignResultNewTausq", tauText, "UrbRurFALSEClusterTRUE", strictPriorText, testText, ".RData"))
      if(sampling == "SRS") {
        designRes$overSampDat = NULL
        designRes$overSampDatPar = NULL
        if(resultType == "pixel")
          designRes[[1]] = designRes$SRSdatPixelUrban
        else if(resultType == "EA")
          designRes[[1]] = designRes$SRSdatClusterUrban
      }
      else {
        designRes$SRSdat = NULL
        designRes$SRSdatPar = NULL
        if(resultType == "pixel")
          designRes[[1]] = designRes$overSampDatPixelUrban
        else if(resultType == "EA")
          designRes[[1]] = designRes$overSampDatClusterUrban
      }
      designResNoUrbClust = designRes
    }
    if("BYM2 uCA" %in% models) {
      includeUrbanRural = FALSE
      includeCluster = TRUE
      aggregateByPopulation = FALSE
      load(paste0('bym2Beta-1.75margVar', round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 'UrbRur',
                  includeUrbanRural, 'Cluster', includeCluster, "aggByPop", aggregateByPopulation, "maxDataSets", 100, strictPriorText, testText, '.RData'))
      if(sampling == "SRS") {
        designRes$overSampDat = NULL
        designRes$overSampDatPar = NULL
        if(resultType == "pixel")
          designRes[[1]] = designRes$SRSdatPixelUrban
        else if(resultType == "EA")
          designRes[[1]] = designRes$SRSdatClusterUrban
      }
      else {
        designRes$SRSdat = NULL
        designRes$SRSdatPar = NULL
        if(resultType == "pixel")
          designRes[[1]] = designRes$overSampDatPixelUrban
        else if(resultType == "EA")
          designRes[[1]] = designRes$overSampDatClusterUrban
      }
      designResNoUrb = designRes
    }
    if("BYM2 uCa'" %in% models) {
      includeUrbanRural = FALSE
      includeCluster = TRUE
      aggregateByPopulation = FALSE
      load(paste0('bym2Beta-1.75margVar', round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 'UrbRur',
                  includeUrbanRural, 'Cluster', includeCluster, "aggByPop", aggregateByPopulation, "debiasedMaxDataSets", 100, strictPriorText, testText, '.RData'))
      if(sampling == "SRS") {
        designRes$overSampDat = NULL
        designRes$overSampDatPar = NULL
        if(resultType == "pixel")
          designRes[[1]] = designRes$SRSdatPixelUrban
        else if(resultType == "EA")
          designRes[[1]] = designRes$SRSdatClusterUrban
      }
      else {
        designRes$SRSdat = NULL
        designRes$SRSdatPar = NULL
        if(resultType == "pixel")
          designRes[[1]] = designRes$overSampDatPixelUrban
        else if(resultType == "EA")
          designRes[[1]] = designRes$overSampDatClusterUrban
      }
      designResNoUrbMod = designRes
    }
    if("BYM2 UcA" %in% models) {
      includeUrbanRural = TRUE
      includeCluster = FALSE
      aggregateByPopulation = FALSE
      load(paste0('bym2Beta-1.75margVar', round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 'UrbRur',
                  includeUrbanRural, 'Cluster', includeCluster, "aggByPop", aggregateByPopulation, "maxDataSets", 100, strictPriorText, testText, '.RData'))
      if(sampling == "SRS") {
        designRes$overSampDat = NULL
        designRes$overSampDatPar = NULL
        if(resultType == "pixel") {
          designRes[[2]] = designRes$SRSdatPixelUrban
          designRes[[1]] = designRes$SRSdatPixelRural
        }
        else if(resultType == "EA") {
          designRes[[2]] = designRes$SRSdatClusterUrban
          designRes[[1]] = designRes$SRSdatClusterRural
        }
      }
      else {
        designRes$SRSdat = NULL
        designRes$SRSdatPar = NULL
        if(resultType == "pixel") {
          designRes[[2]] = designRes$overSampDatPixelUrban
          designRes[[1]] = designRes$overSampDatPixelRural
        }
        else if(resultType == "EA") {
          designRes[[2]] = designRes$overSampDatClusterUrban
          designRes[[1]] = designRes$overSampDatClusterRural
        }
      }
      designResNoClust = designRes
    }
    if("BYM2 UCA'" %in% models) {
      includeUrbanRural = TRUE
      includeCluster = TRUE
      aggregateByPopulation = FALSE
      load(paste0('bym2Beta-1.75margVar', round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 'UrbRur',
                  includeUrbanRural, 'Cluster', includeCluster, "aggByPop", aggregateByPopulation, "debiasedMaxDataSets", 100, strictPriorText, testText, '.RData'))
      if(sampling == "SRS") {
        designRes$overSampDat = NULL
        designRes$overSampDatPar = NULL
        if(resultType == "pixel") {
          designRes[[2]] = designRes$SRSdatPixelUrban
          designRes[[1]] = designRes$SRSdatPixelRural
        }
        else if(resultType == "EA") {
          designRes[[2]] = designRes$SRSdatClusterUrban
          designRes[[1]] = designRes$SRSdatClusterRural
        }
      }
      else {
        designRes$SRSdat = NULL
        designRes$SRSdatPar = NULL
        if(resultType == "pixel") {
          designRes[[2]] = designRes$overSampDatPixelUrban
          designRes[[1]] = designRes$overSampDatPixelRural
        }
        else if(resultType == "EA") {
          designRes[[2]] = designRes$overSampDatClusterUrban
          designRes[[1]] = designRes$overSampDatClusterRural
        }
      }
      designResMod = designRes
    }
    if("BYM2 UCA" %in% models) {
      includeUrbanRural = TRUE
      includeCluster = TRUE
      aggregateByPopulation = FALSE
      load(paste0('bym2Beta-1.75margVar', round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 'UrbRur',
                  includeUrbanRural, 'Cluster', includeCluster, "aggByPop", aggregateByPopulation, "maxDataSets", 100, strictPriorText, testText, '.RData'))
      if(sampling == "SRS") {
        designRes$overSampDat = NULL
        designRes$overSampDatPar = NULL
        if(resultType == "pixel") {
          designRes[[2]] = designRes$SRSdatPixelUrban
          designRes[[1]] = designRes$SRSdatPixelRural
        }
        else if(resultType == "EA") {
          designRes[[2]] = designRes$SRSdatClusterUrban
          designRes[[1]] = designRes$SRSdatClusterRural
        }
      }
      else {
        designRes$SRSdat = NULL
        designRes$SRSdatPar = NULL
        if(resultType == "pixel") {
          designRes[[2]] = designRes$overSampDatPixelUrban
          designRes[[1]] = designRes$overSampDatPixelRural
        }
        else if(resultType == "EA") {
          designRes[[2]] = designRes$overSampDatClusterUrban
          designRes[[1]] = designRes$overSampDatClusterRural
        }
      }
      designResTemp = designRes
    }
    if("BYM2 uca" %in% models) {
      includeUrbanRural = FALSE
      includeCluster = FALSE
      aggregateByPopulation = TRUE
      load(paste0('bym2Beta-1.75margVar', round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 'UrbRur',
                  includeUrbanRural, 'Cluster', includeCluster, "aggByPop", aggregateByPopulation, "maxDataSets", 100, strictPriorText, testText, '.RData'))
      # out = load(paste0("kenyaSpatialDesignResultNewTausq", tauText, "UrbRurFALSEClusterTRUE", strictPriorText, testText, ".RData"))
      if(sampling == "SRS") {
        designRes$overSampDat = NULL
        designRes$overSampDatPar = NULL
        if(resultType == "pixel")
          designRes[[1]] = designRes$SRSdatPixelUrban
        else if(resultType == "EA")
          designRes[[1]] = designRes$SRSdatClusterUrban
      }
      else {
        designRes$SRSdat = NULL
        designRes$SRSdatPar = NULL
        if(resultType == "pixel")
          designRes[[1]] = designRes$overSampDatPixelUrban
        else if(resultType == "EA")
          designRes[[1]] = designRes$overSampDatClusterUrban
      }
      designResNoUrbClustPopAgg = designRes
    }
    if("BYM2 uCa" %in% models) {
      includeUrbanRural = FALSE
      includeCluster = TRUE
      aggregateByPopulation = TRUE
      load(paste0('bym2Beta-1.75margVar', round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 'UrbRur',
                  includeUrbanRural, 'Cluster', includeCluster, "aggByPop", aggregateByPopulation, "maxDataSets", 100, strictPriorText, testText, '.RData'))
      if(sampling == "SRS") {
        designRes$overSampDat = NULL
        designRes$overSampDatPar = NULL
        if(resultType == "pixel")
          designRes[[1]] = designRes$SRSdatPixelUrban
        else if(resultType == "EA")
          designRes[[1]] = designRes$SRSdatClusterUrban
      }
      else {
        designRes$SRSdat = NULL
        designRes$SRSdatPar = NULL
        if(resultType == "pixel")
          designRes[[1]] = designRes$overSampDatPixelUrban
        else if(resultType == "EA")
          designRes[[1]] = designRes$overSampDatClusterUrban
      }
      designResNoUrbPopAgg = designRes
    }
    if("BYM2 uCa'" %in% models) {
      includeUrbanRural = FALSE
      includeCluster = TRUE
      aggregateByPopulation = TRUE
      load(paste0('bym2Beta-1.75margVar', round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 'UrbRur',
                  includeUrbanRural, 'Cluster', includeCluster, "aggByPop", aggregateByPopulation, "debiasedMaxDataSets", 100, strictPriorText, testText, '.RData'))
      if(sampling == "SRS") {
        designRes$overSampDat = NULL
        designRes$overSampDatPar = NULL
        if(resultType == "pixel")
          designRes[[1]] = designRes$SRSdatPixelUrban
        else if(resultType == "EA")
          designRes[[1]] = designRes$SRSdatClusterUrban
      }
      else {
        designRes$SRSdat = NULL
        designRes$SRSdatPar = NULL
        if(resultType == "pixel")
          designRes[[1]] = designRes$overSampDatPixelUrban
        else if(resultType == "EA")
          designRes[[1]] = designRes$overSampDatClusterUrban
      }
      designResNoUrbModPopAgg = designRes
    }
    if("BYM2 Uca" %in% models) {
      includeUrbanRural = TRUE
      includeCluster = FALSE
      aggregateByPopulation = TRUE
      load(paste0('bym2Beta-1.75margVar', round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 'UrbRur',
                  includeUrbanRural, 'Cluster', includeCluster, "aggByPop", aggregateByPopulation, "maxDataSets", 100, strictPriorText, testText, '.RData'))
      if(sampling == "SRS") {
        designRes$overSampDat = NULL
        designRes$overSampDatPar = NULL
        if(resultType == "pixel") {
          designRes[[2]] = designRes$SRSdatPixelUrban
          designRes[[1]] = designRes$SRSdatPixelRural
        }
        else if(resultType == "EA") {
          designRes[[2]] = designRes$SRSdatClusterUrban
          designRes[[1]] = designRes$SRSdatClusterRural
        }
      }
      else {
        designRes$SRSdat = NULL
        designRes$SRSdatPar = NULL
        if(resultType == "pixel") {
          designRes[[2]] = designRes$overSampDatPixelUrban
          designRes[[1]] = designRes$overSampDatPixelRural
        }
        else if(resultType == "EA") {
          designRes[[2]] = designRes$overSampDatClusterUrban
          designRes[[1]] = designRes$overSampDatClusterRural
        }
      }
      designResNoClustPopAgg = designRes
    }
    if("BYM2 UCa'" %in% models) {
      includeUrbanRural = TRUE
      includeCluster = TRUE
      aggregateByPopulation = TRUE
      load(paste0('bym2Beta-1.75margVar', round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 'UrbRur',
                  includeUrbanRural, 'Cluster', includeCluster, "aggByPop", aggregateByPopulation, "debiasedMaxDataSets", 100, strictPriorText, testText, '.RData'))
      if(sampling == "SRS") {
        designRes$overSampDat = NULL
        designRes$overSampDatPar = NULL
        if(resultType == "pixel") {
          designRes[[2]] = designRes$SRSdatPixelUrban
          designRes[[1]] = designRes$SRSdatPixelRural
        }
        else if(resultType == "EA") {
          designRes[[2]] = designRes$SRSdatClusterUrban
          designRes[[1]] = designRes$SRSdatClusterRural
        }
      }
      else {
        designRes$SRSdat = NULL
        designRes$SRSdatPar = NULL
        if(resultType == "pixel") {
          designRes[[2]] = designRes$overSampDatPixelUrban
          designRes[[1]] = designRes$overSampDatPixelRural
        }
        else if(resultType == "EA") {
          designRes[[2]] = designRes$overSampDatClusterUrban
          designRes[[1]] = designRes$overSampDatClusterRural
        }
      }
      designResModPopAgg = designRes
    }
    if("BYM2 UCa" %in% models) {
      includeUrbanRural = TRUE
      includeCluster = TRUE
      aggregateByPopulation = TRUE
      load(paste0('bym2Beta-1.75margVar', round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 'UrbRur',
                  includeUrbanRural, 'Cluster', includeCluster, "aggByPop", aggregateByPopulation, "maxDataSets", 100, strictPriorText, testText, '.RData'))
      if(sampling == "SRS") {
        designRes$overSampDat = NULL
        designRes$overSampDatPar = NULL
        if(resultType == "pixel") {
          designRes[[2]] = designRes$SRSdatPixelUrban
          designRes[[1]] = designRes$SRSdatPixelRural
        }
        else if(resultType == "EA") {
          designRes[[2]] = designRes$SRSdatClusterUrban
          designRes[[1]] = designRes$SRSdatClusterRural
        }
      }
      else {
        designRes$SRSdat = NULL
        designRes$SRSdatPar = NULL
        if(resultType == "pixel") {
          designRes[[2]] = designRes$overSampDatPixelUrban
          designRes[[1]] = designRes$overSampDatPixelRural
        }
        else if(resultType == "EA") {
          designRes[[2]] = designRes$overSampDatClusterUrban
          designRes[[1]] = designRes$overSampDatClusterRural
        }
      }
      designResPopAgg = designRes
    }
    if("BYM2 UCA" %in% models)
      designRes = designResTemp
    if("SPDE uc" %in% models) {
      urbanEffect = FALSE
      includeClustEffect = FALSE
      testText = ifelse(test, "Test", "")
      fileName = paste0("resultsSPDEBeta", round(beta0, 4), "margVar", round(margVar, 4), "tausq", 
                        round(tausq, 4), "gamma", round(gamma, 4), "HHoldVar0urbanOverSamplefrac0", 
                        "urbanEffect", urbanEffect, "clustEffect", includeClustEffect, strictPriorText, testText, ".RData")
      out = load(fileName)
      if(sampling == "SRS")
        spdeNoUrbClust = spdeSRS
      else
        spdeNoUrbClust = spdeOverSamp
    }
    if("SPDE uC" %in% models) {
      urbanEffect = FALSE
      includeClustEffect = TRUE
      testText = ifelse(test, "Test", "")
      fileName = paste0("resultsSPDEBeta", round(beta0, 4), "margVar", round(margVar, 4), "tausq", 
                        round(tausq, 4), "gamma", round(gamma, 4), "HHoldVar0urbanOverSamplefrac0", 
                        "urbanEffect", urbanEffect, "clustEffect", includeClustEffect, strictPriorText, testText, ".RData")
      out = load(fileName)
      if(sampling == "SRS")
        spdeNoUrb = spdeSRS
      else
        spdeNoUrb = spdeOverSamp
    }
    if("SPDE Uc" %in% models) {
      urbanEffect = TRUE
      includeClustEffect = FALSE
      testText = ifelse(test, "Test", "")
      fileName = paste0("resultsSPDEBeta", round(beta0, 4), "margVar", round(margVar, 4), "tausq", 
                        round(tausq, 4), "gamma", round(gamma, 4), "HHoldVar0urbanOverSamplefrac0", 
                        "urbanEffect", urbanEffect, "clustEffect", includeClustEffect, strictPriorText, testText, ".RData")
      out = load(fileName)
      if(sampling == "SRS")
        spdeNoClust = spdeSRS
      else
        spdeNoClust = spdeOverSamp
    }
    if("SPDE UC" %in% models) {
      urbanEffect = TRUE
      includeClustEffect = TRUE
      testText = ifelse(test, "Test", "")
      fileName = paste0("resultsSPDEBeta", round(beta0, 4), "margVar", round(margVar, 4), "tausq", 
                        round(tausq, 4), "gamma", round(gamma, 4), "HHoldVar0urbanOverSamplefrac0", 
                        "urbanEffect", urbanEffect, "clustEffect", includeClustEffect, strictPriorText, testText, ".RData")
      out = load(fileName)
      if(sampling == "SRS")
        spde = spdeSRS
      else
        spde = spdeOverSamp
    }
    # if("SPDE uC'" %in% models) {
    #   urbanEffect = FALSE
    #   includeClustEffect = TRUE
    #   testText = ifelse(test, "Test", "")
    #   fileName = paste0("resultsSPDEBeta", round(beta0, 4), "margVar", round(margVar, 4), "tausq", 
    #                     round(tausq, 4), "gamma", round(gamma, 4), "HHoldVar0urbanOverSamplefrac0", 
    #                     "urbanEffect", urbanEffect, "clustEffect", includeClustEffect, strictPriorText, testText, 
    #                     ".RData")
    #   out = load(fileName)
    #   if(sampling == "SRS")
    #     spdeNoUrbMod = spdeSRS
    #   else
    #     spdeNoUrbMod = spdeOverSamp
    # }
    # if("SPDE UC'" %in% models) {
    #   urbanEffect = TRUE
    #   includeClustEffect = TRUE
    #   testText = ifelse(test, "Test", "")
    #   fileName = paste0("resultsSPDEBeta", round(beta0, 4), "margVar", round(margVar, 4), "tausq", 
    #                     round(tausq, 4), "gamma", round(gamma, 4), "HHoldVar0urbanOverSamplefrac0", 
    #                     "urbanEffect", urbanEffect, "clustEffect", includeClustEffect, strictPriorText, testText, 
    #                     ".RData")
    #   out = load(fileName)
    #   if(sampling == "SRS")
    #     spdeMod = spdeSRS
    #   else
    #     spdeMod = spdeOverSamp
    # }
    
    # relabel direct, naive, and mercer county names
    counties=sort(unique(poppc$admin1))
    for(i in 1:maxDataSets) {
      if("Smoothed Direct" %in% models)
        mercer[[i]]$admin1 = counties
      if("Direct" %in% models)
        directEst[[i]]$admin1 = counties
      if("Naive" %in% models)
        naive[[i]]$admin1 = counties
    }
    
    # commented out the below code, since it's not clear within strata predictions are necessary
    # # modify discrete estimates stratified by urban/rural if resultType is at finer scale than county
    # # (currently, further resolving discrete estimates is only supported for the direct estimates)
    # if(resultType == "EA") {
    #   filterByUrban = function(i) {
    #     urbanI = 1:nrow(directEstSRS[[i]])
    #   }
    #   directEstSRS = lapply(1:100, )
    # }
    
    # function for generating discrete model pixel and EA level predictions 
    # given the county level predictions. In the case that two result tables are given, 
    # the first is assumed to be for rural areas, and the second is for urban areas. 
    # Also, urbanVec is a vector of urban logical labels for the desired aggregation level
    getSubLevelResults = function(resultTable, resultTableUrban=NULL, urbanVec=NULL) {
      if(!is.null(resultTableUrban) && !is.null(urbanVec)) {
        tabRural = getSubLevelResults(resultTableRural)
        tabUrban = getSubLevelResults(resultTableUrban)
        tabRural[urbanVec,] = tabUrban[UrbanVec,]
        return(tabRural)
      }
      
      # get the order of the counties that resultTable is in
      counties = sort(unique(eaDat$admin1))
      
      # convert results to the desired aggregation level if necessary for discrete models:
      if(resultType == "pixel") {
        # compute pixel results
        pixelToAdmin = match(popGrid$admin1, counties)
        resultTable = resultTable[pixelToAdmin,]
      }
      else if(resultType == "EA") {
        # compute EA level results
        eaToAdmin = match(eaDat$admin1, counties)
        resultTable = resultTable[eaToAdmin,]
      }
      
      # modify result row and column table names according to aggregation level
      whichName = which(names(resultTable) == "admin1")
      names(resultTable)[whichName] = resultType
      if((resultType == "EA") && (is.data.frame(resultTable))) {
        # resultTable[[resultType]] = resultTable$eaI
        resultTable[[resultType]] = 1:nrow(resultTable)
      } else if((resultType == "pixel") && (is.data.frame(resultTable))) {
        # resultTable[[resultType]] = resultTable$eaI
        resultTable[[resultType]] = 1:nrow(resultTable)
      }
      
      if(resultType == "pixel")
        resultTable = resultTable[as.numeric(as.character(truth$pixel)),]
      
      if(is.data.frame(resultTable) && resultType == "county" && is.null(resultTable[[resultType]]))
        resultTable[[resultType]] = counties
      
      resultTable
    }
    
    # convert the truth to the desired aggregation level
    if(resultType == "county")
      truth = getSubLevelResults(truth)
    
    # calculate the binomial n (number of children) for each prediction unit
    numChildren = truth$numChildren
    
    # compute scores
    scoresDirect = scoresNaive = scoresMercer = scoresBYMNoUrb = scoresBYM = scoresBYMNoUrbMod = scoresBYMMod = scoresBYMNoUrbClust = scoresBYMNoClust = 
      scoresBYMNoUrbPopAgg = scoresBYMPopAgg = scoresBYMNoUrbModPopAgg = scoresBYMModPopAgg = scoresBYMNoUrbClustPopAgg = scoresBYMNoClustPopAgg = 
      scoresSPDENoUrbClust = scoresSPDENoUrb = scoresSPDENoClust = scoresSPDE = data.frame()
    
    # convert results to the desired aggregation level 
    # not including urban or cluster effect
    if("BYM2 ucA" %in% models) {
      designResNoUrbClust[[1]]$Q10 = getSubLevelResults(designResNoUrbClust[[1]]$Q10)
      designResNoUrbClust[[1]]$Q50 = getSubLevelResults(designResNoUrbClust[[1]]$Q50)
      designResNoUrbClust[[1]]$Q90 = getSubLevelResults(designResNoUrbClust[[1]]$Q90)
      designResNoUrbClust[[1]]$mean = getSubLevelResults(designResNoUrbClust[[1]]$mean)
      designResNoUrbClust[[1]]$stddev = getSubLevelResults(designResNoUrbClust[[1]]$stddev)
    }
    
    # not including urban effect
    if("BYM2 uCA" %in% models) {
      designResNoUrb[[1]]$Q10 = getSubLevelResults(designResNoUrb[[1]]$Q10)
      designResNoUrb[[1]]$Q50 = getSubLevelResults(designResNoUrb[[1]]$Q50)
      designResNoUrb[[1]]$Q90 = getSubLevelResults(designResNoUrb[[1]]$Q90)
      designResNoUrb[[1]]$mean = getSubLevelResults(designResNoUrb[[1]]$mean)
      designResNoUrb[[1]]$stddev = getSubLevelResults(designResNoUrb[[1]]$stddev)
    }
    
    # including urban effect
    if("BYM2 UcA" %in% models) {
      if(resultType != "county") {
        designResNoClust[[1]]$Q10 = getSubLevelResults(designResNoClust[[1]]$Q10, designResNoClust[[2]]$Q10, truth$urban)
        designResNoClust[[1]]$Q50 = getSubLevelResults(designResNoClust[[1]]$Q50, designResNoClust[[2]]$Q50, truth$urban)
        designResNoClust[[1]]$Q90 = getSubLevelResults(designResNoClust[[1]]$Q90, designResNoClust[[2]]$Q90, truth$urban)
        designResNoClust[[1]]$mean = getSubLevelResults(designResNoClust[[1]]$mean, designResNoClust[[2]]$mean, truth$urban)
        designResNoClust[[1]]$stddev = getSubLevelResults(designResNoClust[[1]]$stddev, designResNoClust[[2]]$stddev, truth$urban)
      }
      else {
        designResNoClust[[1]]$Q10 = getSubLevelResults(designResNoClust[[1]]$Q10)
        designResNoClust[[1]]$Q50 = getSubLevelResults(designResNoClust[[1]]$Q50)
        designResNoClust[[1]]$Q90 = getSubLevelResults(designResNoClust[[1]]$Q90)
        designResNoClust[[1]]$mean = getSubLevelResults(designResNoClust[[1]]$mean)
        designResNoClust[[1]]$stddev = getSubLevelResults(designResNoClust[[1]]$stddev)
      }
    }
    
    # including urban effect
    if("BYM2 UCA" %in% models) {
      if(resultType != "county") {
        designRes[[1]]$Q10 = getSubLevelResults(designRes[[1]]$Q10, designRes[[2]]$Q10, truth$urban)
        designRes[[1]]$Q50 = getSubLevelResults(designRes[[1]]$Q50, designRes[[2]]$Q50, truth$urban)
        designRes[[1]]$Q90 = getSubLevelResults(designRes[[1]]$Q90, designRes[[2]]$Q90, truth$urban)
        designRes[[1]]$mean = getSubLevelResults(designRes[[1]]$mean, designRes[[2]]$mean, truth$urban)
        designRes[[1]]$stddev = getSubLevelResults(designRes[[1]]$stddev, designRes[[2]]$stddev, truth$urban)
      }
      else {
        designRes[[1]]$Q10 = getSubLevelResults(designRes[[1]]$Q10)
        designRes[[1]]$Q50 = getSubLevelResults(designRes[[1]]$Q50)
        designRes[[1]]$Q90 = getSubLevelResults(designRes[[1]]$Q90)
        designRes[[1]]$mean = getSubLevelResults(designRes[[1]]$mean)
        designRes[[1]]$stddev = getSubLevelResults(designRes[[1]]$stddev)
      }
    }
    
    # not including urban effect, modified to be debiased using marginal rather than conditional effect as prediction
    if("BYM2 uCA'" %in% models) {
      designResNoUrbMod[[1]]$Q10 = getSubLevelResults(designResNoUrbMod[[1]]$Q10)
      designResNoUrbMod[[1]]$Q50 = getSubLevelResults(designResNoUrbMod[[1]]$Q50)
      designResNoUrbMod[[1]]$Q90 = getSubLevelResults(designResNoUrbMod[[1]]$Q90)
      designResNoUrbMod[[1]]$mean = getSubLevelResults(designResNoUrbMod[[1]]$mean)
      designResNoUrbMod[[1]]$stddev = getSubLevelResults(designResNoUrbMod[[1]]$stddev)
    }
    
    # including urban effect, modified to be debiased using marginal rather than conditional effect as prediction
    if("BYM2 UCA'" %in% models) {
      if(resultType != "county") {
        designResMod[[1]]$Q10 = getSubLevelResults(designResMod[[1]]$Q10, designResMod[[2]]$Q10, truth$urban)
        designResMod[[1]]$Q50 = getSubLevelResults(designResMod[[1]]$Q50, designResMod[[2]]$Q50, truth$urban)
        designResMod[[1]]$Q90 = getSubLevelResults(designResMod[[1]]$Q90, designResMod[[2]]$Q90, truth$urban)
        designResMod[[1]]$mean = getSubLevelResults(designResMod[[1]]$mean, designResMod[[2]]$mean, truth$urban)
        designResMod[[1]]$stddev = getSubLevelResults(designResMod[[1]]$stddev, designResMod[[2]]$stddev, truth$urban)
      }
      else {
        designResMod[[1]]$Q10 = getSubLevelResults(designResMod[[1]]$Q10)
        designResMod[[1]]$Q50 = getSubLevelResults(designResMod[[1]]$Q50)
        designResMod[[1]]$Q90 = getSubLevelResults(designResMod[[1]]$Q90)
        designResMod[[1]]$mean = getSubLevelResults(designResMod[[1]]$mean)
        designResMod[[1]]$stddev = getSubLevelResults(designResMod[[1]]$stddev)
      }
    }
    
    if("BYM2 uca" %in% models) {
      designResNoUrbClustPopAgg[[1]]$Q10 = getSubLevelResults(designResNoUrbClustPopAgg[[1]]$Q10)
      designResNoUrbClustPopAgg[[1]]$Q50 = getSubLevelResults(designResNoUrbClustPopAgg[[1]]$Q50)
      designResNoUrbClustPopAgg[[1]]$Q90 = getSubLevelResults(designResNoUrbClustPopAgg[[1]]$Q90)
      designResNoUrbClustPopAgg[[1]]$mean = getSubLevelResults(designResNoUrbClustPopAgg[[1]]$mean)
      designResNoUrbClustPopAgg[[1]]$stddev = getSubLevelResults(designResNoUrbClustPopAgg[[1]]$stddev)
    }
    
    # not including urban effect
    if("BYM2 uCa" %in% models) {
      designResNoUrbPopAgg[[1]]$Q10 = getSubLevelResults(designResNoUrbPopAgg[[1]]$Q10)
      designResNoUrbPopAgg[[1]]$Q50 = getSubLevelResults(designResNoUrbPopAgg[[1]]$Q50)
      designResNoUrbPopAgg[[1]]$Q90 = getSubLevelResults(designResNoUrbPopAgg[[1]]$Q90)
      designResNoUrbPopAgg[[1]]$mean = getSubLevelResults(designResNoUrbPopAgg[[1]]$mean)
      designResNoUrbPopAgg[[1]]$stddev = getSubLevelResults(designResNoUrbPopAgg[[1]]$stddev)
    }
    
    # including urban effect
    if("BYM2 Uca" %in% models) {
      if(resultType != "county") {
        designResNoClustPopAgg[[1]]$Q10 = getSubLevelResults(designResNoClustPopAgg[[1]]$Q10, designResNoClustPopAgg[[2]]$Q10, truth$urban)
        designResNoClustPopAgg[[1]]$Q50 = getSubLevelResults(designResNoClustPopAgg[[1]]$Q50, designResNoClustPopAgg[[2]]$Q50, truth$urban)
        designResNoClustPopAgg[[1]]$Q90 = getSubLevelResults(designResNoClustPopAgg[[1]]$Q90, designResNoClustPopAgg[[2]]$Q90, truth$urban)
        designResNoClustPopAgg[[1]]$mean = getSubLevelResults(designResNoClustPopAgg[[1]]$mean, designResNoClustPopAgg[[2]]$mean, truth$urban)
        designResNoClustPopAgg[[1]]$stddev = getSubLevelResults(designResNoClustPopAgg[[1]]$stddev, designResNoClustPopAgg[[2]]$stddev, truth$urban)
      }
      else {
        designResNoClustPopAgg[[1]]$Q10 = getSubLevelResults(designResNoClustPopAgg[[1]]$Q10)
        designResNoClustPopAgg[[1]]$Q50 = getSubLevelResults(designResNoClustPopAgg[[1]]$Q50)
        designResNoClustPopAgg[[1]]$Q90 = getSubLevelResults(designResNoClustPopAgg[[1]]$Q90)
        designResNoClustPopAgg[[1]]$mean = getSubLevelResults(designResNoClustPopAgg[[1]]$mean)
        designResNoClustPopAgg[[1]]$stddev = getSubLevelResults(designResNoClustPopAgg[[1]]$stddev)
      }
    }
    
    # including urban effect
    if("BYM2 UCa" %in% models) {
      if(resultType != "county") {
        designResPopAgg[[1]]$Q10 = getSubLevelResults(designResPopAgg[[1]]$Q10, designResPopAgg[[2]]$Q10, truth$urban)
        designResPopAgg[[1]]$Q50 = getSubLevelResults(designResPopAgg[[1]]$Q50, designResPopAgg[[2]]$Q50, truth$urban)
        designResPopAgg[[1]]$Q90 = getSubLevelResults(designResPopAgg[[1]]$Q90, designResPopAgg[[2]]$Q90, truth$urban)
        designResPopAgg[[1]]$mean = getSubLevelResults(designResPopAgg[[1]]$mean, designResPopAgg[[2]]$mean, truth$urban)
        designResPopAgg[[1]]$stddev = getSubLevelResults(designResPopAgg[[1]]$stddev, designResPopAgg[[2]]$stddev, truth$urban)
      }
      else {
        designResPopAgg[[1]]$Q10 = getSubLevelResults(designResPopAgg[[1]]$Q10)
        designResPopAgg[[1]]$Q50 = getSubLevelResults(designResPopAgg[[1]]$Q50)
        designResPopAgg[[1]]$Q90 = getSubLevelResults(designResPopAgg[[1]]$Q90)
        designResPopAgg[[1]]$mean = getSubLevelResults(designResPopAgg[[1]]$mean)
        designResPopAgg[[1]]$stddev = getSubLevelResults(designResPopAgg[[1]]$stddev)
      }
    }
    
    # not including urban effect, modified to be debiased using marginal rather than conditional effect as prediction
    if("BYM2 uCa'" %in% models) {
      designResNoUrbModPopAgg[[1]]$Q10 = getSubLevelResults(designResNoUrbModPopAgg[[1]]$Q10)
      designResNoUrbModPopAgg[[1]]$Q50 = getSubLevelResults(designResNoUrbModPopAgg[[1]]$Q50)
      designResNoUrbModPopAgg[[1]]$Q90 = getSubLevelResults(designResNoUrbModPopAgg[[1]]$Q90)
      designResNoUrbModPopAgg[[1]]$mean = getSubLevelResults(designResNoUrbModPopAgg[[1]]$mean)
      designResNoUrbModPopAgg[[1]]$stddev = getSubLevelResults(designResNoUrbModPopAgg[[1]]$stddev)
    }
    
    # including urban effect, modified to be debiased using marginal rather than conditional effect as prediction
    if("BYM2 UCa'" %in% models) {
      if(resultType != "county") {
        designResModPopAgg[[1]]$Q10 = getSubLevelResults(designResModPopAgg[[1]]$Q10, designResModPopAgg[[2]]$Q10, truth$urban)
        designResModPopAgg[[1]]$Q50 = getSubLevelResults(designResModPopAgg[[1]]$Q50, designResModPopAgg[[2]]$Q50, truth$urban)
        designResModPopAgg[[1]]$Q90 = getSubLevelResults(designResModPopAgg[[1]]$Q90, designResModPopAgg[[2]]$Q90, truth$urban)
        designResModPopAgg[[1]]$mean = getSubLevelResults(designResModPopAgg[[1]]$mean, designResModPopAgg[[2]]$mean, truth$urban)
        designResModPopAgg[[1]]$stddev = getSubLevelResults(designResModPopAgg[[1]]$stddev, designResModPopAgg[[2]]$stddev, truth$urban)
      }
      else {
        designResModPopAgg[[1]]$Q10 = getSubLevelResults(designResModPopAgg[[1]]$Q10)
        designResModPopAgg[[1]]$Q50 = getSubLevelResults(designResModPopAgg[[1]]$Q50)
        designResModPopAgg[[1]]$Q90 = getSubLevelResults(designResModPopAgg[[1]]$Q90)
        designResModPopAgg[[1]]$mean = getSubLevelResults(designResModPopAgg[[1]]$mean)
        designResModPopAgg[[1]]$stddev = getSubLevelResults(designResModPopAgg[[1]]$stddev)
      }
    }
    
    for(i in c(1:maxDataSets)) { # for problem fitting mercerSRS for SRS sampling, tausq=0
      # for(i in 1:100) {
      if((i %% printIEvery == 0) || (i == 1))
        print(i)
      resultName = paste0(resultType, "Results")
      if(resultType == "EA")
        resultName = "eaResults"
      
      # convert results to the desired aggregation level
      if("Direct" %in% models)
        directEsti = getSubLevelResults(directEst[[i]])
      if("Naive" %in% models)
        naivei = getSubLevelResults(naive[[i]])
      if("Smoothed Direct" %in% models)
        merceri = getSubLevelResults(mercer[[i]])
      if(resultType != "county") {
        # if("SPDE uc" %in% models)
        #   spdeNoUrbClusti = spdeNoUrbClust[[resultName]][[i]][as.numeric(as.character(truth[[resultType]])),]
        # if("SPDE uC" %in% models)
        #   spdeNoUrbi = spdeNoUrb[[resultName]][[i]][as.numeric(as.character(truth[[resultType]])),]
        # if("SPDE Uc" %in% models)
        #   spdeNoClusti = spdeNoClust[[resultName]][[i]][as.numeric(as.character(truth[[resultType]])),]
        # if("SPDE UC" %in% models)
        #   spdei = spde[[resultName]][[i]][as.numeric(as.character(truth[[resultType]])),]
        if("Direct" %in% models)
          directEsti = getSubLevelResults(directEst[[i]])
      } else {
        # if("SPDE uc" %in% models)
        #   spdeNoUrbClusti = spdeNoUrbClust[[resultName]][[i]]
        # if("SPDE uC" %in% models)
        #   spdeNoUrbi = spdeNoUrb[[resultName]][[i]]
        # if("SPDE Uc" %in% models)
        #   spdeNoClusti = spdeNoClust[[resultName]][[i]]
        # if("SPDE UC" %in% models)
        #   spdei = spde[[resultName]][[i]]
      }
      
      if(resultType == "EA") {
        # set first row of spde results to be the EA index
        # if("SPDE uc" %in% models) {
        #   spdeNoUrbClusti[[resultType]] = 1:nrow(spdeNoUrbClusti)
        #   
        #   whichName = which(names(spdeNoUrbClusti) == "EA")
        #   spdeNoUrbClusti = cbind(spdeNoUrbClusti[,whichName], spdeNoUrbClusti[,-whichName])
        #   
        #   names(spdeNoUrbClusti)[1] = "EA"
        # }
        # if("SPDE uC" %in% models) {
        #   spdeNoUrbi[[resultType]] = 1:nrow(spdeNoUrbi)
        #   
        #   whichName = which(names(spdeNoUrbi) == "EA")
        #   spdeNoUrbi = cbind(spdeNoUrbi[,whichName], spdeNoUrbi[,-whichName])
        #   
        #   names(spdeNoUrbi)[1] = "EA"
        # }
        # if("SPDE Uc" %in% models) {
        #   spdeNoClusti[[resultType]] = 1:nrow(spdeNoClusti)
        #   
        #   whichName = which(names(spdeNoClusti) == "EA")
        #   spdeNoClusti = cbind(spdeNoClusti[,whichName], spdeNoClusti[,-whichName])
        #   
        #   names(spdeNoClusti)[1] = "EA"
        # }
        # if("SPDE UC" %in% models) {
        #   spdei[[resultType]] = 1:nrow(spdei)
        #   
        #   whichName = which(names(spdei) == "EA")
        #   spdei = cbind(spdei[,whichName], spdei[,-whichName])
        #   
        #   names(spdei)[1] = "EA"
        # }
      }
      
      # for spde results, modify the name of the results
      # modify result row and column table names according to aggregation level
      if(resultType == "county") {
        # if("SPDE uc" %in% models) {
        #   whichName = which(names(spdeNoUrbClusti) == "admin1")
        #   names(spdeNoUrbClusti)[whichName] = resultType
        # }
        # if("SPDE uC" %in% models) {
        #   whichName = which(names(spdeNoUrbi) == "admin1")
        #   names(spdeNoUrbi)[whichName] = resultType
        # }
        # if("SPDE Uc" %in% models) {
        #   whichName = which(names(spdeNoClusti) == "admin1")
        #   names(spdeNoClusti)[whichName] = resultType
        # }
        # if("SPDE UC" %in% models) {
        #   # with urban effect:
        #   whichName = which(names(spdei) == "admin1")
        #   names(spdei)[whichName] = resultType
        # }
      }
      
      if(resultType == "pixel") {
        # set first row of spde results to be the pixel index
        # if("SPDE uc" %in% models) {
        #   spdeNoUrbClusti[[resultType]] = truth$pixel
        #   
        #   whichName = which(names(spdeNoUrbClusti) == "pixel")
        #   spdeNoUrbClusti = cbind(spdeNoUrbClusti[,whichName], spdeNoUrbClusti[,-whichName])
        #   
        #   names(spdeNoUrbClusti)[1] = "pixel"
        # }
        # if("SPDE uC" %in% models) {
        #   spdeNoUrbi[[resultType]] = truth$pixel
        #   
        #   whichName = which(names(spdeNoUrbi) == "pixel")
        #   spdeNoUrbi = cbind(spdeNoUrbi[,whichName], spdeNoUrbi[,-whichName])
        #   
        #   names(spdeNoUrbi)[1] = "pixel"
        # }
        # if("SPDE Uc" %in% models) {
        #   spdeNoClusti[[resultType]] = truth$pixel
        #   
        #   whichName = which(names(spdeNoClusti) == "pixel")
        #   spdeNoClusti = cbind(spdeNoClusti[,whichName], spdeNoClusti[,-whichName])
        #   
        #   names(spdeNoClusti)[1] = "pixel"
        # }
        # if("SPDE UC" %in% models) {
        #   spdei[[resultType]] = truth$pixel
        #   
        #   whichName = which(names(spdei) == "pixel")
        #   spdei = cbind(spdei[,whichName], spdei[,-whichName])
        #   
        #   names(spdei)[1] = "pixel"
        # }
      }
      
      # TODO: fix the below code block
      # change names of table variables in spde model with no urban effect to reflect that
      # if("SPDE uc" %in% models)
      #   names(spdeNoUrbClusti)[2:6] = paste0(names(spdeNoUrbClusti)[2:6], " uc")
      # if("SPDE uC" %in% models)
      #   names(spdeNoUrbi)[2:6] = paste0(names(spdeNoUrbi)[2:6], "NoUrb")
      # if("SPDE uC" %in% models)
      #   names(spdeNoUrbi)[2:6] = paste0(names(spdeNoUrbi)[2:6], "NoUrb")
      # if("SPDE uC" %in% models)
      #   names(spdeNoUrbi)[2:6] = paste0(names(spdeNoUrbi)[2:6], "NoUrb")
      
      if("Direct" %in% models) {
        allres = merge(truth, directEsti, by=resultType)
        colnames(allres) = c(resultType, "truth", paste(colnames(allres)[3:8], "Direct", sep=""))
      } else {
        stop("direct estimates must be included at this point in order to name the estimate table columns")
      }
      if("Naive" %in% models)
        allres = merge(allres, naivei, by=resultType)
      if("Smoothed Direct" %in% models)
        allres = merge(allres, merceri, by=resultType)
      # if("SPDE uc" %in% models)
      #   allres = merge(allres, spdeNoUrbClusti, by=resultType)
      # if("SPDE uC" %in% models)
      #   allres = merge(allres, spdeNoUrbi, by=resultType)
      # if("SPDE Uc" %in% models)
      #   allres = merge(allres, spdeNoClusti, by=resultType)
      # if("SPDE UC" %in% models)
      #   allres = merge(allres, spdei, by=resultType)
      
      # set whether or not to calculate scores on logit scale depending on result type
      thisTruth = allres$truth
      
      # if(is.null(useLogit)) {
      #   useLogit = FALSE
      #   if(resultType != "EA" && resultType != "pixel") {
      #     thisTruth = logit(thisTruth)
      #     useLogit=TRUE
      #   }
      # } else if(useLogit == TRUE) {
      #   thisTruth = logit(thisTruth)
      # }
      
      if("Direct" %in% models) {
        my.scoresdirect = getScores(thisTruth, numChildren, allres$logit.estDirect, allres$var.estDirect, nsim=nsim)
        scoresDirect <- rbind(scoresDirect,
                              cbind(data.frame(dataset=i, region=allres[[resultType]]), my.scoresdirect))
      }
      if("Naive" %in% models) {
        my.scoresnaive = getScores(thisTruth, numChildren, allres$logit.est, allres$var.est, nsim=nsim)
        scoresNaive <- rbind(scoresNaive,
                             cbind(data.frame(dataset=i, region=allres[[resultType]]), my.scoresnaive))
      }
      if("Smoothed Direct" %in% models) {
        my.scoresmercer = getScores(thisTruth, numChildren, allres$logit.est.mercer, allres$var.est.mercer, nsim=nsim)
        scoresMercer <- rbind(scoresMercer,
                              cbind(data.frame(dataset=i, region=allres[[resultType]]), my.scoresmercer))
      }
      if("BYM2 ucA" %in% models) {
        my.scoresbymNoUrbClust = getScores(thisTruth, numChildren, designResNoUrbClust[[1]]$mean[,i], (designResNoUrbClust[[1]]$stddev[,i])^2, nsim=nsim)
        scoresBYMNoUrbClust <- rbind(scoresBYMNoUrbClust,
                                     cbind(data.frame(dataset=i, region=allres[[resultType]]), my.scoresbymNoUrbClust))
      }
      if("BYM2 uCA" %in% models) {
        my.scoresbymNoUrb = getScores(thisTruth, numChildren, designResNoUrb[[1]]$mean[,i], (designResNoUrb[[1]]$stddev[,i])^2, nsim=nsim)
        scoresBYMNoUrb <- rbind(scoresBYMNoUrb,
                                cbind(data.frame(dataset=i, region=allres[[resultType]]), my.scoresbymNoUrb))
      }
      if("BYM2 uCa'" %in% models) {
        my.scoresbymNoUrbMod = getScores(thisTruth, numChildren, designResNoUrbMod[[1]]$mean[,i], (designResNoUrbMod[[1]]$stddev[,i])^2, nsim=nsim)
        scoresBYMNoUrbMod <- rbind(scoresBYMNoUrbMod,
                                   cbind(data.frame(dataset=i, region=allres[[resultType]]), my.scoresbymNoUrbMod))
      }
      if("BYM2 UcA" %in% models) {
        my.scoresbymNoClust = getScores(thisTruth, numChildren, designResNoClust[[1]]$mean[,i], (designResNoClust[[1]]$stddev[,i])^2, nsim=nsim)
        scoresBYMNoClust <- rbind(scoresBYMNoClust,
                                  cbind(data.frame(dataset=i, region=allres[[resultType]]), my.scoresbymNoClust))
      }
      if("BYM2 UCA" %in% models) {
        my.scoresbym = getScores(thisTruth, numChildren, designRes[[1]]$mean[,i], (designRes[[1]]$stddev[,i])^2, nsim=nsim)
        scoresBYM <- rbind(scoresBYM,
                           cbind(data.frame(dataset=i, region=allres[[resultType]]), my.scoresbym))
      }
      if("BYM2 UCA'" %in% models) {
        my.scoresbymMod = getScores(thisTruth, numChildren, designResMod[[1]]$mean[,i], (designResMod[[1]]$stddev[,i])^2, nsim=nsim)
        scoresBYMMod <- rbind(scoresBYMMod,
                              cbind(data.frame(dataset=i, region=allres[[resultType]]), my.scoresbymMod))
      }
      if("BYM2 uca" %in% models) {
        my.scoresbymNoUrbClustPopAgg = getScores(thisTruth, numChildren, designResNoUrbClustPopAgg[[1]]$mean[,i], (designResNoUrbClustPopAgg[[1]]$stddev[,i])^2, nsim=nsim)
        scoresBYMNoUrbClustPopAgg <- rbind(scoresBYMNoUrbClustPopAgg,
                                           cbind(data.frame(dataset=i, region=allres[[resultType]]), my.scoresbymNoUrbClustPopAgg))
      }
      if("BYM2 uCa" %in% models) {
        my.scoresbymNoUrbPopAgg = getScores(thisTruth, numChildren, designResNoUrbPopAgg[[1]]$mean[,i], (designResNoUrbPopAgg[[1]]$stddev[,i])^2, nsim=nsim)
        scoresBYMNoUrbPopAgg <- rbind(scoresBYMNoUrbPopAgg,
                                      cbind(data.frame(dataset=i, region=allres[[resultType]]), my.scoresbymNoUrbPopAgg))
      }
      if("BYM2 uCa'" %in% models) {
        my.scoresbymNoUrbModPopAgg = getScores(thisTruth, numChildren, designResNoUrbModPopAgg[[1]]$mean[,i], (designResNoUrbModPopAgg[[1]]$stddev[,i])^2, nsim=nsim)
        scoresBYMNoUrbModPopAgg <- rbind(scoresBYMNoUrbModPopAgg,
                                         cbind(data.frame(dataset=i, region=allres[[resultType]]), my.scoresbymNoUrbModPopAgg))
      }
      if("BYM2 Uca" %in% models) {
        my.scoresbymNoClustPopAgg = getScores(thisTruth, numChildren, designResNoClustPopAgg[[1]]$mean[,i], (designResNoClustPopAgg[[1]]$stddev[,i])^2, nsim=nsim)
        scoresBYMNoClustPopAgg <- rbind(scoresBYMNoClustPopAgg,
                                        cbind(data.frame(dataset=i, region=allres[[resultType]]), my.scoresbymNoClustPopAgg))
      }
      if("BYM2 UCa" %in% models) {
        my.scoresbymPopAgg = getScores(thisTruth, numChildren, designResPopAgg[[1]]$mean[,i], (designResPopAgg[[1]]$stddev[,i])^2, nsim=nsim)
        scoresBYMPopAgg <- rbind(scoresBYMPopAgg,
                                 cbind(data.frame(dataset=i, region=allres[[resultType]]), my.scoresbymPopAgg))
      }
      if("BYM2 UCa'" %in% models) {
        my.scoresbymModPopAgg = getScores(thisTruth, numChildren, designResModPopAgg[[1]]$mean[,i], (designResModPopAgg[[1]]$stddev[,i])^2, nsim=nsim)
        scoresBYMModPopAgg <- rbind(scoresBYMModPopAgg,
                                    cbind(data.frame(dataset=i, region=allres[[resultType]]), my.scoresbymModPopAgg))
      }
      if("SPDE uC" %in% models) {
        # stop("determine if the spde code should compute all of these directly")
        # my.biasspdeNoUrb = bias(thisTruth, allres$logit.est.spdeNoUrb, logit=useLogit, my.var=allres$var.est.spdeNoUrb)
        # my.msespdeNoUrb = mse(thisTruth, allres$logit.est.spdeNoUrb, logit=useLogit, my.var=allres$var.est.spdeNoUrb)
        # my.dssspdeNoUrb = dss(thisTruth, allres$logit.est.spdeNoUrb, allres$var.est.spdeNoUrb)
        # # the below line used to be commented for some reason
        # my.crpsspdeNoUrb = crpsNormal(thisTruth, allres$logit.est.spdeNoUrb, allres$var.est.spdeNoUrb, logit=useLogit, n=numChildren)
        # my.crpsspdeNoUrb = allres$crps.spdeNoUrb
        # my.coveragespdeNoUrb = coverage(thisTruth, allres$lower.spdeNoUrb, allres$upper.spdeNoUrb, logit=useLogit)
        # my.lengthspdeNoUrb = intervalWidth(allres$lower.spdeNoUrb, allres$upper.spdeNoUrb, logit=useLogit)
        # scoresSPDENoUrb <- rbind(scoresSPDENoUrb,
        #                          cbind(data.frame(dataset=i, region=allres[[resultType]]), my.scoresspdeNoUrb))
        # scoresSPDENoUrb = rbind(scoresSPDENoUrb # TODO: fix this
      }
      if("SPDE UC" %in% models) {
        # stop("determine if the spde code should compute all of these directly")
        # my.biasspde = bias(thisTruth, allres$logit.est.spde, logit=useLogit, my.var=allres$var.est.spde)
        # my.msespde = mse(thisTruth, allres$logit.est.spde, logit=useLogit, my.var=allres$var.est.spde)
        # my.dssspde = dss(thisTruth, allres$logit.est.spde, allres$var.est.spde)
        # # the below line needs to be commented for some reason
        # my.crpsspde = crpsNormal(thisTruth, allres$logit.est.spde, allres$var.est.spde, logit=useLogit, n=numChildren)
        # my.crpsspde = allres$crps.spde
        # my.coveragespde = coverage(thisTruth, allres$lower.spde, allres$upper.spde, logit=useLogit)
        # my.lengthspde = intervalWidth(allres$lower.spde, allres$upper.spde, logit=useLogit)
        # scoresSPDE <- rbind(scoresSPDE,
        #                     cbind(data.frame(dataset=i, region=allres[[resultType]]), my.scoresspde))
      }
    }
    
    # save progress
    runId = paste0("Beta-1.75margVar", round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 
                   "HHoldVar0urbanOverSamplefrac0", strictPriorText, testText, bigText, sampling, 
                   "models", do.call("paste0", as.list(modelsI)), "nsim", nsim, "MaxDataSetI", maxDataSets)
    
    # first collect all the results. Save everything except for the postprocessing arguments: 
    # produceFigures, digits
    objectNames = ls()
    objectNames = objectNames[-match(c("produceFigures", "xtable.args", "tableFormat", "colScale", 
                                       "colUnits", "colDigits"), objectNames)]
    save(list=objectNames, file=paste0("scoresTemp", runId, ".RData"))
  }
  else {
    # in this case, we have already computed the results so just load them into the environment
    print("Loading results...")
    temp = doFancyTables
    if(loadResults)
      load(paste0("scores", runId, ".RData"))
    else if(loadTempProgress) {
      temp = loadTempProgress
      load(paste0("scoresTemp", runId, ".RData"))
      loadTempProgress = temp
    }
    
    doFancyTables = temp
    allNames = c("Naive", "Direct", "Smoothed Direct", "BYM2 ucA", "BYM2 uCA", "BYM2 uCA'", "BYM2 UcA", "BYM2 UCA", "BYM2 UCA'", 
                 "BYM2 uca", "BYM2 uCa", "BYM2 uCa'", "BYM2 Uca", "BYM2 UCa", "BYM2 UCa'", 
                 "SPDE uc", "SPDE uC", "SPDE Uc", "SPDE UC", "SPDE uC'", "SPDE UC'")
    allNamesBinomial = paste0(allNames, " Bin.")
    allModels = allNames
    models = allModels[modelsI]
    
    # also, don't resave these results that we've already saved
    if(!loadTempProgress)
      saveResults = FALSE
    else
      saveResults = TRUE
  }
  
  # final table: 
  if("Naive" %in% models)
    naive = apply(scoresNaive[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("Direct" %in% models)
    direct = apply(scoresDirect[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("Smoothed Direct" %in% models)
    mercer = apply(scoresMercer[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("BYM2 ucA" %in% models)
    bymNoUrbClust = apply(scoresBYMNoUrbClust[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("BYM2 uCA" %in% models)
    bymNoUrb = apply(scoresBYMNoUrb[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("BYM2 uCa'" %in% models)
    bymNoUrbMod = apply(scoresBYMNoUrbMod[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("BYM2 UcA" %in% models)
    bymNoClust = apply(scoresBYMNoClust[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("BYM2 UCA" %in% models)
    bym = apply(scoresBYM[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("BYM2 UCA'" %in% models)
    bymMod = apply(scoresBYMMod[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("BYM2 uca" %in% models)
    bymNoUrbClustPopAgg = apply(scoresBYMNoUrbClustPopAgg[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("BYM2 uCa" %in% models)
    bymNoUrbPopAgg = apply(scoresBYMNoUrbPopAgg[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("BYM2 uCa'" %in% models)
    bymNoUrbModPopAgg = apply(scoresBYMNoUrbModPopAgg[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("BYM2 Uca" %in% models)
    bymNoClustPopAgg = apply(scoresBYMNoClustPopAgg[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("BYM2 UCa" %in% models)
    bymPopAgg = apply(scoresBYMPopAgg[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("BYM2 UCa'" %in% models)
    bymModPopAgg = apply(scoresBYMModPopAgg[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("SPDE uc" %in% models) {
    theseNames = names(spdeNoUrbClust)
    namesI = grepl(tolower(resultType), tolower(theseNames))
    first = TRUE
    for(i in 1:length(theseNames)) {
      if(namesI[i]) {
        if(first) {
          spdeNoUrbClustScores = apply(spdeNoUrbClust[[i]], 2, mean)
          first = FALSE
        }
        else {
          spdeNoUrbClustScores = rbind(spdeNoUrbClustScores, apply(spdeNoUrbClust[[i]], 2, mean))
        }
      }
    }
  }
  if("SPDE uC" %in% models || "SPDE uC'" %in% models) {
    theseNames = names(spdeNoUrb)
    namesI = grepl(tolower(resultType), tolower(theseNames))
    first = TRUE
    for(i in 1:length(theseNames)) {
      if(namesI[i]) {
        if(first) {
          spdeNoUrbScores = apply(spdeNoUrb[[i]], 2, mean)
          first = FALSE
        }
        else {
          spdeNoUrbScores = rbind(spdeNoUrbScores, apply(spdeNoUrb[[i]], 2, mean))
        }
      }
    }
  }
  if("SPDE Uc" %in% models) {
    theseNames = names(spdeNoClust)
    namesI = grepl(tolower(resultType), tolower(theseNames))
    first = TRUE
    for(i in 1:length(theseNames)) {
      if(namesI[i]) {
        if(first) {
          spdeNoClustScores = apply(spdeNoClust[[i]], 2, mean)
          first = FALSE
        }
        else {
          spdeNoClustScores = rbind(spdeNoClustScores, apply(spdeNoClust[[i]], 2, mean))
        }
      }
    }
  }
  if("SPDE UC" %in% models || "SPDE UC'" %in% models) {
    theseNames = names(spde)
    namesI = grepl(tolower(resultType), tolower(theseNames))
    first = TRUE
    for(i in 1:length(theseNames)) {
      if(namesI[i]) {
        if(first) {
          spdeScores = apply(spde[[i]], 2, mean)
          first = FALSE
        }
        else {
          spdeScores = rbind(spdeScores, apply(spde[[i]], 2, mean))
        }
      }
    }
  }
  idx = 1:9
  tab = c()
  if("Naive" %in% models)
    tab = rbind(tab, c(naive[idx]))
  if("Direct" %in% models)
    tab = rbind(tab, c(direct[idx]))
  if("Smoothed Direct" %in% models)
    tab = rbind(tab, c(mercer[idx]))
  if("BYM2 ucA" %in% models)
    tab = rbind(tab, c(bymNoUrbClust[idx]))
  if("BYM2 uCA" %in% models)
    tab = rbind(tab, c(bymNoUrb[idx]))
  if("BYM2 uCa'" %in% models)
    tab = rbind(tab, c(bymNoUrbMod[idx]))
  if("BYM2 UcA" %in% models)
    tab = rbind(tab, c(bymNoClust[idx]))
  if("BYM2 UCA" %in% models)
    tab = rbind(tab, c(bym[idx]))
  if("BYM2 UCA'" %in% models)
    tab = rbind(tab, c(bymMod[idx]))
  if("BYM2 uca" %in% models)
    tab = rbind(tab, c(bymNoUrbClustPopAgg[idx]))
  if("BYM2 uCa" %in% models)
    tab = rbind(tab, c(bymNoUrbPopAgg[idx]))
  if("BYM2 uCa'" %in% models)
    tab = rbind(tab, c(bymNoUrbModPopAgg[idx]))
  if("BYM2 Uca" %in% models)
    tab = rbind(tab, c(bymNoClustPopAgg[idx]))
  if("BYM2 UCa" %in% models)
    tab = rbind(tab, c(bymPopAgg[idx]))
  if("BYM2 UCa'" %in% models)
    tab = rbind(tab, c(bymModPopAgg[idx]))
  # SPDE models are in format 2
  # if("SPDE uc" %in% models)
  #   tab = rbind(tab, c(spdeNoUrbClust[idx]))
  # if("SPDE uC" %in% models)
  #   tab = rbind(tab, c(spdeNoUrb[idx]))
  # if("SPDE Uc" %in% models)
  #   tab = rbind(tab, c(spdeNoClust[idx]))
  # if("SPDE UC" %in% models)
  #   tab = rbind(tab, c(spde[idx]))
  colnames(tab) = c("Bias", "Var", "MSE", "CRPS", "CRPS Bin.", "80\\% Cvg", "80\\% Cvg Bin.", "CI Width", "CI Width Bin.")
  finalNames = allNames[modelsI]
  finalNamesBinomial = allNamesBinomial[modelsI]
  
  if(tableFormat == "1")
    rownames(tab) = finalNames
  else {
    # separate out the binomial and non-binomial models into separate rows
    binomialColumns = c(1:3, 5, 7, 9)
    otherColumns = c(1:3, 4, 6, 8)
    binomialTable = tab[, binomialColumns]
    otherTable = tab[, otherColumns]
    tab = otherTable[numeric(0),]
    thisFinalNames = c()
    
    for(i in 1:nrow(otherTable)) {
      if(includeBVarResults) {
        tab = rbind(tab, 
                    otherTable[i,], 
                    binomialTable[i,])
        thisFinalNames = c(thisFinalNames, finalNames[i], finalNamesBinomial[i])
      }
      else {
        tab = rbind(tab, 
                    otherTable[i,])
        thisFinalNames = c(thisFinalNames, finalNames[i])
      }
    }
    
    # add in SPDE models if necessary
    if("SPDE uc" %in% models) {
      if(continuousSPDEonly) {
        thisFinalNames = c(thisFinalNames, "SPDE uc")
        tab = rbind(tab, spdeNoUrbClustScores[1,])
      }
      else {
        thisFinalNames = c(thisFinalNames, paste0("SPDE uc", c(" Cts.", " Discrete", " Exact")))
        tab = rbind(tab, spdeNoUrbClustScores)
      }
    }
    if("SPDE uC" %in% models || "SPDE uC'" %in% models) {
      if(continuousSPDEonly) {
        thisFinalNames = c(thisFinalNames, "SPDE uC", "SPDE uC'")
        tab = rbind(tab, spdeNoUrbScores[2:1,])
      }
      else {
        thisFinalNames = c(thisFinalNames, paste0("SPDE uC", c(" Cts.", " Discrete", " Exact")))
        tab = rbind(tab, spdeNoUrbScores)
      }
    }
    if("SPDE Uc" %in% models) {
      if(continuousSPDEonly) {
        thisFinalNames = c(thisFinalNames, "SPDE Uc")
        tab = rbind(tab, spdeNoClustScores[1,])
      }
      else {
        thisFinalNames = c(thisFinalNames, paste0("SPDE Uc", c(" Cts.", " Discrete", " Exact")))
        tab = rbind(tab, spdeNoClustScores)
      }
    }
    if("SPDE UC" %in% models || "SPDE UC'" %in% models) {
      if(continuousSPDEonly) {
        thisFinalNames = c(thisFinalNames, "SPDE UC", "SPDE UC'")
        tab = rbind(tab, spdeScores[2:1,])
      }
      else {
        thisFinalNames = c(thisFinalNames, paste0("SPDE UC", c(" Cts.", " Discrete", " Exact")))
        tab = rbind(tab, spdeScores)
      }
    }
    # if("SPDE uC'" %in% models) {
    #   if(continuousSPDEonly) {
    #     thisFinalNames = c(thisFinalNames, "SPDE uC'")
    #     tab = rbind(tab, spdeNoUrbModScores[1,])
    #   }
    #   else {
    #     thisFinalNames = c(thisFinalNames, paste0("SPDE uC'", c(" Cts.", " Discrete", " Exact")))
    #     tab = rbind(tab, spdeNoUrbModScores)
    #   }
    # }
    # if("SPDE UC" %in% models) {
    #   if(continuousSPDEonly) {
    #     thisFinalNames = c(thisFinalNames, "SPDE UC'")
    #     tab = rbind(tab, spdeModScores[1,])
    #   }
    #   else {
    #     thisFinalNames = c(thisFinalNames, paste0("SPDE UC'", c(" Cts.", " Discrete", " Exact")))
    #     tab = rbind(tab, spdeModScores)
    #   }
    # }
    
    rownames(tab) = thisFinalNames
  }
  
  # round the columns of tab and modify the column names to include the scale
  unroundedTab = tab
  for(i in 1:ncol(tab)) {
    tab[,i] = as.numeric(round(tab[,i] * colScale[i], digits=colDigits[i]))
    colnames(tab)[i] = paste0(colnames(tab)[i], colUnits[i])
    unroundedTab[,i] = as.numeric(unroundedTab[,i] * colScale[i])
    colnames(unroundedTab)[i] = paste0(colnames(unroundedTab)[i], colUnits[i])
  }
  
  # remove "Direct Bin." model, since direct estimates already account for Binomial variation
  # if("Direct" %and% models) {
  #   rowI = thisFinalNames == "Direct Bin."
  #   tab = tab[!rowI,]
  # }
  
  # remove redundant rows of BYM2 model
  require(stringr)
  modelTypes = word(models, 1)
  modelTypes[modelTypes == "Smoothed"] = "Smoothed Direct"
  uniqueModelTypes = unique(modelTypes)
  modelTypeGroups = lapply(uniqueModelTypes, function(x) {(1:length(modelTypes))[modelTypes == x]})
  modelVariations = word(models, 2)
  modelVariations[is.na(modelVariations)] = ""
  modelVariations[modelVariations == "Direct"] = ""
  if("BYM2 ucA" %in% models)
    models[models == "BYM2 ucA"] = "BYM2 uc"
  if("BYM2 uCA" %in% models)
    models[models == "BYM2 uCA"] = "BYM2 uC"
  if("BYM2 uCA'" %in% models)
    models[models == "BYM2 uCA'"] = "BYM2 uC'"
  if("BYM2 uca" %in% models) {
    tab = tab[models != "BYM2 uca",]
    unroundedTab = unroundedTab[models != "BYM2 uca",]
    models = models[models != "BYM2 uca"]
  }
  if("BYM2 uCa" %in% models) {
    tab = tab[models != "BYM2 uCa",]
    unroundedTab = unroundedTab[models != "BYM2 uCa",]
    models = models[models != "BYM2 uCa"]
  }
  if("BYM2 uCa'" %in% models) {
    tab = tab[models != "BYM2 uCa'",]
    unroundedTab = unroundedTab[models != "BYM2 uCa'",]
    models = models[models != "BYM2 uCa'"]
  }
  if(identical(modelsI, 1:21)) {
    models = models[c(1:12, 13, 14, 17, 15, 16, 18)]
  }
  rownames(tab) = models
  rownames(unroundedTab) = models
  
  # recalculate model types and variations
  modelTypes = word(models, 1)
  modelTypes[modelTypes == "Smoothed"] = "Smoothed Direct"
  uniqueModelTypes = unique(modelTypes)
  modelTypeGroups = lapply(uniqueModelTypes, function(x) {(1:length(modelTypes))[modelTypes == x]})
  modelVariations = word(models, 2)
  modelVariations[is.na(modelVariations)] = ""
  modelVariations[modelVariations == "Direct"] = ""
  
  if(!doFancyTables && printScoreTable)
    print(do.call("xtable", c(list(tab), xtable.args)), 
          include.colnames=TRUE,
          hline.after=0, 
          math.style.exponents=TRUE, 
          sanitize.text.function=function(x){x})
  
  if(doFancyTables && printScoreTable) {
    require(stringr)
    require(dplyr)
    require(kableExtra)
    
    options(knitr.table.format = "latex")
    
    # bold the best entries of each column, italicize worst entries of each column
    centers = c(rep(0, 4), 80, 0)
    columnBest = apply(cbind(abs(tab[,1]), tab[,2:4], abs(tab[,5]-80), tab[,6]), 2, min)
    columnWorst = apply(cbind(abs(tab[,1]), tab[,2:4], abs(tab[,5]-80), tab[,6]), 2, max)
    dat = data.table(tab)
    test = dat %>% mutate(Bias = cell_spec(tab[,1], "latex", bold=abs(tab[,1] - centers[1]) <= columnBest[1], italic = abs(tab[,1] - centers[1]) >= columnWorst[1], 
                                           monospace=FALSE, underline=FALSE, strikeout=FALSE), 
                          Var = cell_spec(tab[,2], "latex", bold=abs(tab[,2] - centers[2]) <= columnBest[2], italic = abs(tab[,2] - centers[2]) >= columnWorst[2], 
                                          monospace=FALSE, underline=FALSE, strikeout=FALSE), 
                          MSE = cell_spec(tab[,3], "latex", bold=abs(tab[,3] - centers[3]) <= columnBest[3], italic = abs(tab[,3] - centers[3]) >= columnWorst[3], 
                                          monospace=FALSE, underline=FALSE, strikeout=FALSE), 
                          CRPS = cell_spec(tab[,4], "latex", bold=abs(tab[,4] - centers[4]) <= columnBest[4], italic = abs(tab[,4] - centers[4]) >= columnWorst[4], 
                                           monospace=FALSE, underline=FALSE, strikeout=FALSE), 
                          CVG = cell_spec(tab[,5], "latex", bold=abs(tab[,5] - centers[5]) <= columnBest[5], italic = abs(tab[,5] - centers[5]) >= columnWorst[5], 
                                          monospace=FALSE, underline=FALSE, strikeout=FALSE), 
                          Width = cell_spec(tab[,6], "latex", bold=abs(tab[,6] - centers[6]) <= columnBest[6], italic = abs(tab[,6] - centers[6]) >= columnWorst[6], 
                                            monospace=FALSE, underline=FALSE, strikeout=FALSE)) %>%
      select(Bias, Var, MSE, CRPS, CVG, Width)
    
    # revert the column names to their true values, set the model variations to be the values in the first column
    colnames(test) = colnames(tab)
    test = cbind(" "=modelVariations, test)
    rownames(test)=NULL
    
    # group the rows by the type of model
    fullTab = test %>%
      kable("latex", escape = F, booktabs = T) %>% kable_styling()
    for(i in 1:length(uniqueModelTypes)) {
      startR = min(modelTypeGroups[[i]])
      endR = max(modelTypeGroups[[i]])
      fullTab = fullTab %>% pack_rows(uniqueModelTypes[i], startR, endR, latex_gap_space = "2em")
    }
    if(printScoreTable)
      print(fullTab)
  }
  
  ## append parameter tables from each smoothing model if necessary
  anySmoothingModels = as.numeric("Naive" %in% models) + as.numeric("Direct" %in% models)
  anySmoothingModels = length(models) > anySmoothingModels
  parTab = c()
  if(anySmoothingModels) {
    parRowNames = c()
    if("Smoothed Direct" %in% models) {
      parTab = rbind(parTab, mercerPar)
      parRowNames = c(parRowNames, rep("Smoothed Direct", nrow(mercerPar)))
    }
    if("BYM2 ucA" %in% models || "BYM2 uca" %in% models || "BYM2 uc" %in% models) {
      parTab = rbind(parTab, designResNoUrbClust[[length(designResNoUrbClust)]])
      parRowNames = c(parRowNames, rep("BYM2 uc", nrow(designResNoUrbClust[[length(designResNoUrbClust)]])))
    }
    if("BYM2 uCA" %in% models || "BYM2 uCA'" %in% models || "BYM2 uCa" %in% models || "BYM2 uCa'" %in% models  || "BYM2 uC" %in% models || "BYM2 uC'" %in% models) {
      parTab = rbind(parTab, designResNoUrb[[length(designResNoUrb)]])
      parRowNames = c(parRowNames, rep("BYM2 uC", nrow(designResNoUrb[[length(designResNoUrb)]])))
    }
    if("BYM2 UcA" %in% models || "BYM2 Uca" %in% models || "BYM2 Uca'" %in% models || "BYM2 Ucb'" %in% models) {
      parTab = rbind(parTab, designResNoClust[[length(designResNoClust)]])
      parRowNames = c(parRowNames, rep("BYM2 Uc", nrow(designResNoClust[[length(designResNoClust)]])))
    }
    if("BYM2 UCA" %in% models || "BYM2 UCa" %in% models || "BYM2 UCA'" %in% models || "BYM2 UCa'" %in% models) {
      parTab = rbind(parTab, designRes[[length(designRes)]])
      parRowNames = c(parRowNames, rep("BYM2 UC", nrow(designRes[[length(designRes)]])))
    }
    spdeParIndices = c(1:2, 4:6) # leave out variance and width
    if("SPDE uc" %in% models) {
      thisParTab = matrix(spdeNoUrbClust$interceptSummary[spdeParIndices], nrow=1)
      theseRowNames = "Intercept"
      if(!is.null(spdeNoUrbClust$urbanSummary[spdeParIndices])) {
        thisParTab = rbind(thisParTab, Urban=spdeNoUrbClust$urbanSummary[spdeParIndices])
        theseRowNames = c(theseRowNames, "Urban")
      }
      if(!is.null(spdeNoUrbClust$rangeSummary[spdeParIndices])) {
        thisParTab = rbind(thisParTab, spdeNoUrbClust$rangeSummary[spdeParIndices])
        theseRowNames = c(theseRowNames, "Range")
      }
      if(!is.null(spdeNoUrbClust$varSummary[spdeParIndices])) {
        thisParTab = rbind(thisParTab, spdeNoUrbClust$varSummary[spdeParIndices])
        theseRowNames = c(theseRowNames, "Spatial Var")
      }
      if(!is.null(spdeNoUrbClust$sdSummary[spdeParIndices])) {
        thisParTab = rbind(thisParTab, spdeNoUrbClust$sdSummary[spdeParIndices])
        theseRowNames = c(theseRowNames, "Spatial SD")
      }
      if(!is.null(spdeNoUrbClust$nuggetVarSummary[spdeParIndices])) {
        thisParTab = rbind(thisParTab, spdeNoUrbClust$nuggetVarSummary[spdeParIndices])
        theseRowNames = c(theseRowNames, "Cluster Var")
      }
      if(!is.null(spdeNoUrbClust$nuggetSDSummary[spdeParIndices])) {
        thisParTab = rbind(thisParTab, spdeNoUrbClust$nuggetSDSummary[spdeParIndices])
        theseRowNames = c(theseRowNames, "Cluster SD")
      }
      colnames(thisParTab) = names(parTab)
      rownames(thisParTab) = theseRowNames
      parTab = rbind(parTab, thisParTab)
      parRowNames = c(parRowNames, rep("SPDE uc", nrow(thisParTab)))
    }
    if("SPDE uC" %in% models || "SPDE uC'" %in% models) {
      thisParTab = matrix(spdeNoUrb$interceptSummary[spdeParIndices], nrow=1)
      theseRowNames = "Intercept"
      if(!is.null(spdeNoUrb$urbanSummary[spdeParIndices])) {
        thisParTab = rbind(thisParTab, Urban=spdeNoUrb$urbanSummary[spdeParIndices])
        theseRowNames = c(theseRowNames, "Urban")
      }
      if(!is.null(spdeNoUrb$rangeSummary[spdeParIndices])) {
        thisParTab = rbind(thisParTab, spdeNoUrb$rangeSummary[spdeParIndices])
        theseRowNames = c(theseRowNames, "Range")
      }
      if(!is.null(spdeNoUrb$varSummary[spdeParIndices])) {
        thisParTab = rbind(thisParTab, spdeNoUrb$varSummary[spdeParIndices])
        theseRowNames = c(theseRowNames, "Spatial Var")
      }
      if(!is.null(spdeNoUrb$sdSummary[spdeParIndices])) {
        thisParTab = rbind(thisParTab, spdeNoUrb$sdSummary[spdeParIndices])
        theseRowNames = c(theseRowNames, "Spatial SD")
      }
      if(!is.null(spdeNoUrb$nuggetVarSummary[spdeParIndices])) {
        thisParTab = rbind(thisParTab, spdeNoUrb$nuggetVarSummary[spdeParIndices])
        theseRowNames = c(theseRowNames, "Cluster Var")
      }
      if(!is.null(spdeNoUrb$nuggetSDSummary[spdeParIndices])) {
        thisParTab = rbind(thisParTab, spdeNoUrb$nuggetSDSummary[spdeParIndices])
        theseRowNames = c(theseRowNames, "Cluster SD")
      }
      colnames(thisParTab) = names(parTab)
      rownames(thisParTab) = theseRowNames
      parTab = rbind(parTab, thisParTab)
      parRowNames = c(parRowNames, rep("SPDE uC", nrow(thisParTab)))
    }
    if("SPDE Uc" %in% models) {
      thisParTab = matrix(spdeNoClust$interceptSummary[spdeParIndices], nrow=1)
      theseRowNames = "Intercept"
      if(!is.null(spdeNoClust$urbanSummary[spdeParIndices])) {
        thisParTab = rbind(thisParTab, Urban=spdeNoClust$urbanSummary[spdeParIndices])
        theseRowNames = c(theseRowNames, "Urban")
      }
      if(!is.null(spdeNoClust$rangeSummary[spdeParIndices])) {
        thisParTab = rbind(thisParTab, spdeNoClust$rangeSummary[spdeParIndices])
        theseRowNames = c(theseRowNames, "Range")
      }
      if(!is.null(spdeNoClust$varSummary[spdeParIndices])) {
        thisParTab = rbind(thisParTab, spdeNoClust$varSummary[spdeParIndices])
        theseRowNames = c(theseRowNames, "Spatial Var")
      }
      if(!is.null(spdeNoClust$sdSummary[spdeParIndices])) {
        thisParTab = rbind(thisParTab, spdeNoClust$sdSummary[spdeParIndices])
        theseRowNames = c(theseRowNames, "Spatial SD")
      }
      if(!is.null(spdeNoClust$nuggetVarSummary[spdeParIndices])) {
        thisParTab = rbind(thisParTab, spdeNoClust$nuggetVarSummary[spdeParIndices])
        theseRowNames = c(theseRowNames, "Cluster Var")
      }
      if(!is.null(spdeNoClust$nuggetSDSummary[spdeParIndices])) {
        thisParTab = rbind(thisParTab, spdeNoClust$nuggetSDSummary[spdeParIndices])
        theseRowNames = c(theseRowNames, "Cluster SD")
      }
      colnames(thisParTab) = names(parTab)
      rownames(thisParTab) = theseRowNames
      parTab = rbind(parTab, thisParTab)
      parRowNames = c(parRowNames, rep("SPDE Uc", nrow(thisParTab)))
    }
    if("SPDE UC" %in% models || "SPDE UC'" %in% models) {
      thisParTab = matrix(spde$interceptSummary[spdeParIndices], nrow=1)
      theseRowNames = "Intercept"
      if(!is.null(spde$urbanSummary[spdeParIndices])) {
        thisParTab = rbind(thisParTab, Urban=spde$urbanSummary[spdeParIndices])
        theseRowNames = c(theseRowNames, "Urban")
      }
      if(!is.null(spde$rangeSummary[spdeParIndices])) {
        thisParTab = rbind(thisParTab, spde$rangeSummary[spdeParIndices])
        theseRowNames = c(theseRowNames, "Range")
      }
      if(!is.null(spde$varSummary[spdeParIndices])) {
        thisParTab = rbind(thisParTab, spde$varSummary[spdeParIndices])
        theseRowNames = c(theseRowNames, "Spatial Var")
      }
      if(!is.null(spde$sdSummary[spdeParIndices])) {
        thisParTab = rbind(thisParTab, spde$sdSummary[spdeParIndices])
        theseRowNames = c(theseRowNames, "Spatial SD")
      }
      if(!is.null(spde$nuggetVarSummary[spdeParIndices])) {
        thisParTab = rbind(thisParTab, spde$nuggetVarSummary[spdeParIndices])
        theseRowNames = c(theseRowNames, "Cluster Var")
      }
      if(!is.null(spde$nuggetSDSummary[spdeParIndices])) {
        thisParTab = rbind(thisParTab, spde$nuggetSDSummary[spdeParIndices])
        theseRowNames = c(theseRowNames, "Cluster SD")
      }
      colnames(thisParTab) = names(parTab)
      rownames(thisParTab) = theseRowNames
      parTab = rbind(parTab, thisParTab)
      parRowNames = c(parRowNames, rep("SPDE UC", nrow(thisParTab)))
    }
    # if("SPDE uC'" %in% models) {
    #   thisParTab = matrix(spdeNoUrbMod$interceptSummary[spdeParIndices], nrow=1)
    #   theseRowNames = "Intercept"
    #   if(!is.null(spdeNoUrbMod$urbanSummary[spdeParIndices])) {
    #     thisParTab = rbind(thisParTab, Urban=spdeNoUrbMod$urbanSummary[spdeParIndices])
    #     theseRowNames = c(theseRowNames, "Urban")
    #   }
    #   if(!is.null(spdeNoUrbMod$rangeSummary[spdeParIndices])) {
    #     thisParTab = rbind(thisParTab, spdeNoUrbMod$rangeSummary[spdeParIndices])
    #     theseRowNames = c(theseRowNames, "Range")
    #   }
    #   if(!is.null(spdeNoUrbMod$varSummary[spdeParIndices])) {
    #     thisParTab = rbind(thisParTab, spdeNoUrbMod$varSummary[spdeParIndices])
    #     theseRowNames = c(theseRowNames, "Spatial Var")
    #   }
    #   if(!is.null(spdeNoUrbMod$sdSummary[spdeParIndices])) {
    #     thisParTab = rbind(thisParTab, spdeNoUrbMod$sdSummary[spdeParIndices])
    #     theseRowNames = c(theseRowNames, "Spatial SD")
    #   }
    #   if(!is.null(spdeNoUrbMod$nuggetVarSummary[spdeParIndices])) {
    #     thisParTab = rbind(thisParTab, spdeNoUrbMod$nuggetVarSummary[spdeParIndices])
    #     theseRowNames = c(theseRowNames, "Cluster Var")
    #   }
    #   if(!is.null(spdeNoUrbMod$nuggetSDSummary[spdeParIndices])) {
    #     thisParTab = rbind(thisParTab, spdeNoUrbMod$nuggetSDSummary[spdeParIndices])
    #     theseRowNames = c(theseRowNames, "Cluster SD")
    #   }
    #   colnames(thisParTab) = names(parTab)
    #   rownames(thisParTab) = theseRowNames
    #   parTab = rbind(parTab, thisParTab)
    #   parRowNames = c(parRowNames, rep("SPDE uC'", nrow(thisParTab)))
    # }
    # if("SPDE UC" %in% models) {
    #   thisParTab = matrix(spdeMod$interceptSummary[spdeParIndices], nrow=1)
    #   theseRowNames = "Intercept"
    #   if(!is.null(spdeMod$urbanSummary[spdeParIndices])) {
    #     thisParTab = rbind(thisParTab, Urban=spdeMod$urbanSummary[spdeParIndices])
    #     theseRowNames = c(theseRowNames, "Urban")
    #   }
    #   if(!is.null(spdeMod$rangeSummary[spdeParIndices])) {
    #     thisParTab = rbind(thisParTab, spdeMod$rangeSummary[spdeParIndices])
    #     theseRowNames = c(theseRowNames, "Range")
    #   }
    #   if(!is.null(spdeMod$varSummary[spdeParIndices])) {
    #     thisParTab = rbind(thisParTab, spdeMod$varSummary[spdeParIndices])
    #     theseRowNames = c(theseRowNames, "Spatial Var")
    #   }
    #   if(!is.null(spdeMod$sdSummary[spdeParIndices])) {
    #     thisParTab = rbind(thisParTab, spdeMod$sdSummary[spdeParIndices])
    #     theseRowNames = c(theseRowNames, "Spatial SD")
    #   }
    #   if(!is.null(spdeMod$nuggetVarSummary[spdeParIndices])) {
    #     thisParTab = rbind(thisParTab, spdeMod$nuggetVarSummary[spdeParIndices])
    #     theseRowNames = c(theseRowNames, "Cluster Var")
    #   }
    #   if(!is.null(spdeMod$nuggetSDSummary[spdeParIndices])) {
    #     thisParTab = rbind(thisParTab, spdeMod$nuggetSDSummary[spdeParIndices])
    #     theseRowNames = c(theseRowNames, "Cluster SD")
    #   }
    #   colnames(thisParTab) = names(parTab)
    #   rownames(thisParTab) = theseRowNames
    #   parTab = rbind(parTab, thisParTab)
    #   parRowNames = c(parRowNames, rep("SPDE UC'", nrow(thisParTab)))
    # }
    
    # add the model to the row names, remove the numbers at the end of the duplicated row names, print out the aggregated parameter table
    rownames(parTab) = paste(parRowNames, rownames(parTab))
    for(i in 1:nrow(parTab)) {
      lastCharacter = substr(rownames(parTab)[i], nchar(rownames(parTab)[i]), nchar(rownames(parTab)[i]))
      if(grepl("\\d", lastCharacter))
        rownames(parTab)[i] = substr(rownames(parTab)[i], 1, nchar(rownames(parTab)[i])-1)
    }
    
    if(!doFancyTables && printParTable)
      print(xtable(parTab, digits=2, display=c("s", rep("fg", ncol(parTab)))))
    
    if(doFancyTables && printParTable) {
      # now make the fancy table using kableExtra by grouping the rows by models and model variations
      # fancyParTable = xtable2kable(xtable(parTab, digits=2, display=c("s", rep("fg", ncol(parTab)))))
      # fancyParTable = xtable2kable(xtable(parTab, digits=2, display=c("s", rep("e", ncol(parTab)))))
      
      # determine the model type and model variation groupings
      SmoothDirectI = (1:nrow(parTab))[grepl("Smoothed Direct", rownames(parTab))]
      BYM2I = (1:nrow(parTab))[grepl("BYM2", rownames(parTab))]
      BYM2I = setdiff(BYM2I, SmoothDirectI)
      SPDEI = (1:nrow(parTab))[grepl("SPDE", rownames(parTab))]
      modelTypes = rep("Smoothed Direct", nrow(parTab))
      modelTypes[BYM2I] = "BYM2"
      modelTypes[SPDEI] = "SPDE"
      uniqueModelTypes = unique(modelTypes)
      modelTypeGroups = lapply(uniqueModelTypes, function(x) {(1:length(modelTypes))[modelTypes == x]})
      modelVariations = rep("", nrow(parTab))
      modelVariations[BYM2I] = word(rownames(parTab)[BYM2I], 2)
      modelVariations[SPDEI] = word(rownames(parTab)[SPDEI], 2)
      
      # determine the parameter names
      require("tm")
      parNames = trimws(removeWords(rownames(parTab), c(uniqueModelTypes, unique(modelVariations))))
      
      # round intercept, phi, and range parameters to 3, 3, and 0 digits respectively
      parTab[parNames == "Intercept",] = round(parTab[parNames == "Intercept",], digits=3)
      parTab[parNames == "Phi",] = round(parTab[parNames == "Phi",], digits=3)
      parTab[parNames == "Range",] = round(parTab[parNames == "Range",], digits=0)
      parTab[parNames == "Urban",] = round(parTab[parNames == "Urban",], digits=3)
      
      # round everything else to approximately three significant figures, paying special note 
      # to the non-SD quantities, which should be rounded out to the same decimal
      otherPar = (parNames != "Intercept") & (parNames != "Phi") & (parNames != "Range") & (parNames != "Urban")
      parTab[otherPar,2] = signif(parTab[otherPar,2], digits=3)
      parTab[otherPar,-2] = t(apply(parTab[otherPar,-2], 1, roundToFirstSigFigs, digits=3))
      
      # add the grouping variables and parameter names as new columns
      parTab = cbind(" "=modelTypes, " "=modelVariations, " "=parNames, parTab)
      rownames(parTab) = NULL
      
      # group the rows by model type and model variation
      row_group_label_fonts <-list(list(bold = T, italic = T), list(bold = F, italic = F))
      print(kable(parTab, "latex", booktabs = T, escape=FALSE, format.args=list(drop0trailing=TRUE, scientific=FALSE), 
                  longtable=TRUE, caption = "Longtable") %>%
              collapse_rows(1:2, row_group_label_position ='stack', latex_hline ='custom', custom_latex_hline = 1:2, 
                            row_group_label_fonts = row_group_label_fonts) %>%
              kable_styling(latex_options =c("repeat_header")))
    }
  }
  
  runId = paste0("Beta-1.75margVar", round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 
                 "HHoldVar0urbanOverSamplefrac0", strictPriorText, testText, bigText, sampling, 
                 "models", do.call("paste0", as.list(modelsI)), "nsim", nsim, "MaxDataSetI", maxDataSets)
  if(saveResults) {
    # first collect all the results. Save everything except for the postprocessing arguments: 
    # produceFigures, digits
    objectNames = ls()
    objectNames = objectNames[-match(c("produceFigures", "xtable.args", "tableFormat", "colScale", 
                                       "colUnits", "colDigits"), objectNames)]
    save(list=objectNames, file=paste0("scores", runId, ".RData"))
  }
  
  if(produceFigures) {
    # compare all five
    pdf(paste0("figures/biasbyregion", runId, "_nSamples", nSamples, ".pdf"), width=20, height=12)
    par(mar=c(10,4,2,1), las=2, cex.lab=2, cex.axis=1.4, cex.main=2)
    boxplot(bias~region, data=scoresDirect, at=seq(-1, 275, by=6), 
            col="yellow", xlim=c(0,279), names=FALSE, xaxt="n")
    boxplot(bias~region, data=scoresNaive, at=seq(0, 276, by=6), col="orange", xlim=c(0,230), add=TRUE)
    boxplot(bias~region, data=scoresMercer, at=seq(1, 277, by=6), col="green", xlim=c(0,230), add=TRUE, xaxt="n")
    boxplot(bias~region, data=scoresBYM, at=seq(2, 278, by=6), col="lightblue", xlim=c(0,230), add=TRUE, xaxt="n")
    boxplot(bias~region, data=scoresSPDE, at=seq(3, 279, by=6), col="purple", xlim=c(0,230), add=TRUE, xaxt="n")
    # axis(2, at=seq(0.5, 276.5, by=6), labels=scoresDirect$region[1:47])
    # axis(1, at=seq(0.5, 276.5, by=6), labels=scoresDirect$region[1:47])
    legend("top", c("Direct estimates", "Naive", "Mercer", "BYM", "SPDE"),
           fill = c("yellow", "orange", "green", "lightblue", "purple"), ncol=4, cex=2)
    abline(h=0, lwd=2, col=2)
    dev.off()
    
    pdf(paste0("figures/crpsbyregion", runId, "_nSamples", nSamples, ".pdf"), width=20, height=12)
    par(mar=c(10,4,2,1), las=2, cex.lab=2, cex.axis=1.4, cex.main=2)
    boxplot(crps~region, data=scoresDirect, at=seq(-1, 275, by=6), 
            col="yellow", xlim=c(0,279), names=FALSE, xaxt="n")
    if("Naive" %in% models)
      boxplot(crps~region, data=scoresNaive, at=seq(0, 276, by=6), col="orange", xlim=c(0,230), add=TRUE)
    if("Smoothed Direct" %in% models)
      boxplot(crps~region, data=scoresMercer, at=seq(1, 277, by=6), col="green", xlim=c(0,230), add=TRUE, xaxt="n")
    if("BYM" %in% models)
      boxplot(crps~region, data=scoresBYM, at=seq(2, 278, by=6), col="lightblue", xlim=c(0,230), add=TRUE, xaxt="n")
    if("SPDE" %in% models)
      boxplot(crps~region, data=scoresSPDE, at=seq(3, 279, by=6), col="purple", xlim=c(0,230), add=TRUE, xaxt="n")
    # axis(2, at=seq(0.5, 276.5, by=6), labels=scoresDirect$region[1:47])
    legend("top", c("Direct estimates", "Naive", "Mercer", "BYM", "SPDE"),
           fill = c("yellow", "orange", "green", "lightblue", "purple"), ncol=4, cex=2)
    dev.off()
  }
  
  list(tab=tab, parTab=parTab, unroundedTab=unroundedTab)
}

compareMixtureModeling = function(sigma2=.1^2, n=900, seed=1, nSamples=100, work=FALSE) {
  nu = 1
  thetas = c(0.08, 0.8) / 2.3
  nTest = n * 0.1
  rho = 1
  
  # set random seeds for each simulation, and get the first simulation
  set.seed(seed)
  allSeeds = sample(1:1000000, nSamples, replace = FALSE)
  set.seed(allSeeds[1])
  
  # load and the generator dataset and grid used for predictions
  mixtureCorFun = function(x) {0.5 * stationary.cov(x, theta=thetas[1], Covariance="Matern", smoothness=nu) + 
      0.5 * stationary.cov(x, theta=thetas[2], Covariance="Matern", smoothness=nu)}
  simulationData = getSimulationDataSetsGivenCovarianceTest(mixtureCorFun, nTotal=n, nTest=nTest, marginalVar=rho, errorVar=sigma2, 
                                                            nDataSets=2, plotNameRoot=paste0("(0.5*Matern(", thetas[1], ") + 0.5*Matern(", thetas[2], "))"), fileNameRoot="mix", 
                                                            saveDataSetPlot=FALSE, doPredGrid=TRUE)
  ysTest = c(simulationData$zTest[,1], simulationData$zTestRural[,1], simulationData$zTestUrban[,1], simulationData$zGrid[,1])
  
  # Plot example simulation
  pdf(paste0("Figures/finalMixture/exampleSimulation", nSamples, ".pdf"), width=5, height=5)
  quilt.plot(cbind(simulationData$xGrid, simulationData$yGrid), simulationData$zGrid[,1], 
             nx=70, ny=70, main="")
  points(simulationData$xTrain[,1], simulationData$yTrain[,1], cex=.2, pch=19)
  abline(h=c(-1 / 3, 1 / 3), lty=2)
  abline(v=c(-1 / 3, 1 / 3), lty=2)
  dev.off()
  
  # Plot true covariance
  pdf(paste0("Figures/finalMixture/trueCorrelation", nSamples, ".pdf"), width=5, height=5)
  ds = seq(0, 1, l=200)
  spatialCorFun = function(x) {0.5 * stationary.cov(x, theta=thetas[1], Covariance="Matern", smoothness=nu, distMat=x) + 
      0.5 * stationary.cov(x, theta=thetas[2], Covariance="Matern", smoothness=nu, distMat=x)}
  plot(ds, spatialCorFun(ds) / (1+0.1^2), type="l", main="", 
       xlab="Distance", ylab="Correlation", col="blue", ylim=c(0, 1))
  points(0, 1, pch=19, col="blue", cex=.2)
  dev.off()
  
  plotNamePrefix = "Figures/finalMixture/finalMixture"
  plotNameSuffix = paste0("_nugV", round(sigma2, 2), "_n", n, "_nSamples", nSamples, ".pdf")
  
  # get the file names of the results from the simulations
  if(!work) {
    SPDEname = paste0("savedOutput/simulations/mixtureSPDEAll_nsim", nSamples, "_n", n, "_nu", nu, "_nugV", round(sigma2, 2), "_Kenya", FALSE, 
                      "_noInt", TRUE, "_urbOversamp", round(0, 4), ".RData")
    LKname = paste0("savedOutput/simulations/mixtureLKAll_nsim", nSamples, ".RData")
    LKINLA3name = paste0("savedOutput/simulations/mixtureLKINLAAll_nsim", nSamples, "_L", 3, "_NC14", "_sepRange", FALSE, "_n", n, "_nu", nu, "_nugV", 
                         round(sigma2, 2), "_Kenya", FALSE, "_noInt", TRUE, "_urbOversamp", round(0, 4), ".RData")
    LKINLA2name = paste0("savedOutput/simulations/mixtureLKINLAAll_nsim", nSamples, "_L", 2, "_NC14_126", "_sepRange", TRUE, "_n", n, "_nu", nu, "_nugV", 
                         round(sigma2, 2), "_Kenya", FALSE, "_noInt", TRUE, "_urbOversamp", round(0, 4), ".RData")
  } else {
    SPDEname = paste0("/work/johnpai/simulations/mixtureSPDEAll_nsim", nSamples, "_n", n, "_nu", nu, "_nugV", round(sigma2, 2), "_Kenya", FALSE, 
                      "_noInt", TRUE, "_urbOversamp", round(0, 4), ".RData")
    LKname = paste0("/work/johnpai/simulations/mixtureLKAll_nsim", nSamples, ".RData")
    LKINLA3name = paste0("/work/johnpai/simulations/mixtureLKINLAAll_nsim", nSamples, "_L", 3, "_NC14", "_sepRange", FALSE, "_n", n, "_nu", nu, "_nugV", 
                         round(sigma2, 2), "_Kenya", FALSE, "_noInt", TRUE, "_urbOversamp", round(0, 4), ".RData")
    LKINLA2name = paste0("/work/johnpai/simulations/mixtureLKINLAAll_nsim", nSamples, "_L", 2, "_NC14_126", "_sepRange", TRUE, "_n", n, "_nu", nu, "_nugV", 
                         round(sigma2, 2), "_Kenya", FALSE, "_noInt", TRUE, "_urbOversamp", round(0, 4), ".RData")
  }
  
  # load the simulation results
  # allScoringRules
  # allFits
  ## allCovInfo
  # allPredictionMatrices
  # allAggregatedScoringRules
  # binnedScoringRulesGrid
  # pooledScoringRulesGrid
  # fullPooledScoringRulesLeftOut
  # pooledScoringRulesLeftOut
  ## covMean
  ## upperCov
  ## lowerCov
  ## corMean
  ## upperCor
  ## lowerCor
  ## fullPredictionMatrix
  # leftOutPredictionMatrix
  # leftInPredictionMatrix
  ## leftOutScores
  ## leftInScores
  ## aggregatedScores
  out = load(paste0(SPDEname))
  allPredsSPDE = do.call("c", lapply(allFits, function(x) {x$preds}))
  allSigmasSPDE = do.call("c", lapply(allFits, function(x) {x$sigmas}))
  timingsSPDE = sapply(allFits, function(x) {x$timings$totalTime})
  binnedScoringRulesGridSPDE = binnedScoringRulesGrid
  pooledScoringRulesGridSPDE = pooledScoringRulesGrid
  pooledScoringRulesLeftOutSPDE = pooledScoringRulesLeftOut
  covInfoSPDE = list(d = allCovInfo[[1]]$d, covMean=covMean, 
                     upperCov=upperCov, 
                     lowerCov=lowerCov, 
                     corMean=corMean, 
                     upperCor=upperCor, 
                     lowerCor=lowerCor)
  fullPredictionMatrixSPDE = fullPredictionMatrix
  leftOutScoresSPDE = leftOutScores
  leftInScoresSPDE = leftInScores
  aggregatedScoresSPDE = aggregatedScores
  
  out = load(paste0(LKname))
  allPredsLK = do.call("c", lapply(allFits, function(x) {x$preds}))
  allSigmasLK = do.call("c", lapply(allFits, function(x) {x$sigmas}))
  timingsLK = sapply(allFits, function(x) {x$timings$totalTime})
  binnedScoringRulesGridLK = binnedScoringRulesGrid
  pooledScoringRulesGridLK = pooledScoringRulesGrid
  pooledScoringRulesLeftOutLK = pooledScoringRulesLeftOut
  covInfoLK = list(d = allCovInfo[[1]]$d, covMean=covMean, 
                     upperCov=upperCov, 
                     lowerCov=lowerCov, 
                     corMean=corMean, 
                     upperCor=upperCor, 
                     lowerCor=lowerCor)
  fullPredictionMatrixLK = fullPredictionMatrix
  leftOutScoresLK = leftOutScores
  leftInScoresLK = leftInScores
  aggregatedScoresLK = aggregatedScores
  
  out = load(paste0(LKINLA3name))
  allPredsLKINLA3 = do.call("c", lapply(allFits, function(x) {x$preds}))
  allSigmasLKINLA3 = do.call("c", lapply(allFits, function(x) {x$sigmas}))
  timingsLKINLA3 = sapply(allFits, function(x) {x$timings$totalTime})
  binnedScoringRulesGridLKINLA3 = binnedScoringRulesGrid
  pooledScoringRulesGridLKINLA3 = pooledScoringRulesGrid
  pooledScoringRulesLeftOutLKINLA3 = pooledScoringRulesLeftOut
  covInfoLKINLA3 = list(d = allCovInfo[[1]]$d, covMean=covMean, 
                     upperCov=upperCov, 
                     lowerCov=lowerCov, 
                     corMean=corMean, 
                     upperCor=upperCor, 
                     lowerCor=lowerCor)
  fullPredictionMatrixLKINLA3 = fullPredictionMatrix
  leftOutScoresLKINLA3 = leftOutScores
  leftInScoresLKINLA3 = leftInScores
  aggregatedScoresLKINLA3 = aggregatedScores
  
  out = load(paste0(LKINLA2name))
  allPredsLKINLA2 = do.call("c", lapply(allFits, function(x) {x$preds}))
  allSigmasLKINLA2 = do.call("c", lapply(allFits, function(x) {x$sigmas}))
  timingsLKINLA2 = sapply(allFits, function(x) {x$timings$totalTime})
  binnedScoringRulesGridLKINLA2 = binnedScoringRulesGrid
  pooledScoringRulesGridLKINLA2 = pooledScoringRulesGrid
  pooledScoringRulesLeftOutLKINLA2 = pooledScoringRulesLeftOut
  covInfoLKINLA2 = list(d = allCovInfo[[1]]$d, covMean=covMean, 
                     upperCov=upperCov, 
                     lowerCov=lowerCov, 
                     corMean=corMean, 
                     upperCor=upperCor, 
                     lowerCor=lowerCor)
  fullPredictionMatrixLKINLA2 = fullPredictionMatrix
  leftOutScoresLKINLA2 = leftOutScores
  leftInScoresLKINLA2 = leftInScores
  aggregatedScoresLKINLA2 = aggregatedScores
  
  browser()
  # modelNames = c("SPDE", "LK", "LK-INLA (3 Layer)", "LK-INLA (2 Layer)")
  # modelNames = c("SPDE", "LK-INLA (3 Layer)", "LK-INLA (2 Layer)")
  modelNames = c("SPDE", "LK", "ELK (3 Layer)", "ELK (2 Layer)")
  
  totalLength = length(allFits[[1]]$preds)
  testIndices = (length(allFits[[1]]$preds) - length(ysTest) + 1):length(allFits[[1]]$preds)
  testIndicesAll = c(outer(testIndices, totalLength * (0:9), "+"))
  leftOutIndices = (length(allFits[[1]]$preds) - length(ysTest) + 1):(length(allFits[[1]]$preds) - length(ysTest) + length(simulationData$zTest[,1]))
  leftOutIndicesAll = c(outer(leftOutIndices, totalLength * (0:9), "+"))
  gridIndices = (length(allFits[[1]]$preds) - length(ysTest) + length(simulationData$zTest[,1]) + 1):length(allFits[[1]]$preds)
  gridIndicesAll = c(outer(gridIndices, totalLength * (0:9), "+"))
  leftOutIndicesTest = match(leftOutIndices, testIndices)
  gridIndicesTest = match(gridIndices, testIndices)
  
  # Plot predictions together
  gridCoords = cbind(simulationData$xGrid, simulationData$yGrid)
  ysGrid = ysTest[gridIndicesTest]
  predsGridSPDE = allPredsSPDE[gridIndices]
  predsGridLK = allPredsLK[gridIndices]
  predsGridLKINLA3 = allPredsLKINLA3[gridIndices]
  predsGridLKINLA2 = allPredsLKINLA2[gridIndices]
  zlim = range(c(predsGridSPDE, predsGridLK, predsGridLKINLA3, predsGridLKINLA2))
  
  png(paste0("Figures/finalMixture/exampleMixturePredictions", nSamples, ".png"), width=1000, height=1000)
  par(mfrow=c(2,2))
  
  quilt.plot(gridCoords, predsGridSPDE, main="SPDE Predictions", nx=70, ny=70, zlim=zlim, cex.main=2)
  abline(h=c(-1 / 3, 1 / 3), lty=2)
  abline(v=c(-1 / 3, 1 / 3), lty=2)
  points(simulationData$xTrain[,1], simulationData$yTrain[,1], cex=.2, pch=19)
  
  quilt.plot(gridCoords, predsGridLK, main="LatticeKrig Predictions", nx=70, ny=70, zlim=zlim, cex.main=2)
  abline(h=c(-1 / 3, 1 / 3), lty=2)
  abline(v=c(-1 / 3, 1 / 3), lty=2)
  points(simulationData$xTrain[,1], simulationData$yTrain[,1], cex=.2, pch=19)
  
  quilt.plot(gridCoords, predsGridLKINLA3, main="LK-INLA (L=3) Predictions", nx=70, ny=70, zlim=zlim, cex.main=2)
  abline(h=c(-1 / 3, 1 / 3), lty=2)
  abline(v=c(-1 / 3, 1 / 3), lty=2)
  points(simulationData$xTrain[,1], simulationData$yTrain[,1], cex=.2, pch=19)
  
  quilt.plot(gridCoords, predsGridLKINLA2, main="LK-INLA (L=2) Predictions", nx=70, ny=70, zlim=zlim, cex.main=2)
  abline(h=c(-1 / 3, 1 / 3), lty=2)
  abline(v=c(-1 / 3, 1 / 3), lty=2)
  points(simulationData$xTrain[,1], simulationData$yTrain[,1], cex=.2, pch=19)
  dev.off()
  
  # Plot predictive standard deviations together
  sigmasGridSPDE = allSigmasSPDE[gridIndices]
  sigmasGridLK = allSigmasLK[gridIndices]
  sigmasGridLKINLA3 = allSigmasLKINLA3[gridIndices]
  sigmasGridLKINLA2 = allSigmasLKINLA2[gridIndices]
  zlim = range(c(sigmasGridSPDE, sigmasGridLK, sigmasGridLKINLA3, sigmasGridLKINLA2))
  
  png(paste0("Figures/finalMixture/exampleMixturePredictiveSDs", nSamples, ".png"), width=1000, height=1000)
  par(mfrow=c(2,2))
  
  quilt.plot(gridCoords, sigmasGridSPDE, main="SPDE Predictive SDs", nx=70, ny=70, zlim=zlim, cex.main=2)
  abline(h=c(-1 / 3, 1 / 3), lty=2)
  abline(v=c(-1 / 3, 1 / 3), lty=2)
  points(simulationData$xTrain[,1], simulationData$yTrain[,1], cex=.2, pch=19)
  
  quilt.plot(gridCoords, sigmasGridLK, main="LatticeKrig Predictive SDs", nx=70, ny=70, zlim=zlim, cex.main=2)
  abline(h=c(-1 / 3, 1 / 3), lty=2)
  abline(v=c(-1 / 3, 1 / 3), lty=2)
  points(simulationData$xTrain[,1], simulationData$yTrain[,1], cex=.2, pch=19)
  
  quilt.plot(gridCoords, sigmasGridLKINLA3, main="LK-INLA (L=3) Predictive SDs", nx=70, ny=70, zlim=zlim, cex.main=2)
  abline(h=c(-1 / 3, 1 / 3), lty=2)
  abline(v=c(-1 / 3, 1 / 3), lty=2)
  points(simulationData$xTrain[,1], simulationData$yTrain[,1], cex=.2, pch=19)
  
  quilt.plot(gridCoords, sigmasGridLKINLA2, main="LK-INLA (L=2) Predictive SDs", nx=70, ny=70, zlim=zlim, cex.main=2)
  abline(h=c(-1 / 3, 1 / 3), lty=2)
  abline(v=c(-1 / 3, 1 / 3), lty=2)
  points(simulationData$xTrain[,1], simulationData$yTrain[,1], cex=.2, pch=19)
  dev.off()
  
  # plot correlation functions together
  spatialCovFun = function(x) {0.5 * stationary.cov(x, theta=thetas[1], Covariance="Matern", distMat=x, smoothness=nu) + 
      0.5 * stationary.cov(x, theta=thetas[2], Covariance="Matern", smoothness=nu, distMat=x)}
  mixtureCovFun = function(x) {
    out = spatialCovFun(x)
    out[x == 0] = 1 + sigma2
    out
  }
  mixtureCorFun = function(x) { mixtureCovFun(x) * (1 / (1 + sigma2)) }
  
  pdf(paste0(plotNamePrefix, "Correlation", plotNameSuffix), width=5, height=5)
  plot(covInfoSPDE$d, covInfoSPDE$corMean, col="blue", 
       main="Correlation functions (and 80% CIs)", ylab="Correlation", xlab="Distance", type="l")
  lines(covInfoSPDE$d, covInfoSPDE$upperCor, col="blue", lty=2)
  lines(covInfoSPDE$d, covInfoSPDE$lowerCor, col="blue", lty=2)
  
  lines(covInfoLK$d, covInfoLK$corMean, col="black")
  # lines(covInfoLK$d, covInfoLK$upperCor, col="black", lty=2)
  # lines(covInfoLK$d, covInfoLK$lowerCor, col="black", lty=2)
  
  lines(covInfoLKINLA3$d, covInfoLKINLA3$corMean, col="purple")
  lines(covInfoLKINLA3$d, covInfoLKINLA3$upperCor, col="purple", lty=2)
  lines(covInfoLKINLA3$d, covInfoLKINLA3$lowerCor, col="purple", lty=2)
  
  lines(covInfoLKINLA2$d, covInfoLKINLA2$corMean, col="red")
  lines(covInfoLKINLA2$d, covInfoLKINLA2$upperCor, col="red", lty=2)
  lines(covInfoLKINLA2$d, covInfoLKINLA2$lowerCor, col="red", lty=2)
  
  lines(covInfoLKINLA2$d, mixtureCorFun(covInfoLKINLA2$d), col="green")
  
  legend("topright", c("SPDE", "LK", "LK-INLA (L=3)", "LK-INLA (L=2)", "Truth"), lty=1, col=c("blue", "black", "purple", "red", "green"))
  # legend("topright", c("SPDE", "LK-INLA (3)", "LK-INLA (2)", "Truth"), lty=1, col=c("blue", "purple", "red", "green"))
  dev.off()
  
  pdf(paste0(plotNamePrefix, "CorrelationNoCIs", plotNameSuffix), width=5, height=5)
  plot(covInfoSPDE$d, covInfoSPDE$corMean, col="blue", 
       main="", ylab="Correlation", xlab="Distance", type="l")
  
  lines(covInfoLK$d, covInfoLK$corMean, col="black")
  
  lines(covInfoLKINLA3$d, covInfoLKINLA3$corMean, col="purple")
  
  lines(covInfoLKINLA2$d, covInfoLKINLA2$corMean, col="red")
  
  lines(covInfoLKINLA2$d, mixtureCorFun(covInfoLKINLA2$d), col="green")
  
  legend("topright", c("SPDE", "LK", "LK-INLA (L=3)", "LK-INLA (L=2)", "Truth"), lty=1, col=c("blue", "black", "purple", "red", "green"))
  # legend("topright", c("SPDE", "LK-INLA (3)", "LK-INLA (2)", "Truth"), lty=1, col=c("blue", "purple", "red", "green"))
  dev.off()
  
  # construct information for separating out prediction types
  totalLength = length(allFits[[1]]$preds)
  testIndices = (length(allFits[[1]]$preds) - length(ysTest) + 1):length(allFits[[1]]$preds)
  testIndicesAll = c(outer(testIndices, totalLength * (0:9), "+"))
  leftOutIndices = (length(allFits[[1]]$preds) - length(ysTest) + 1):(length(allFits[[1]]$preds) - length(ysTest) + length(simulationData$zTest[,1]))
  leftOutIndicesAll = c(outer(leftOutIndices, totalLength * (0:9), "+"))
  gridIndices = (length(allFits[[1]]$preds) - length(ysTest) + length(simulationData$zTest[,1]) + 1):length(allFits[[1]]$preds)
  gridIndicesAll = c(outer(gridIndices, totalLength * (0:9), "+"))
  leftOutIndicesTest = match(leftOutIndices, testIndices)
  gridIndicesTest = match(gridIndices, testIndices)
  
  set.seed(1)
  maxPoints = 7500
  gridIndicesSample = sample(gridIndicesAll, maxPoints, replace=FALSE)
  
  # plot prediction standard deviations
  pdf(paste0("Figures/finalMixture/pairSigmas", nSamples, ".pdf"), width=8, height=8)
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
    points(x[gridIndicesSample],y[gridIndicesSample],..., col="black", cex=.1, pch=".")
    # points(x[leftOutIndicesAll],y[leftOutIndicesAll],..., col="blue", cex=.5, pch=19)
  }
  
  values = data.frame(allSigmasSPDE, allSigmasLK, allSigmasLKINLA3, allSigmasLKINLA2)
  zlim = range(c(as.matrix(values)))
  lims = rep(list(zlim), length(values))
  myPairs(values, 
          labels=modelNames, 
          lower.panel=my_line, upper.panel = my_line, 
          main=paste0("Predictive SDs"), 
          lims=lims, oma=c(3,3,6,7))
  dev.off()
  
  # aggregated scoring rules
  fullAggregatedScoresSPDE = getScores(fullPredictionMatrixSPDE$Truth, fullPredictionMatrixSPDE$Est, fullPredictionMatrixSPDE$SDs^2, fullPredictionMatrixSPDE$Lower, fullPredictionMatrixSPDE$Upper, getAverage=FALSE)
  fullAggregatedScoresLK = getScores(fullPredictionMatrixLK$Truth, fullPredictionMatrixLK$Est, fullPredictionMatrixLK$SDs^2, fullPredictionMatrixLK$Lower, fullPredictionMatrixLK$Upper, getAverage=FALSE)
  fullAggregatedScoresLKINLA3 = getScores(fullPredictionMatrixLKINLA3$Truth, fullPredictionMatrixLKINLA3$Est, fullPredictionMatrixLKINLA3$SDs^2, fullPredictionMatrixLKINLA3$Lower, fullPredictionMatrixLKINLA3$Upper, getAverage=FALSE)
  fullAggregatedScoresLKINLA2 = getScores(fullPredictionMatrixLKINLA2$Truth, fullPredictionMatrixLKINLA2$Est, fullPredictionMatrixLKINLA2$SDs^2, fullPredictionMatrixLKINLA2$Lower, fullPredictionMatrixLKINLA2$Upper, getAverage=FALSE)
  
  leftOutIndices = seq(5, nrow(fullPredictionMatrixSPDE), by=9)
  leftInIndices = (1:nrow(fullPredictionMatrixSPDE))[-leftOutIndices]
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
    points(x[leftInIndices],y[leftInIndices],..., col="black", cex=.5, pch=19)
    points(x[leftOutIndices],y[leftOutIndices],..., col="blue", cex=.5, pch=19)
  }
  
  values = data.frame(fullPredictionMatrixSPDE$SDs, fullPredictionMatrixLK$SDs, fullPredictionMatrixLKINLA2$SDs, fullPredictionMatrixLKINLA3$SDs)
  zlim = range(c(as.matrix(values)))
  lims = rep(list(zlim), length(values))
  myPairs(values, 
          labels=modelNames, 
          lower.panel=my_line, upper.panel = my_line, 
          main=paste0("Aggregated Predictive SDs"), 
          lims=lims, oma=c(3,3,6,7))
  
  
  values = data.frame(SPDE=fullAggregatedScoresSPDE$RMSE, LK=fullAggregatedScoresLK$RMSE, LKINLA3=fullAggregatedScoresLKINLA3$RMSE, LKINLA2=fullAggregatedScoresLKINLA2$RMSE)[leftOutIndices,]
  boxplot(values, names=modelNames, col=rainbow(4), main="RMSE (Central Block)")
  
  values = data.frame(SPDE=fullAggregatedScoresSPDE$RMSE, LK=fullAggregatedScoresLK$RMSE, LKINLA3=fullAggregatedScoresLKINLA3$RMSE, LKINLA2=fullAggregatedScoresLKINLA2$RMSE)[leftInIndices,]
  boxplot(values, names=modelNames, col=rainbow(4), main="RMSE (Outer Blocks)")
  
  values = data.frame(SPDE=fullAggregatedScoresSPDE$CRPS, LK=fullAggregatedScoresLK$CRPS, LKINLA3=fullAggregatedScoresLKINLA3$CRPS, LKINLA2=fullAggregatedScoresLKINLA2$CRPS)[leftOutIndices,]
  boxplot(values, names=modelNames, col=rainbow(4), main="CRPS (Central Block)")
  
  values = data.frame(SPDE=fullAggregatedScoresSPDE$CRPS, LK=fullAggregatedScoresLK$CRPS, LKINLA3=fullAggregatedScoresLKINLA3$CRPS, LKINLA2=fullAggregatedScoresLKINLA2$CRPS)[leftInIndices,]
  boxplot(values, names=modelNames, col=rainbow(4), main="CRPS (Outer Blocks)")
  
  
  allLeftOutAggregatedScores = rbind(leftOutScoresSPDE, leftOutScoresLK, leftOutScoresLKINLA3, leftOutScoresLKINLA2)
  allLeftInAggregatedScores = rbind(leftInScoresSPDE, leftInScoresLK, leftInScoresLKINLA3, leftInScoresLKINLA2)
  allAggregatedScores = rbind(aggregatedScoresSPDE, aggregatedScoresLK, aggregatedScoresLKINLA3, aggregatedScoresLKINLA2)
  rownames(allLeftOutAggregatedScores) = modelNames
  rownames(allLeftInAggregatedScores) = modelNames
  rownames(allAggregatedScores) = modelNames
  
  allLeftOutAggregatedScores$Coverage = allLeftOutAggregatedScores$Coverage * 100
  allLeftInAggregatedScores$Coverage = allLeftInAggregatedScores$Coverage * 100
  allAggregatedScores$Coverage = allAggregatedScores$Coverage * 100
  
  print("Left out aggregated scores:")
  print(xtable(format(allLeftOutAggregatedScores[,-3], digits=3)))
  # temp = allLeftOutAggregatedScores[,-3]
  # theseNames = dimnames(temp)
  # temp = matrix(unlist(temp), ncol=ncol(temp))
  # dimnames(temp) = theseNames
  # temp = cbind(round(temp[,1:2], digits=4), round(temp[,3:4], digits=4), Coverage=round(temp[,5], digits=0), Width=round(temp[,6], digits=3))
  # mode(temp) = "character"
  # print(xtable(temp))
  
  print("Left in aggregated scores:")
  print(xtable(format(allLeftInAggregatedScores[,-3], digits=2)))
  
  print("All aggregated scores:")
  print(xtable(format(allAggregatedScores[,-3], digits=2)))
  
  # computation times for the ENTIRE simulation study, not just fitting the models.  This includes model fitting, predictions, uncertainty calculations, aggregation, and covariance calculations:
  # -rw------- 1 johnpai posixgroup 123956328 Jan 23 03:27 mixtureLKINLAsim1_L2_NC14_126_sepRangeTRUE_n900_nu1_nugV0.01_KenyaFALSE_noIntTRUE_urbOversamp0.RData
  # -rw------- 1 johnpai posixgroup 123047142 Jan 23 01:45 mixtureLKINLAsim1_L3_NC14_sepRangeFALSE_n900_nu1_nugV0.01_KenyaFALSE_noIntTRUE_urbOversamp0.RData
  # -rw------- 1 johnpai posixgroup 124691837 Jan 23 00:31 mixtureSPDEsim1_n900_nu1_nugV0.01_KenyaFALSE_noIntTRUE_urbOversamp0.RData
  # -rw------- 1 johnpai posixgroup   116011406 Jan 29 13:33 mixtureLKsim1.RData
  
  # -rw------- 1 johnpai posixgroup   123880261 Jan 28 15:57 mixtureLKINLAsim100_L2_NC14_126_sepRangeTRUE_n900_nu1_nugV0.01_KenyaFALSE_noIntTRUE_urbOversamp0.RData
  # -rw------- 1 johnpai posixgroup   122967029 Jan 27 17:36 mixtureLKINLAsim100_L3_NC14_sepRangeFALSE_n900_nu1_nugV0.01_KenyaFALSE_noIntTRUE_urbOversamp0.RData
  # -rw------- 1 johnpai posixgroup   124621949 Jan 24 08:31 mixtureSPDEsim100_n900_nu1_nugV0.01_KenyaFALSE_noIntTRUE_urbOversamp0.RData
  # -rw------- 1 johnpai posixgroup   115944856 Feb  3 12:11 mixtureLKsim100.RData
  totalTimeSPDE = 60 * 24 * (24 - 23) + 60 * (8 - 0) + 1 * (31 - 31)
  totalTimeLK = 60 * 24 * (5) + 60 * (12 - 13) + 1 * (11 - 33)
  totalTimeLKINLA3 = 60 * 24 * (27 - 23) + 60 * (17 - 1) + 1 * (36 - 45)
  totalTimeLKINLA2 = 60 * 24 * (28 - 23) + 60 * (15 - 3) + 1 * (57 - 27)
  averageTimeSPDE = totalTimeSPDE / 99
  averageTimeLK = totalTimeLK / 99
  averageTimeLKINLA3 = totalTimeLKINLA3 / 99
  averageTimeLKINLA2 = totalTimeLKINLA2 / 99
  pooledScores = rbind(pooledScoringRulesGridSPDE, pooledScoringRulesGridLK, 
                       pooledScoringRulesGridLKINLA3, pooledScoringRulesGridLKINLA2)
  pooledScores = cbind(pooledScores[,-c(1, 2, 5)], "Time elapsed (min.)"=c(averageTimeSPDE, averageTimeLK, averageTimeLKINLA3, averageTimeLKINLA2))
  rownames(pooledScores) = modelNames
  pooledScores$Coverage = pooledScores$Coverage * 100
  
  print("All pointwise scores:")
  print(xtable(format(round(pooledScores, digits=3), scientific=FALSE)))
  print(xtable(format(round(pooledScores[,-c(1, 2)], digits=3), scientific=FALSE)))
  # binned scoring rules
  
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
    points(x[gridIndicesTest],y[gridIndicesTest],..., col="black", cex=.1, pch=".")
    points(x[leftOutIndicesTest],y[leftOutIndicesTest],..., col="blue", cex=.5, pch=19)
  }
  
  # calculate CRPS
  values = data.frame(crps(ysTest, fitSPDE$preds[testIndices], fitSPDE$sigmas[testIndices]^2, getAverage=FALSE), 
                      crps(ysTest, fitLKINLA3$preds[testIndices], fitLKINLA3$sigmas[testIndices]^2, getAverage=FALSE), 
                      crps(ysTest, fitLKINLA2$preds[testIndices], fitLKINLA2$sigmas[testIndices]^2, getAverage=FALSE))
  zlim = range(c(as.matrix(values)))
  lims = rep(list(zlim), length(values))
  myPairs(values, 
          labels=modelNames, 
          lower.panel=my_line, upper.panel = my_line, 
          main=paste0("CRPS"), 
          lims=lims, oma=c(3,3,6,7), log="xy")
  
  ##### Aggregated results
  gridCoords = cbind(simulationData$xGrid, simulationData$yGrid)
  inSampleIndices = (1:9)[-5]
  leftOutIndices = 5
  
  A = makeNumericalIntegralMat(gridCoords, mx=3, my=3)
  aggregatedTruth = A %*% ysTest[gridIndicesTest]
  
  sdsSPDE = apply(A %*% fitSPDE$predMat[gridIndices,], 1, sd)
  # sdsLK = apply(A %*% fitLK$predMat[gridIndices,], 1, sd)
  sdsLKINLA3 = apply(A %*% fitLKINLA3$predMat[gridIndices,], 1, sd)
  sdsLKINLA2 = apply(A %*% fitLKINLA2$predMat[gridIndices,], 1, sd)
  
  # est = A %*% preds[gridIndices]
  # aggregatedPredMatSPDlat = A %*% predictionMatrixSPDE[gridIndices,]
  # vars = apply(aggregatedPredMat, 1, var)
  # sds = apply(aggregatedPredMat, 1, sd)
  
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
    points(x[inSampleIndices],y[inSampleIndices],..., col="black", cex=1, pch=19)
    points(x[leftOutIndices],y[leftOutIndices],..., col="blue", cex=1, pch=19)
  }
  
  values = data.frame((A %*% fitSPDE$preds[gridIndices]-aggregatedTruth)^2, (A %*% fitLKINLA3$preds[gridIndices]-aggregatedTruth)^2, (A %*% fitLKINLA2$preds[gridIndices]-aggregatedTruth)^2)
  zlim = range(c(as.matrix(values)))
  lims = rep(list(zlim), length(values))
  myPairs(values, 
          labels=modelNames, 
          lower.panel=my_line, upper.panel = my_line, 
          main=paste0("Square Error"), 
          lims=lims, oma=c(3,3,6,7), log="xy")
  
  values = data.frame((A %*% fitSPDE$preds[gridIndices]-aggregatedTruth)/sdsSPDE, 
                      (A %*% fitLKINLA3$preds[gridIndices]-aggregatedTruth)/sdsLKINLA3, 
                      (A %*% fitLKINLA2$preds[gridIndices]-aggregatedTruth)/sdsLKINLA2)
  zlim = range(c(as.matrix(values)))
  lims = rep(list(zlim), length(values))
  myPairs(values, 
          labels=modelNames, 
          lower.panel=my_line, upper.panel = my_line, 
          main=paste0("Studentized Error"), 
          lims=lims, oma=c(3,3,6,7))
  
  values = data.frame(sdsSPDE, sdsLKINLA3, sdsLKINLA2)
  zlim = range(c(as.matrix(values)))
  lims = rep(list(zlim), length(values))
  myPairs(values, 
          labels=modelNames, 
          lower.panel=my_line, upper.panel = my_line, 
          main=paste0("Predictive SDs"), 
          lims=lims, oma=c(3,3,6,7), log="xy")
  
  # calculate CRPS
  values = data.frame(crps(aggregatedTruth, A %*% fitSPDE$preds[gridIndices], sdsSPDE^2, getAverage=FALSE), 
                      crps(aggregatedTruth, A %*% fitLKINLA3$preds[gridIndices], sdsLKINLA3^2, getAverage=FALSE), 
                      crps(aggregatedTruth, A %*% fitLKINLA2$preds[gridIndices], sdsLKINLA2^2, getAverage=FALSE))
  zlim = range(c(as.matrix(values)))
  lims = rep(list(zlim), length(values))
  myPairs(values, 
          labels=modelNames, 
          lower.panel=my_line, upper.panel = my_line, 
          main=paste0("CRPS"), 
          lims=lims, oma=c(3,3,6,7), log="xy")
  
  ##### binned results versus distance
  ns = binnedScoringRulesGridSPDE$nPerBin
  includeBins = 1:(length(ns)-2)
  cex = 5/sqrt(ns)
  pdf(paste0("Figures/finalMixture/binnedMSE", nSamples, ".pdf"), width=5, height=5)
  values = data.frame(binnedScoringRulesGridSPDE$MSE, binnedScoringRulesGridLK$MSE, binnedScoringRulesGridLKINLA3$MSE, binnedScoringRulesGridLKINLA2$MSE)[includeBins,]
  zlim = range(c(as.matrix(values)))
  plot(binnedScoringRulesGridSPDE$NNDist[includedBins], values[[1]][includedBins], ylim=zlim, col="blue", 
       main="MSE vs Distance to Observation", xlab="Distance to Observation", ylab="MSE")
  lines(binnedScoringRulesGridSPDE$NNDist[includeBins], values[[1]], col="blue")
  points(binnedScoringRulesGridSPDE$NNDist[includeBins], values[[2]], col="black")
  lines(binnedScoringRulesGridSPDE$NNDist[includeBins], values[[2]], col="black")
  points(binnedScoringRulesGridSPDE$NNDist[includeBins], values[[3]], col="purple")
  lines(binnedScoringRulesGridSPDE$NNDist[includeBins], values[[3]], col="purple")
  points(binnedScoringRulesGridSPDE$NNDist[includeBins], values[[4]], col="red")
  lines(binnedScoringRulesGridSPDE$NNDist[includeBins], values[[4]], col="red")
  legend("bottomright", c("SPDE", "LK", "ELK (3 Layers)", "ELK (2 Layers)"), lty=1, pch=1, col=c("blue", "black", "purple", "red"))
  dev.off()
  
  pdf(paste0("Figures/finalMixture/binnedRMSE", nSamples, ".pdf"), width=5, height=5)
  values = data.frame(binnedScoringRulesGridSPDE$RMSE, binnedScoringRulesGridLK$RMSE, binnedScoringRulesGridLKINLA3$RMSE, binnedScoringRulesGridLKINLA2$RMSE)[includeBins,]
  zlim = range(c(as.matrix(values)))
  plot(binnedScoringRulesGridSPDE$NNDist[includedBins], values[[1]][includedBins], ylim=zlim, col="blue", 
       main="RMSE vs Distance to Observation", xlab="Distance to Observation", ylab="RMSE")
  lines(binnedScoringRulesGridSPDE$NNDist[includeBins], values[[1]], col="blue")
  points(binnedScoringRulesGridSPDE$NNDist[includeBins], values[[2]], col="black")
  lines(binnedScoringRulesGridSPDE$NNDist[includeBins], values[[2]], col="black")
  points(binnedScoringRulesGridSPDE$NNDist[includeBins], values[[3]], col="purple")
  lines(binnedScoringRulesGridSPDE$NNDist[includeBins], values[[3]], col="purple")
  points(binnedScoringRulesGridSPDE$NNDist[includeBins], values[[4]], col="red")
  lines(binnedScoringRulesGridSPDE$NNDist[includeBins], values[[4]], col="red")
  legend("bottomright", c("SPDE", "LK", "ELK (3 Layers)", "ELK (2 Layers)"), lty=1, pch=1, col=c("blue", "black", "purple", "red"))
  dev.off()
  
  ##### TODO: fix plot titles?
  pdf(paste0("Figures/finalMixture/binnedCRPS", nSamples, ".pdf"), width=5, height=5)
  values = data.frame(binnedScoringRulesGridSPDE$CRPS, binnedScoringRulesGridLK$CRPS, binnedScoringRulesGridLKINLA3$CRPS, binnedScoringRulesGridLKINLA2$CRPS)[includeBins,]
  zlim = range(c(as.matrix(values)))
  plot(binnedScoringRulesGridSPDE$NNDist[includeBins], values[[1]], ylim=zlim, col="blue", 
       main="CRPS vs Distance to Observation", xlab="Distance to Observation", ylab="CRPS")
  lines(binnedScoringRulesGridSPDE$NNDist[includeBins], values[[1]], col="blue")
  points(binnedScoringRulesGridSPDE$NNDist[includeBins], values[[2]], col="black")
  lines(binnedScoringRulesGridSPDE$NNDist[includeBins], values[[2]], col="black")
  points(binnedScoringRulesGridSPDE$NNDist[includeBins], values[[3]], col="purple")
  lines(binnedScoringRulesGridSPDE$NNDist[includeBins], values[[3]], col="purple")
  points(binnedScoringRulesGridSPDE$NNDist[includeBins], values[[4]], col="red")
  lines(binnedScoringRulesGridSPDE$NNDist[includeBins], values[[4]], col="red")
  legend("bottomright", c("SPDE", "LK", "ELK (3 Layers)", "ELK (2 Layers)"), pch=1, lty=1, col=c("blue", "black", "purple", "red"))
  dev.off()
  
  pdf(paste0("Figures/finalMixture/binnedCoverage", nSamples, ".pdf"), width=5, height=5)
  values = data.frame(binnedScoringRulesGridSPDE$Coverage, binnedScoringRulesGridLK$Coverage, binnedScoringRulesGridLKINLA3$Coverage, binnedScoringRulesGridLKINLA2$Coverage)[includeBins,]
  zlim = range(c(as.matrix(values)))
  plot(binnedScoringRulesGridSPDE$NNDist[includeBins], values[[1]], ylim=zlim, col="blue", 
       main="80% Coverage vs Distance to Observation", xlab="Distance to Observation", ylab="Coverage")
  abline(h=.8, lty=2)
  lines(binnedScoringRulesGridSPDE$NNDist[includeBins], values[[1]], col="blue")
  points(binnedScoringRulesGridSPDE$NNDist[includeBins], values[[2]], col="black")
  lines(binnedScoringRulesGridSPDE$NNDist[includeBins], values[[2]], col="black")
  points(binnedScoringRulesGridSPDE$NNDist[includeBins], values[[3]], col="purple")
  lines(binnedScoringRulesGridSPDE$NNDist[includeBins], values[[3]], col="purple")
  points(binnedScoringRulesGridSPDE$NNDist[includeBins], values[[4]], col="red")
  lines(binnedScoringRulesGridSPDE$NNDist[includeBins], values[[4]], col="red")
  legend("bottomright", c("SPDE", "LK", "ELK (3 Layers)", "ELK (2 Layers)"), pch=1, lty=1, col=c("blue", "black", "purple", "red"))
  dev.off()
  
  pdf(paste0("Figures/finalMixture/binnedWidth", nSamples, ".pdf"), width=5, height=5)
  values = data.frame(binnedScoringRulesGridSPDE$Width, binnedScoringRulesGridLK$Width, binnedScoringRulesGridLKINLA3$Width, binnedScoringRulesGridLKINLA2$Width)[includeBins,]
  zlim = range(c(as.matrix(values)))
  plot(binnedScoringRulesGridSPDE$NNDist[includeBins], values[[1]], ylim=zlim, col="blue", 
       main="80% CI Width vs Distance to Observation", xlab="Distance to Observation", ylab="80% CI Width")
  lines(binnedScoringRulesGridSPDE$NNDist[includeBins], values[[1]], col="blue")
  points(binnedScoringRulesGridSPDE$NNDist[includeBins], values[[2]], col="black")
  lines(binnedScoringRulesGridSPDE$NNDist[includeBins], values[[2]], col="black")
  points(binnedScoringRulesGridSPDE$NNDist[includeBins], values[[3]], col="purple")
  lines(binnedScoringRulesGridSPDE$NNDist[includeBins], values[[3]], col="purple")
  points(binnedScoringRulesGridSPDE$NNDist[includeBins], values[[4]], col="red")
  lines(binnedScoringRulesGridSPDE$NNDist[includeBins], values[[4]], col="red")
  legend("bottomright", c("SPDE", "LK", "ELK (3 Layers)", "ELK (2 Layers)"), pch=1, lty=1, col=c("blue", "black", "purple", "red"))
  dev.off()
}


