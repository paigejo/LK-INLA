# script for getting predictions from the classical Gaussian process to data

# Function for computing conditional mean and variance for normal distribution given data.
# Xp | Xd has (conditional) mean:
# muc = muP + SigmaPtoD %*% SigmaD^(-1) %*% (Xd - muD)
# and (conditional) variance:
# Sigmac = SigmaP - SigmaPtoD %*% SigmaD^(-1) %*% SigmaDtoP
conditionalNormal = function(Xd, muP, muD, SigmaP, SigmaD, SigmaPtoD) {
  
  # SigmaDInv = solve(SigmaD) # NOTE: luckily we only need to do this once.
  # muc = muP + SigmaPtoD %*% SigmaDTildeInv %*% (Xd - muD)
  # Sigmac = SigmaP - SigmaPtoD %*% SigmaDInv %*% SigmaDtoP
  
  # compute conditional mean and variance of zeta
  muc = muP + SigmaPtoD %*% solve(SigmaD, Xd - muD)
  Sigmac = SigmaP - SigmaPtoD %*% solve(SigmaD, t(SigmaPtoD))
  
  return(list(muc=muc, Sigmac=Sigmac))
}

# NOTE: cov.fun should be the TRUE underlying covariance
# NOTE: betas should be the TRUE covariate coefficients for the underlying model
# A: an observation/aggregation matrix transforming the predictions to another vector of interest
GPpreds = function(obsCoords, obsValues, xObs=matrix(rep(1, length(obsValues)), ncol=1), predCoords, 
                   xPred = matrix(rep(1, nrow(predCoords)), ncol=1), betas=0, 
                   cov.fun = "stationary.cov", significanceCI=.8, A=NULL) {
  
  # First construct relevant covariance and cross covariance matrices
  SigmaObs = do.call(cov.fun, list(obsCoords))
  SigmaPredsObs = do.call(cov.fun, list(predCoords, obsCoords))
  SigmaPred = do.call(cov.fun, list(predCoords))
  
  # Now construct relevant mean vectors
  if(length(xObs) != 0) {
    muObs = xObs * betas
    muPred = xPred * betas
  } else {
    muObs = rep(0, length(obsValues))
    muPred = rep(0, nrow(predCoords))
  }
  
  # calculate the predictive distribution
  predictions = conditionalNormal(obsValues, muPred, muObs, SigmaPred, SigmaObs, SigmaPredsObs)
  muc = predictions$muc
  Sigmac = predictions$Sigmac
  sigmac = sqrt(diag(Sigmac))
  
  # transform the predictive distribution if necessary
  if(!is.null(A)) {
    mucAggregated = A %*% muc[1:ncol(A)]
    SigmacAggregated = A %*% Sigmac[1:ncol(A),1:ncol(A)] %*% t(A)
    sigmacAggregated = sqrt(diag(SigmacAggregated))
    
    lowerAggregated = mucAggregated + qnorm((1 - significanceCI) / 2, sd=sigmacAggregated)
    medsAggregated = mucAggregated + qnorm(.5, sd=sigmacAggregated)
    upperAggregated = mucAggregated + qnorm(1 - (1 - significanceCI) / 2, sd=sigmacAggregated)
  }
  
  # calculate confidence intervals
  lower = muc + qnorm((1 - significanceCI) / 2, sd=sigmac)
  meds = muc + qnorm(.5, sd=sigmac)
  upper = muc + qnorm(1 - (1 - significanceCI) / 2, sd=sigmac)
  
  ## preds
  ## sigmas
  ## lower
  ## medians
  ## upper
  ## interceptSummary
  ## rangeSummary
  ## sdSummary
  ## varSummary
  
  interceptSummary = c(Est=betas[1], SD=0, Qlower=betas[1], Q50=betas[1], Qupper=betas[1])
  rangeSummary = c(Est=NA, SD=NA, Qlower=NA, Q50=NA, Qupper=NA)
  sdSummary = c(Est=1, SD=1, Qlower=1, Q50=1, Qupper=1)
  varSummary = c(Est=1, SD=1, Qlower=1, Q50=1, Qupper=1)
  
  # return results
  list(preds=muc, sigmas=sigmac, lower=lower, medians=meds, upper=upper, 
       predsAggregated=mucAggregated, sigmasAggregated=sigmacAggregated, 
       lowerAggregated=lowerAggregated, mediansAggregated=medsAggregated, upperAggregated=upperAggregated, 
       interceptSummary=interceptSummary, rangeSummary=rangeSummary, 
       sdSummary=sdSummary, varSummary=varSummary)
}

# this function generates results for the simulation study for the GP model
# input arguments:
#   argument specifying the dataset type
resultsGP = function(randomSeeds=NULL, covType=c("exponential", "matern", "mixture"), rangeText=c("01", "05", "1", "mix"), 
                     maxDataSets=NULL) {
  # determine the type of covariance for the data set
  covType = match.arg(covType)
  
  # determine the spatial range for the data set. No range text means it's a mixture
  rangeText = match.arg(rangeText)
  
  # construct the file name for the desired data set and load it
  if(rangeText == "mix")
    dataText = paste0(covType, "DataSet.RData")
  else
    dataText = paste0(covType, rangeText, "DataSet.RData")
  out = load(dataText)
  dataSets = simulationData
  
  # generate results for all data sets and return results
  resultsGP = fitModelToDataSets(GPpreds, dataSets, randomSeeds=randomSeeds, maxDataSets=maxDataSets)
  
  # save results
  fileName = paste0("resultsGP_cov", covType, "Range", rangeText, "maxDataSets", maxDataSets, ".RData")
  save(resultsGP, file=fileName)
  
  resultsGP
}









