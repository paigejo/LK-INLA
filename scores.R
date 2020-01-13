# this script contains all used scoring rules used to evaluate the models for a given data set

# this function computes all scoring rules
# truth: the true values
# est: the estimates, of the same length as truth. By default calculated from estMat
# var: the estimates, of the same length as truth. By default calculated from estMat
# lower: the lower end of the credible interval. By default calculated from estMat
# upper: the upper end of the credible interval. By default calculated from estMat
# estMat: a matrix of joint estimate draws, with number of rows equal to the length of truth, a number of 
#          columns equal to the number of draws. If not included, a gaussian distribution is assumed.
# significance: the significance level of the credible interval. By default 80%
# distances: the distances to the nearest observation if not NULL, scores are broken up 
#            as a function of nearest neighbor distances
# breaks: the number of equal spaced bins to break the scores into as a function of distance
# NOTE: Discrete, count level credible intervals are estimated based on the input estMat along with coverage and CRPS
getScores = function(truth, est=NULL, var=NULL, lower=NULL, upper=NULL, estMat=NULL, significance=.8, 
                     distances=NULL, breaks=30, doRandomReject=FALSE, doFuzzyReject=TRUE) {
  
  # if distances is included, must also break down scoring rules by distance bins
  if(!is.null(distances)) {
    # construct the distance bins with which to group the data and compute scores within
    if(length(breaks) == 1)
      breaks = seq(0, max(distances), l=breaks)
    binsI = cut(distances, breaks, labels=1:(length(breaks)-1), include.lowest=TRUE)
    centers = breaks[1:(length(breaks)-1)] + diff(breaks)/2
    uniqueBinsI = sort(unique(binsI))
    
    # determine the number of observations per bin
    nPerBin = as.numeric(table(binsI))
    
    # helper function to compute the scoring rules for a given bin
    getSubScores = function(uniqueBinI, truth, est, var, lower, upper, estMat, significance) {
      thisDatI = binsI == uniqueBinI
      
      getScores(truth[thisDatI], est[thisDatI], var[thisDatI], lower[thisDatI], upper[thisDatI], 
                estMat[thisDatI,], significance, doRandomReject=doRandomReject, doFuzzyReject=doFuzzyReject)
    }
    
    # calculate scores for each bin individually
    binnedScores = t(sapply(uniqueBinsI, getSubScores, truth=truth, est=est, var=var, lower=lower, upper=upper, 
                          estMat=estMat, significance=significance))
    
    # make sure each variable in binnedScores is a numeric, not a list...
    temp = matrix(unlist(binnedScores), nrow=length(uniqueBinsI))
    theseNames = colnames(binnedScores)
    binnedScores = data.frame(temp)
    names(binnedScores) = theseNames
    binnedScores = as.data.frame(cbind(NNDist=centers[uniqueBinsI], nPerBin=nPerBin[uniqueBinsI], binnedScores))
  }
  
  # compute central estimates if estMat is not null
  if(!is.null(estMat)) {
    if(is.null(est))
      est = rowMeans(estMat)
  }
  
  # first calculate bias, variance, and MSE
  out = mse(truth, est)
  thisMSE = out$MSE
  thisBias = out$bias
  thisVar = out$var
  
  # calculate coverage and credible interval width with and without binomial variation
  coverage = coverage(truth, est, var, lower, upper, estMat=estMat, 
                      significance=significance, returnIntervalWidth=TRUE, 
                      doRandomReject=doRandomReject, doFuzzyReject=doFuzzyReject)
  thisCoverage = coverage[1]
  thisWidth = coverage[2]
  
  # calculate CRPS
  thisCRPS = crps(truth, est, var, estMat=estMat)
  
  # collect the results in a data frame
  results = matrix(c(thisBias, thisVar, thisMSE, sqrt(thisMSE), thisCRPS, thisCoverage, 
                     thisWidth), nrow=1)
  colnames(results) = c("Bias", "Var", "MSE", "RMSE", "CRPS", "Coverage", "Width")
  results = as.data.frame(results)
  
  # include both binned and pooled results in one final table
  if(!is.null(distances)) {
    results = list(pooledResults=results, binnedResults=binnedScores)
  }
  
  results
}

# calculate bias, variance, and MSE
mse <- function(truth, est, weights=NULL){
  if(!is.null(weights))
    weights = weights / sum(weights)
  
  res = est - truth
  thisBias = est - truth
  
  if(!is.null(weights)) {
    MSE = sum(res^2 * weights)
    bias=sum(res * weights)
    thisVar = (res - sum(res*weights))^2
    thisVar = sum(thisVar * weights)
    out = list(MSE=MSE, bias=bias, var=thisVar)
  }
  else {
    MSE = mean(res^2)
    bias=mean(res)
    thisVar = (res - mean(res))^2
    thisVar = mean(thisVar)
    out = list(MSE=MSE, bias=bias, var=thisVar)
  }
  
  out
}

# either include both lower and upper, or include either: 
#    - the joint estimate draw matrix
#    - estimates and variances (assumes gaussian)
# truth: the true empirical proportions of mortality rates within the regions or enumeration areas of interest
# lower: the lower end of the credible interval
# upper: the upper end of the credible interval
# estMat: a matrix of joint draws of estimates, with number of rows equal to the length of truth, a number of 
#         columns equal to the number of draws. If not included, a gaussian distribution is assumed.
# significance: the significance level of the credible interval. By default 80%
# doRandomReject, doFuzzyReject: based on https://www.jstor.org/stable/pdf/20061193.pdf
coverage = function(truth, est=NULL, var=NULL, lower=NULL, upper=NULL, 
                    estMat=NULL, significance=.8, returnIntervalWidth=FALSE, doRandomReject=FALSE, doFuzzyReject=TRUE){
  
  if(any(is.null(lower)) || any(is.null(upper))) {
    # if the user did not supply their own credible intervals, we must get them ourselves given the other information
    
    if(is.null(estMat) && (is.null(est) || is.null(var)))
    stop("either include both lower and upper, est and var, or estMat")
    
    if(!is.null(est) && !is.null(var) && is.null(estMat)) {
      # in this case, we must calculate lower and upper assuming gaussianity
      lower = qnorm((1 - significance) / 2, est, sqrt(var))
      upper = qnorm(1 - (1 - significance) / 2, est, sqrt(var))
    }
    else {
      # we don't have information about the predictive distribution, and don't assume normality. 
      # Instead, use the user supplied to probability matrix estMat
      
      # take the quantiles of the probability draws
      CIs = apply(estMat, 1, function(ps) {quantile(ps, probs=c((1 - significance) / 2, 1 - (1 - significance) / 2))})
      lower = CIs[1,]
      upper = CIs[2,]
    }
  }
  
  if(any(lower > upper)) {
    warning("lower > upper, reordering")
    tmp = lower
    wrongOrder = lower > upper
    lower[wrongOrder] = upper[wrongOrder]
    upper[wrongOrder] = tmp[wrongOrder]
  }
  
  res = lower <= truth & upper >= truth
  
  if(returnIntervalWidth)
    width = upper - lower
  
  if(doRandomReject) {
    # in this case, we sometimes randomly reject if the truth is at the edge of the coverage interval. First 
    # determine what values are at the edge of the intervals, then determine the probability of rejection 
    # for each, then randomly reject
    atLowerEdge = which(lower == truth)
    atUpperEdge = which(upper == truth)
    
    probRejectLower = sapply(atLowerEdge, function(i) {((1 - significance) / 2 - mean(estMat[i,] < lower[i])) / mean(estMat[i,] == lower[i])})
    probRejectUpper = sapply(atUpperEdge, function(i) {((1 - significance) / 2 - mean(estMat[i,] > upper[i])) / mean(estMat[i,] == upper[i])})
    
    if(doFuzzyReject) {
      rejectLower = probRejectLower
      rejectUpper = probRejectUpper
    } else {
      rejectLower = runif(length(atLowerEdge)) <= probRejectLower
      rejectUpper = runif(length(atUpperEdge)) <= probRejectUpper
    }
    
    res[atLowerEdge] = sapply(1:length(atLowerEdge), function(i) {min(res[atLowerEdge][i], (1-rejectLower[i]))})
    res[atUpperEdge] = sapply(1:length(atUpperEdge), function(i) {min(res[atUpperEdge][i], (1-rejectUpper[i]))})
    # res[atLowerEdge] = res[atLowerEdge] & (!rejectLower)
    # res[atUpperEdge] = res[atUpperEdge] & (!rejectUpper)
  }
  
  allResults = c(coverage=mean(res))
  if(returnIntervalWidth)
    allResults = c(allResults, width=mean(width))
  
  allResults
}

# truth: a vector of observations on the desired scale
# est: a vector of logit-scale predictions of the same length as truth 
# my.var: a vector of logit-scale predictive variances of the same length as truth
# estMat: if available, use these probability draws in the integration. Use this argument 
#         when a gaussian approximation to the (possibly transformed) posterior is unreasonable
crps <- function(truth, est=NULL, my.var=NULL, estMat=NULL){
  if(!is.null(est) && !is.null(my.var) && is.null(estMat)) {
    sig = sqrt(my.var)
    x0 <- (truth - est) / sig
    res <- sig * (1 / sqrt(pi) -  2 * dnorm(x0) - x0 * (2 * pnorm(x0) - 1))
    
    ## sign as in Held (2008)
    res <- -res
  }
  else {
    # Integrate numerically using estMat
    if(is.null(estMat))
      stop("must include either or both est and my.var, or estMat")
    
    # the following code was commented out since it has been modified to be computationally efficient
    # # compute the crps for this row of truth
    # crpsRow = function(rowI) {
    #   thisTruth = truth[rowI]
    #   
    #   # either build the predictive cdf assuming normality on the logit scale or from the 
    #   # empirical distribution given by estMat if the user supplies it
    #   if(is.null(estMat)) {
    #     thisEst = est[rowI]
    #     thisVar = my.var[rowI]
    #     thisCdf = function(ws) {pnorm(logit(ws), thisEst, sqrt(thisVar))}
    #   } else {
    #     thisCdf = ecdf(estMat[rowI,])
    #   }
    #   
    #   intFun = function(ws) {
    #     (thisCdf(ws) - (ws >= thisTruth))^2
    #   }
    #   
    #   if(is.null(estMat)) {
    #     # when integrating we will set bounds on the integral to be reasonable to avoid 
    #     # faulty "divergent integral" error. The bounds will be 20 standard errors out 
    #     # of the estimate, making sure to include the truth.
    #     lowerBound = max(0, min(thisTruth - .01, expit(thisEst - 20 * sqrt(thisVar))))
    #     upperBound = min(1, max(thisTruth + .01, expit(thisEst + 20 * sqrt(thisVar))))
    #     integrate(intFun, lowerBound, upperBound)$value
    #   }
    #   else {
    #     # since we are using the empirical distribution, there is a closed form for the integral
    #     ps = estMat[rowI,]
    #     allPoints = sort(c(ps, thisTruth))
    #     deltas = diff(allPoints)
    #     sum(deltas * intFun(allPoints[1:length(ps)]))
    #   }
    # }
    # 
    # crpsRow2 = function(rowI) {
    #   thisTruth = truth[rowI]
    #   
    #   # either build the predictive cdf assuming normality on the logit scale or from the 
    #   # empirical distribution given by estMat if the user supplies it
    #   if(is.null(estMat)) {
    #     thisEst = est[rowI]
    #     thisVar = my.var[rowI]
    #     thisCdf = function(ws) {pnorm(logit(ws), thisEst, sqrt(thisVar))}
    #   } else {
    #     # thisCdf = ecdf(estMat[rowI,])
    #     sorted = sort(estMat[rowI,])
    #     thisCdf = approxfun(sorted, (1:length(sorted))/length(sorted), 
    #                         method = "constant", yleft = 0, yright = 1, f = 0, ties = "ordered")
    #   }
    #   
    #   intFun = function(ws) {
    #     (thisCdf(ws) - (ws >= thisTruth))^2
    #   }
    #   
    #   if(is.null(estMat)) {
    #     # when integrating we will set bounds on the integral to be reasonable to avoid 
    #     # faulty "divergent integral" error. The bounds will be 20 standard errors out 
    #     # of the estimate, making sure to include the truth.
    #     lowerBound = max(0, min(thisTruth - .01, expit(thisEst - 20 * sqrt(thisVar))))
    #     upperBound = min(1, max(thisTruth + .01, expit(thisEst + 20 * sqrt(thisVar))))
    #     integrate(intFun, lowerBound, upperBound)$value
    #   }
    #   else {
    #     # since we are using the empirical distribution, there is a closed form for the integral
    #     allPoints = sort(c(sorted, thisTruth))
    #     firstGreater = match(TRUE, sorted >= thisTruth)
    #     if(is.na(firstGreater))
    #       allPoints = c(sorted, thisTruth)
    #     else if(firstGreater == 1)
    #       allPoints = c(thisTruth, sorted)
    #     else
    #       allPoints = c(sorted[1:(firstGreater - 1)], thisTruth, sorted[firstGreater:length(sorted)])
    #     
    #     deltas = diff(allPoints)
    #     sum(deltas * intFun(allPoints[1:length(sorted)]))
    #   }
    # }
    
    crpsRow = function(rowI) {
      thisTruth = truth[rowI]
      
      # build the predictive cdf assuming from the empirical distribution given by 
      # estMat
      
      # thisCdf = ecdf(estMat[rowI,])
      sorted = estMat[rowI,] # already sorted
      
      # since we are using the empirical distribution, there is a closed form for the integral
      allPoints = sort(c(sorted, thisTruth))
      deltas = diff(allPoints)
      firstGreater = match(TRUE, sorted >= thisTruth)
      vals = (1:length(sorted))/length(sorted)
      if(is.na(firstGreater))
        return(sum((vals)^2 * deltas))
      else if(firstGreater == 1)
        return(deltas[1] + sum((1-vals[1:(length(sorted)-1)])^2 * deltas[2:length(deltas)]))
      else {
        left = sum(vals[1:(firstGreater-1)]^2 * deltas[1:(firstGreater-1)])
        mid = sum((1 - vals[firstGreater-1])^2 * deltas[firstGreater])
        right = ifelse(firstGreater == length(vals), 0, sum((1 - vals[firstGreater:(length(vals)-1)])^2 * deltas[(firstGreater+1):length(deltas)]))
        return(left+mid+right)
      }
      
      # intFun = function(ws) {
      #   (thisCdf(ws) - (ws >= thisTruth))^2
      # }
    }
    
    if(!is.null(estMat))
      estMat = t(apply(estMat, 1, sort))
    res = sapply(1:length(truth), crpsRow)
  }
  
  mean(res)
}