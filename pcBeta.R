# PC prior for the beta distribution variance conditional on p from line 188 of https://static-content.springer.com/esm/art%3A10.1038%2Fncomms12145/MediaObjects/41467_2016_BFncomms12145_MOESM379_ESM.pdf
# NOTE: v in (0, p*(1-p))
dpcBeta = function(v, p=0.5, lambda, doLog=FALSE) {
  if(length(p) > 1)
    stop("p must be a single number, not a vector")
  # if(length(lambda) > 1)
  #   stop("lambda must be a single number, not a vector")
  v[v == 0] = 1e-6
  v[v == p*(1-p)] = p*(1-p) - 1e-6
  
  # reparameterization from Eq. (5) on line 116
  a = (p*(1-p)/v - 1)*p
  b = (p*(1-p)/v - 1)*(1-p)
  
  # there is a typo in the supplementary material: the square root sign should be applied to the distance within the exponential
  # term1 = log(lambda) + sqrt(-lambda*(digamma(a + b) - p*digamma(a) - (1-p)*digamma(b)))
  term1 = log(lambda) + -lambda*sqrt((digamma(a + b) - p*digamma(a) - (1-p)*digamma(b)))
  numer = log(p) + log(1-p) + log(abs(p^2 * trigamma(a) + (1-p)^2*trigamma(b) - trigamma(a+b)))
  denom = log(2) + 2*log(v) + 0.5 * log(digamma(a+b) - p*digamma(a) - (1-p)*digamma(b))
  starVal = log(abs(exp(numer)/exp(denom)))
  
  logDensity = term1 + starVal
  
  # normalization constant must be modified since there is a maximum KL distance possible: 
  # exponential distribution is truncated to account for this.
  # NOTE: in the parameterization from the supplementary materials, the distance is not (1/lambda) 
  vmax = p*(1-p) - 1e-6
  vmin = 1e-6
  amax = (p*(1-p)/vmax - 1)*p
  bmax = (p*(1-p)/vmax - 1)*(1-p)
  amin = (p*(1-p)/vmin - 1)*p
  bmin = (p*(1-p)/vmin - 1)*(1-p)
  kldmin = kldDistance(vmin, p)
  kldmax = kldDistance(vmax, p)
  distanceRange = range(c(kldmin, kldmax))
  # norm = -log(pexp(distanceRange[2], rate=lambda) - pexp(distanceRange[1], rate=lambda))
  norm = -(-lambda*distanceRange[1] + log(1-exp(-lambda*(distanceRange[2]-distanceRange[1]))))
  
  if(doLog) {
    out = logDensity + norm 
    out[v < 0] = -Inf
    out[v > p*(1-p)] = -Inf
  }
  else {
    out = exp(logDensity) * exp(norm)
    out[v <= 0] = 0
    out[v >= p*(1-p)] = 0
  }
  
  out
}

testdpcBeta = function(p=.5, lambda=c(.1, 0.5, 1, 10, 100), l = 1e3, subdivisions=100, rel.tol = .Machine$double.eps^0.25) {
  
  vs = seq(0, .25, l=l)
  densityMat = outer(vs, lambda, FUN=function(v, l) {dpcBeta(v, l, p=p)})
  
  for(i in 1:length(lambda)) {
    plot(vs, densityMat[,i], type="l", col="blue", main=paste0("Lambda: ", lambda[i], ", p: ", p))
  }
  
  print(paste0("Integrated probability mass for each lambda:"))
  print(colSums(densityMat)*(.25/l))
  
  print(paste0("Integrated probability mass for each lambda (using integrate()):"))
  
  print(sapply(lambda, 
         function(l) {
           integrand = function(v) {dpcBeta(v, lambda=l, p=p)}
           integrate(integrand, 0, p*(1-p), subdivisions=subdivisions, rel.tol=rel.tol)$value
         }
  ))
}

# PC prior for the beta distribution variance conditional on p from line 188 of https://static-content.springer.com/esm/art%3A10.1038%2Fncomms12145/MediaObjects/41467_2016_BFncomms12145_MOESM379_ESM.pdf
# NOTE: this is a reparameterization using overdispersion rho = v/[p(1-p)]
# NOTE: rho in (0, 1)
dpcBeta2 = function(rho, p=0.5, lambda, doLog=FALSE, normalize=FALSE) {
  v = rho * p * (1 - p)
  out = dpcBeta(v, p, lambda, doLog)
  if(doLog) {
    out = out + log(p*(1-p))
    if(normalize) {
      integrand = function(r) {dpcBeta2(r, lambda=lambda, p=p)}
      out = out - log(integrate(integrand, 0, 1)$value)
    }
  }
  else {
    out = out * p*(1-p)
    if(normalize) {
      integrand = function(r) {dpcBeta2(r, lambda=lambda, p=p)}
      out = out / integrate(integrand, 0, 1)$value
    }
  }
  
  out
}

testdpcBeta2 = function(p=.5, lambda=c(.25, 0.5, 1, 10, 100), l = 1e3, subdivisions=100, rel.tol = .Machine$double.eps^0.25, normalize=FALSE, saveResults=FALSE) {
  
  rho = seq(0, 1, l=l)
  # densityMat = outer(rho, lambda, FUN=function(r, l) {dpcBeta2(r, l, p=p, normalize=normalize)})
  densityMat = sapply(lambda, function(l) {dpcBeta2(rho, lambda=l, p=p, normalize=normalize)})
  cols = rainbow(length(lambda))
  
  if(saveResults)
    pdf("Figures/PCBetaTest.pdf", width=5, height=5)
  for(i in 1:length(lambda)) {
    if(i == 1)
      plot(rho, densityMat[,i], type="l", col=cols[i], main="PC prior for rho", ylim=c(0, max(densityMat)))
    else
      lines(rho, densityMat[,i], col=cols[i])
  }
  legend("top", paste0("lambda=", lambda), lty=1, col=cols)
  if(saveResults)
    dev.off()
  
  print(paste0("Integrated probability mass for each lambda:"))
  print(colSums(densityMat)*(1/l))
  
  print(paste0("Integrated probability mass for each lambda (using integrate()):"))
  
  print(sapply(lambda, 
               function(l) {
                 integrand = function(r) {dpcBeta2(r, lambda=l, p=p, normalize=normalize)}
                 integrate(integrand, 0, 1, subdivisions=subdivisions, rel.tol=rel.tol)$value
               }
  ))
}

# PC prior cdf for the beta distribution overdispersion from line 188 of https://static-content.springer.com/esm/art%3A10.1038%2Fncomms12145/MediaObjects/41467_2016_BFncomms12145_MOESM379_ESM.pdf
# NOTE: this is a reparameterization using overdispersion rho = v/[p(1-p)]
# NOTE: rho in (0, 1)
# NOTE: Since it's unclear what a closed form expression for the KLD-based distance is as a function of rho, 
#       the cdf is computed numerically using the integrate function
ppcBeta2 = function(q, p=0.5, lambda, doLog=FALSE, normalize=FALSE) {
  # define integrand and integrate
  integrand = function(r) {dpcBeta2(r, p=p, lambda=lambda)}
  out = integrate(integrand, 0, q)$value
  
  # make sure total integrated probability density is exactly one if desired by user (makes little difference)
  if(normalize)
    out = out * (1 / integrate(integrand, 0, 1)$value)
  
  # compute log if necessary
  if(doLog)
    log(out)
  else
    out
}

# PC prior inverse cdf for the beta distribution overdispersion from line 188 of https://static-content.springer.com/esm/art%3A10.1038%2Fncomms12145/MediaObjects/41467_2016_BFncomms12145_MOESM379_ESM.pdf
# NOTE: this is a reparameterization using overdispersion rho = v/[p(1-p)]
# NOTE: rho in (0, 1)
# NOTE: Since it's unclear what a closed form expression for the KLD-based distance is as a function of rho, 
#       the inverse cdf is computed numerically using the integrate function and uniroot algorithm
qpcBeta2 = function(prob, p=0.5, lambda, normalize=FALSE) {
  # define cdf
  cdf = function(q) {ppcBeta2(q, p=p, lambda=lambda, normalize=normalize)}
  # Use a root finding function to invert the cdf
  invcdf <- function(q){
    uniroot(function(x){cdf(x) - q}, c(0, 1))$root
  }
  
  sapply(prob, invcdf)
}

# Computes lambda based on rho (or logit rho by default) threshold U and exceedance probability alpha
getLambdapcBeta = function(U=1, logitU=TRUE, alpha=0.01, p=0.5, normalize=FALSE) {
  if(logitU)
    U = expit(U)
  exp(optim(0, function(llambda){(qpcBeta2(1-alpha, p=p, lambda=exp(llambda), normalize=normalize) - U)^2}, c(0, 1))$par)
}



# NOTE: (1/lambda) term not included in the distance. This matches the supplementary materials 
#       distance metric with its rate parameter lambda as used in the actual exponential 
#       distribution rather than the distance metric in the exact KLD distance.
kldDistance = function(v, p=0.5) {
  
  # reparameterization from Eq. (5) on line 116
  a = (p*(1-p)/v - 1)*p
  b = (p*(1-p)/v - 1)*(1-p)
  
  # NOTE: (1/lambda) term not included below to match the supplementary materials distance metric with its rate parameter lambda
  sqrt((digamma(a + b) - p*digamma(a) - (1-p)*digamma(b)))
}

# PC prior for the beta distribution logit overdispersion conditional on p from line 188 of https://static-content.springer.com/esm/art%3A10.1038%2Fncomms12145/MediaObjects/41467_2016_BFncomms12145_MOESM379_ESM.pdf
# NOTE: this is a reparameterization using overdispersion rho = v/[p(1-p)]
# NOTE: rho in (0, 1), but lrho = logit(rho) is in (-infty, infty)
dpcBetaLogit = function(lrho, p=0.5, lambda, doLog=FALSE, normalize=FALSE) {
  rho = expit(lrho)
  out = dpcBeta2(rho, p, lambda, doLog, normalize)
  if(doLog)
    out = out + c(sapply(lrho, function(z) {log(multivariateExpitJacobian(z))}))
  else
    out = out * c(sapply(lrho, function(z) {multivariateExpitJacobian(z)}))
  
  out
}

# instructions for defining a log prior density for a parameters on INLA's 
# internal scale is available at: http://www.r-inla.org/faq#TOC-I-want-to-use-a-prior-for-my-hyperparameter-that-seems-not-to-be-implemented-what-should-I-do-
# This function returns a table of logit overdispersion values and log 
# densities for the pc prior for the logit overdispersion parameter of the 
# Beta distribution. To use this in a call to INLA, just use something like 
# control.family = list(hyper = list(prec = list(prior = prior.table))), 
# where prior.table is the results of this function.
getpcBetaLogitTableForINLA = function(lambda, p=0.5, tailProb=1e-3, 
                                      rhoRange=qpcBeta2(c(tailProb, 1-tailProb), lambda=lambda, p=p), 
                                      n=100) {
  logitRhoRange = logit(rhoRange)
  logitRhoVals = seq(logitRhoRange[1], logitRhoRange[2], l=n)
  tabVals = dpcBetaLogit(logitRhoVals, p, lambda, doLog=TRUE, normalize=TRUE)
  
  ## use suitable support points x
  lprec = seq(-10, 10, len=100)
  ## link the x and corresponding y values into a string which begins with "table:""
  INLA:::inla.paste(c("table:", cbind(logitRhoVals, tabVals)))
}

# # on alternative prior for the beta overdispersion parameter
# 
# # set median at .04 and upper 97.5th pctile at 0.2
# mu = logit(0.04)
# prec = 1/((logit(.2)-logit(.04))/qnorm(.975))^2
dlogitNormal = function(x, mu=logit(0.04), prec=1/((logit(.2)-logit(.04))/qnorm(.975))^2) {
  z = logit(x)
  J = sapply(z, function(y) {abs(multivariateExpitJacobian(y))})
  out = dnorm(z, mu, sqrt(1/prec)) * (1/J)
  out[x <= 0] = 0
  out[x >= 1] = 0
  out
}








