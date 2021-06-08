# simple rgeneric example for testing

inla.rgeneric.ar1.model = function(
  cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
          "log.prior", "quit"),
  theta = NULL) {
  
  # for reference and potential storage for objects to
  # cache, this is the environment of this function
  # which holds arguments passed as `...` in
  # `inla.rgeneric.define()`.
  envir = parent.env(environment())
  interpret.theta = function() { 
    return(list(prec = exp(theta[1L]),
                rho = 2 * exp(theta[2L])/(1 + exp(theta[2L])) - 1))
  }
  
  graph = function() {
    return(Q())
  }
  
  Q = function() {
    p = interpret.theta()
    i = c(1L, n, 2L:(n - 1L), 1L:(n - 1L))
    j = c(1L, n, 2L:(n - 1L), 2L:n)
    x = p$prec/(1 - p$rho^2) *
      c(1L, 1L, rep(1 + p$rho^2, n - 2L), rep(-p$rho, n - 1L))
    return(sparseMatrix(i = i, j = j, x = x, giveCsparse = FALSE))
  }
  
  mu = function() { return(numeric(0)) }
  
  log.norm.const = function() { return (numeric(0)) }
  
  log.prior = function() {
    p = interpret.theta()
    val = dgamma(p$prec, shape = 1, rate = 1, log=TRUE) + theta[1L] +
      dnorm(theta[2L], mean = 0, sd = 1, log=TRUE)
    return(val)
  }
  
  initial = function() { return (rep(1, 2)) }
  
  quit = function() { return (invisible()) }
  # sometimes this is useful, as argument 'graph' and 'quit'
  # will pass theta=NULL as the values of theta are not
  # required for defining the graph. however, this statement
  # will ensure that theta is always defined.
  if (is.null(theta)) theta = initial()
  
  val = do.call(match.arg(cmd), args = list())
  return(val)
}

n = 100
rho=0.9
set.seed(123)
x = arima.sim(n, model = list(ar = rho)) * sqrt(1-rho^2)
y = x + rnorm(n, sd = 0.1)
model = inla.rgeneric.define(inla.rgeneric.ar1.model, n=n)
formula = y ~ -1 + f(idx, model=model)
r = inla(formula, data = data.frame(y, idx = 1:n), verbose=TRUE)


fformula = y ~ -1 + f(idx, model = "ar1",
                      hyper = list(prec = list(prior = "loggamma", param = c(1,1)),
                                   rho = list(prior = "normal", param = c(0,1))))
rr = inla(fformula, data = data.frame(y, idx = 1:n))
