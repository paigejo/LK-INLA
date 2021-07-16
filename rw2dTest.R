# Based on example obtained from: 
# https://inla.r-inla-download.org/r-inla.org/doc/latent/rw2d-example.R
nrow=50
ncol=25
set.seed(123)
n = nrow*ncol
s.mat=matrix(NA,nrow=nrow,ncol=ncol)
j=1:ncol
for(i in 1:nrow)
  s.mat[i,j] = 0.1*(i+2*j)

## a covariate
z.mat=matrix(runif(nrow*ncol),nrow,ncol)

## noise
noise.mat=matrix(rnorm(nrow*ncol, sd=0.3),nrow,ncol)

## make simulated data
y.mat = s.mat  + 0.5*z.mat + noise.mat

## convert matrices to the internal representation in INLA
y = inla.matrix2vector(y.mat)
z = inla.matrix2vector(z.mat)
node = 1:n
formula= y ~ z + f(node, model="rw2d", nrow=nrow, ncol=ncol)
data=data.frame(y=y,z=z,node=node)
## fit the model
result=inla(formula, family="gaussian", data=data)

# now try with a stack
matchWith2dKnots = function(coords, xlim=range(coords[,1]), ylim=range(coords[,2]), 
                            nrow=30, ncol=30) {
  # generate breaks and knots (cell centers) in x and y directions
  xBreaks = seq(xlim[1], xlim[2]+.000001, l=ncol+1)
  yBreaks = seq(ylim[1], ylim[2]+.000001, l=nrow+1)
  xKnots = xBreaks[1:ncol] + diff(xBreaks) / 2
  yKnots = yBreaks[1:nrow] + diff(yBreaks) / 2
  
  # match coordinates with the cell
  xInd = sapply(coords[,1], function (x) {(length(xBreaks) - match(TRUE, x >= rev(xBreaks))) + 1})
  yInd = sapply(coords[,2], function (y) {(length(yBreaks) - match(TRUE, y >= rev(yBreaks))) + 1})
  ind = (xInd - 1)*nrow + yInd
  list(ind=ind, xInd=xInd, yInd=yInd)
}
coords = cbind(c(col(z.mat)), c(row(z.mat)))
out = matchWith2dKnots(coords, nrow=nrow, ncol=ncol, 
                       xlim=c(.5, ncol+.5), ylim=c(.5, nrow+.5))
inds = out$ind
out = matchWith2dKnots(coords, nrow=nrow, ncol=ncol)
inds2 = out$ind # inds and inds 2 are the same: identical(inds, inds2)
image(matrix(inds, nrow=nrow, ncol=ncol))
stack.est = inla.stack(A =list(1, 1), 
                       effects =list(z=z, rw2d=inds), 
                       data =list(y=y, link=1), 
                       tag ="est", 
                       remove.unused=FALSE)
dat = inla.stack.data(stack.est, remove.unused=FALSE)

node = 1:n
formula= y ~ -1 + z + f(rw2d, model="rw2d", nrow=nrow, ncol=ncol)
## fit the model
result2=inla(formula, family="gaussian", data=dat, 
            control.predictor=list(A=inla.stack.A(stack.est), compute=FALSE))

summary(result)
summary(result2)

head(result$summary.random$node)
head(result2$summary.random$rw2d)

#plot the posterior mean for `node' with the truth
dev.new()
INLA:::inla.display.matrix(s.mat)
dev.new()
INLA:::inla.display.matrix(INLA:::inla.vector2matrix(result$summary.random$node$mean,nrow,ncol))
INLA:::inla.display.matrix(INLA:::inla.vector2matrix(result2$summary.random$rw2d$mean,nrow,ncol))


