# script for converting between LatticeKrig and ELK in terms of:
# lattice basis
# parameters

ELKBasis2LKinfo = function(latInfo, nu=1.5, lambdaStart=.1, kappa.start=1, a.wghtStart=4+kappa.start^2, 
                           normalize=TRUE, fixedFunctionArgs=list(m = 1)) {
  
  domainCoords = cbind(latInfo[[1]]$xRangeDat, latInfo[[1]]$yRangeDat)
  NC = latInfo[[1]]$NC
  nBuffer = length(unique(latInfo[[1]]$latCoords[latInfo[[1]]$latCoords[,2] > latInfo[[1]]$yRangeDat[2],2]))
  nLayer = length(latInfo)
  
  # set up the lattice, the arguments to LatticeKrig
  LKinfo = LKrigSetup(domainCoords, nlevel=nLayer, nu=nu, NC=NC, normalize=normalize, NC.buffer=nBuffer, 
                      lambda=lambdaStart, a.wght=a.wghtStart, fixedFunctionArgs=fixedFunctionArgs, alpha=rep(NA, nLayer))
  
  LKinfo
}

LKinfo2ELKBasis = function(LKinfo) {
  nlevel = LKinfo$nlevel
  nBuffer = LKinfo$latticeInfo$NC.buffer
  NC = max(LKinfo$latticeInfo$mx[1,]) - 2*nBuffer
  nlevel = LKinfo$nlevel
  xRange = LKinfo$latticeInfo$rangeLocations[,1]
  yRange = LKinfo$latticeInfo$rangeLocations[,2]
  
  latInfo = makeLatGrids(xRange, yRange, NC=NC, nLayer=nlevel, nBuffer=nBuffer)
  latInfo
}