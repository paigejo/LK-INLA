# this script is for generating the data sets for the simulation study

# generate the data sets that will be used in the simulation study
getSimulationDataSets = function(nTotal=1000, nTest=100, marginalVar=1, errorVar=sqrt(.1), nDataSets=100, seed=123) {
  if(is.null(seed))
    set.seed(seed)
  
  # set the correlation functions for the simulation study
  exponentialCorFuns = list(function(x) {Exp.cov(x, theta=0.1)}, 
                            function(x) {Exp.cov(x, theta=0.5)}, 
                            function(x) {Exp.cov(x, theta=1)})
  maternCorFuns = list(function(x) {stationary.cov(x, theta=0.1, Covariance="Matern", nu=1)}, 
                       function(x) {stationary.cov(x, theta=0.5, Covariance="Matern", nu=1)}, 
                       function(x) {stationary.cov(x, theta=1, Covariance="Matern", nu=1)})
  mixtureCorFun = function(x) {0.4 * Exp.cov(x, theta=0.1) + 0.6 * Exp.cov(x, theta=3)}
  
  ##### simulate the corresponding data sets
  ## simulate the exponential correlation function data sets
  print("Simulating dataset with first exponential correlation function...")
  simulationData = getSimulationDataSetsGivenCovariance(exponentialCorFuns[[1]], nTotal=nTotal, nTest=nTest, marginalVar=marginalVar, errorVar=errorVar, 
                                                        nDataSets=nDataSets, plotNameRoot="(Exp(0.1))", fileNameRoot="exp01")
  print("Saving dataset...")
  save(simulationData, file="exponential01DataSet.RData")
  
  print("Simulating data set with second exponential correlation function...")
  simulationData = getSimulationDataSetsGivenCovariance(exponentialCorFuns[[2]], nTotal=nTotal, nTest=nTest, marginalVar=marginalVar, errorVar=errorVar, 
                                                        nDataSets=nDataSets, plotNameRoot="(Exp(0.5))", fileNameRoot="exp05")
  print("Saving dataset...")
  save(simulationData, file="exponential05DataSet.RData")
  
  print("Simulating data set with third exponential correlation function...")
  simulationData = getSimulationDataSetsGivenCovariance(exponentialCorFuns[[3]], nTotal=nTotal, nTest=nTest, marginalVar=marginalVar, errorVar=errorVar, 
                                                        nDataSets=nDataSets, plotNameRoot="(Exp(1))", fileNameRoot="exp1")
  print("Saving dataset...")
  save(simulationData, file="exponential1DataSet.RData")
  
  ## simulate the Matern correlation function data sets
  print("Simulating dataset with first matern correlation function...")
  simulationData = getSimulationDataSetsGivenCovariance(maternCorFuns[[1]], nTotal=nTotal, nTest=nTest, marginalVar=marginalVar, errorVar=errorVar, 
                                                        nDataSets=nDataSets, plotNameRoot="(Matern(1,0.1))", fileNameRoot="mat01")
  print("Saving dataset...")
  save(simulationData, file="matern01DataSet.RData")
  
  print("Simulating data set with second matern correlation function...")
  simulationData = getSimulationDataSetsGivenCovariance(maternCorFuns[[2]], nTotal=nTotal, nTest=nTest, marginalVar=marginalVar, errorVar=errorVar, 
                                                        nDataSets=nDataSets, plotNameRoot="(Matern(1,.5))", fileNameRoot="mat05")
  print("Saving dataset...")
  save(simulationData, file="matern05DataSet.RData")
  
  print("Simulating data set with third matern correlation function...")
  simulationData = getSimulationDataSetsGivenCovariance(maternCorFuns[[3]], nTotal=nTotal, nTest=nTest, marginalVar=marginalVar, errorVar=errorVar, 
                                                        nDataSets=nDataSets, plotNameRoot="(Matern(1,1))", fileNameRoot="mat1")
  print("Saving dataset...")
  save(simulationData, file="matern1DataSet.RData")
  
  # simulate the mixture correlation function data set
  print("Simulating data set with mixture correlation function...")
  simulationData = getSimulationDataSetsGivenCovariance(mixtureCorFun, nTotal=nTotal, nTest=nTest, marginalVar=marginalVar, errorVar=errorVar, 
                                                        nDataSets=nDataSets, plotNameRoot="(0.4*Exp(0.1) + 0.6*Exp(3))", fileNameRoot="mix")
  print("Saving dataset...")
  save(simulationData, file="mixtureDataSet.RData")
}

getSimulationDataSetsGivenCovariance = function(corFun, nTotal=1000, nTest=100, marginalVar=1, errorVar=sqrt(.1), nDataSets=100, 
                                                printEvery=10, saveDataSetPlot=TRUE, fileNameRoot="", plotNameRoot="") {
  
  # set the spatial domain
  xRange = c(-1, 1)
  yRange = c(-1, 1)
  
  # simulate observation spatial locations
  xValues = matrix(runif(nTotal * nDataSets, xRange[1], xRange[2]), nrow=nTotal, ncol=nDataSets)
  yValues = matrix(runif(nTotal * nDataSets, yRange[1], yRange[2]), nrow=nTotal, ncol=nDataSets)
  
  # preallocate observation matrix, and pregenerate standard normal draws
  observations = matrix(nrow=nTotal, ncol=nDataSets)
  zsims = matrix(rnorm(nTotal * nDataSets), nrow=nTotal, ncol=nDataSets)
  
  # generate spatial component of observation values
  for(i in 1:nDataSets) {
    if(i %% printEvery == 0 || i == 1)
      print(paste0("Simulating data set ", i, "/", nDataSets))
    
    thisx = xValues[,i]
    thisy = yValues[,i]
    L = t(chol(corFun(cbind(thisx, thisy))))
    observations[,i] = L %*% zsims[,i]
  }
  
  # scale by marginal standard deviation and add in error variance
  observations = observations * sqrt(marginalVar) + matrix(rnorm(nTotal * nDataSets, sd=sqrt(errorVar)), nrow=nTotal, ncol=nDataSets)
  
  # separate out test and train results
  trainI = 1:(nTotal - nTest)
  testI = (nTotal - nTest + 1):nTotal
  xTest = xValues[testI,]
  xTrain = xValues[trainI,]
  yTest = yValues[testI,]
  yTrain = yValues[trainI,]
  zTest = observations[testI,]
  zTrain = observations[trainI,]
  
  # put relevant values into a list
  out = list(xTrain=xTrain, yTrain=yTrain, zTrain=zTrain, xTest=xTest, yTest=yTest, zTest=zTest, 
             corFun=corFun, marginalVar=marginalVar, errorVar=errorVar, xRange=xRange, yRange=yRange)
  
  # plot the results
  plotExampleDataSets(out, saveDataSetPlot=saveDataSetPlot, plotNameRoot=plotNameRoot, fileNameRoot=fileNameRoot)
  
  # return the results
  out
}

plotExampleDataSets = function(simulatedDataSets, saveDataSetPlot=TRUE, plotNameRoot="", fileNameRoot="", nx=100, ny=100) {
  # load and variable names
  xTrain = simulatedDataSets$xTrain
  yTrain = simulatedDataSets$yTrain
  zTrain = simulatedDataSets$zTrain
  xRange = simulatedDataSets$xRange
  yRange = simulatedDataSets$yRange
  
  if(saveDataSetPlot)
    pdf(file=paste0("Figures/sampleDataSet_", fileNameRoot, ".pdf"))
  
  quilt.plot(xTrain[,1], yTrain[,1], zTrain[,1], xlim=xRange, ylim=yRange, 
             main=paste("Example simulated dataset", plotNameRoot), nx=nx, ny=ny)
  
  if(saveDataSetPlot)
    dev.off()
}
