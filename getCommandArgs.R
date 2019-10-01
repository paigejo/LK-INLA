# this file contains functions for generating lists of command arguments for each model results function, 
# saving the files to be used later

# make the command arguments file for runAllSPDE.R
getSpdeCommandArgs = function(covType=c("exponential", "matern", "mixture"), rangeText=c("01", "05", "1", "")) {
  spdeCommandArgs = list()
  i = 1
  for(i1 in 1:length(covType)) {
    thisCovType=covType[i1]
    
    for(i2 in 1:length(rangeText)) {
      thisRangeText = rangeText[i2]
      
      spdeCommandArgs[[i]] = list(covType=thisCovType, rangeText=thisRangeText)
    }
  }
  
  save(spdeCommandArgs, file="spdeCommandArgs.RData")
}