# this file contains functions for generating lists of command arguments for each model results function, 
# saving the files to be used later

# make the command arguments file for resultsSPDE
getGPCommandArgs = function(covType=c("exponential", "matern", "mixture"), rangeText=c("01", "05", "1", "")) {
  gpCommandArgs = list()
  i = 1
  for(i1 in 1:length(covType)) {
    thisCovType=covType[i1]
    
    for(i2 in 1:length(rangeText)) {
      thisRangeText = rangeText[i2]
      
      if(covType == "mixture" && rangeText != "")
        next
      else if(covType != "mixture" && rangeText == "")
        next
      
      gpCommandArgs[[i]] = list(covType=thisCovType, rangeText=thisRangeText)
    }
  }
  
  save(gpCommandArgs, file="gpCommandArgs.RData")
}

# make the command arguments file for resultsSPDE
getSpdeCommandArgs = function(covType=c("exponential", "matern", "mixture"), rangeText=c("01", "05", "1", "")) {
  spdeCommandArgs = list()
  i = 1
  for(i1 in 1:length(covType)) {
    thisCovType=covType[i1]
    
    for(i2 in 1:length(rangeText)) {
      thisRangeText = rangeText[i2]
      
      if(covType == "mixture" && rangeText != "")
        next
      else if(covType != "mixture" && rangeText == "")
        next
      
      spdeCommandArgs[[i]] = list(covType=thisCovType, rangeText=thisRangeText)
    }
  }
  
  save(spdeCommandArgs, file="spdeCommandArgs.RData")
}

# make the command arguments file for resultsLK
getLKCommandArgs = function(covType=c("exponential", "matern", "mixture"), rangeText=c("01", "05", "1", ""), 
                            nLayer=c(3, 4)) {
  LKCommandArgs = list()
  i = 1
  for(i1 in 1:length(covType)) {
    thisCovType=covType[i1]
    
    for(i2 in 1:length(rangeText)) {
      thisRangeText = rangeText[i2]
      
      if(covType == "mixture" && rangeText != "")
        next
      else if(covType != "mixture" && rangeText == "")
        next
      
      for(i3 in 1:length(nLayer)) {
        thisnLayer = nLayer[i3]
        
        LKCommandArgs[[i]] = list(covType=thisCovType, rangeText=thisRangeText, nLayer=thisnLayer)
      }
    }
  }
  
  save(LKCommandArgs, file="LKCommandArgs.RData")
}

# make the command arguments file for resultsLKINLA
getLKINLACommandArgs = function(covType=c("exponential", "matern", "mixture"), rangeText=c("01", "05", "1", ""), 
                            nLayer=c(3, 4)) {
  LKINLACommandArgs = list()
  i = 1
  for(i1 in 1:length(covType)) {
    thisCovType=covType[i1]
    
    for(i2 in 1:length(rangeText)) {
      thisRangeText = rangeText[i2]
      
      if(covType == "mixture" && rangeText != "")
        next
      else if(covType != "mixture" && rangeText == "")
        next
      
      for(i3 in 1:length(nLayer)) {
        thisnLayer = nLayer[i3]
        
        LKINLACommandArgs[[i]] = list(covType=thisCovType, rangeText=thisRangeText, nLayer=thisnLayer)
      }
    }
  }
  
  save(LKINLACommandArgs, file="LKINLACommandArgs.RData")
}









