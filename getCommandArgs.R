# this file contains functions for generating lists of command arguments for each model results function, 
# saving the files to be used later

# make the command arguments file for resultsSPDE
getGPCommandArgs = function(covType=c("exponential", "matern", "mixture"), rangeText=c("01", "05", "1", "mix")) {
  gpCommandArgs = list()
  i = 1
  for(i1 in 1:length(covType)) {
    thisCovType=covType[i1]
    
    for(i2 in 1:length(rangeText)) {
      thisRangeText = rangeText[i2]
      
      if(thisCovType == "mixture" && thisRangeText != "mix")
        next
      else if(thisCovType != "mixture" && thisRangeText == "mix")
        next
      
      gpCommandArgs[[i]] = list(covType=thisCovType, rangeText=thisRangeText)
      i=i+1
    }
  }
  
  save(gpCommandArgs, file="gpCommandArgs.RData")
}

# make the command arguments file for resultsSPDE
getSpdeCommandArgs = function(covType=c("exponential", "matern", "mixture"), rangeText=c("01", "05", "1", "mix")) {
  spdeCommandArgs = list()
  i = 1
  for(i1 in 1:length(covType)) {
    thisCovType=covType[i1]
    
    for(i2 in 1:length(rangeText)) {
      thisRangeText = rangeText[i2]
      
      if(thisCovType == "mixture" && thisRangeText != "mix")
        next
      else if(thisCovType != "mixture" && thisRangeText == "mix")
        next
      
      spdeCommandArgs[[i]] = list(covType=thisCovType, rangeText=thisRangeText)
      i=i+1
    }
  }
  
  save(spdeCommandArgs, file="spdeCommandArgs.RData")
}

# make the command arguments file for resultsLK
getLKCommandArgs = function(covType=c("exponential", "matern", "mixture"), rangeText=c("01", "05", "1", "mix"), 
                            nLayer=c(3, 4)) {
  LKCommandArgs = list()
  i = 1
  for(i1 in 1:length(covType)) {
    thisCovType=covType[i1]
    
    for(i2 in 1:length(rangeText)) {
      thisRangeText = rangeText[i2]
      
      if(thisCovType == "mixture" && thisRangeText != "mix")
        next
      else if(thisCovType != "mixture" && thisRangeText == "mix")
        next
      
      for(i3 in 1:length(nLayer)) {
        thisnLayer = nLayer[i3]
        
        LKCommandArgs[[i]] = list(covType=thisCovType, rangeText=thisRangeText, nLayer=thisnLayer)
        i=i+1
      }
    }
  }
  
  save(LKCommandArgs, file="LKCommandArgs.RData")
}

# make the command arguments file for resultsLKINLA
getLKINLACommandArgs = function(covType=c("exponential", "matern", "mixture"), rangeText=c("01", "05", "1", "mix"), 
                            nLayer=c(3, 4)) {
  LKINLACommandArgs = list()
  i = 1
  for(i1 in 1:length(covType)) {
    thisCovType=covType[i1]
    
    for(i2 in 1:length(rangeText)) {
      thisRangeText = rangeText[i2]
      
      if(thisCovType == "mixture" && thisRangeText != "mix")
        next
      else if(thisCovType != "mixture" && thisRangeText == "mix")
        next
      
      for(i3 in 1:length(nLayer)) {
        thisnLayer = nLayer[i3]
        
        LKINLACommandArgs[[i]] = list(covType=thisCovType, rangeText=thisRangeText, nLayer=thisnLayer)
        i=i+1
      }
    }
  }
  
  save(LKINLACommandArgs, file="LKINLACommandArgs.RData")
}

# getGPCommandArgs()
# getSpdeCommandArgs()
# getLKCommandArgs()
# getLKINLACommandArgs()







