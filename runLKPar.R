source("setupParallel.R")
index = as.numeric(commandArgs(trailingOnly = TRUE))
load("LKCommandArgs.RData")
argList = LKCommandArgs[[index]]
do.call("resultsLK", argList)
stopCluster(cl)