source("setupParallel.R")
index = as.numeric(commandArgs(trailingOnly = TRUE))
load("gpCommandArgs.RData")
argList = gpCommandArgs[[index]]
do.call("resultsGP", argList)
stopCluster(cl)