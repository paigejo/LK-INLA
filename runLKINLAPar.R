source("setupParallel.R")
index = as.numeric(commandArgs(trailingOnly = TRUE))
load("LKINLACommandArgs.RData")
argList = LKINLACommandArgs[[index]]
do.call("resultsLKINLA", argList)
stopCluster(cl)