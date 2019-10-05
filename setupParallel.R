# sets up global variables and libraries to run code in parallel
assign("doParallel", TRUE, envir=.GlobalEnv)

library(parallel)
library(doParallel)
library(foreach)
setwd("~/git/LK-INLA/")
source("setup.R")

# assign("cl", detectCores() - 1, envir=.GlobalEnv)
assign("cl", makeCluster(7, outfile="log.txt"), envir=.GlobalEnv)
inf = sessionInfo()
if(inf$platform != "x86_64-w64-mingw32/x64 (64-bit)") {
  clusterEvalQ(cl, {setwd("~/git/LK-INLA/"); source("setup.R")})
} else {
  clusterEvalQ(cl, {setwd("U:/git/LK-INLA/"); source("setup.R")})
}
registerDoParallel(cl)
options(error=recover)
# stopCluster(cl)