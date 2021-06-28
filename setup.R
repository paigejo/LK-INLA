# load required packages and R scripts
library(Matrix)
library(spam)
library(fields)
library(LatticeKrig)
library(invgamma)
library(latex2exp)
library(xtable)
library(profvis)
library(geosphere)
library(viridis)
library(sp)
library(raster)
library(MCMCpack)
library(numDeriv)
library(INLA)
library(spam64)


inf = sessionInfo()
if(inf$platform != "x86_64-apple-darwin15.6.0 (64-bit)" && inf$platform != "x86_64-w64-mingw32/x64 (64-bit)" && inf$platform != "x86_64-pc-linux-gnu (64-bit)") {
  INLA:::inla.dynload.workaround()
  # avoid setting too many threads and thereby using too much memory
  inla.setOption(num.threads=1)
  options(error=traceback)
  setwd("~/git/LK-INLA/")
} else if(inf$platform != "x86_64-w64-mingw32/x64 (64-bit)" && inf$platform != "x86_64-pc-linux-gnu (64-bit)") {
  setwd("~/git/LK-INLA/")
  options(error=recover)
} else if(inf$platform == "x86_64-pc-linux-gnu (64-bit)") {
  setwd("~/git/LK-INLA/")
  inla.setOption(num.threads=1) # consider raising
  options(error=recover)
} else {
  setwd("U:/git/LK-INLA/")
  inla.setOption(num.threads=1)
  options(error=recover)
}

source('LKinla.R')
source('LKinla_rgeneric.R')
source('modGP.R')
source('modSPDE.R')
source('modLKinla.R')
source('modLK.R')
source('modelResults.R')
source("scores.R")
source("test.R")
source('compareModels.R')
source("utilityFuns.R")
source('getCommandArgs.R')
# source('mapOver.R')
source('utilityFuns.R')
# source('~/git/M9/exploratoryAnalysisFuns.R')
source('getSimulationDataSets.R')
source('generateExampleResults.R')
source('plotGenerator.R')
source('validation.R')
source('pcBeta.R')
source('modBCEF.R')
source('LK2ELK.R')

# parallelization
if(!exists("doParallel") || (exists("doParallel") && doParallel == FALSE)) {
  assign("doParallel", FALSE, envir=.GlobalEnv)
  assign("cores", NULL, envir=.GlobalEnv)
  assign("cl", NULL, envir=.GlobalEnv)
}