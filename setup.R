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

setwd("~/git/LK-INLA/")
source('compareModels.R')
source('LKinla.R')
source('LKinla_rgeneric.R')
source('modGP.R')
source('modSPDE.R')
source('modLKinla.R')
source('modLK.R')
source('modelResults.R')
source('getCommandArgs.R')
# source('mapOver.R')
source('utilityFuns.R')
# source('~/git/M9/exploratoryAnalysisFuns.R')

inf = sessionInfo()
if(inf$platform != "x86_64-apple-darwin15.6.0 (64-bit)") {
  INLA:::inla.dynload.workaround()
  # avoid setting too many threads and thereby using too much memory
  inla.setOption(num.threads=1)
  options(error=traceback)
} else {
  options(error=recover)
}

# parallelization
if(!exists("doParallel") || (exists("doParallel") && doParallel == FALSE)) {
  assign("doParallel", FALSE, envir=.GlobalEnv)
  assign("cores", NULL, envir=.GlobalEnv)
  assign("cl", NULL, envir=.GlobalEnv)
}