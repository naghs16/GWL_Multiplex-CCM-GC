#'!/usr/bin/ Rscript

args <- commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("Usage: script.R <layer.type [FOR, INV]>", call.=FALSE)
}

layer.type <- as.character(args[1])

library('rEDM')
library('compiler')

SURR <- 50
layer.name <- "MON"
method.name <- "CCM"
sys.path <- "~/PostDoc/GWL/"
data.dir <- paste0(sys.path,"Data/MON/")

run.dir <- paste0(sys.path,"modeling/")
dir.create(run.dir, showWarnings = F)
run.dir <- paste0(sys.path,"modeling/", layer.name, "/")
dir.create(run.dir, showWarnings = F)
run.dir <- paste0(sys.path,"modeling/", layer.name, "/", method.name, "/")
dir.create(run.dir, showWarnings = F)

########Define Functions
E.function <- function(ts){
  E <- vector()
  
  lapply(1:N, function(i) {
    ts.used <- ts[,i]
    simplex_output <- rEDM::simplex(ts.used, lib, pred)
    E[i] <<- which.max(simplex_output$rho) #this is the bestE
  })
  
  return(E)
}
E.fun <- compiler::cmpfun(E.function)

Get_CCM_Mat <- function(ts){
  M_rho <- matrix(0, nrow = N, ncol = N)
  
  i <- 1
  runAgain <- TRUE
  while(runAgain)
  {
    runAgain <- FALSE
    lapply(1:N, function(j) {
      tss <- cbind(ts[,i], ts[,j])
      ccm.tmp <- rEDM::ccm(tss, lib, pred, E = E[j], lib_column = 1, target_column = 2, lib_sizes = M,  random_libs = FALSE, replace = FALSE, silent = TRUE)
      M_rho[i, j] <<- ccm.tmp$rho
    }
    )
    
    if(i < N){
      runAgain <- TRUE
      i <- i+1
    }
    
  }
  return(M_rho)
}
Get_CCM_Matrix <- compiler::cmpfun(Get_CCM_Mat)
########

for (sur in 0:SURR) {
  name.data.sur <- paste0(data.dir, "GWL-Tabriz_sur", sur, "_Layer-", layer.type, ".RDS")
  TS.MAT <- readRDS(name.data.sur)
  #TS.MAT is the time series and its format is matrix, i.e. the columns are the normalised data (mean and variance) for each station  
  N <- dim(TS.MAT)[2] #number of nodes
  M <- dim(TS.MAT)[1] #length of time series
  lib <- c(1, as.integer(2*M/3)) #this part of data will be used to construct the model (2/3), # portion of data to train
  pred <- c(as.integer(2*(M/3))+1, M) #this part of data will be used to test the model (1/3) # portion of data to predict
  
  E <- E.fun(TS.MAT)
  
  cormat.name <- paste0(run.dir, "GWL-Tabriz_", method.name, "-CORMAT-sur", sur, "_Layer-", layer.type, ".RDS")
  M_rho <- Get_CCM_Matrix(TS.MAT)
  saveRDS(M_rho, cormat.name)
}
