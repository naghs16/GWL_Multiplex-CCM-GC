#'!/usr/bin/ Rscript

# args <- commandArgs(trailingOnly=TRUE)
# 
# if (length(args)==0) {
#   stop("Usage: script.R <layer.type [FOR, INV]>", call.=FALSE)
# }
# 
# layer.type <- as.character(args[1])
layer.type <- "INV"#FOR or INV

library('grangers')
library('compiler')

SURR <- 50
layer.name <- "MON"
method.name <- "GC"
sys.path <- "D:/PostDoc/GWL/"
data.dir <- paste0(sys.path,"Data/MON/")

run.dir <- paste0(sys.path,"modeling/")
dir.create(run.dir, showWarnings = F)
run.dir <- paste0(sys.path,"modeling/", layer.name, "/")
dir.create(run.dir, showWarnings = F)
run.dir <- paste0(sys.path,"modeling/", layer.name, "/", method.name, "/")
dir.create(run.dir, showWarnings = F)

########Define Functions
Get_GC_Mat <- function(TS.MAT){
  
  N <- dim(TS.MAT)[2] #number of nodes
  M <- dim(TS.MAT)[1] #length of time series
  M_rho <- matrix(0, nrow = N, ncol = N)
  
  i <- 1
  runAgain <- TRUE
  while(runAgain)
  {
    runAgain <- FALSE
    for (j in 1:N) {
      
      if (i!=j) {
        
        te.tmp <- grangers::bc_test_uncond(TS.MAT[,i], TS.MAT[,j], ic.chosen = "AIC", max.lag = min(4, M - 1))
        M_rho[i,j] <- max(te.tmp$`F-test`)
        
      } else {
        next
      }
      
    } 
      
    if(i < N){
      runAgain <- TRUE
      i <- i+1
    } 
  }
  
  return(M_rho)
}
Get_GC_Matrix <- compiler::cmpfun(Get_GC_Mat)
########

for (sur in 0:SURR) {
  name.data.sur <- paste0(data.dir, "GWL-Tabriz-WFFT_sur", sur, "_Layer-", layer.type, ".RDS")
  TS.MAT <- readRDS(name.data.sur) #TS.MAT is the time series and its format is matrix, i.e. the columns are the normalised data (mean and variance) for each station  
 
  cormat.name <- paste0(run.dir, "GWL-Tabriz-WFFT_", method.name, "-CORMAT-sur", sur, "_Layer-", layer.type, ".RDS")
  M_rho <- Get_GC_Matrix(TS.MAT)
  saveRDS(M_rho, cormat.name)
}
