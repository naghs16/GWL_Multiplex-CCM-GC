#'!/usr/bin/ Rscript

# args <- commandArgs(trailingOnly=TRUE)
# 
# if (length(args)==0) {
#   stop("Usage: script.R <method.name [GC,CCM]>", call.=FALSE)
# }
# 
# method.name <- as.character(args[1])

method.name <- "CCM"

#Z_Score <- c(2.576, 2.326, 1.96, 1.645, 1.28) #z:1,2,5,10,20%
layer.type <- "FOR"
layer.name <- "MON"
sys.path <- "D:/PostDoc/GWL/"
z.dir <- paste0(sys.path, "Analysis/")
dir.create(z.dir, showWarnings = F)
z.dir <- paste0(sys.path, "Analysis/ZSCORE/")
dir.create(z.dir, showWarnings = F)
run.dir <- paste0(sys.path,"modeling/", layer.name, "/", method.name, "/")
sur <- 0
cormat.name <- paste0(run.dir, "GWL-Tabriz-WFFT_", method.name, "-CORMAT-sur", sur, "_Layer-", layer.type, ".RDS")
cor.mat <- readRDS(cormat.name)

N <- dim(cor.mat)[2]
Surr <- 50
mu.mat <- matrix(0, N, N)
sig.mat <- matrix(0, N, N)
z_score <- matrix(0, N, N)
cnt.surr <- 0

for (i in 1:N) {
  for (j in 1:N) {
    
    rho.mean <- cor.mat[i,j]
    
    #keep the informative entries
    if(rho.mean < 0 || is.na(rho.mean) || i==j){
      cor.mat[i,j] <- 0
    }
    
  }
}

for (sur in 1:Surr) {
  cormat.name.surr <- paste0(run.dir, "GWL-Tabriz-WFFT_", method.name, "-CORMAT-sur", sur, "_Layer-", layer.type, ".RDS")
  cor.mat.surr <- readRDS(cormat.name.surr)
  
  ### For the CCM method, this part needs to be active
  for (i in 1:N) {
    for (j in 1:N) {
      
      rho.mean <- cor.mat.surr[i,j]
      
      #keep the informative entries
      if(rho.mean < 0 || is.na(rho.mean) || i==j){
        cor.mat.surr[i,j] <- 0
      }
      
    }
  }
  
  mu.mat <- mu.mat + cor.mat.surr
  sig.mat <- sig.mat + cor.mat.surr^2
  
  #this counter allows to normalize correctly (accounting for possible missing files)
  cnt.surr <- cnt.surr + 1
  cat(paste("Surr", Surr, "/", cnt.surr, "\n"))
}

c.ave <- mu.mat/cnt.surr
c.std <- sqrt( (sig.mat/cnt.surr) - (c.ave)^2 )

z_score <- (cor.mat - c.ave)/c.std

for (i in 1:N) {
  for (j in 1:N) {
    
    if(i==j || is.na(z_score[i,j]) ) {
      z_score[i,j] <- 0
    }
    
  }
}

#Save results
saveRDS(z_score, paste0(z.dir, "GWL-Tabriz-WFFT_", method.name, "_Layer-", layer.name, "_ZSCORE-Surr", cnt.surr, "_", layer.type, ".RDS"))