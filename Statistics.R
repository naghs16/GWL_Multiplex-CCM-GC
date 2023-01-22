# args <- commandArgs(trailingOnly=TRUE)
# 
# if (length(args)==0) {
#   stop("Usage: script.R <group [I, II, III, IV]> <layer.name [MON]> <layer.type [FOR, INV]>", call.=FALSE)
# }
rm(list=ls(all=TRUE))
Sys.setlocale(locale = "persian")

group <- c("I", "II", "III","IV")
layer.name <- "MON"
layer.type <- "FOR"

library('rEDM')
library('compiler')
library(e1071)  

II <- 1 #No of the files in subDir, CHANGE


stat.data <- data.frame()
for (i in 1:length(group)) {
  mainDir <- paste0("D:/Thesis_Ph.D/GWL_Tabriz/Data/Selection/", group[i], "/missing/")
  data.files <- list.files(path = mainDir)
  name.file <- paste0(mainDir, data.files[II])
  ts <- read.csv(name.file)
a1 <- mean(ts[,4])
a2 <- min(ts[,4])
a3 <- max(ts[,4])
a4 <- sd(ts[,4])
a5 <- var(ts[,4])
a6 <- skewness(ts[,4])
a7 <- kurtosis(ts[,4])
stat.data <- rbind(stat.data, data.frame(a1=a1,a2=a2,a3=a3,a4=a4,a5=a5,a6=a6,a7=a7))
}