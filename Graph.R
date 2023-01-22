#'!/usr/bin/ Rscript

# args <- commandArgs(trailingOnly=TRUE)
# 
# if (length(args)==0) {
#   stop("Usage: script.R <method.name [GC,CCM]> <layer.type [FOR,INV,MUX]> <z [1:5]>", call.=FALSE)
# }
# 
# method.name <- as.character(args[1])
# layer.type <- as.character(args[2])
# z <- as.numeric(args[3])

method.name <- "GC"#GC,CCM
layer.type <- "MUX"#FOR, MUX
z <- 3

library('igraph')

layer.name <- "MON"
sys.path <- "D:/PostDoc/GWL/"
g.dir <- paste0(sys.path, "Analysis/GRAPH/")
dir.create(g.dir, showWarnings = F)
z.dir <- paste0(sys.path, "Analysis/ZSCORE/")
cnt <- 50 
if (layer.type=="MUX"){
  z.name <- paste0(z.dir, "GWL-Tabriz-WFFT_", method.name, "_Layer-", layer.name, "_ZSCORE-Surr", cnt, "_Z", z, "_", layer.type, ".RDS")
} else {
  z.name <-  paste0(z.dir, "GWL-Tabriz-WFFT_", method.name, "_Layer-", layer.name, "_ZSCORE-Surr", cnt, "_", layer.type, ".RDS")
}

z.score <- readRDS(z.name)
Z_Score <- c(2.576, 2.326, 1.96, 1.645, 1.28) #z:1,2,5,10,20%

g.name <- paste0(g.dir, "GWL-Tabriz-WFFT_", method.name, "_Layer-", layer.name, "_undir-GRAPH_Surr", cnt, "_Z", z, "_", layer.type, ".RDS")
graph_z_score <- igraph::graph.adjacency(z.score, weighted=T, mode="directed")
graph_z_scoreMo <- igraph::delete.edges(graph_z_score, E(graph_z_score)[weight < Z_Score[z]])

saveRDS(graph_z_scoreMo, g.name)