#' Local measures

# args <- commandArgs(trailingOnly=TRUE)
# 
# if (length(args)==0) {
#   stop("Usage: script.R <method.name [GC,TE,CCM]> <layer.type [FOR,INV,MUX]> <z [1:5]>", call.=FALSE)
# }

# method.name <- as.character(args[1])
# layer.type <- as.character(args[2])
# z <- as.numeric(args[3])

method.name <- "GC"
layer.type <- "MUX"
z <- 3

library('igraph')

cnt <- 50
layer.name <- "MON"
sys.path <- "D:/PostDoc/GWL/"
g.dir <- paste0(sys.path, "Analysis/GRAPH/")
g.name <- paste0(g.dir, "GWL-Tabriz-WFFT_", method.name, "_Layer-", layer.name, "_undir-GRAPH_Surr", cnt, "_Z", z, "_", layer.type, ".RDS")
com.graph <- readRDS(g.name)
ana.dir <- paste0(sys.path, "Analysis/Measures/") 
dir.create(ana.dir, showWarnings = F)
cen.count <- 11
cen.mat <- matrix(0,vcount(com.graph), cen.count)

cen.mat[,1] <- igraph::degree(com.graph, v = V(com.graph), mode = c("all"), loops = FALSE, normalized = TRUE)
cen.mat[,2] <- igraph::graph.strength(com.graph, vids = V(com.graph), mode = c("all"), loops = FALSE) #weighted degree
cen.mat[,3] <- igraph::betweenness(com.graph, v = V(com.graph), directed = FALSE, nobigint = TRUE, normalized = TRUE)
#cen.mat[,4] <- igraph::closeness(com.graph, vids = V(com.graph), mode = c("all"), normalized = TRUE)
cen.mat[,5] <- igraph::eigen_centrality(com.graph, directed = FALSE, scale = TRUE)$vector#scale:to have a maximum score of one
cen.mat[,6] <- igraph::alpha_centrality(com.graph, nodes = V(com.graph), alpha = 1)
#cen.mat[,7] <- igraph::power_centrality(com.graph, nodes = V(com.graph), loops = FALSE, rescale = TRUE, sparse = TRUE)
cen.mat[,8] <- igraph::page_rank(com.graph, vids = V(com.graph), directed = FALSE)$vector
cen.mat[,9] <- 1/igraph::eccentricity(com.graph, vids = V(com.graph), mode = c("all"))
cen.mat[,10] <- igraph::authority_score(com.graph, scale = TRUE)$vector #for undirected graphs this score is same as hub score
cen.mat[,11] <- igraph::subgraph_centrality(com.graph, diag = FALSE)

nam.mat <- paste0(ana.dir, "GWL-Tabriz-WFFT_", method.name,"_undir-GRAPH-Local-measures-Surr", cnt, "_Z", z, "_", layer.type, ".RDS")
saveRDS(cen.mat,nam.mat)