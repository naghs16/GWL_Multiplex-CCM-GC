#' Computation of edge directionality as a local measure for the networks
#' Reference: Wolf,2019.Edge directionality (PHYSICAL REVIEW E 99, 012301 (2019))

# args <- commandArgs(trailingOnly=TRUE)
# 
# if (length(args)==0) {
#   stop("Usage: script.R <method.name [GC,CCM]> <layer.type [FOR,INV,MUX]> <z [1:5]>", call.=FALSE)
# }
# 
# method.name <- as.character(args[1])
# layer.type <- as.character(args[2])
# z <- as.numeric(args[3])

method.name <- "GC"#CCM, GC
layer.type <- "FOR"#FOR,MUX
z <- 3

library('igraph')

###########################################
####### Define parameters, directions #####
###########################################

cnt <- 50
layer.name <- "MON"
sys.path <- "D:/PostDoc/GWL/"
ana.dir <- paste0(sys.path, "Analysis/Measures/") 
dir.create(ana.dir, showWarnings = F)
g.dir <- paste0(sys.path, "Analysis/GRAPH/")
g.name <- paste0(g.dir, "GWL-Tabriz-WFFT_", method.name, "_Layer-", layer.name, "_undir-GRAPH_Surr", cnt, "_Z", z, "_", layer.type, ".RDS")
com.graph <- readRDS(g.name)
wei.vec <- E(com.graph)$weight
adj.mat <- igraph::as_adjacency_matrix(com.graph, type = c("both"))

###########################################
####### Import latitude and longitude #####
###########################################
gwl.shp.dir <- "D:/Thesis_Ph.D/GWL_Tabriz/Data/shape_files/"
shp.dir <- paste(gwl.shp.dir, "shps/", sep="")
shape.name <- paste(shp.dir, "Azarbayjan_basin.shp", sep="")
layer <- ogrListLayers(shape.name)
EA.map <- readOGR(shape.name, layer=layer, use_iconv = TRUE, encoding="UTF-8")

df.station <- readRDS("D:/PostDoc/GWL/Data/df.station.RDS")

df.station1 <- data.frame()

for (nc in 1:N) { #nodes
  utmx <- df.station$utmx[nc]
  utmy <- df.station$utmy[nc]
  df.station1 <- rbind(df.station1, data.frame(utmx = utmx, utmy = utmy))
}

sp.all <- SpatialPoints(coords=df.station1, proj4string = crs(EA.map))
sp.all_wgs84 <- spTransform(sp.all, CRS("+proj=longlat +datum=WGS84")) #Transform
df.station1 <- as.data.frame(sp.all_wgs84)

lat.v <- df.station1$utmy*pi/180
nlat <- length(lat.v)
lon.v <- df.station1$utmx*pi/180
nlon <- length(lon.v)
N <- nlat 
cen.mat <- vector() #Local anisotropy
mean.dir <- vector("list", N) #mean edge direction

#Function to compute norm of a vector
norm_vec <- function(x) sqrt(sum(x^2))

###########################################
####### corrected weighted local anisotropy
####### Mean edge direction ###############
###########################################

cnt.one <- 0
m <- 1
runAgain <- TRUE
while (runAgain)
{
  runAgain <- FALSE
  
  loc.one <- which(adj.mat[m,]==1)
  K <- 0
  r.m <- c(0,0)
  
  lapply(loc.one, function(n) {
    
      cnt.one <<- cnt.one + 1
     
      c.mn1 <- cos(lat.v[m])*sin(lat.v[n]) - cos((lon.v[m]-lon.v[n]))*cos(lat.v[n])*sin(lat.v[m])
      c.mn2 <- cos((lon.v[m]-lon.v[n]))*cos(lat.v[m])*cos(lat.v[n])+sin(lat.v[m])*sin(lat.v[n])
      c.mn3 <- sqrt(1-(c.mn2^2))
      c.mn <- c.mn1/c.mn3
      c.mn <- pmin(pmax(c.mn,-1.0), 1.0)
      
      if(is.infinite(c.mn) || is.na(c.mn)){
        c.mn <- 0
      }
      
      ###Correction for westward/eastward component
      if (lon.v[m] <= 0 && lon.v[m] <= lon.v[n] && lon.v[n] <= (lon.v[m] + pi)) {
        beta.mn <- acos(c.mn)
        
      } else if (lon.v[m] <= 0) {
        beta.mn <- 2 * pi - acos(c.mn)
        
      } else if (lon.v[m] > 0 && (lon.v[m] - pi) <= lon.v[n] && lon.v[n] <= lon.v[m]) {
        beta.mn <- 2 * pi - acos(c.mn)
        
      } else if (lon.v[m] > 0) {
        beta.mn <- acos(c.mn)
        
      } 
      
      #cat(paste("m", m, "-n", n, "-beta", beta.mn, "\n"))
      
      e.mn <- c(sin(beta.mn), cos(beta.mn))
      r.m <<- r.m + cos(lat.v[n])*wei.vec[cnt.one]*e.mn
      #cat(paste("m", m, "-n", n,"r.m", r.m, "\n"))
      K <<- K + cos(lat.v[n])*wei.vec[cnt.one]
  
  })
    
  mean.dir[[m]] <- r.m/K
  cen.mat[m] <- norm_vec(r.m)/K
  cat(paste("m", m, "\n"))
  
  if (m < N) {
    runAgain <- TRUE
    m <- m + 1
  }
  
}

nam.cen.mat <- paste0(ana.dir, "GWL-Tabriz-WFFT_", method.name,"_Layer-", layer.name, "_Local-anisotropy_Surr", cnt, "_Z", z, "_", layer.type, ".RDS")
saveRDS(cen.mat, nam.cen.mat)

nam.mean.dir <- paste0(ana.dir, "GWL-Tabriz-WFFT_", method.name,"_Layer-", layer.name, "_Mean-direction_Surr", cnt, "_Z", z, "_", layer.type, ".RDS")
saveRDS(mean.dir, nam.mean.dir)