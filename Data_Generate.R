
rm(list=ls(all=TRUE))
library('raster')
library('rgdal')
library('rnaturalearth')
library('rnaturalearthdata')
library('sf')
library('ggspatial')
library('rasterVis')
library('sp')

Sys.setlocale(locale = "persian")

gwl.shp.dir <- "D:/Thesis_Ph.D/GWL_Tabriz/Data/shape_files/"
shp.dir <- paste(gwl.shp.dir, "shps/", sep="")
shape.name <- paste(shp.dir, "Azarbayjan_basin.shp", sep="")
layer <- ogrListLayers(shape.name)
EA.map <- readOGR(shape.name, layer=layer, use_iconv = TRUE, encoding="UTF-8")

######################################################
#######Generate data of GWL from Tabriz aquifer#######
######################################################

subDir <- "D:/Thesis_Ph.D/GWL_Tabriz/Data/Selection/"
nam.data <- paste0(subDir, "GWL_LOC.RDS")
station.df <- readRDS(nam.data)
station.df.used <- data.frame(utmx=station.df$utmx,utmy=station.df$utmy)

df.station1 <- data.frame()
data.dir <- paste0("D:/Thesis_Ph.D/GWL_Tabriz/Data/Selection/I/raw/")
data.files <- list.files(path = data.dir)
num.data <- length(data.files)

for (nc in 1:num.data) {#MON and MUX
  loc <- which(station.df$NAME==data.files[nc])
  
  utmx <- station.df$utmx[loc]
  utmy <- station.df$utmy[loc]
  df.station1 <- rbind(df.station1, data.frame(utmx = utmx, utmy = utmy))
  
}

sp.all <- SpatialPoints(coords=df.station1, proj4string = crs(EA.map))
sp.all_wgs84 <- spTransform(sp.all, CRS("+proj=longlat +datum=WGS84")) #Transform
sp.all_df <- data.frame(sp.all_wgs84)
sp.all_df2 <- cbind(sp.all_df, ID=c(1:40))
sp.all_df <- sp.all_df[-c(3:9,40),]
sp.all_df2 <- cbind(sp.all_df, ID=c(1:32))

df.station <- data.frame()
for (nc in 1:num.data) {#MON and MUX
  loc <- which(station.df$NAME==data.files[nc])
  
  utmx <- station.df$utmx[loc]
  utmy <- station.df$utmy[loc]
  df.station1 <- data.frame(utmx = utmx, utmy = utmy)
  sp.all <- SpatialPoints(coords=df.station1, proj4string = crs(EA.map))
  sp.all_wgs84 <- spTransform(sp.all, CRS("+proj=longlat +datum=WGS84"))
  sp.all_df <- data.frame(sp.all_wgs84)
  df.station <- rbind(df.station, data.frame(name=station.df$NAME[loc],ID.name=station.df$ID[loc],utmx=utmx,m.utmx = sp.all_df$utmx, utmy=utmy,m.utmy = sp.all_df$utmy))
  
}
df.station <- df.station[-c(3:9,40),]
df.station <- cbind(df.station, ID=c(1:32))
saveRDS(df.station,"D:/PostDoc/GWL/Data/df.station.RDS")
