## Creates wave transects 100 km offshore from the station locations and 1.2 deg alongshore
## This has now been changed to 100 km (the wavewathc resolution is ~ 40 km)
##  offshore (line 207) from the station as the 1000 km csv files produced
##  are ~ 200 MB per station 11 July 2017

## The script is adapted from the wind transects and first creates an along shore 
##  transect from which a perpendicular line is estimated

# created by Robert Williamson 18 May 2017


library(raster)
library(rasterVis)
library(maptools)
library(sp)
library(rgeos)
library(polyclip)
library(parallel)


doGraphic <- T

# Local functions ------------------(--------------------------------------------------------------------
# check if any part of the coastline subset falls onto land
rm.coast <- function(coast.lines){
  test.line <- as.matrix(coast.lines@coords)
  test.height <- mean(extract(test.bathy, test.line), na.rm=T)
  if (nrow(coast.lines@coords) < 50){
    return(1)
  } else if (is.nan(test.height) | is.na(test.height) | test.height > 10) {
    return(1)
  } else {
    return(0)
  }
}

earth.bear <- function (long1, lat1, long2, lat2) 
{
  rad <- pi/180
  a1 <- lat1 * rad
  a2 <- long1 * rad
  b1 <- lat2 * rad
  b2 <- long2 * rad
  dlon <- b2 - a2
  bear <- atan2(sin(dlon) * cos(b1), cos(a1) * sin(b1) - sin(a1) * 
                  cos(b1) * cos(dlon))
  deg <- (bear%%(2 * pi)) * (180/pi)
  return(deg)
}

new.lon.lat <-
  function (lon, lat, bearing, distance) 
  {
    rad <- pi/180
    a1 <- lat * rad
    a2 <- lon * rad
    tc <- bearing * rad
    d <- distance/6378.145
    nlat <- asin(sin(a1) * cos(d) + cos(a1) * sin(d) * cos(tc))
    dlon <- atan2(sin(tc) * sin(d) * cos(a1), cos(d) - sin(a1) * 
                    sin(nlat))
    nlon <- ((a2 + dlon + pi)%%(2 * pi)) - pi
    npts <- cbind(nlon/rad, nlat/rad)
    return(npts)
  }

cl <- makeCluster(4)
clusterExport(cl, c("rm.coast", "extract"))
#-----------------------------------------------------------------------------------------------------
setwd("~/R/myFunctions/")

source("func_csv2pts.R")

# --------------------------- for mapping
crs <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
bathy <- raster("/home/robert/R/x86_64-pc-linux-gnu-library/3.3/bathy-GEBCO/gebco-southernAfrica_1min.nc")
mapThemeBathy <- colorRampPalette(c("blue","white","brown"),alpha=F)
at <- c(5000,2000,1000,500,200,100,50,0,-50,-100,-200,-500,-1000,-2000,-5000)

# use local R map (for coarse coastline)
gshhs.file <- system.file("share/gshhs_c.b", package="maptools")
xlims <- c(10,40)
ylims <- c(-36,-18)
world <- Rgshhs(gshhs.file, xlim=xlims, ylim=ylims, level=1, checkPolygons=TRUE)
border.coll <- as(world[[3]],"SpatialLines")

setwd("~/R_projects-CS/ame-temporalbabe/ExtractWaves/")

# ------------------------- Surface Wind
baseURL <- "/media/robert/KELP-HDD-Portable/"
product <- "CCMP"
# 
# listFiles <- list.files(paste(baseURL,type,product,"csv/",sep = "/"),full.names = TRUE)
# pts <- csv2pts(listFiles[1])
# r.list <- lapply(pts, function(f) {rasterFromXYZ(f)})
# r.stack <- stack(r.list)

# open existing raster data
# open.file1 <- paste0("~/R_projects-CS/ame-temporalbabe/ExtractWindFields/Output/",product,"_300kmMask.grd")
# bathy.mask <- raster(open.file1)
# rm(open.file1)

# get station locations and createn spatial points
station.loc <- read.table("~/R/myFunctions/station_locations.Rdata",stringsAsFactors=F)
station.pts <- SpatialPointsDataFrame(station.loc[,1:2],station.loc,proj4string = crs)

# -------------------------------------------------- Create Transects
# create empty list for all transects of spatialLines class with station names and coords
# transWind <- list()

for (i in 1:nrow(station.loc)){

  stat.name <- station.loc$station[i]
  cat(paste0("Starting transects for ",stat.name,'\n'))
  
  a <- station.pts@data[i,1:2]
  aa <- a
  
  # create a bounding box of 1.2 deg N/S and E/W to clip the buffer, This sets the limits of the 
  #   transect. Use longer length of coastline from SA border to create buffer
  bbox <- matrix(c(a[,1]-1.2,a[,2]-1.2,a[,1]+1.2,a[,2]+1.2),nrow=2,ncol=2)  # for 1000 km
  x <- extent(bbox)
  coast.loc <- crop(border.coll, x)
  #test.bathy <- crop(bathy, x)
  rm(x, bbox)
  
  # ---------------------------------- USE WITH HI RES COASTLINE
  # check if any of the hi res line objects in coast.loc extract heights > 10 m above sea level
  #   or which lines have length < 50
  #   these lines will be inland borders or fragments and will be removed
  # clusterExport(cl, "test.bathy")
  # test.coast <- parLapply(cl, slot(coast.loc@lines[[1]], "Lines"), function(r) rm.coast(r))
  # slot(coast.loc@lines[[1]], "Lines")[which(test.coast==1)] <- NULL
  # -------------------------------------------------------------
  
  bbox <- matrix(c(a[,1]-5,a[,2]-5,a[,1]+5,a[,2]+5),nrow=2,ncol=2)
  x <- extent(bbox)
  bathy.loc <- crop(bathy, x)
  rm(bbox, x)

  # Use a larger bathy for the plot
  #bbox <- matrix(c(a[,1]-10,a[,2]-10,a[,1]+10,a[,2]+10),nrow=2,ncol=2) # for 1000 km
  bbox <- matrix(c(a[,1]-2,a[,2]-2,a[,1]+2,a[,2]+2),nrow=2,ncol=2) # for 100 km
  x <- extent(bbox)
  bathy.plot <- crop(bathy, x)
  rm(bbox,x)

  # get the coordinates of the coastal line sections from the SpatialLines object 
  #   and create an offset polyline parallel to the coastline 0.25 deg / 25 km away for 1000 km
  
  res <- coordinates(coast.loc)
  res2 <- lapply(res[[1]], function(x) as.list(as.data.frame(x)))
  
  trans.along <- polylineoffset(res2, 0.3, jointype = "round", endtype = "openround")
  
  trans.along <- SpatialLines(list(Lines(Line(cbind(trans.along[[1]]$x,trans.along[[1]]$y)),
                                         ID="a")),proj4string = crs)
  trans.along <- gSimplify(trans.along,tol=0.1)
  bbox <- matrix(c(aa[,1]-1,aa[,2]-1,aa[,1]+1,aa[,2]+1),nrow=2,ncol=2)
  x <- extent(bbox)
  trans.along <- crop(trans.along, x)
 
  line.stor <- list()
  depth.stor <- list()
  for (j in 1:length(trans.along@lines[[1]]@Lines)){ 
    x <- slot(trans.along@lines[[1]],"Lines")[[j]]@coords
    line.stor[[j]] <- SpatialLines(list(Lines(Line(cbind(x[,1],x[,2])), ID="a")),proj4string = crs)
    depth.stor[[j]] <- mean(extract(bathy.loc, line.stor[[j]])[[1]],na.rm=T)
  }
  
  # for each line segment select those lines with a mean depth < 0
  idx <- which(depth.stor < 0)
  
  ll0 <- lapply(line.stor[idx], function(x) `@`(x, "lines"))
  ll1 <- lapply(unlist(ll0), function(y) `@`(y,"Lines"))
  trans.along <- SpatialLines(list(Lines(unlist(ll1), ID = 1)))
  
  # create a transect line perpendicular to average coastal orientation --------------------------------
  # if there is only a single 'Line' use the coords else if there a multiple lines combine and 
  # reorder the lines
  if (length(slot(trans.along@lines[[1]],"Lines")) == 1){
    line.coords <- slot(trans.along@lines[[1]],"Lines")[[1]]@coords
  } else {
    line.coords <- do.call(rbind, lapply(unlist(ll1), function(z) '@' (z, "coords")))
    line.coords <- line.coords[order(line.coords[,1]),]
    #trans.along <- SpatialLines(list(Lines(Line(ll2), ID = 1)))
  }
  rm(ll0, ll1)
  
  loc1 <- line.coords[1,]
  loc2 <- line.coords[nrow(line.coords),]
  
  # make sure most westerly lon is first, then transect always +90 degrees
  if (loc1[1] > loc2[1]){
    loc2 <- line.coords[1,]
    loc1 <- line.coords[nrow(line.coords),]
  }
  rm(line.coords)
  
  # get the bearing using 'earth.bear' function
  bear.coast <- earth.bear(loc1[1],loc1[2],loc2[1],loc2[2])
  bear.perp <- bear.coast+90
  if (bear.perp > 360){
    bear.perp <- bear.perp-360
  }
  # create a line perpendicular from the station 100 km long
  loc.new <- new.lon.lat(a$lon, a$lat, bear.perp, 100)
  lons <- rbind(a[1],loc.new[,1])
  lats <- rbind(a[2],loc.new[,2])
  rm(bear.perp, bear.coast)
  
  # coerce 'a','trans.offshore' to spatial points
  coordinates(a) <- ~lon+lat
  trans.offshore <- SpatialLines(list(Lines(Line(cbind(lons,lats)),
                                            ID="a")),proj4string = crs)

  if (doGraphic){
  filename1 <- paste0("Output/",stat.name,"_WaveTransects.png")
  png(filename = filename1,width = 4.5,height = 4.8,units = 'in',res = 300)
  
  p <- contourplot(bathy.plot, main=list(paste0(stat.name," Sample Transects"), x=.6, just="center"),
              region=T,col.regions=mapThemeBathy, alpha.regions=.3,
              at=at, margin = F, labels = F, ylab="Latitude", xlab="Longitude") +
    layer(sp.lines(trans.offshore,col="red",lwd=10)) +
    layer(sp.points(a, col="blue", pch=19, cex=1.5, fill=F, lwd=2)) +
    layer(sp.points(a, col="blue", pch=1, cex=3, fill=F, lwd=2)) +
    layer(sp.points(a, col="blue", pch=1, cex=6, fill=F, lwd=2))
    # vectorplot(along.off.plot, par.settings=mapTheme,region=F,scaleSlope=0.99,isField=TRUE,
    #            unit="degrees",colorkey=FALSE,alpha.regions=.3,lwd.arrows=4,
    #            col.arrow="blue") +
    # vectorplot(subset(wind.local.plot,4:5), par.settings=mapTheme,region=F,scaleSlope=0.99,isField=TRUE,
    #            unit="degrees",colorkey=FALSE,alpha.regions=.3,lwd.arrows=6,
    #            col.arrow="yellow")
  print(p)
  dev.off()
  }
  
  station.list <- list(trans.offshore)
  station.list <- setNames(station.list, "offshore")
  
  filename1 <- paste0("~/R_projects-CS/ame-temporalbabe/ExtractWaves/Output/",
                      stat.name,"_WaveTransects.Rdata")
  save(station.list, file = filename1)
  cat(paste0("Wave transects for ",stat.name," saved\n"))
  rm(station.list, trans.along, trans.offshore, region.wind)
}
stopCluster(cl)

