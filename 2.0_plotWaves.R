# Plot a Hovmoller for each station consisting of the transect from the station offshore to 100 km
# Each figure is the transect on a map and the Hovmoller for a particular year and variable


library(raster)
library(rasterVis)
library(data.table)
library(maptools)
library(grid)
library(gridBase)
library(gridExtra)
library(ggplot2)
library(magrittr)
library(multipanelfigure)
library(mosaic)

# ---------------------------------------------------------------------- DIRECTORIES
# directory containing functions and additional required files
func.Dir <- "~/R/myFunctions/"

# directory to store output files
out.Dir <- "/media/robert/KELP-HDD-Portable/WaveWatchIII/extractCSV/"

# directory for transect files
trans.Dir <- "/media/robert/KELP-HDD-Portable/WaveWatchIII/Transects/"

# directory for images
img.Dir <- "/media/robert/KELP-HDD-Portable/WaveWatchIII/Images/"
# ----------------------------------------------------------------------------------
# get SA borders
crs <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
coastline <- readShapePoly(paste0(func.Dir,"SAfrica-coast.shp"), proj4string = crs)
coastline <- fortify(coastline)

mapTheme <- colorRampPalette(c("blue","lightgoldenrodyellow","brown"),alpha=F)
# minLat <- -40;minLon <- 10;maxLat <- -20;maxLon <- 40
# locbb <- matrix(c(minLon,minLat,maxLon,maxLat),nrow=2,ncol=2)
# x <- extent(locbb)
# border.sa <- crop(borders,x)

# extract data for each station
station.loc <- read.table(paste0(func.Dir,"station_locations.Rdata"),stringsAsFactors=F)
station.pts <- SpatialPointsDataFrame(station.loc[,1:2],station.loc,proj4string = crs)

# list of the wave variables
var.list <- c("HTSGW","DIRPW","UGRD","PERPW")

# for each year collate the monthly files data and save in a single file
year.seq <- seq(1978, 2015, 1)

# for each station
#for (i in 1:nrow(station.loc)){
for (i in 1){
  
  stat.name <- station.loc$station[i]
  a <- station.pts@data[i,1:2]
  #coordinates(a) <- ~lon+lat
  
  # load the transect for plotting
  file.load <- paste0(trans.Dir,stat.name,"_WaveTransects.Rdata")
  load(file.load)
  trans.offshore <- station.list$offshore
  trans <- SpatialLinesDataFrame(trans.offshore, data.frame(ID=1),match.ID=F)
  rm(file.load, trans.offshore)
  
  # create the inset map
  map.inset <- ggplot() + coord_equal() + theme_bw() + 
    geom_polygon(aes(x=long, y=lat, group=group), colour="black", fill="grey80", data=coastline) +
    labs(x = "Longitude", y = "Latitude") +
    coord_map(xlim = c(8,40), ylim = c(-40,-25), projection = "mercator") +
    geom_path(data=trans, aes(x=long,y=lat),colour="red",lwd=1.5) + 
    geom_point(data=a, aes(x=lon,y=lat), colour="blue", size=4)
  
  # for each wave variable
  #for (h in 1:length(var.list)){
  for (h in 1){
    
    var.name <- var.list[h]
    
    #for (j in 1:length(year.seq)){
    for (j in 2){
      
      year.name <- year.seq[j]
      
      # file to load
      file.name1 <- paste0(out.Dir,stat.name,"_",var.name,"_surface_",year.name,".csv")
      data.wave <- as.data.frame(fread(file.name1), stringAsFactors=F)
      
      # file to save
      file.name2 <- paste0(img.Dir,stat.name,"_",var.name,"_surface_",year.name,".png")
      
      # extract details on the dates
      data.var <- colnames(data.wave)[ncol(data.wave)]
      data.year <- substr(data.wave["date"][1,1],1,4)
      data.date <- as.Date(as.character(unique(data.wave$date)),"%Y%m%d")
      data.day <- format(data.date,"%d-%b")
      
      data.hour <- unique(data.wave$time)
      data.hour <- sprintf("%04d", data.hour)
      data.hour <- format(strptime(data.hour, format="%H%M"), format = "%H:%M")
      
      date.names <- expand.grid(data.hour, data.day, stringsAsFactors = F)
      date.names <- date.names[c("Var2","Var1")]
      date.names$labels <- with(date.names, paste(Var2,Var1))
      
      # find the first day of the month for the tick marks
      data.day1 <- data.day[grep("01",data.day)]
      data.day2 <- split(data.day1, seq(length(data.day1)))
      
      date.ticks <- lapply(data.day2, function(x) which(date.names[,1] == x)[1])
      date.ticks <- unlist(date.ticks)
      
      #obs.names <- unique(data.wave[c("date","time")])
      #obs.data <- split(obs.data, seq(nrow(obs.data)))
      obs.data <- unique(data.wave["date"])
      obs.data <- split(obs.data, seq(nrow(obs.data)))
      
      # get the set of lon lat for the transect pixels
      obs.loc <- unique(data.wave[c("lon","lat")])
      
      # plot the Hovmoller ----------------------------------------------------------------------------------
      # get the daily transects (8 x 3-hourly records perday)
      data.sub <- lapply(obs.data, function(x) matrix(data.wave[data.var][which(data.wave["date"]==x[1,1]),],
                                                      ncol = length(data.hour)))
      names(data.sub) <- data.day # may have to check the last entry for consistent matrix dims
      
      #================================================================================
      
      data.sub <- data.sub[-length(data.sub)] # may only apply to test run
      
      # ===============================================================================
      
      #sp <- SpatialPointsDataFrame(obs.loc, as.data.frame(data.sub[[1]]))
      temp <- as.data.frame(data.sub)
      sp <- SpatialPointsDataFrame(obs.loc, temp)
      #sp <- lapply(data.sub, function(x) SpatialPointsDataFrame(obs.loc, as.data.frame(x)))
      
      gridded(sp) <- TRUE
      s <- brick(sp)
      # (an alternative route would be rasterFromXYZ)
      
      #idx <- seq(as.Date('1970-01-01'), as.Date('2003-03-01'), by='month')
      #s <- setZ(s, idx)
      s <- setZ(s,date.names$labels[1:length(temp)])
      y.scale <- list(at=date.ticks,labels=data.day1)
      x.scale <- list(at=seq(s@extent@xmin,s@extent@xmax),1)
      
      main.title <- paste0(stat.name," ",year.name,"\n",data.var)
      mytext <- textGrob(main.title, just=c("right","top"), gp=gpar(cex=1.2,font=2), hjust=.65, vjust=1)
      
      myplot <- hovmoller(s, dirXY=x, xlab="Longitude", ylab="Time (3-hour intervals)",
                          scales=list(y=y.scale),
                          colorkey=list(at=seq(0,6,.2)), col.regions=mapTheme,
                          page = function(page) grid.text("Sign. Wave Height (m)",
                                                          x=unit(0.9,"npc") ,y=unit(0.55,"npc"),rot=90,
                                                          gp=gpar(cex=.8)),
                          par.settings=list(layout.heights=list(top.padding=0),
                                            layout.widths=list(right.padding=8)))
      
      # using multipanelfigure
      fig1 <- multi_panel_figure(width=c(50,60), height=c(30,120), panel_clip="off",
                                 row_spacing=0, column_spacing=0, panel_label_type="none")
      fig1 %<>% fill_panel(map.inset, row=1, col=1)
      fig1 %<>% fill_panel(mytext, row=1, column=c(2))
      fig1 %<>% fill_panel(myplot, row=2, column=c(1,2))
      
      save_multi_panel_figure(fig1,file.name2, dpi = 300)
    }
  }
}
