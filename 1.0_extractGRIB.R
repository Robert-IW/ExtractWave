# This script opens each monthly file containing three hourly estimates as a brick
#   and extracts the wave data along a 1000 km transect created by '0.0_createWaveTransects.R'

# Extracts significan height (hs), peak direction (dp), peak period (tp) and wind (u and v components)

# Data is collected and collated for each year to be consistent with the wind transects
# There are three choices of wave archive files:
#   1. climatology - uses NCEP Reanalysis hourly winds 1979 - 2009 e.g. 'multi_reanal.glo_30m_ext...'
#     NOTE the filename containing "_30m _ext" the 'ext' indicates phase 2 data which uses a more recent method
#   2. non-climatolgy - operational nww3 early version 1999-2009 e.g. 'nww3...'
#   3. non-climatolgoy - operational multi 2005 - 2017 e.g 'multi_1.glo_30m...'

# To convert from GRIB2 to netcdf use the bash script 'grib2ncdf' from the grib directory
#   -  this script produces an subset domain 0 - 50 WE and -50 - -15 SN
# To convert from GRIB1 to netcdf using ncl (NCAR language)
#   - these files need to be subset


library(ncdf4)
library(raster)
library(rasterVis)
library(reshape2)
library(data.table)

# ---------------------------------------------------------------------- DIRECTORIES
# directory containing functions and additional required files
func.Dir <- "~/R/myFunctions/"

# directory to store output files
out.Dir <- "/media/robert/KELP-HDD-Portable/WaveWatchIII/extractCSV/"

# directory to store output images
img.Dir <- NULL

# directory for transect files
trans.Dir <- "/media/robert/KELP-HDD-Portable/WaveWatchIII/Transects/"
# -----------------------------------------------------------------------------------

# set domain
maxLat <- -15
minLat <- -50
maxLon <- 50
minLon <- 0

# which version?
ver.pat <- "multi_reanal.glo_30m_ext*"

# mapping data
mapThemeBathy <- colorRampPalette(c("blue","white","brown"),alpha=F)
my.at <- seq(0,18,length.out = 10)

dom.extent <- extent(c(minLon,maxLon,minLat,maxLat))
crs <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

# extract data for each station
station.loc <- read.table(paste0(func.Dir,"station_locations.Rdata"),stringsAsFactors=F)
station.pts <- SpatialPointsDataFrame(station.loc[,1:2],station.loc,proj4string = crs)

# for each variable hs dp tp wind (NOTE wind has two variables!!)
var.list <- c("*hs*","*dp*","*tp*","*wind*")

# for each year collate the monthly files data and save in a single file
year.seq <- seq(1978, 2015, 1)

for (z in 1:4){
  #for (z in 1){
  var.pat <- var.list[z]
  
  for (y in 1:length(year.seq)){
    #for (y in 2){
    year.pat <- as.character(year.seq[y])
    
    # get a list of all the nc files in the directory containing a particular year
    filelist_un <- list.files(".", pattern = glob2rx(paste0(ver.pat,var.pat,year.pat,"*.nc")))
    
    # get a list of all the nc files in the directory containing a particular year and with 'subset'
    filelist_sub <- list.files(".", pattern = glob2rx(paste0(ver.pat,var.pat,year.pat,"*subset.nc")))
    
    # if the list is empty move on
    if (length(filelist_un)==0 && length(filelist_sub)==0){
      print(paste0("No data for ",year.pat))
      next
    }
    
    # if there are subset files use those preferentially 
    if (length(filelist_sub)!=0){
      filelist <- filelist_sub
      isSub <- T
    } else {
      filelist <- filelist_un
      isSub <- F
    }
    
    # ind_subset <- grep("subset",filelist_un, value = F)
    # filelist_sub <- filelist_un[ind_subset]
    # filelist_un <- filelist_un[-ind_subset]
    # rm(ind_subset)
    
    # create list to store the year data
    data.list <- list()
    # create list to store station month data
    stat.list <- vector("list",nrow(station.loc))
    names(stat.list) <- station.loc$station
    
    # for each file not subsetted
    for (i in 1:length(filelist)){
      
      print(paste0("Working on ",var.pat," monthly data for year ",year.pat))
      
      # if wind data there are two variables
      nc.stor <- nc_open(filelist[i])
      
      lon <- ncvar_get(nc.stor, "longitude")
      lat <- ncvar_get(nc.stor, "latitude")
      
      # time label for graph
      time <- ncvar_get(nc.stor, "time")
      tunits <- ncatt_get(nc.stor, "time", "units")
      torigin <- substr(tunits[2], 15, 100)
      time_store <- as.Date(as.POSIXct(time, origin = torigin, format="%Y-%m-%d %H:%M:%S", tz="UTC"))
      time_label <- as.Date(as.POSIXct(time, origin = torigin, format="%Y-%m-%d", tz="UTC"))
      time_hour <- strftime(as.POSIXct(time, origin = torigin), "%H%M", tz="UTC")
      date_file <- gsub("-", "", strftime(time_label, format = "%Y%m%d"))
      
      var.name <- as.list(names(nc.stor[['var']]))
      
      # graph label
      #main.title <- ncatt_get(nc.stor, var.name, "long_name")[[2]]
      main.title <- lapply(var.name, function(x) ncatt_get(nc.stor, x, "long_name")[[2]])
      
      # # check if file exists
      # filename1 <- paste(sourceURL,saveImage,sourceTitle,"_",date_file,"_",don,".png",sep = '')
      # filename2 <- paste(sourceURL,saveData,sourceTitle,"_",date_file,"_",don,".csv",sep = '')
      # 
      # if (file.exists(filename2)){
      #   cat("file exists\n")
      #   nc_close(nc.stor)
      #   next
      # }
      
      if (isSub){       # subsetting has already adjusted the image
        #var.brick <- brick(ncvar_get(nc.stor, var.name),xmn=min(lon),xmx=max(lon),ymn=min(lat),ymx=max(lat),crs=crs, transpose = T)
        #var.brick <- flip(var.brick, 2)
        var.brick <- lapply(var.name, function(x) brick(ncvar_get(nc.stor, x),
                                                        xmn=min(lon),xmx=max(lon),ymn=min(lat),ymx=max(lat),
                                                        crs=crs, transpose = T))
        var.brick <- lapply(var.brick, function(x) flip(x, 2))
        
      } else {
        #var.brick <- brick(ncvar_get(nc.stor, var.name),xmn=0,xmx=360,ymn=min(lat),ymx=max(lat),crs=crs, transpose = T)
        #var.brick <- flip(var.brick, 2)
        #var.brick <- rotate(var.brick)
        var.brick <- lapply(var.name, function(x) brick(ncvar_get(nc.stor, x),
                                                        xmn=0,xmx=360,ymn=min(lat),ymx=max(lat),
                                                        crs=crs, transpose = T))
        var.brick <- lapply(var.brick, function(x) flip(x, 2))
        var.brick <- lapply(var.brick, function(x) rotate(x))
      }
      
      nc_close(nc.stor)
      
      #hs <- brick(hs.temp,xmn=min(lon),xmx=max(lon),ymn=min(lat),ymx=max(lat),crs=crs)
      #var.brick <- brick(var.temp,xmn=0,xmx=360,ymn=min(lat),ymx=max(lat),crs=crs, transpose = T)
      
      #var.brick.loc <- crop(var.brick, dom.extent)
      var.brick.loc <- lapply(var.brick, function(x) crop(x, dom.extent))
      
      for (h in 1:nrow(station.loc)){
        
        stat.name <- station.loc$station[h]
        
        file.load <- paste0(trans.Dir,stat.name,"_WaveTransects.Rdata")
        load(file.load)
        rm(file.load)
        
        trans.offshore <- station.list$offshore
        rm(station.list)
        
        a <- station.pts@data[h,1:2]
        coordinates(a) <- ~lon+lat
        
        # create matrix
        #wave.trans <- extract(var.brick.loc[[1]], trans.offshore, cellnumbers=TRUE)
        wave.trans <- lapply(var.brick.loc, function(x) extract(x, trans.offshore, cellnumbers=TRUE))
        
        # get the coordinates of the extracted points  (locations are the same for each U and V wind)
        wave.trans.loc <- xyFromCell(var.brick.loc[[1]][[1]], wave.trans[[1]][[1]][,'cell'])   # a matrix
        
        # remove first column 'cell' and combine with lon lat columns
        wave.trans <- lapply(wave.trans, function(x) x[[1]][,-1])
        
        #wave.trans.data <- cbind(wave.trans.loc, wave.trans)
        #wave.trans.data <- cbind(wave.trans.loc, do.call("cbind", lapply(wave.trans, function(x) as.matrix(x))))
        wave.trans.data <- lapply(wave.trans, function(x) cbind(wave.trans.loc, x))
        
        #colnames(wave.trans.data) <- c("lon","lat",date_file)
        
        #data.wave <- as.data.frame(wave.trans.data, stringsAsFactors=F)
        data.wave <- lapply(wave.trans.data, function(x) as.data.frame(x, stringAsFactors=F))
        
        rm(wave.trans.data, wave.trans)
        
        # data.wave <- cbind(melt(data.wave,id=1:2,measure=3:ncol(data.wave))[,1:2],
        #                    melt(data.wave[3:ncol(data.wave)],id=NULL),
        #                    rep(time_hour, each=nrow(wave.trans.loc)))
        data.wave <- lapply(data.wave, function(x) cbind(melt(x,id=1:2,measure=3:ncol(x))[,1:2],
                                                         melt(x[3:ncol(x)],id=NULL),
                                                         rep(time_hour, each=nrow(wave.trans.loc))))
        
        if (var.pat != "*wind*"){
          data.wave <- data.wave[[1]]
          names(data.wave) <- c("lon","lat","date",var.name,"time")
          data.wave <- data.wave[c("date","time","lon","lat",var.name[[1]])]
          data.wave$date <- rep(date_file, each=nrow(wave.trans.loc))
        } else {
          data.wave1 <- data.wave[[1]]
          names(data.wave1) <- c("lon","lat","date",var.name[[1]],"time")
          data.wave1 <- data.wave1[c("date","time","lon","lat",var.name[[1]])]
          data.wave2 <- data.wave[[2]]
          names(data.wave2) <- c("lon","lat","date",var.name[[2]],"time")
          data.wave <- cbind(data.wave1, data.wave2[var.name[[2]]])
        }
        
        # convert to characters
        data.wave <- data.frame(lapply(data.wave, as.character), stringsAsFactors=FALSE)
        
        rm(wave.trans.loc)
        
        stat.list[[h]][[length(stat.list[[h]]) +1]] <- list(data.wave)
        rm(data.wave)
        
      } # for h in station
      gc()
      # 
      # vectorplot(var.brick.loc[[10]], region=F,
      #            isField=F, narrows=1000,unit="degrees",lwd.arrows=2, col="black",
      #            alpha.regions=.99,length=.1,
      #            colorkey=list(at=my.at), reverse=F,
      #            #scales=list(tick.number=4, tck=1),
      #            #aspX=.05,aspY=.05,
      #            #xlim=c(extent.tmp@xmin,extent.tmp@xmax),
      #            #ylim=c(extent.tmp@ymin,extent.tmp@ymax),
      #            main=list(label=main.title, cex=1.4),
      #            ylab=list("Latitude",cex=1.2), xlab=list("Longitude",cex=1.2))
      
      # extract data along transects
      
    } # for i in file list
    
    # save the station data as a single year
    print("Beginning writing data files")
    
    for (k in 1:nrow(station.loc)){
      stat.name <- names(stat.list[k])
      temp1 <- lapply(stat.list[[k]], function(x) do.call("rbind", x))
      temp2 <- do.call("rbind", temp1)
      
      # remove duplicates
      temp3 <- temp2[!duplicated(temp2),]
      rm(temp1,temp2)
      
      # replace '0' values with NA
      temp3[which(temp3[,5]==0),5] <- NA
      
      # save the file
      print(paste0("Saving ",var.name, "data for ",stat.name))
      filename1 <- paste0(out.Dir,stat.name,"_",var.name,"_",year.pat,".csv")
      fwrite(temp3, filename1, col.names = T)
    }
    
  } # for y in year
} # for z in variable
