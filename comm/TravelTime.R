#===========================================================
#    Objectives
#===========================================================
# Objective: given a target point, and a set of reference points, 
# compute travel time accounting for topography (elevation, rivers) and man-made lanscape (roads),
# prescribing diferent types of roads different travel speeds.

#===========================================================
#    To Do
#===========================================================
# Assign travel speeds based on the socio-economic class/assets distribution.

#===========================================================
#    Folders structure
#===========================================================
# Make sure that there are 2 folders on the same level:
# - comm (for scripts)
# - data (for data)

#===========================================================
#    Libraries
#===========================================================
rm(list=ls()) 

# set directory to current. works given the RStudio IDE is used
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# check directory
getwd()
require(raster)
#install.packages('gdistance')
require(gdistance)

#===========================================================
#    Read all relevant rasters
#===========================================================
elev_raster <- raster("../data/elev_raster.tif")
lc_raster <- raster("../data/lc_raster.tif")
motorways_raster <- raster("../data/motorways_raster.tif")
primroads_raster <- raster("../data/primroads_raster.tif")
tertroads_raster <- raster("../data/tertroads_raster.tif")
rivers_raster <- raster("../data/rivers_raster.tif")

#===========================================================
#    Take care of extents and resolutions, if needed
#===========================================================
# e <- extent(31.0, 31.2, -26.2, -26.0)
# rc <- crop(elev_raster, e)
# plot(rc)
# writeRaster(rc, "../data/elev_raster.tif")

 # Check whether rasters have same resolution. If not, resample as shown below.

 if (sum(res(elev_raster) == res(lc_raster))<2){
    # Resample Landcover raster to fit elevation raster
    lc_raster<- resample(lc_raster,elev_raster,method="ngb")
 }

identical( res(elev_raster), 
           res(lc_raster), 
           res(motorways_raster),  
           res(primroads_raster),
           res(tertroads_raster),
           res(rivers_raster))
 
 #===========================================================
 #    Compute
 #===========================================================

 motorways_raster[is.na(motorways_raster[])]<-0
 primroads_raster[is.na(primroads_raster[])]<-0
 tertroads_raster[is.na(tertroads_raster[])]<-0
 lc_raster[is.na(lc_raster[])]<-0
 rivers_raster[is.na(rivers_raster[])]<-0
 
 # Calculate transition across different directions from each pixel
 transitionFunction <- function(x){x[2] - x[1]}
 tr <- transition(elev_raster, transitionFunction, directions=8, symm=FALSE)
 slope <- geoCorrection(tr, scl=FALSE)
 
 # restrict estimates of speed to adjacent cells
 adj <- adjacent(elev_raster, cells=1:ncell(elev_raster), pairs=TRUE, directions=8)
 
 # Use Tobler's hiking function to define speed (km/h) as a function of slope
 speed <- slope
 speed[adj] <- 6*exp(-3.5 * abs(slope[adj] + 0.05))
 
 
 #===========================================================
 #    set travel speeds here
 #===========================================================
 WaterAdj <- adj[lc_raster[adj[,2]]==210,] # waterbodies
 TreeOther<-adj[lc_raster[adj[,2]]>=100 & lc_raster[adj[,2]]<=120,] # Other tree cover
 Herbaceous<-adj[lc_raster[adj[,2]]==140,] # Herbaceous
 SparseHerbaceous<-adj[lc_raster[adj[,2]]==150,] # Sparse herbaceous
 BareAreas<-adj[lc_raster[adj[,2]]==200,] # Sparse herbaceous
 speed[WaterAdj]<-0
 speed[TreeOther]<-speed[TreeOther]*0.4 # (i.e. 2 km/h on flat (2/5ths of normal 5 km/h))
 speed[Herbaceous]<-speed[Herbaceous]*0.6
 speed[SparseHerbaceous]<-speed[SparseHerbaceous]*0.8
 speed[BareAreas]<-speed[BareAreas]*0.4
 
 
 # Block movement acros major rivers (before adding roads)
 RiversAdj<-adj[lc_raster[adj[,2]]==1,]
 speed[RiversAdj]<-0
 
 # Add speed of roads. 
 # Restrict to adjacent 'road' cells. If going to a road cell then assign as road movement, coming from 
 # road to off road assign as walk
 TertiaryRoadAdj<-adj[tertroads_raster[adj[,2]]==1,]
 speed[TertiaryRoadAdj]<-10 # assumes 10km/h for roads
 
 PrimaryRoadAdj<-adj[primroads_raster[adj[,2]]==1,]
 speed[PrimaryRoadAdj]<-60 # assumes 60km/h for roads
 
 MotorwayAdj<-adj[motorways_raster[adj[,2]]==1,]
 speed[MotorwayAdj]<-80 # assumes 80km/h for roads
 
 # Now geocorrect speed to get conductance values (required for next step)
 x <- geoCorrection(speed, scl=FALSE)
 

 #===========================================================
 #    Compute
 #===========================================================
 set.seed(123)
 fromCoords <- data.frame(x=runif(1, 31.01, 31.199), y=runif(1, -26.19, -26.01))
 fromCoords
 
 toCoords <- data.frame(x=runif(10, 31.01, 31.199), y=runif(10, -26.19, -26.01))
 toCoords
 
 plot(extent(elev_raster))
 points(fromCoords, pch=16, col=2, cex=2)
 points(toCoords, pch=16, col=4)
 
 DistToPoints<-costDistance(x, fromCoords=as.matrix(fromCoords), toCoords=as.matrix(toCoords))
 dim(DistToPoints)
 
 png(filename = 'png1.png')
 plot(extent(elev_raster))
 points(fromCoords, pch=16, col=2, cex=2)
 points(toCoords, pch=16, col=4)
 points(toCoords[which(DistToPoints[1,]==min(DistToPoints[1,])),], pch=8, col=4, cex=3)
 dev.off()

 png(filename = paste('../
 /elev.png', sep=''))
 plot(elev_raster, main='Elevation')
 points(fromCoords, pch=16, col=2, cex=2)
 points(toCoords, pch=16, col=4)
 dev.off()
 
 png(filename = 'lc.png')
 plot(lc_raster, main='Land Cover')
 points(fromCoords, pch=16, col=2, cex=2)
 points(toCoords, pch=16, col=4)
 dev.off()

png(filename = 'prim_road.png')
plot(primroads_raster, main='Primary roads')
points(fromCoords, pch=16, col=2, cex=2)
points(toCoords, pch=16, col=4)
dev.off()

dim(toCoords)
r <- elev_raster
r[] <- DistToPoints
plot(r)
points(fromCoords, pch=16, col=2, cex=2)
