########################################
### Get Past PRISM Climate Estimates ###
########################################

# Clear lists and graphics.
rm(list=ls())
graphics.off()

# Set working directory.
setwd("/Users/kenengoodwin/Desktop/USU/Graduate/Projects/Butterflies/Analyses/Data/Climate/Past")

# Read in butterfly population locations.
loc<-read.csv(file="../Butterfly_Locations_v2.csv",stringsAsFactors=FALSE)

# Load the sp package.
library(sp)

# Create a spatial data frame of butterfly population locations.
loc_sp<-SpatialPointsDataFrame(coords=loc[,c("Longitude.E","Latitude.N")],
                               data=loc,
                               proj4string=CRS("+proj=longlat +datum=WGS84"))

# # Check spatial data frame.
# coordinates(loc_sp)
# proj4string(loc_sp)
# plot(loc_sp)

# Get PRISM climate data.
## Load the prism package.
library(prism)
## Set PRISM path.
prism_set_dl_dir(path="/Users/kenengoodwin/Desktop/USU/Graduate/Projects/Butterflies/Analyses/Data/Climate/Past/Rasters",create=FALSE)
## Uncomment the lines below to download PRISM data to the PRISM download directory.
## PRISM climate normals are for the period 1991-2020. 800 m grid cell resolution is
## not available for monthly or daily estimates. Monthly and daily estimates
## are only available in 4km grid cell resolution.
##get_prism_normals(type="tmean",resolution="800m",annual=T,keepZip=F)
##get_prism_normals(type="ppt",resolution="800m",annual=T,keepZip=F)
# get_prism_monthlys(type="tmean",years=1981:2020,mon=1:12,keepZip=FALSE)
# get_prism_monthlys(type="ppt",years=1981:2020,mon=1:12,keepZip=FALSE)

### Resume from here.

# Read in PRISM data as a raster stack.
library(raster)
climate_raster<-stack(pd_to_file(prism_archive_ls()))

# Check raster stack coordinate system.
proj4string(climate_raster)

# Transform spatial points data frame to same coordinate system as raster stack.
library(sf)
loc_sp<-as_Spatial(st_transform(st_as_sf(loc_sp),proj4string(climate_raster)))

# Extract raster vales at points.
climate<-as.data.frame(extract(x=climate_raster,y=loc_sp))

# Rename climate variables.
## Get relevant portions of the climate information field names.
colnames(climate)<-sapply(X=strsplit(x=colnames(climate),split="_"),FUN=function(x) paste0(x[2],"_",x[5]))
## Temperature (C).
#colnames(climate)[colnames(climate)=="PRISM_tmean_30yr_normal_800mM4_annual_bil"]<-"Temp.C"
colnames(climate)<-sub(pattern="^tmean_",replacement="Temp.C_",x=colnames(climate))
## Precipitation (mm).
#colnames(climate)[colnames(climate)=="PRISM_ppt_30yr_normal_800mM4_annual_bil"]<-"Precip.mm"
colnames(climate)<-sub(pattern="^ppt_",replacement="Precip.mm_",x=colnames(climate))

# Add climate variables to butterfly data frame.
climate<-cbind(loc[,"Locality",drop=FALSE],climate)

# Reshape climate information from wide to long format.
climate<-stats::reshape(data=climate,direction="long",
                        idvar="Locality",timevar="Time",
                        varying=colnames(climate)[2:ncol(climate)],
                        sep="_")

# Create a year field.
climate$Year<-as.numeric(substr(x=climate$Time,start=1,stop=4))

# Create a month field.
climate$Month<-as.numeric(substr(x=climate$Time,start=5,stop=6))

# Get fields of interest.
climate<-climate[,c("Locality","Year","Month","Temp.C","Precip.mm")]

# Re-order climate records.
climate<-climate[order(climate$Locality,climate$Year,climate$Month),]

# Rename locality field to population.
colnames(climate)[colnames(climate)=="Locality"]<-"Population"

# Reset row names.
row.names(climate)<-1:nrow(climate)

# Convert climate data to fully long format.
climate<-stats::reshape(data=climate,direction="long",
                        idvar=c("Population","Year","Month"),
                        varying=c("Temp.C","Precip.mm"),
                        v.names="Value",
                        timevar="Variable",times=c("Temp.C","Precip.mm"))

# Re-order fields.
climate<-climate[,c("Population","Variable","Year","Month","Value")]

# Re-order climate records.
climate<-climate[order(climate$Population,climate$Variable,climate$Year,climate$Month),]

# Reset row names.
row.names(climate)<-1:nrow(climate)

# Write out the butterfly past climate estimates.
write.csv(x=climate,file="Butterfly_PastMonths_v1.csv",row.names=FALSE)

# Done!
print("Done!")
