#################################
### Get PRISM Climate Normals ###
#################################

# Clear lists and graphics.
rm(list=ls())
graphics.off()

# Set working directory.
setwd("/Users/kenengoodwin/Desktop/USU/Graduate/Projects/Butterflies/Analyses/Data/Climate")

# Read in butterfly population locations.
loc<-read.csv(file="Butterfly_Locations_v2.csv",stringsAsFactors=F)

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
prism_set_dl_dir(path="/Users/kenengoodwin/Desktop/USU/Graduate/Projects/Butterflies/Analyses/Data/Climate/PRISM",create=F)
## Uncomment the lines below to download PRISM data to the PRISM download directory.
## PRISM climate normals are for the period 1991-2020. Use 800 m grid cell resolution.
#get_prism_normals(type="tmean",resolution="800m",annual=T,keepZip=F)
#get_prism_normals(type="ppt",resolution="800m",annual=T,keepZip=F)

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
## Temperature (C).
colnames(climate)[colnames(climate)=="PRISM_tmean_30yr_normal_800mM4_annual_bil"]<-"Temp.C"
## Precipitation (mm).
colnames(climate)[colnames(climate)=="PRISM_ppt_30yr_normal_800mM4_annual_bil"]<-"Precip.mm"

# Add climate variables to butterfly data frame.
loc<-cbind(loc,climate)

# Load ggplot2 and stringr.
library(ggplot2)
library(stringr)

# Create a scatterplot of temperature and precipitation data.
climate_plot<-ggplot()+
  geom_text(aes(x=Temp.C,y=Precip.mm,label=stringr::str_wrap(string=Locality,width=15),
                angle=33.75,color=Taxon),data=loc,size=1)+
  theme_light()+
  theme(plot.title=element_text(hjust=0.5,face="bold"),
        legend.title.align=0.5)+
  ggtitle("30-Year Climate Normals (1991-2020)")+
  xlab("Mean Temperature (\u00B0C)")+
  ylab("Cumulative Annual Precipitation (mm)")

# Display climate plot.
print(climate_plot)

# Save climate plot.
ggsave(filename="Climate_Plot_v2.pdf",plot=climate_plot,width=7,height=5,units="in")

# Write out the butterfly information with climate data.
write.csv(x=loc,file="Butterfly_Climate_v2.csv",row.names=F)

# Done!
print("Done!")
