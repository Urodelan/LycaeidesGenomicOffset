##########################################
### Extract RegCM3 Climate Projections ###
##########################################

# Clear lists and graphics.
rm(list=ls())
graphics.off()

# Set working directory.
setwd("/Users/kenengoodwin/Desktop/USU/Graduate/Projects/Butterflies/Analyses/Data/Climate/RegCM3/Extract_Values")

# Specify directory with folders containing raster layers.
file_dir<-"/Users/kenengoodwin/Desktop/USU/Graduate/Projects/Butterflies/Analyses/Data/Climate/RegCM3/Files"

# Read in butterfly population locations.
loc<-read.csv(file="/Users/kenengoodwin/Desktop/USU/Graduate/Projects/Butterflies/Analyses/Data/Climate/Butterfly_Locations_v2.csv",stringsAsFactors=F)

# Load the sp, raster, ncdf4, and sf packages.
library(sp)
library(raster)
library(ncdf4)
library(sf)

# Create a spatial data frame of butterfly population locations.
loc_sp<-SpatialPointsDataFrame(coords=loc[,c("Longitude.E","Latitude.N")],
                               data=loc,
                               proj4string=CRS("+proj=longlat +datum=WGS84"))

# # Check spatial data frame.
# coordinates(loc_sp)
# proj4string(loc_sp)
# plot(loc_sp)

# Get vector of all raster layers to extract from.
files<-list.files(path=file_dir,pattern=".ncml.nc$",full.names=T,recursive=T)

# Create storage data frame for climate data.
climate<-data.frame(NULL)

# For both temperature and precipitation.
for(i in c("TA","RT")){
  
  # Loop through each raster layer.
  for(j in files){
    
    # Read in the NetCDF file as a raster stack.
    climate_raster<-brick(x=j,varname=i)
    
    # The raster package thinks the units for this raster is meters,
    # but it's actually kilometers. Let's fix this.
    ## Get the raster's original project 4 string.
    orig_proj4string<-proj4string(climate_raster)
    ## Change the untis from meters to kilometers.
    new_proj4string<-sub(pattern="\\+units=m",replacement="\\+units=km",x=orig_proj4string)
    ## Update the raster's project 4 string to indicate kilometers as the units.
    proj4string(climate_raster)<-new_proj4string
    
    # Transform spatial points data frame to same coordinate system as raster stack.
    loc_sp_reproj<-as_Spatial(st_transform(st_as_sf(loc_sp),proj4string(climate_raster)))
    
    # Extract raster vales at points.
    climate_sub<-as.data.frame(extract(x=climate_raster,y=loc_sp_reproj))
    
    # Add a population field to the extracted values data frame.
    climate_sub$Population<-loc$Locality
    
    # Reshape the extracted values data frame from wide to long format.
    climate_sub<-reshape(data=climate_sub,
                         varying=colnames(climate_sub)[-ncol(climate_sub)],
                         v.names="Value",
                         timevar="Time",
                         times=colnames(climate_sub)[-ncol(climate_sub)],
                         idvar="Population",
                         direction="long")
    
    # Remove NA values.
    climate_sub<-climate_sub[!is.na(climate_sub$Value),]
    
    # Get the domain name from the spatial layer.
    domain<-strsplit(x=basename(j),split="_")[[1]][3]
    
    # Add a field denoted the climate variable.
    climate_sub<-cbind(data.frame(Variable=i,Domain=domain),climate_sub)
    
    # Rename records.
    row.names(climate_sub)<-1:nrow(climate_sub)
    
    # Add information to the storage data frame.
    climate<-rbind(climate,climate_sub)
    
    # Print progress.
    print(paste0("Finished variable ",i,", file ",basename(j)))
    
  }
  
}

# Some observations are represented multiple times across the files.
# E.g., a download starting at the year 2000 will sometimes have observations from 1999.
## Create a field which is unique to each combination of variable,
## domain, population, and time.
climate$Combn<-paste0(climate$Variable,"-",climate$Domain,"-",climate$Population,"-",climate$Time)
## Remove duplicated observations.
climate<-climate[!duplicated(climate$Combn),]
## Remove combination field.
climate<-climate[,colnames(climate)!="Combn"]

# Format the time field.
## Remove the leading X.
climate$Time<-sub(pattern="^X",replacement="",x=climate$Time)
## Replace periods with hyphens.
climate$Time<-gsub(pattern="\\.",replacement="-",x=climate$Time)

# Add a field for state.
climate$State<-sapply(X=strsplit(x=climate$Population,split=", "),FUN=function(x) x[length(x)])

# Define priority domains for each state.
# Example: Point in Idaho.
## First check if point falls within the Northern Rocky Mountain raster domain.
## If not, then check if point falls within the Southern Rocky Mountain raster domain.
## If not, then check if point falls within the Pacific Northwest raster domain.
## If not, then check if point falls within the Pacific Southwest raster domain.
## Use the first domain the point falls within as the domain for the point.
## An NA in the domain priority data frame means that the previous domain
## priorities cover all of the state, and there is no need for additional
## domain priorities for the state.
#domain_prior<-data.frame(State=c("SD","MT","WY","ID","OR","CA","CO","NV","UT"),
#           Primary_Domain=c("ena","nrm","nrm","nrm","pnw","psw","srm","psw","srm"),
#           Secondary_Domain=c("nrm",NA,"srm","srm","psw",NA,NA,NA,NA),
#           Tertiary_Domain=c("srm",NA,NA,"pnw",NA,NA,NA,NA,NA),
#           Quaternary_Domain=c(NA,NA,NA,"psw",NA,NA,NA,NA,NA))
# South Dakota is the only state which falls in the Eastern North America domain,
# which has a shorter climate forecast than other domains (e.g., reaches 2080 but not 2090).
# Let's prioritize South Dakota's domain to be Northern Rocky Mountain (primary) and
# Southern Rocky Mountain (secondary). These two domains cover all of the state.
domain_prior<-data.frame(State=c("SD","MT","WY","ID","OR","CA","CO","NV","UT"),
                         Primary_Domain=c("nrm","nrm","nrm","nrm","pnw","psw","srm","psw","srm"),
                         Secondary_Domain=c("srm",NA,"srm","srm","psw",NA,NA,NA,NA),
                         Tertiary_Domain=c(NA,NA,NA,"pnw",NA,NA,NA,NA,NA),
                         Quaternary_Domain=c(NA,NA,NA,"psw",NA,NA,NA,NA,NA))

# Create an empty storage data frame.
pop_df<-data.frame(NULL)

# Loop through each population.
for(i in unique(climate$Population)){
  
  # Get the population of interest.
  pop<-climate[climate$Population==i,]
  
  # Get names of domains for which climate values are available.
  domain_avail<-unique(pop$Domain)
  
  # Get the domain priority for the state the population resides in.
  domain_order<-unlist(domain_prior[domain_prior$State==pop$State[1],-1])
  
  # Identify the highest priority domain for which there are values.
  domain_used<-domain_order[min(match(domain_avail,domain_order),na.rm=T)]
  
  # Subset climate values for the population to just the highest priority domain.
  pop<-pop[pop$Domain==domain_used,]
  
  # Add information to the storage data frame.
  pop_df<-rbind(pop_df,pop)
  
}

# Rename the storage data frame to climate.
climate<-pop_df

# Add a field for year.
climate$Year<-as.numeric(sapply(X=strsplit(x=climate$Time,split="-"),FUN=function(x) x[1]))

# What years are present in the dataset?
#sort(unique(climate$Year)) # 1968-1999 & 2011-2099.
# How many populations are in the dataset?
#nrow(loc) # 75
# How many populations have data for the given year?
## All have data for 1968-1999 & 2011-2099.
#length(unique(climate[climate$Year==2090,]$Population))

# Create an empty storage data frame for 10-year climate normals.
normals_df<-data.frame(NULL)

# Loop through each decade.
for(i in seq(from=2020,to=2090,by=10)){
  # Subset climate data to monthly values from the decade.
  sub<-climate[climate$Year %in% (i-9):i,]
  # Get mean monthly values for the decade.
  sub<-aggregate(Value~Population+Variable,data=sub,FUN=mean)
  # Add a field denoting the decade.
  sub<-cbind(Decade=paste0((i-9),"-",i),sub)
  # Add 10-year climate normals to the storage data frame.
  normals_df<-rbind(normals_df,sub)
}

# Convert total precipitation in mm/day to mm/year.
normals_df$Value[normals_df$Variable=="RT"]<-normals_df$Value[normals_df$Variable=="RT"]*365

# Give variables more meaningful names.
normals_df$Variable[normals_df$Variable=="TA"]<-"Temp.C"
normals_df$Variable[normals_df$Variable=="RT"]<-"Precip.mm"

# Create empty storage data frame for differenced observations.
norm_diffs<-data.frame(NULL)

# Loop through each population.
for(i in unique(normals_df$Population)){
  # Loop through each climate variable.
  for(j in unique(normals_df$Variable)){
    # Get the normal data for the population and climate variable.
    sub<-normals_df[normals_df$Population==i & normals_df$Variable==j,]
    # Get the initial decade value.
    init_value<-sub$Value[sub$Decade=="2011-2020"]
    # Calculate the change in value between the initial decade and all later decades.
    sub$Change<-sub$Value-init_value
    # Store the differences in the storage data frame.
    norm_diffs<-rbind(norm_diffs,sub)
  }
}

# Specify that the values in the normal differences data frame
# are from the RegCM3 climate projections.
## Values.
colnames(norm_diffs)[colnames(norm_diffs)=="Value"]<-"Value.RegCM3"
## Differences.
colnames(norm_diffs)[colnames(norm_diffs)=="Change"]<-"Change.RegCM3"

# Change the field name "Decade" to "Year".
colnames(norm_diffs)[colnames(norm_diffs)=="Decade"]<-"Year"

# Denote the decade by its final year.
norm_diffs$Year<-sapply(X=strsplit(x=norm_diffs$Year,split="-"),FUN=function(x) x[2])

# Read in butterfly PRISM data.
butterflies<-read.csv(file="/Users/kenengoodwin/Desktop/USU/Graduate/Projects/Butterflies/Analyses/Data/Climate/Butterfly_Climate_v2.csv")

# Rename the locality field to population on the butterfly PRISM data.
colnames(butterflies)[colnames(butterflies)=="Locality"]<-"Population"

# Subset butterfly PRISM data to just the climate variables of interest.
butterflies<-butterflies[,c("Population","Temp.C","Precip.mm")]

# Create a storage data frame for climate projections.
projections<-data.frame(NULL)

# Loop through each population.
for(i in unique(norm_diffs$Population)){
  # Loop through each climate variable.
  for(j in unique(norm_diffs$Variable)){
    # Get the data for the population and climate variable.
    sub<-norm_diffs[norm_diffs$Population==i & norm_diffs$Variable==j,]
    # Get the PRISM normal for the population.
    prism_value<-butterflies[butterflies$Population==i,j]
    # Add the PRISM normal value to the subsetted data frame.
    sub$PRISM.normal<-prism_value
    # Calculate the projected value of the climate variable as the PRISM normal
    # plus the difference observed in the RegCM3 climate projections.
    sub$Value.projected<-sub$PRISM.normal+sub$Change.RegCM3
    # Store information in the storage data frame.
    projections<-rbind(projections,sub)
  }
}

# Re-order fields of the climate projections data frame.
projections<-projections[,c("Population","Variable","Year","Value.RegCM3","Change.RegCM3","PRISM.normal","Value.projected")]

# Write out the butterfly climate projections.
write.csv(x=projections,file="Butterfly_Projections_v2.csv",row.names=F)

# Done!
print("Done!")
