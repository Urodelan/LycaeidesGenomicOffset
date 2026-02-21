########################################################
### Process Butterfly Past Monthly Climate Estimates ###
########################################################

# Clear R environment.
rm(list=ls())

# Clear graphics.
graphics.off()

# Set working directory.
setwd("/Users/kenengoodwin/Desktop/USU/Graduate/Projects/Butterflies/Analyses/Data/Climate/Past")

# Read in past climate monthlys.
climate<-read.csv(file="/Users/kenengoodwin/Desktop/USU/Graduate/Projects/Butterflies/Analyses/Data/Climate/Past/Butterfly_PastMonths_v1.csv")

# Create an empty storage data frame for 10-year climate normals.
normals_df<-data.frame(NULL)

# Loop through each decade.
for(i in seq(from=1990,to=2020,by=10)){
  # Subset climate data to monthly values from the decade.
  sub<-climate[climate$Year %in% (i-9):i,]
  # Get mean monthly values for the decade.
  sub<-aggregate(Value~Population+Variable,data=sub,FUN=mean)
  # Add a field denoting the decade.
  sub<-cbind(Decade=paste0((i-9),"-",i),sub)
  # Add 10-year climate normals to the storage data frame.
  normals_df<-rbind(normals_df,sub)
}

# Convert total precipitation in mm/month to mm/year.
normals_df$Value[normals_df$Variable=="Precip.mm"]<-normals_df$Value[normals_df$Variable=="Precip.mm"]*12

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

# Denote the decade by its final year.
norm_diffs$Year<-sapply(X=strsplit(x=norm_diffs$Decade,split="-"),FUN=function(x) x[2])

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
    sub$Value.projected<-sub$PRISM.normal+sub$Change
    # Store information in the storage data frame.
    projections<-rbind(projections,sub)
  }
}

# Re-order fields of the climate projections data frame.
projections<-projections[,c("Population","Variable","Year","Value","Change","PRISM.normal","Value.projected")]

# Write out the past butterfly climate projections.
write.csv(x=projections,file="Butterfly_Past_v1.csv",row.names=FALSE)

# Done!
print("Done!")
