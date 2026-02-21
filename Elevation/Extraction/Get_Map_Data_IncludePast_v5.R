####################
### Get Map Data ###
####################

# Includes past projections.

# Clear R enivornment.
rm(list=ls())

# Clear graphics.
graphics.off()

# Set working directory.
setwd("/Users/kenengoodwin/Desktop/USU/Graduate/Projects/Butterflies/Analyses/Data/Climate/GO_Map/Include_Past")

# Read in butterfly data.
df<-read.csv(file="/Users/kenengoodwin/Desktop/USU/Graduate/Projects/Butterflies/Analyses/Data/Product/butterfly_climate_v3.csv")

# Get just desired fields from the butterfly data.
df<-df[,c("Population","PopID","Longitude.E","Latitude.N")]

# Read in GO projections.
go<-read.csv(file="/Users/kenengoodwin/Desktop/USU/Graduate/Projects/Butterflies/Analyses/Models/Product/Output/Process/IncludePast/Genomic_Offset_Projections_IncludesPast_v1.csv")

# Get just desired projection fields.
go<-go[,c("Population","Year","Genomic_Offset")]

# Get the original GO values in 2020.
go_orig<-go[go$Year==2020,]

# Remove the year field from the original GO values.
go_orig<-go_orig[,-which(colnames(go_orig)=="Year")]

# Rename the GO field in the original GO values to GO_Initial.
colnames(go_orig)[colnames(go_orig)=="Genomic_Offset"]<-"GO_Initial"

# Add the original GO values to the GO data.
go<-merge(x=go,y=go_orig)

# Calculate the change in GO values.
go$Delta_GO<-go$Genomic_Offset-go$GO_Initial

# Calcalate the change in GO values as a percentage.
go$Delta_GO_pct<-go$Delta_GO/go$GO_Initial

# Add genomic offset projections to the butterfly climate data frame.
df<-merge(x=df,y=go)

# Write out the butterfly climate data with genomic offset projections.
write.csv(x=df,file="butterfly_map_IncludesPast_v1.csv",row.names=FALSE)

# Load the sp package.
library(sp)

# Create a spatial data frame of butterfly population locations.
df_sp<-SpatialPointsDataFrame(coords=df[,c("Longitude.E","Latitude.N")],
                              data=df,
                              proj4string=CRS("+proj=longlat +datum=WGS84"))

# Load ggplot.
library(ggplot2)

# Load package for saving unicode characters in pdfs.
library(Cairo)

# Load sf.
library(sf)

# Load ggpspatial.
library(ggspatial)

# Load the scales package.
library(scales)

# Load the raster package.
library(raster)

# Load the ggnewscale package.
library(ggnewscale)

# Load cowplot.
library(cowplot)

# Convert the butterfly data from sp to sf.
df_sf<-st_as_sf(df_sp)

# Transform the butterfly data to the pseudo-Mercator projection.
df_sf<-st_transform(df_sf,crs=3857)

# Read in US boundary.
nation<-read_sf(dsn="/Users/kenengoodwin/Desktop/USU/Graduate/Projects/Butterflies/Analyses/Data/Climate/State_Boundaries/cb_2021_us_all_20m/cb_2021_us_nation_20m/cb_2021_us_nation_20m.shp")

# Read in US state boundaries.
states<-read_sf(dsn="/Users/kenengoodwin/Desktop/USU/Graduate/Projects/Butterflies/Analyses/Data/Climate/State_Boundaries/cb_2021_us_all_20m/cb_2021_us_state_20m/cb_2021_us_state_20m.shp")

# Read in elevation data.
dem<-raster(x="/Users/kenengoodwin/Desktop/USU/Graduate/Projects/Butterflies/Analyses/Data/Climate/Elevation/Elevation_PseudoMercator_Clipped.tif")

# Throw an error if the coordinate system between the elevation raster
# and the butterflylocations are not the same.
if(!identical(proj4string(dem),proj4string(as_Spatial(df_sf)))) stop("The coordinate system between the elevation raster and the butterfly locations are not the same!")

# Define the cube root sign function.
cube_root_sign<-function(x) sign(x) * abs(x)^(1/3)

# Get the final genomic offset values.
go_elev<-df_sf[df_sf$Year==2090,]

# Extract elevation at butterfly locations.
go_elev$Elevation.m<-raster::extract(x=dem,y=go_elev)

# Transform the change in genomic offset values.
go_elev$Delta_GO_trans<-cube_root_sign(go_elev$Delta_GO)

# Get elevation plus genomic offset data as a data frame.
go_elev<-as.data.frame(go_elev)

# Get desired fields from the elevation plus genomic offset data.
go_elev<-go_elev[,c("Population","PopID","Elevation.m","GO_Initial","Delta_GO_trans")]

# Fit a linear model between elevation and the 2020 genomic offset values.
lm.initial<-lm(log(GO_Initial)~Elevation.m,data=go_elev)

# Get 5% of the elevation range.
buffer<-diff(range(go_elev$Elevation.m))*0.05

# Create elevation data to predict on.
pred_data<-data.frame(Elevation.m=seq(from=min(go_elev$Elevation.m)-buffer,
                                      to=max(go_elev$Elevation.m)+buffer,
                                      by=1))

# Predict the initial genomic offset by elevation.
pred.initial<-cbind(pred_data,
                    as.data.frame(
                      predict(object=lm.initial,newdata=pred_data,interval="confidence")))

# Display summary of the linar model between elevation and the 2020 genomic offset values.
#summary(lm.initial)
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  3.593e+00  1.587e-01  22.647  < 2e-16 ***
# Elevation.m -2.778e-04  7.931e-05  -3.503  0.00115 ** 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2979 on 40 degrees of freedom
# Multiple R-squared:  0.2348,	Adjusted R-squared:  0.2156 
# F-statistic: 12.27 on 1 and 40 DF,  p-value: 0.001147

# Fit a linear model between elevation and the transformed difference
# in 2090-2020 genomic offset values.
lm.delta<-lm(Delta_GO_trans~Elevation.m,data=go_elev)

# Predict the change in genomic offset by elevation.
pred.delta<-cbind(pred_data,
                  as.data.frame(
                    predict(object=lm.delta,newdata=pred_data,interval="confidence")))

# Display summary of the linear model between elevation and
# the transformed difference in 2090-2020 genomic offset values.
#summary(lm.delta)
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)
# (Intercept) 6.924e-01  6.402e-01   1.082    0.286
# Elevation.m 8.076e-05  3.200e-04   0.252    0.802
# 
# Residual standard error: 1.202 on 40 degrees of freedom
# Multiple R-squared:  0.00159,	Adjusted R-squared:  -0.02337 
# F-statistic: 0.06369 on 1 and 40 DF,  p-value: 0.802

# Create first GO vs elevation plot.
elev_plot1<-ggplot(data=go_elev)+
  geom_point(aes(x=Elevation.m,y=log(GO_Initial)))+
  theme_light()+
  ylab("ln(2020 Genomic Offset)")+
  xlab(NULL)+
  theme_light()+
  ggtitle("Genomic Offset by Elevation",subtitle="Point: Observation \u2013 Solid line: Predicted value \u2013 Dashed line: 95% confidence interval")+
  geom_line(data=pred.initial,aes(x=Elevation.m,y=fit),color="black",linetype="solid")+
  geom_line(data=pred.initial,aes(x=Elevation.m,y=lwr),color="black",linetype="dashed")+
  geom_line(data=pred.initial,aes(x=Elevation.m,y=upr),color="black",linetype="dashed")+
  theme(plot.margin=margin(t=5.5,r=5.5,b=1,l=8.8,unit="pt"),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        plot.title=element_text(face="bold",hjust=0.5),
        plot.subtitle=element_text(hjust=0.5))+
  scale_x_continuous(expand=expansion(mult=c(0,0)))

# Create second GO vs elevation plot.
elev_plot2<-ggplot(data=go_elev)+
  geom_point(aes(x=Elevation.m,y=Delta_GO_trans))+
  theme_light()+
  ylab(expression("(\u0394 2090-2020 Genomic Offset)"^~frac(1,3)))+
  xlab("Elevation (m)")+
  geom_line(data=pred.delta,aes(x=Elevation.m,y=fit),color="black",linetype="solid")+
  geom_line(data=pred.delta,aes(x=Elevation.m,y=lwr),color="black",linetype="dashed")+
  geom_line(data=pred.delta,aes(x=Elevation.m,y=upr),color="black",linetype="dashed")+
  theme(plot.margin=margin(t=1,r=5.5,b=5.5,l=5.5,unit="pt"))+
  scale_x_continuous(expand=expansion(mult=c(0,0)))

# Combine elevation plots.
elev_plots<-plot_grid(elev_plot1,elev_plot2,ncol=1,rel_heights=c(1,1))

# Save elevation vs genomic offset plot.
ggsave(filename="Elevation_vs_GO_v2.pdf",device=cairo_pdf,plot=elev_plots,units="in",width=7,height=7)

# Convert elevation data from raster to matrix array.
dem<-rasterToPoints(x=dem)

# Convert elevation data into a data frame.
dem<-as.data.frame(dem)

# Rename the elevation field.
colnames(dem)[which(colnames(dem)=="Elevation_PseudoMercator_Clipped")]<-"Elevation"

# # Apply a 3-standard deviation stretch to the elevation data.
# ## Get mean elevation.
# mean_elevation<-mean(dem$Elevation)
# ## Get standard deviation of elevation.
# sd_elevation<-sd(dem$Elevation)
# ## Get the mean plus three standard deviations.
# upper<-mean_elevation+3*sd_elevation
# ## Get the mean minus three standard deviations.
# lower<-mean_elevation-3*sd_elevation
# ## If elevation values exceed the upper limit, set them to the upper limit.
# dem$Elevation<-ifelse(dem$Elevation > upper,upper,dem$Elevation)
# ## If elevation values exceed the lower limit, set them to the lower limit.
# dem$Elevation<-ifelse(dem$Elevation < lower,lower,dem$Elevation)

# Get the maximum absolute value of the transformed genomic offsets.
color_limit<-max(abs(cube_root_sign(df$Delta_GO)))

# Define plotting breaks for the initial genomic offset values.
brks<-seq(from=min(log(df$GO_Initial)),to=max(log(df$GO_Initial)),length.out=5)

# Copy the dem raster data frame to create a masking layer.
mask<-dem

# Create genomic offset map.
go_map<-ggplot()+
  facet_wrap(facets="Year",
             ncol=3, # 4: horizontal setup. 3: vertical setup.
             nrow=4, # 3: horizontal setup. 4: vertical setup.
             scales="fixed",
             strip.position="top",
             dir="h")+
  geom_rect(data=data.frame(xmin=min(dem$x),
                            xmax=max(dem$x),
                            ymin=min(dem$y),
                            ymax=max(dem$y)),
            aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),color=NA,fill="lightblue")+
  geom_raster(data=mask,aes(x=x,y=y),fill="white")+
  geom_raster(data=dem,aes(x=x,y=y,fill=Elevation),alpha=0.5)+
  scale_fill_gradientn(name="Base Layer: Elevation (m)",colors=c("darkgreen","green","yellow","blue","white"))+
  new_scale_fill()+
  geom_sf(data=nation,color="black",fill=NA,linewidth=1)+
  geom_sf(data=states,color="black",fill=NA,linewidth=0.5)+
  geom_sf(data=df_sf,aes(size=log(GO_Initial),fill=cube_root_sign(Delta_GO)),color="black",shape=21,alpha=0.9)+
  scale_size(name="Points: ln(2020 Genomic Offset)",breaks=brks,labels=format(round(x=brks,digits=1),nsmall=1))+
  scale_fill_gradientn(name=expression("Points: (\u0394 Genomic Offset Since 2020)"^~frac(1,3)),colors=c("blue","white","red"),limits=c(-color_limit,color_limit))+
  coord_sf(crs=3857,xlim=c(-13800000,-10750000),ylim=c(3900000,6200000))+
  theme(panel.background=element_rect(fill=NA),
        panel.border=element_rect(colour="black",fill=NA,linewidth=0.5),
        legend.key=element_rect(color=NA,fill=NA),
        plot.title=element_text(face="bold",hjust=0.5,size=24))+
  annotation_north_arrow(location="bl",which_north="true",
                         height=unit(0.25,"in"),
                         pad_x=unit(-0.15,"in"),style=north_arrow_minimal)+
  guides(size=guide_legend(order=1))+
  ggtitle("Genomic Offset Projections")+
  xlab(NULL)+
  ylab(NULL)

# Save map.
## Horizontal setup.
#ggsave(filename="GO_Map_IncludesPast_v1.pdf",device=cairo_pdf,plot=go_map,units="in",width=20,height=10.5)
## Vertical setup.
ggsave(filename="GO_Map_IncludesPast_Vertical_v1.pdf",device=cairo_pdf,plot=go_map,units="in",width=15,height=13.5)

# Done!
print("Done!")
