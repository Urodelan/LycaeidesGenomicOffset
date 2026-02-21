#########################################################################
### Plot Effective Population Size Against 2020 Genomic Offset Values ###
#########################################################################

# Clear R environment.
rm(list=ls())

# Clear graphics.
graphics.off()

# Set working directory.
setwd("/Users/kenengoodwin/Desktop/USU/Graduate/Projects/Butterflies/Analyses/Models/Product/Output/Process")

# Read in effective population size.
Ne<-read.delim(file="/Users/kenengoodwin/Desktop/USU/Graduate/Projects/Butterflies/Analyses/Data/From_Zach/EffectivePopSize/ldNe.txt",header=TRUE,sep=" ")

# Rename the effective population size fields.
## Population ID field.
colnames(Ne)[colnames(Ne)=="pops"]<-"PopID"
## Quantiles fields.
colnames(Ne)[colnames(Ne) %in% c("median","q25","q75")]<-paste0("Ne.",colnames(Ne)[colnames(Ne) %in% c("median","q25","q75")])

# Rename the population SV to SVY.
Ne$PopID[Ne$PopID=="SV"]<-"SVY"

# Read in genomic offset values.
go<-read.csv(file="Genomic_Offset_Projections.csv")

# Get genomic offset values for 2020.
go<-go[go$Year==2020,]

# Combine the effective population size and genomic offset information.
df<-merge(x=go,y=Ne)

# Fit a linear model between 2020 genomic offset values and effective population size.
lm.Ne<-lm(log(Ne.median)~log(Genomic_Offset),data=df)

# Get 5% of the elevation range.
buffer<-diff(range(log(df$Genomic_Offset)))*0.05

# Create elevation data to predict on.
pred_data<-data.frame(Genomic_Offset=exp(seq(from=min(log(df$Genomic_Offset))-buffer,
                                             to=max(log(df$Genomic_Offset))+buffer,
                                             by=buffer/1e2)))

# Predict the effective population size by 2020 genomic offset.
pred.Ne<-cbind(pred_data,
               as.data.frame(
                 predict(object=lm.Ne,newdata=pred_data,interval="confidence")))

# Display summary of the linear model between 2020 genomic offset
# values and effective population size.
#summary(lm.Ne)
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           6.6100     1.2206   5.415 3.13e-06 ***
# log(Genomic_Offset)  -0.5501     0.3964  -1.388    0.173    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.8538 on 40 degrees of freedom
# Multiple R-squared:  0.04593,	Adjusted R-squared:  0.02208 
# F-statistic: 1.926 on 1 and 40 DF,  p-value: 0.1729

# Load ggplot.
library(ggplot2)

# Load package for saving unicode characters in pdfs.
library(Cairo)

# Create Ne vs GO plot.
Ne_GO_plot<-ggplot(data=df)+
  geom_point(aes(x=log(Genomic_Offset),y=log(Ne.median)))+
  theme_light()+
  ylab("ln(Effective Population Size)")+
  xlab("ln(2020 Genomic Offset)")+
  theme_light()+
  ggtitle("Effective Population Size by Genomic Offset",subtitle="Point: Observation \u2013 Solid line: Predicted value \u2013 Dashed line: 95% confidence interval")+
  geom_line(data=pred.Ne,aes(x=log(Genomic_Offset),y=fit),color="black",linetype="solid")+
  geom_line(data=pred.Ne,aes(x=log(Genomic_Offset),y=lwr),color="black",linetype="dashed")+
  geom_line(data=pred.Ne,aes(x=log(Genomic_Offset),y=upr),color="black",linetype="dashed")+
  theme(plot.margin=margin(t=5.5,r=5.5,b=1,l=8.8,unit="pt"),
        plot.title=element_text(face="bold",hjust=0.5),
        plot.subtitle=element_text(hjust=0.5))+
  scale_x_continuous(expand=expansion(mult=c(0,0)))

# Save plot.
ggsave(filename="Ne_vs_GO.pdf",device=cairo_pdf,plot=Ne_GO_plot,units="in",width=7,height=5)

# Done!
print("Done!")
