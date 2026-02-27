##############################
### Temporal GO Regression ###
##############################

# Clear environment.
rm(list=ls())

# Clear graphics.
graphics.off()

# Set working directory.
setwd("/Users/kenengoodwin/Desktop/USU/Graduate/Projects/Butterflies/Analyses/Models/GO_Correlations/Temporal")

# Load mixed effects package.
suppressPackageStartupMessages(library(lmerTest))

# Read in data.
df<-read.csv(file="/Users/kenengoodwin/Desktop/USU/Graduate/Projects/Butterflies/Analyses/Models/Product/Output/Process/IncludePast/Genomic_Offset_Projections_IncludesPast_v1.csv")

# Fit mixed effects model.
mod<-lmer(log(Genomic_Offset)~Year+(1|PopID),data=df,REML=TRUE)

# Produce residual quantile plot.
pdf(file="Output/Plots/qqplot.pdf",width=7,height=5)
qqnorm(y=residuals(mod)/sigma(mod),xlab="Theoretical Quantiles",
       ylab="Standardized residuals",main="Q-Q Residuals")
qqline(y=residuals(mod)/sigma(mod),col="red",lty="dashed")
dev.off()

# Produce residuals vs leverage plot.
pdf(file="Output/Plots/leverage.pdf",width=7,height=5)
plot(x=hatvalues(mod),y=residuals(mod)/sigma(mod),xlab="Leverage",
     ylab="Standardized residuals",main="Residuals vs Leverage")
abline(h=0,col="red",lty="dashed")
dev.off()

# Retrieve random effects.
## Extract estimates.
R<-coef(mod)$PopID
## Rename intercept field.
colnames(R)[colnames(R)=="(Intercept)"]<-"Intercept"
## Rename slope field.
colnames(R)[colnames(R)=="Year"]<-"Slope"
## Create population field.
R$PopID<-as.factor(row.names(R))

# Retrieve fixed effects.
## Extract estimates.
B<-mod@beta
## Create data frame.
B<-data.frame(Intercept=B[1],Slope=B[2])

# Define prediction function for mixed-effects model.
predict.lmer<-function(X,mod,conf=0.95){
  # Extract regression coefficients.
  B<-mod@beta
  # Extract covariance matrix.
  Sigma<-stats::vcov(mod)
  # Retrieve degrees of freedom.
  df<-stats::df.residual(mod)
  # Create storage vectors.
  ## Fitted values.
  fit<-c()
  ## Upper limits.
  upr<-c()
  ## Lower limits.
  lwr<-c()
  # Loop through predictor record.
  for(i in 1:nrow(X)){
    # Retrieve predictor values.
    x<-unlist(X[i,])
    # Calculate expected value.
    mu<-B %*% x
    # Calculate standard error.
    se<-as.vector(sqrt(t(x) %*% Sigma %*% x))
    # Derive critical value.
    t<-stats::qt(p=(conf+1)/2,df=df)
    # Store calculations.
    ## Fitted value.
    fit<-c(fit,mu)
    ## Upper limit.
    upr<-c(upr,mu+t*se)
    ## Lower limit.
    lwr<-c(lwr,mu-t*se)
  }
  # Collect results.
  results<-data.frame(mu=fit,lwr=lwr,upr=upr)
  # Return results.
  return(results)
}

# Generate predictions.
## Get 5% of the year range.
buffer<-diff(range(df$Year))*0.05
## Create data to predict on.
X<-data.frame(Intercept=1,Year=seq(from=min(df$Year)-buffer,
                                   to=max(df$Year)+buffer,
                                   by=1e-2))
## Predict genomic offset by year.
## 95% confindence interval.
pred.go95<-cbind(X,predict.lmer(X=X,mod=mod,conf=0.95))
## 50% confidence interval.
pred.go50<-cbind(X,predict.lmer(X=X,mod=mod,conf=0.50))

# Load color functions.
library(scales)

# Generate colors.
## Get number of random effect colors.
N<-nrow(R)
## Store random effect population names.
name<-R$PopID
## Generate one color per population.
cols<-hue_pal()(N)
## Supply site names to color vector.
names(cols)<-name

# Generate x-axis breaks.
bks<-seq(from=min(df$Year),to=max(df$Year),by=10)

# Generate y-axis breaks.
## Get the minimum and maximum discharge values.
lim<-range(log10(df$Genomic_Offset))
## Compute the range of the discharge values.
rng<-diff(x=lim)
## Round down the minimum discharge value.
lim[1]<-floor(lim[1])
## Round up the maximum discharge value.
lim[2]<-ceiling(lim[2])
## Create a vector containing the range of discharge values.
lim<-lim[1]:lim[2]
## Back-transform the range of discharge values.
lim<-10^lim
## Create an empty storage vector.
y<-c()
## Loop through each log10 increment.
for(i in 1:(length(lim)-1)){
  ## Retrieve the lower end of the range.
  lwr<-lim[i]
  ## Retrieve the upper end of the range.
  upr<-lim[i+1]
  ## Generate a range of values on the original scale.
  vals<-seq(from=lwr,to=(upr-lwr),by=lwr)
  ## Store values in storage vector.
  y<-c(y,vals)
}

# Load ggplot2.
suppressPackageStartupMessages(library(ggplot2))

# Load Cairo.
library(Cairo)

# Create plot.
plt<-ggplot()+
  geom_point(aes(x=Year,y=log(Genomic_Offset),color=PopID),shape=1,alpha=0.25,data=df)+
  geom_abline(aes(intercept=Intercept,slope=Slope,color=PopID),data=R,
              linewidth=0.25,alpha=0.25,linetype="dashed")+
  geom_abline(aes(intercept=Intercept,slope=Slope),data=B,
              linewidth=0.75,linetype="solid")+
  geom_line(aes(x=Year,y=upr),data=pred.go50,
            linewidth=0.5,linetype="dashed")+
  geom_line(aes(x=Year,y=lwr),data=pred.go50,
            linewidth=0.5,linetype="dashed")+
  geom_line(aes(x=Year,y=upr),data=pred.go95,
            linewidth=0.375,linetype="dotted")+
  geom_line(aes(x=Year,y=lwr),data=pred.go95,
            linewidth=0.375,linetype="dotted")+
  scale_x_continuous(breaks=bks,expand=expansion(mult=0))+
  scale_y_continuous(name="Genomic Offset",breaks=log(y),labels=y)+
  xlab(NULL)+
  ggtitle("Genomic Offset Projections",subtitle="Colored dashed lines: Population predictions \u2013 Black solid line: Expected values\nBlack dashed lines: 50% confidence limits \u2013 Black dotted lines: 95% confidence limits")+
  theme_light()+
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1),
        legend.title=element_text(hjust=0.5),
        legend.text=element_text(size=10),
        plot.title=element_text(hjust=0.5,face="bold"),
        plot.subtitle=element_text(hjust=0.5),
        panel.grid.minor=element_blank(),
        legend.position="none")

# Display plot.
print(plt)

# Save plot.
ggsave(filename=paste0("Output/Plots/Temporal_GO_v1.pdf"),
       plot=plt,device=cairo_pdf,
       width=7,height=5,units="in")

# Get model summaries.
## Fixed effects.
output.fixed<-summary(mod)$coefficients
## Random effects.
output.random<-summary(mod)$varcor

# Format random effects summary.
output.random<-data.frame(Groups=c(names(output.random),"Residual"),
                          SD=c(unname(attr(output.random$PopID,"stddev")),
                               attr(output.random,"sc")))

# Write out model summaries.
## Fixed effects.
write.csv(x=output.fixed,file="Output/Tables/fixed_effects.csv",row.names=TRUE)
## Random effects.
write.csv(x=output.random,file="Output/Tables/random_effects.csv",row.names=FALSE)

# Done!
print("Done!")
