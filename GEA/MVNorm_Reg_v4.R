####################################################################################
### Multivariate Normal Regression with Variable Selection by Adjusted R-Squared ###
####################################################################################

# Script for use on the cluster.
# Residual covariance and correlation matrices not
# computed to avoid maxing out RAM on cluster.

# Below is an example command to run this script from the command line.
## R --vanilla --quiet --slave --file=MVNorm_Reg_v4.R > Log_MVNorm_Reg_v4.Rout

# Admixture proportion between melissa and idas as predictor,
# representing population structure. Interactions between climate
# variables and population structure.

# Includes partial R-squared calculations using the formula from Wiki:
# https://en.wikipedia.org/wiki/Coefficient_of_determination#Coefficient_of_partial_determination

# Reduced model in partial R-squared calculations is the full model
# without the predictor. Partial R-squared can be interpreted as the
# percent reduction in unexplained variance when the predictor is
# included in the model, compared to a model without the predictor.

# R-squared is calculated using the canonical R-squared formula from
# "Numerical Ecology"'s RDA section. This calculation of R-squared and
# adjusted R-squared results in the same values from using the eigen values
# approach in vegan. This is because PCA rotation does not affect the total
# variance of the variables. PCAs are now completely avoided, and the RDA is
# simply a multivariate linear regression model (which is itself just a
# series of univariate linear regression models).

# References for PCA rotation having no effect on total variance of the variables:
## 1) https://stats.stackexchange.com/questions/22569/pca-and-proportion-of-variance-explained
## 2) https://ro-che.info/articles/2017-12-11-pca-explained-variance

##############################
### User-Defined Variables ###
##############################

# Declare working directory.
working_directory<-"/uufs/chpc.utah.edu/common/home/gompert-group2/data/kg_butterflies/Models/Product"

# Declare path to butterfly climate data.
path_to_climate_data<-"/uufs/chpc.utah.edu/common/home/gompert-group2/data/kg_butterflies/Data/Product/butterfly_climate_v3.csv"

# Declare path to butterfly allele data.
path_to_allele_data<-"/uufs/chpc.utah.edu/common/home/gompert-group2/data/kg_butterflies/Data/Product/butterfly_alleles_v3.csv"

####################
### Begin Script ###
####################

# Report starting script.
print("Starting script...")

# Set working directory.
setwd(working_directory)

# Define function for displaying elasped time.
## Arguments:
### result_of_initial_proc.time: Result of an earlier call of the proc.time function.
## Returns: A character string of formatted elasped time.
elapsed.time<-function(result_of_initial_proc.time){
  # Get number of seconds since the initial proc.time result.
  tot_num_sec<-c(proc.time()-result_of_initial_proc.time)[3]
  # Get number of full hours since initial proc.time result.
  num_hrs<-floor(tot_num_sec/3600)
  # Get number of full minutes minus full hours since initial proc.time result.
  num_min<-floor((tot_num_sec-3600*num_hrs)/60)
  # Get number of full seconds minus full hours and full minutes since initial proc.time result.
  num_sec<-floor(tot_num_sec-3600*num_hrs-60*num_min)
  # Get formatted elapsed time.
  elapsed_time<-paste0(num_hrs," hrs, ",num_min," min, ",num_sec," sec")
  # Return character string of formatted elasped time.
  return(elapsed_time)
}

# Store results of initial proc.time for timing purposes.
ptm<-proc.time()

### Create predictor matrix.

# Begin creating predictor matrix.
print(paste0("- Creating predictor matrix: ",elapsed.time(ptm)))

# Read in butterfly climate data.
climate<-read.csv(file=path_to_climate_data)

# Get the means of predictor values.
## Admixture.
admix_mean<-mean(climate$Admix.prop)
## Temperature.
temp_mean<-mean(climate$Temp.C)
## Precipitation.
precip_mean<-mean(climate$Precip.mm)

# Get the standard deviations of predictor values.
## Admixture.
admix_SD<-sd(climate$Admix.prop)
## Temperature.
temp_SD<-sd(climate$Temp.C)
## Precipitation.
precip_SD<-sd(climate$Precip.mm)

# Create a data frame containing the scaling parameters for the predictor variables.
scaling<-data.frame(Predictor=c("Admixture","Temperature","Precipitation"),
                    Mean=c(admix_mean,temp_mean,precip_mean),
                    SD=c(admix_SD,temp_SD,precip_SD))

# Center and scale predictor variables.
## Admixture.
climate$Admixture<-(climate$Admix.prop-admix_mean)/admix_SD
## Temperature.
climate$Temperature<-(climate$Temp.C-temp_mean)/temp_SD
## Precipitation.
climate$Precipitation<-(climate$Precip.mm-precip_mean)/precip_SD

# Create the predictor matrix with interactions between the predictor variables.
X<-model.matrix(~Admixture*Temperature*Precipitation,data=climate)

# Rename the intercept field.
colnames(X)[colnames(X)=="(Intercept)"]<-"Intercept"

# Add population information as row names.
row.names(X)<-climate$PopID

### Prepare response data.

# Read in response data.
print(paste0("- Preparing response data: ",elapsed.time(ptm)))

# Read in allele frequencies.
Y<-read.csv(file=path_to_allele_data,row.names=1)

# Format the response data as a matrix.
Y<-as.matrix(Y)

# Throw an error if the order of populations does not
# match between the predictor and response matrices.
if(!identical(row.names(X),row.names(Y))) stop("The order of populations does not match between the predictor and response matrices.")

# Throw an error if any allele frequencies show no variation across populations.
if(any(apply(X=Y,MARGIN=2,FUN=stats::sd)==0)) stop("Some allele frequencies have no variation across populations.")

### Fit multivariate normal regression model.

# Generate all combinations of predictors.
## Create an empty storage list to hold combination vectors.
combs<-vector(mode="list",length=ncol(X))
## Set the first combination to just the intercept.
combs[[1]]<-list(1)
## Loop through the remaining elements of the combinations list.
for(i in 2:ncol(X)){
  ## Get all combinations of a given number of predictor variables.
  combs[[i]]<-combn(x=2:ncol(X),m=i-1,simplify=F)
}
## Turn the list of lists into a list of vectors.
combs<-unlist(combs,recursive=F)
## Loop through all but the first combination in the list.
for(i in 2:length(combs)){
  ## Add the intercept to each combination.
  combs[[i]]<-c(1,combs[[i]])
}

# Report beginning to fit models.
print(paste0("- Beginning to fit models for each predictor combination: ",elapsed.time(ptm)))

# Report initial progress.
print(paste0("-- Modeling progress: 0 of ",length(combs)," (",elapsed.time(ptm),")"))

# Define function for calculating adjusted R-squared from R-squared,
# number of observations, and degrees of freedom.
R2Adj<-function(R_squared,n,df) 1-(1-R_squared)*(n-1)/df

# Create a storage data frame for R-squared calculations.
R2_df<-data.frame(NULL)

# Loop through each combination of predictors.
for(k in 1:length(combs)){
  
  # Get the current combination of predictors.
  comb<-combs[[k]]
  
  # Subset the matrix of predictor values to the given predictor combination.
  X_partial<-X[,comb,drop=F]
  
  # Fit multivariate normal regression.
  multi_fit<-lm(Y~X_partial-1)
  
  # Calculate R-squared.
  ## Calculate the total sum of squares of the fitted values matrix.
  total_sum_of_squares_of_Y_hat<-sum(scale(fitted(multi_fit),scale=F)^2)
  ## Calculate the total sum of squares of the response matrix.
  total_sum_of_squares_of_Y<-sum(scale(Y,scale=F)^2)
  ## Calculate R-squared as the ratio between the total sum of squares of the
  ## fitted values matrix and the total sum of squares of the response matrix.
  R2<-total_sum_of_squares_of_Y_hat/total_sum_of_squares_of_Y
  
  # Calculate adjusted R-squared.
  R2.adj<-R2Adj(R_squared=R2,n=nrow(Y),df=df.residual(multi_fit))
  
  # Get predictor combination string.
  ## Get predictor names.
  pred_comb_text<-colnames(X_partial)
  ## Concatenate predictor names with " + ".
  pred_comb_text<-paste(pred_comb_text,collapse=" + ")
  
  # Combine R-squared and adjusted R-squared into a data frame.
  R2.both<-data.frame(Predictor_Combination_ID=k,
                      Best="No",
                      Predictor_Combination=pred_comb_text,
                      Number_of_Predictors=ncol(X_partial),
                      SS.fitted=total_sum_of_squares_of_Y_hat,
                      SS.res=total_sum_of_squares_of_Y-total_sum_of_squares_of_Y_hat,
                      SS.tot=total_sum_of_squares_of_Y,
                      R2=R2,
                      R2.adj=R2.adj)
  
  # Store R-squared information in the storage data frame.
  R2_df<-rbind(R2_df,R2.both)
  
  # Report modeling progress.
  print(paste0("-- Modeling progress: ",k," of ",length(combs)," (",elapsed.time(ptm),")"))
  
}

# Identify best-fitting model by adjusted R-squared.
print(paste0("- Identifying best-fitting model by adjusted R-squared: ",elapsed.time(ptm)))
## Get the maximum adjusted R-squared value.
R2.adj_max_value<-max(R2_df$R2.adj)
## Get the indices of records which have the highest adjusted R-squared.
R2.adj_max_records<-which(R2_df$R2.adj==R2.adj_max_value)
## Get the number of predictors the highest-adjusted R-squared models have.
R2.adj_max_records_num_pred<-R2_df$Number_of_Predictors[R2.adj_max_records]
## Subset highest-adjusted R-squared models to those with the fewest predictors.
R2.adj_max_records<-R2.adj_max_records[R2.adj_max_records_num_pred==min(R2.adj_max_records_num_pred)]
## If multiple highest-adjusted R-squared models have the same number of predictors,
## select the first of these models.
R2.adj_max_record<-R2.adj_max_records[1]
## Denote this model as the best-fitting model.
R2_df$Best[R2.adj_max_record]<-"Yes"

# Subset predictor matrix to predictors included in the best-fitting model.
print(paste0("- Subsetting predictor matrix to predictors included in the best-fitting model: ",elapsed.time(ptm)))
## Get the predictor combination ID of the best-fitting model.
comb_ID_best<-R2_df$Predictor_Combination_ID[R2_df$Best=="Yes"]
## Get the predictor combination of the best-fitting model.
comb_best<-combs[[comb_ID_best]]
## Subset the matrix of predictor values to the given predictor combination.
X_partial<-X[,comb_best,drop=F]

# Report fitting final model.
print(paste0("- Fitting final model: ",elapsed.time(ptm)))

# Fit multivariate normal regression.
multi_fit<-lm(Y~X_partial-1)

# Get regression coefficients from the multivariate normal regression.
## Extract point estimates of model coefficients.
B<-t(coef(multi_fit))
## Remove the leading X_partial string from the coefficient names.
colnames(B)<-sub(pattern="^X_partial",replacement="",x=colnames(B))
## Extract model coefficient summaries.
coefs<-coef(summary(multi_fit))
## Add a field denoting the SNP (response) for each matrix within the list.
coefs<-lapply(1:length(coefs),function(id) cbind(coefs[[id]],SNP=names(coefs)[id]))
## Turn row names into a field for each matrix within the list.
coefs<-lapply(1:length(coefs),function(id) cbind(coefs[[id]],Predictor=row.names(coefs[[id]])))
## Combine all matrices within the list into a single data frame.
coefs<-as.data.frame(do.call(rbind,coefs))
## Reset the row names.
row.names(coefs)<-1:nrow(coefs)
## Re-order fields of coefficients data frame.
coefs<-coefs[,c((ncol(coefs)-1):ncol(coefs),1:(ncol(coefs)-2))]
## Remove the leading "Response " string from the SNP (response) names.
coefs$SNP<-sub(pattern="^Response ",replacement="",x=coefs$SNP)
## Remove the leading X_partial string from the predictor names.
coefs$Predictor<-sub(pattern="^X_partial",replacement="",x=coefs$Predictor)

# Adjust the p-values of the regression coefficients for each predictor
# using the false discovery rate (FDR) method.
## Create an empty field for FDR-adjusted p-values.
coefs$`Pr(>|t|) FDR-Adjusted`<-NA
## Loop through each predictor included in the final model.
for(i in unique(coefs$Predictor)){
  ## Get the records in the coefficients data frame which belong to the predictor.
  predictor_coef_records<-which(coefs$Predictor==i)
  ## Get the p-values of the predictor's regression coefficients.
  raw_p_values<-coefs$`Pr(>|t|)`[predictor_coef_records]
  ## Adjust the p-values using the FDR method.
  adjusted_p_values<-p.adjust(p=raw_p_values,method="fdr")
  ## Add the adjusted p-values to the storage field.
  coefs$`Pr(>|t|) FDR-Adjusted`[predictor_coef_records]<-adjusted_p_values
}

# # Check residual correlation estimates from the multivariate regression.
# ## Get the residual covariance matrix.
# Sigma_hat<-estVar(multi_fit)
# ## Get the residual standard deviations from the residual covariance matrix.
# residual_sd<-sqrt(diag(Sigma_hat))
# ## Get the residual standard deviations in a data frame.
# residual_sd<-data.frame(SNP=names(residual_sd),Residual_SD=residual_sd,row.names=1:length(residual_sd))
# ## Get the residual correlation matrix from the residual covariance matrix.
# R_hat<-cov2cor(Sigma_hat)

# Calculate partial R-squared statistic for each predictor in the best-fitting model.
## Create an empty storage data frame for partial R-squared statistics.
R2.partial_df<-data.frame(NULL)
## Get the residual sum of squares for the full (best-fitting) model.
SS.resid.full<-R2_df$SS.res[R2_df$Best=="Yes"]
## Get the predictor combination for the full (best-fitting) model.
pred_comb_best<-R2_df$Predictor_Combination[R2_df$Best=="Yes"]
## If the best-fitting model includes more than just the intercept.
if(ncol(X_partial) > 1){
  ## Loop through each predictor in the best-fitting model, excluding the intercept.
  for(i in colnames(X_partial)[-1]){
    ## Get the predictor string for the reduced model by removing the predictor
    ## of interest from the predictor combination of the full (best-fitting) model.
    reduced_mod_predictors<-sub(pattern=paste0(" \\+ ",i),replacement="",x=pred_comb_best)
    ## Get the residual sum of squares for the reduced model.
    SS.resid.reduced<-R2_df$SS.res[R2_df$Predictor_Combination==reduced_mod_predictors]
    ## Calculate the partial R-squared for the predictor.
    R2.partial<-(SS.resid.reduced-SS.resid.full)/SS.resid.reduced
    ## Add the partial R-squared information to the storage data frame.
    R2.partial_df<-rbind(R2.partial_df,data.frame(Predictor=i,
                                                  SS.resid.reduced=SS.resid.reduced,
                                                  SS.resid.best=SS.resid.full,
                                                  R2.partial=R2.partial))
  }
} else { ## If the best-fitting model includes just the intercept.
  ## Put a note in the storage data frame.
  R2.partial_df<-data.frame(Note="Best fit model is intercept only.")
}

# Calculate genomic offset.
### Euclidean distance between observed and predicted values.
GO<-apply(X=resid(multi_fit),MARGIN=1,FUN=function(x) sqrt(sum(x^2)))
### Get genomic offset values in a data frame.
GO<-data.frame(PopID=names(GO),
               Genomic_Offset=GO,
               row.names=1:length(GO))

# Write out model information.
print(paste0("- Writing out model information: ",elapsed.time(ptm)))
## Scaling parameters for climate predictor values.
write.csv(x=scaling,file="Output/Scaling.csv",row.names=F)
## Regression coefficient estimates.
write.csv(x=B,file="Output/MVNorm_Reg_Coef_Estimates.csv")
## Regression coefficient summaries.
write.csv(x=coefs,file="Output/MVNorm_Reg_Coef_Summaries.csv",row.names=F)
# ## Residual standard deviation estimates.
# write.csv(x=residual_sd,file="Output/Residual_SDs.csv",row.names=F)
# ## Residual correlation matrix estimate.
# write.csv(x=R_hat,file="Output/R_hat.csv")
## Genomic offset values.
write.csv(x=GO,file="Output/Genomic_Offset.csv",row.names=F)
## R-squared and adjusted R-squared values.
write.csv(x=R2_df,file="Output/R_squared.csv",row.names=F)
## Partial R-squared values.
write.csv(x=R2.partial_df,file="Output/Partial_R_squared.csv",row.names=F)

# Done!
print(paste0("Done! (",elapsed.time(ptm),")"))
