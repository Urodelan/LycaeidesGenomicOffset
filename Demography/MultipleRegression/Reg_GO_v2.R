############################################################################################
### Multiple Linear Regression with Variable Selection by Leave-One-Out Cross Validation ###
############################################################################################

# Response: Natural log of genomic offset
# Predictors: Intercept, elevation, effective population size, migration, and diversity.
# No interactions.

# Variable selection by leave-one-out cross validation log-likelihood.

# v2 writes out regression coefficient summaries for all models.

# Clear R environment.
rm(list=ls())

# Clear graphics.
graphics.off()

##############################
### User-Defined Variables ###
##############################

# Declare working directory.
working_directory<-"/Users/kenengoodwin/Desktop/USU/Graduate/Projects/Butterflies/Analyses/Models/GO_Correlations/GO_Regression"

# Declare path to butterfly genomic offset covariate data.
path_to_covariate_data<-"/Users/kenengoodwin/Desktop/USU/Graduate/Projects/Butterflies/Analyses/Models/GO_Correlations/GO_Covariates.csv"

# Declare response variable.
response<-"ln_GO"

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

# Read in genomic offset covariate data.
covariates<-read.csv(file=path_to_covariate_data)

# Rename the elevation field.
colnames(covariates)[colnames(covariates)=="Elevation.m"]<-"Elevation"

# Get the means of predictor values.
## Elevation.
elevation_mean<-mean(covariates$Elevation)
## Effective population size.
Ne_mean<-mean(covariates$Ne)
## Migration.
migration_mean<-mean(covariates$Migration)
## Diversity.
diversity_mean<-mean(covariates$Diversity)

# Get the standard deviations of predictor values.
## Elevation.
elevation_SD<-sd(covariates$Elevation)
## Effective population size.
Ne_SD<-sd(covariates$Ne)
## Migration.
migration_SD<-sd(covariates$Migration)
## Diversity.
diversity_SD<-sd(covariates$Diversity)

# Create a data frame containing the scaling parameters for the predictor variables.
scaling<-data.frame(Predictor=c("Elevation","Ne","Migration","Diversity"),
                    Mean=c(elevation_mean,Ne_mean,migration_mean,diversity_mean),
                    SD=c(elevation_SD,Ne_SD,migration_SD,diversity_SD))

# Center and scale predictor variables.
## Elevation.
covariates$Elevation<-(covariates$Elevation-elevation_mean)/elevation_SD
## Effective population size.
covariates$Ne<-(covariates$Ne-Ne_mean)/Ne_SD
## Migration.
covariates$Migration<-(covariates$Migration-migration_mean)/migration_SD
## Diversity.
covariates$Diversity<-(covariates$Diversity-diversity_mean)/diversity_SD

# Create the predictor matrix with interactions between the predictor variables.
X<-model.matrix(~Elevation+Ne+Migration+Diversity,data=covariates)

# Rename the intercept field.
colnames(X)[colnames(X)=="(Intercept)"]<-"Intercept"

# Add population information as row names.
row.names(X)<-covariates$PopID

### Prepare response data.

# Read in response data.
print(paste0("- Preparing response data: ",elapsed.time(ptm)))

# Throw an error if the user-defined response variable cannot be found in the input data.
if(!(response %in% colnames(covariates))) stop(paste0("The response variable ",response," could not be found in the input data."))

# Format the response data as a matrix.
Y<-as.matrix(covariates[,response,drop=FALSE])

# Add population information as row names.
row.names(Y)<-covariates$PopID

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

# The least-squares estimate of sigma (which is returned by stats::sigma
# using stats::lm) differs from the maximum likelihood estimate of sigma
# (which is used by stats::logLik in AIC calculations). For consistency
# with stats::logLik and stats::AIC, we'll convert the least-squares estimate
# of sigma to the maximum likelihood estimate of sigma in log-likelihood calculations.
# For details, see: https://stats.stackexchange.com/questions/73196/recalculate-log-likelihood-from-a-simple-r-lm-model

# A correct but old function. No need to convert least-squares sigma to MLE sigma,
# just calculate MLE sigma using the new function below.
# # Define function for converting a linear model's least-squares estimate
# # of sigma to the maximum likelihood estimate of sigma.
# sigma.MLE<-function(model){
#   # Extract the model matrix.
#   X<-model.matrix(model)
#   # Get the number of predictors in the model.
#   p<-ncol(X)
#   # Get the number of observations in the model.
#   n<-nrow(X)
#   # Convert the least-squares estimate of sigma to
#   # the maximum likelihood estimate of sigma.
#   sigma.ML<-sigma(model)*sqrt((n-p)/n)
#   # Return the maximum likelihood estimate of sigma.
#   return(sigma.ML)
# }

# Define function for calculating the maximum likelihood
# estimate of sigma from a linear model object.
sigma.MLE<-function(model){
  # Extract the model residuals.
  r<-resid(model)
  # Get the number of observations in the model.
  n<-length(r)
  # Calculate the maximum likelihood estimate of sigma.
  sigma.ML<-sqrt(sum(r^2)/n)
  # Return the maximum likelihood estimate of sigma.
  return(sigma.ML)
}

# Define function for calculating adjusted R-squared from R-squared,
# number of observations, and degrees of freedom.
R2Adj<-function(R_squared,n,df) 1-(1-R_squared)*(n-1)/df

# Create a storage data frame for model fit information.
mod.fit<-data.frame(NULL)

# Create a storage data frame for all model coefficients.
coefs.all<-data.frame(NULL)

# Loop through each combination of predictors.
for(k in 1:length(combs)){
  
  # Get the current combination of predictors.
  comb<-combs[[k]]
  
  # Subset the matrix of predictor values to the given predictor combination.
  X_partial<-X[,comb,drop=F]
  
  # Fit linear regression model.
  fit<-lm(Y~X_partial-1)
  
  # Create an empty storage vector for leave-one-out (LOO) pointwise log-likelihoods.
  pntwise_loo_loglik<-c()
  
  # Loop through each observation.
  for(i in 1:nrow(Y)){
    
    # Remove the observation from the response matrix.
    Y_sub<-Y[-i,,drop=FALSE]
    
    # Remove the observation from the predictor matrix.
    X_partial_sub<-X_partial[-i,,drop=FALSE]
    
    # Throw an error if the row identifiers of the subsetted data do not match
    # between the response and predictor matrices.
    if(!identical(row.names(Y_sub),row.names(X_partial_sub))) stop("The row identifiers of the subsetted data do not match between the response and predictor matrices.")
    
    # Fit linear regression model without the observation.
    fit_sub<-lm(Y_sub~X_partial_sub-1)
    
    # Predict the response of the missing observation.
    pred<-as.vector(X_partial[i,,drop=FALSE] %*% coef(fit_sub))
    
    # Calculate the pointwise leave-one-out log-likelihood
    # and add it to the storage vector.
    pntwise_loo_loglik<-c(pntwise_loo_loglik,
                          dnorm(x=Y[i,],
                                mean=pred,
                                sd=sigma.MLE(fit_sub),
                                log=TRUE))
    
  }
  
  # Calculate the leave-one-out log-likelihood by summing
  # the pointwise leave-one-out log-likelihoods.
  loo_loglik<-sum(pntwise_loo_loglik)
  
  # Calculate R-squared.
  ## Calculate the total sum of squares of the fitted values matrix.
  total_sum_of_squares_of_Y_hat<-sum(scale(fitted(fit),scale=FALSE)^2)
  ## Calculate the total sum of squares of the response matrix.
  total_sum_of_squares_of_Y<-sum(scale(Y,scale=FALSE)^2)
  ## Calculate R-squared as the ratio between the total sum of squares of the
  ## fitted values matrix and the total sum of squares of the response matrix.
  R2<-total_sum_of_squares_of_Y_hat/total_sum_of_squares_of_Y
  
  # Calculate adjusted R-squared.
  R2.adj<-R2Adj(R_squared=R2,n=nrow(Y),df=df.residual(fit))
  
  # Get regression coefficients from the multiple linear regression.
  ## Extract model coefficient summaries.
  coefs<-as.data.frame(coef(summary(fit)))
  ## Turn row names into a field for each matrix within the list.
  coefs<-cbind(Predictor=row.names(coefs),coefs)
  ## Reset the row names.
  row.names(coefs)<-1:nrow(coefs)
  ## If the only predictor is the intercept.
  if(nrow(coefs)==1){
    ## Fix the intercept's name in the predictor field.
    coefs$Predictor[1]<-paste0(coefs$Predictor[1],"Intercept")
  }
  ## Remove the leading X_partial string from the predictor names.
  coefs$Predictor<-sub(pattern="^X_partial",replacement="",x=coefs$Predictor)
  
  # Get predictor combination string.
  ## Get predictor names.
  pred_comb_text<-colnames(X_partial)
  ## Concatenate predictor names with " + ".
  pred_comb_text<-paste(pred_comb_text,collapse=" + ")
  
  ## Add the predictor combination ID to the coefficients data frame.
  coefs<-cbind(Predictor_Combination_ID=k,
               Predictor_Combination=pred_comb_text,
               coefs)
  ## Store coefficients in the storage data frame.
  coefs.all<-rbind(coefs.all,coefs)
  
  # Combine model fit information into a data frame.
  mod.info<-data.frame(Predictor_Combination_ID=k,
                       Best="No",
                       Predictor_Combination=pred_comb_text,
                       Number_of_Predictors=ncol(X_partial),
                       SS.fitted=total_sum_of_squares_of_Y_hat,
                       SS.res=total_sum_of_squares_of_Y-total_sum_of_squares_of_Y_hat,
                       SS.tot=total_sum_of_squares_of_Y,
                       R2=R2,
                       R2.adj=R2.adj,
                       loo_loglik=loo_loglik)
  
  # Store R-squared information in the storage data frame.
  mod.fit<-rbind(mod.fit,mod.info)
  
  # Report modeling progress.
  print(paste0("-- Modeling progress: ",k," of ",length(combs)," (",elapsed.time(ptm),")"))
  
}

# Multiple the LOO log-likelihood by negative two to put it on a similar scale to AIC.
mod.fit$loo_loglik.neg2<-mod.fit$loo_loglik*-2

# Calculate the difference in negative log-likelihood.
mod.fit$delta.loo_loglik.neg2<-mod.fit$loo_loglik.neg2 - min(mod.fit$loo_loglik.neg2)

# Identify best-fitting model by -2 x LOO log-likelihood.
print(paste0("- Identifying best-fitting model by -2 x LOO log-likelihood: ",elapsed.time(ptm)))
## Get the minimum -2 x LOO log-likelihood.
min_loo_loglik.neg2<-min(mod.fit$loo_loglik.neg2)
## Get the indices of records which have the lowest -2 x LOO log-likelihood.
min_loo_loglik.neg2_records<-which(mod.fit$loo_loglik.neg2==min_loo_loglik.neg2)
## Get the number of predictors the lowest -2 x LOO log-likelihood models have.
min_loo_loglik.neg2_records_num_pred<-mod.fit$Number_of_Predictors[min_loo_loglik.neg2_records]
## Subset lowest -2 x LOO log-likelihood models to those with the fewest predictors.
min_loo_loglik.neg2_records<-min_loo_loglik.neg2_records[min_loo_loglik.neg2_records_num_pred==min(min_loo_loglik.neg2_records_num_pred)]
## If multiple lowest -2 x LOO log-likelihood models have the same number of predictors,
## select the first of these models.
min_loo_loglik.neg2_record<-min_loo_loglik.neg2_records[1]
## Denote this model as the best-fitting model.
mod.fit$Best[min_loo_loglik.neg2_record]<-"Yes"

# Subset predictor matrix to predictors included in the best-fitting model.
print(paste0("- Subsetting predictor matrix to predictors included in the best-fitting model: ",elapsed.time(ptm)))
## Get the predictor combination ID of the best-fitting model.
comb_ID_best<-mod.fit$Predictor_Combination_ID[mod.fit$Best=="Yes"]
## Get the predictor combination of the best-fitting model.
comb_best<-combs[[comb_ID_best]]
## Subset the matrix of predictor values to the given predictor combination.
X_partial<-X[,comb_best,drop=F]

# Report fitting final model.
print(paste0("- Fitting final model: ",elapsed.time(ptm)))

# Fit multiple linear regression.
fit<-lm(Y~X_partial-1)

# Get regression coefficients from the multiple linear regression.
## Extract model coefficient summaries.
coefs<-as.data.frame(coef(summary(fit)))
## Turn row names into a field for each matrix within the list.
coefs<-cbind(Predictor=row.names(coefs),coefs)
## Reset the row names.
row.names(coefs)<-1:nrow(coefs)
## If the only predictor is the intercept.
if(nrow(coefs)==1){
  ## Fix the intercept's name in the predictor field.
  coefs$Predictor[1]<-paste0(coefs$Predictor[1],"Intercept")
}
## Remove the leading X_partial string from the predictor names.
coefs$Predictor<-sub(pattern="^X_partial",replacement="",x=coefs$Predictor)

# Calculate partial R-squared statistic for each predictor in the best-fitting model.
## Create an empty storage data frame for partial R-squared statistics.
R2.partial_df<-data.frame(NULL)
## Get the residual sum of squares for the full (best-fitting) model.
SS.resid.full<-mod.fit$SS.res[mod.fit$Best=="Yes"]
## Get the predictor combination for the full (best-fitting) model.
pred_comb_best<-mod.fit$Predictor_Combination[mod.fit$Best=="Yes"]
## If the best-fitting model includes more than just the intercept.
if(ncol(X_partial) > 1){
  ## Loop through each predictor in the best-fitting model, excluding the intercept.
  for(i in colnames(X_partial)[-1]){
    ## Get the predictor string for the reduced model by removing the predictor
    ## of interest from the predictor combination of the full (best-fitting) model.
    reduced_mod_predictors<-sub(pattern=paste0(" \\+ ",i),replacement="",x=pred_comb_best)
    ## Get the residual sum of squares for the reduced model.
    SS.resid.reduced<-mod.fit$SS.res[mod.fit$Predictor_Combination==reduced_mod_predictors]
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

# Rename fields in primary data frame.
## -2 x LOO log-likehood.
colnames(mod.fit)[colnames(mod.fit)=="loo_loglik.neg2"]<-"neg2xloo_loglik"
## Delta of -2 x LOO log-likelihood.
colnames(mod.fit)[colnames(mod.fit)=="delta.loo_loglik.neg2"]<-"delta.neg2xloo_loglik"

# Write out model information.
print(paste0("- Writing out model information: ",elapsed.time(ptm)))
## Scaling parameters for climate predictor values.
write.csv(x=scaling,file=paste0("Output/Scaling.",response,".csv"),row.names=FALSE)
## Regression coefficient summaries for best-fitting model.
write.csv(x=coefs,file=paste0("Output/Reg_Coef_Summaries_Best.",response,".csv"),
          row.names=FALSE)
## Regression coefficient summaries for best-fitting model.
write.csv(x=coefs.all,file=paste0("Output/Reg_Coef_Summaries_All.",response,".csv"),
          row.names=FALSE)
## Model fit information.
write.csv(x=mod.fit,file=paste0("Output/Mod_Fit_Info.",response,".csv"),row.names=FALSE)
## Partial R-squared values.
write.csv(x=R2.partial_df,file=paste0("Output/Partial_R_squared.",response,".csv"),
          row.names=FALSE)

# Done!
print(paste0("Done! (",elapsed.time(ptm),")"))
