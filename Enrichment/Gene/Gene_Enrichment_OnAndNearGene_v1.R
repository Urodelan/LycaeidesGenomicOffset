#######################################################
### Test for Gene Enrichment Among Significant SNPs ###
#######################################################

# Performs gene enrichment tests for on-and-near-gene SNPs.

# Clear environment.
rm(list=ls())

# Clear graphics.
graphics.off()

##############################
### User-Defined Variables ###
##############################

# Declare working directory.
working_directory<-"/uufs/chpc.utah.edu/common/home/gompert-group2/data/kg_butterflies/Annotation/GeneEnrichment"

# Declare path to file containing SNP p-values.
path_to_SNP_p.values<-"/uufs/chpc.utah.edu/common/home/gompert-group2/data/kg_butterflies/Models/Product/Output/MVNorm_Reg_Coef_Summaries.csv"

# Declare to file containing SNP positions.
path_to_SNP_positions<-"/uufs/chpc.utah.edu/common/home/gompert-group2/data/kg_butterflies/Data/Product/butterfly_SNPs_v3.csv"

# Declare path to SNP structural annotations.
path_to_SNP_structural_annotations<-"/uufs/chpc.utah.edu/common/home/gompert-group2/data/kg_butterflies/Annotation/Output/snp_annotations.txt"

# Declare path to output file.
path_to_output_file<-"/uufs/chpc.utah.edu/common/home/gompert-group2/data/kg_butterflies/Annotation/GeneEnrichment/gene_enrichment_OnAndNearGene_output.csv"

# Set the number of iterations for gene enrichment p-value calculations.
num_iter<-1e6

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

# Read in data.
print(paste0("- Reading in data: ",elapsed.time(ptm)))
## SNP p-values.
p.values<-read.csv(file=path_to_SNP_p.values,check.names=FALSE,stringsAsFactors=FALSE)
## SNP positions.
positions<-read.csv(file=path_to_SNP_positions,stringsAsFactors=FALSE)
## Structural annotations.
structure<-read.delim(file=path_to_SNP_structural_annotations,stringsAsFactors=FALSE)

# Format structural information.
print(paste0("- Formatting structural information: ",elapsed.time(ptm)))
## Add field for on gene or near gene.
structure$Gene<-(structure$neargene==1) | (structure$ongene==1)
## Subset to desired fields.
structure<-structure[,c("scaffold","position","Gene")]
## Capitalize the scaffold field.
colnames(structure)[colnames(structure)=="scaffold"]<-"Scaffold"
## Capitalize the position field.
colnames(structure)[colnames(structure)=="position"]<-"Position"

# Add gene structural information to SNP names.
gene<-merge(x=positions,y=structure)

# Subset gene structural information to desired fields.
gene<-gene[,c("SNP","Gene")]

# Add structural information to the model output information.
data<-merge(x=p.values,y=gene)

# Add a significance field.
data$Significant<-data$`Pr(>|t|) FDR-Adjusted` <= 0.05

# Subset data to desired fields.
data<-data[,c("SNP","Predictor","Gene","Significant")]

# Begin gene enrichment tests.
print(paste0("- Beginning gene enrichment tests: ",elapsed.time(ptm)))

# Get predictor names.
preds<-unique(data$Predictor)

# Remove the intercept from the predictor names.
preds<-preds[preds!="Intercept"]

# Create empty storage data frame.
results<-data.frame(NULL)

# Set seed.
set.seed(1234)

# Loop through each predictor.
for(p in preds){
  
  # Get results for the predictor.
  pred<-data[data$Predictor==p,]
  
  # Count how many SNPs are gene-associated.
  all_gene<-sum(pred$Gene)
  
  # Calculate the proportion of SNPs that are gene-associated.
  prop_all_gene<-all_gene/nrow(pred)
  
  # Count how many SNPs were significantly associated with the predictor.
  sig<-sum(pred$Significant)
  
  # Count how many gene-associated SNPs were significantly associated with the predictor.
  sig_gene<-sum(pred$Significant & pred$Gene)
  
  # Create an empty storage vector for the null distribution.
  null<-c()
  
  # Loop through each number of user-specified iterations.
  for(i in 1:num_iter){
    
    # Get a random subset of SNPs of length significance-SNPs.
    sub<-pred[sample.int(n=nrow(pred),size=sig,replace=FALSE),]
    
    # Count how many of the randomly-subsetted SNPs are gene-associated.
    sub_gene<-sum(sub$Gene)
    
    # Generate the null distribution.
    null<-c(null,sub_gene)
    
  }
  
  # Calculate the upper tail probability.
  prob<-sum(null >= sig_gene)/length(null)
  
  # Summarize the null distribution.
  null.summary<-as.data.frame(t(quantile(x=null/sig,probs=c(0.025,0.25,0.5,0.75,0.975))))
  
  # Rename null summary fields.
  colnames(null.summary)<-paste0("Quantile_",as.numeric(sub(pattern="%$",replacement="",x=colnames(null.summary)))/100)
  
  # Calculate the IQR of the null distribution.
  null.iqr<-null.summary$Quantile_0.75-null.summary$Quantile_0.25
  
  # Calculate the mean of the null distribution.
  null.mean<-mean(null/sig)
  
  # Calculate the SD of the null distribution.
  null.sd<-stats::sd(null/sig)
  
  # Combine results into a data frame.
  result<-data.frame(Predictor=p,
                     Number_SNPs=nrow(pred),
                     Number_Gene.Associated_SNPs=all_gene,
                     Proportion_SNPs_Gene.Associated=prop_all_gene,
                     Number_Significant_SNPs=sig,
                     Number_Significant_Gene.Associated_SNPs=sig_gene,
                     Proportion_Significant_SNPs_Gene.Associated=sig_gene/sig,
                     Proportion_SNPs_Gene.Associated_Null.mean=null.mean,
                     Proportion_SNPs_Gene.Associated_Null.SD=null.sd,
                     Proportion_SNPs_Gene.Associated_Null=null.summary,
                     Proportion_SNPs_Gene.Associated_Null.IQR=null.iqr,
                     Number_Iterations=num_iter,
                     p.value=prob,
                     stringsAsFactors=FALSE)
  
  # Add results to the storage data frame.
  results<-rbind(results,result)
  
  # Report progress.
  print(paste0("-- Completed: ",p," (",elapsed.time(ptm),")"))
  
}

# Write out results.
print(paste0("- Writing out results: ",elapsed.time(ptm)))
write.csv(x=results,file=path_to_output_file,row.names=FALSE)

# Done!
print(paste0("Done! (",elapsed.time(ptm),")"))
