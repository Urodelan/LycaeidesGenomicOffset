#####################################################
### Test for GO Enrichment Among Significant SNPs ###
#####################################################

# Performs GO enrichment tests for on-or-near-gene SNPs.
# Uses traditional resampling approach for enrichment test.

# Clear environment.
rm(list=ls())

# Clear graphics.
graphics.off()

##############################
### User-Defined Variables ###
##############################

# Declare working directory.
working_directory<-"/uufs/chpc.utah.edu/common/home/gompert-group2/data/kg_butterflies/Annotation/GO_Enrichment/Iterative"

# Declare path to file containing SNP p-values.
path_to_SNP_p.values<-"/uufs/chpc.utah.edu/common/home/gompert-group2/data/kg_butterflies/Models/Product/Output/MVNorm_Reg_Coef_Summaries.csv"

# Declare to file containing SNP positions.
path_to_SNP_positions<-"/uufs/chpc.utah.edu/common/home/gompert-group2/data/kg_butterflies/Data/Product/butterfly_SNPs_v3.csv"

# Declare path to SNP functional annotations.
path_to_SNP_functional_annotations<-"/uufs/chpc.utah.edu/common/home/gompert-group2/data/kg_butterflies/Annotation/Output/snp_functions.txt"

# Declare path to output file.
path_to_output_file<-"/uufs/chpc.utah.edu/common/home/gompert-group2/data/kg_butterflies/Annotation/GO_Enrichment/Iterative/go_enrichment_traditional_OnOrNearGene_output.csv"

# Set the number of iterations for gene enrichment p-value calculations.
num_iter<-10000

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
## Functional annotations.
structure<-read.delim(file=path_to_SNP_functional_annotations,stringsAsFactors=FALSE)

# Format structural information.
print(paste0("- Formatting structural information: ",elapsed.time(ptm)))
# ## Add field for on gene or near gene.
# structure$Gene<-(structure$neargene==1) | (structure$ongene==1)
## Subset SNPs to just those on or near genes.
structure<-structure[structure$location=="on-gene" | structure$location=="near-gene",]
## Remove SNPs without GO numbers.
structure<-structure[structure$GO_number!="-",]
## Store GO number descriptions for future use.
descriptions<-structure
## Subset to desired fields.
structure<-structure[,c("scaffold","position","GO_number")]
## Capitalize the scaffold field.
colnames(structure)[colnames(structure)=="scaffold"]<-"Scaffold"
## Capitalize the position field.
colnames(structure)[colnames(structure)=="position"]<-"Position"
# ## Capitalize the location field.
# colnames(structure)[colnames(structure)=="location"]<-"Location"
## Get all GO numbers associated with each SNP.
structure<-stats::aggregate(GO_number~Scaffold+Position,data=structure,FUN=function(x) paste(sort(unique(unlist(strsplit(x=x,split=",")))),collapse=","))

# Add gene structural information to SNP names.
gene<-merge(x=positions,y=structure)

# Subset gene structural information to desired fields.
gene<-gene[,c("SNP","GO_number")]

# Rename GO number to GO.
colnames(gene)[colnames(gene)=="GO_number"]<-"GO"

# Convert functional information to long format.
## Create an empty storage data frame.
go<-data.frame(NULL)
## Loop through each SNP.
for(i in 1:nrow(gene)){
  ## Get the SNP name.
  snp_name<-gene[i,"SNP"]
  ## Get the GO numbers.
  go_numbers<-gene[i,"GO"]
  ## Get each individual GO number.
  go_numbers<-strsplit(x=go_numbers,split=",")[[1]]
  ## Add the information to the storage data frame.
  go<-rbind(go,data.frame(SNP=snp_name,GO=go_numbers))
}

# Add structural information to the model output information.
data<-merge(x=p.values,y=go)

# Add a significance field.
## To functional annotation data.
data$Significant<-data$`Pr(>|t|) FDR-Adjusted` <= 0.05
## To model output.
p.values$Significant<-p.values$`Pr(>|t|) FDR-Adjusted` <= 0.05

# Subset data to desired fields.
data<-data[,c("SNP","Predictor","GO","Significant")]

# Prepare data for enrichment tests.
print(paste0("- Preparing data for enrichment tests: ",elapsed.time(ptm)))

# For each predictor, count the number of significant SNPs.
total.sig<-stats::aggregate(Significant~Predictor,data=p.values,FUN=sum)
colnames(total.sig)[colnames(total.sig)=="Significant"]<-"Sig.Total"

# Begin gene enrichment tests.
print(paste0("- Beginning gene enrichment tests: ",elapsed.time(ptm)))

# Get predictor names.
preds<-unique(p.values$Predictor)

# Remove the intercept from the predictor names.
preds<-preds[preds!="Intercept"]

# Get GO numbers.
go_numbers<-unique(data$GO)

# Get a vector containing all SNPs.
all_SNPs<-unique(p.values$SNP)

# Get the number of total SNPs.
num_SNPs<-length(all_SNPs)

# Create empty storage data frame.
results<-data.frame(NULL)

# Set seed.
set.seed(1234)

# Loop through each predictor.
for(p in preds){
  
  # Get results for the predictor.
  pred<-data[data$Predictor==p,]
  
  # Loop through each GO number.
  for(g in 1:length(go_numbers)){
    
    # Get the GO number.
    go_number<-go_numbers[g]
    
    # Create a field denoting whether the SNP is associated with the GO number.
    pred$IsGO<-pred$GO==go_number
    
    # Count how many SNPs are GO-number-associated.
    all_go<-sum(pred$IsGO)
    
    # Calculate the proportion of SNPs that are GO-number-associated.
    prop_all_go<-all_go/num_SNPs
    
    # Count how many SNPs were significantly associated with the predictor.
    sig<-total.sig[total.sig$Predictor==p,"Sig.Total"]
    
    # Count how many GO-number-associated SNPs were significantly
    # associated with the predictor.
    sig_go<-sum(pred$Significant & pred$IsGO)
    
    # Create an empty storage vector for the null distribution.
    null<-c()
    
    # Loop through each number of user-specified iterations.
    for(i in 1:num_iter){
      
      # Get a random subset of SNPs of length significance-SNPs.
      sub_SNPs<-sample(x=all_SNPs,size=sig,replace=FALSE)
      
      # If any of the randomly-selected SNPs are associated with a GO number.
      if(any(pred$SNP %in% sub_SNPs)){
        
        # Subset SNP-GO information to the randomly selected SNPs.
        sub<-pred[pred$SNP %in% sub_SNPs,,drop=FALSE]
        
        # Count how many of the randomly-subsetted SNPs are GO-associated.
        sub_go<-sum(sub$IsGO)
        
      } else {
        
        # If none of the randomly-selected SNPs are associated with a GO number.
        
        # Set the number of randomly-subsetted SNPs that are GO-associated to zero.
        sub_go<-0
        
      }
      
      # Generate the null distribution.
      null<-c(null,sub_go)
      
    }
    
    # Calculate the upper tail probability.
    prob<-sum(null >= sig_go)/length(null)
    
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
                       GO=go_number,
                       Number_SNPs=num_SNPs,
                       Number_GO.Associated_SNPs=all_go,
                       Proportion_SNPs_GO.Associated=prop_all_go,
                       Number_Significant_SNPs=sig,
                       Number_Significant_GO.Associated_SNPs=sig_go,
                       Proportion_Significant_SNPs_GO.Associated=sig_go/sig,
                       Proportion_SNPs_GO.Associated_Null.mean=null.mean,
                       Proportion_SNPs_GO.Associated_Null.SD=null.sd,
                       Proportion_SNPs_GO.Associated_Null=null.summary,
                       Proportion_SNPs_GO.Associated_Null.IQR=null.iqr,
                       Number_Iterations=num_iter,
                       p.value=prob,
                       stringsAsFactors=FALSE)
    
    # Add results to the storage data frame.
    results<-rbind(results,result)
    
  }
  
  # Report progress.
  print(paste0("-- Completed: ",p," (",elapsed.time(ptm),")"))
  
}

# Format output.
print(paste0("- Formatting output: ",elapsed.time(ptm)))

# Format GO number descriptions.
## Subset to desired fields.
descriptions<-descriptions[,c("GO_number","GO_name","GO_description")]
## Remove duplicate records.
descriptions<-unique(descriptions)
## Split descriptions by delimiters.
### GO numbers.
go_numbers_split<-strsplit(x=descriptions$GO_number,split=",")
### GO names.
go_names_split<-strsplit(x=descriptions$GO_name,split="\\|")
### GO descriptions.
go_descriptions_split<-strsplit(x=descriptions$GO_description,split="\\|")
## Check that lengths match up.
### GO numbers and GO names.
if(!identical(sapply(X=go_numbers_split,FUN=length),sapply(X=go_names_split,FUN=length))){
  stop("The number of GO numbers do not match up with the number of GO names.")
}
### GO numbers and GO descriptions.
if(!identical(sapply(X=go_numbers_split,FUN=length),sapply(X=go_descriptions_split,FUN=length))){
  stop("The number of GO numbers do not match up with the number of GO descriptions.")
}
## Collect GO description information into a data frame.
descriptions<-data.frame(GO=unlist(go_numbers_split),
                         GO_Name=unlist(go_names_split),
                         GO_Description=unlist(go_descriptions_split))
## Remove duplicated GO descriptions.
descriptions<-unique(descriptions)
## Through an error if not all GO numbers have a single name and description.
if(any(duplicated(descriptions$GO))) stop("Not all GO numbers have a single name and description.")

# Add GO description information to the results data frame.
results<-merge(results,descriptions)

# Re-order fields.
results<-results[,c("Predictor","GO","GO_Name","GO_Description","Number_SNPs","Number_GO.Associated_SNPs","Proportion_SNPs_GO.Associated","Number_Significant_SNPs","Number_Significant_GO.Associated_SNPs","Proportion_Significant_SNPs_GO.Associated","Proportion_SNPs_GO.Associated_Null.mean","Proportion_SNPs_GO.Associated_Null.SD","Proportion_SNPs_GO.Associated_Null.Quantile_0.025","Proportion_SNPs_GO.Associated_Null.Quantile_0.25","Proportion_SNPs_GO.Associated_Null.Quantile_0.5","Proportion_SNPs_GO.Associated_Null.Quantile_0.75","Proportion_SNPs_GO.Associated_Null.Quantile_0.975","Proportion_SNPs_GO.Associated_Null.IQR","Number_Iterations","p.value")]

# Rename GO field to GO number.
colnames(results)[colnames(results)=="GO"]<-"GO_Number"

# Perform multiple hypothesis correction.
## Create empty storage data frame.
output<-data.frame(NULL)
## Loop through each predictor (in original order).
for(p in unique(p.values$Predictor)){
  ## Omit the intercept term.
  if(p!="Intercept"){
    ## Subset results to the predictor of interest.
    pred<-results[results$Predictor==p,]
    ## Perform FDR p-value adjustment.
    pred$p.value.fdr<-p.adjust(p=pred$p.value,method="fdr")
    ## Order records by GO number.
    pred<-pred[order(pred$GO_Number),]
    ## Add results to storage data frame.
    output<-rbind(output,pred)
  }
}

# Add a significance field.
output$Significant<-output$p.value.fdr <= 0.05

# Write out results.
print(paste0("- Writing out results: ",elapsed.time(ptm)))
write.csv(x=output,file=path_to_output_file,row.names=FALSE)

# Done!
print(paste0("Done! (",elapsed.time(ptm),")"))
