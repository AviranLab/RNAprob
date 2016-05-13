## Generate data for the cross-validation analysis - from previously simulated data
# Author: Mirko Ledda | Date: 09/16/2015 | Version: v1.0

## Remove all active variables from the workspace
rm(list = ls())

# Set directory/file paths ------------------------------------------------
source_dir <- 'src/' # location of source scripts
root_ct <- 'data/ct_sequence_file/' # ct data directory
fn_name <- 'data/data_pool_name.txt' # RNA names
root_data <- 'data/generated_data/sim_data_from_std_fit_distributions/' # sim shape data directory

# Output directory
root_out <- 'output/sim_data_noise_crossval/'
dir.create(root_out)
root_out <- paste(root_out,'percentile_',toString(n_bin),'/',sep='')
dir.create(root_out)

# Variable inputs --------------------------------------------------------
rna_ix <- 1:23 # Select RNAs based on name_list
err_multiplier <- c(0,1,2,3,4,5,6) # Multiplier for the variance of the normal distribution of random errors
n_bin <- 5 # number of bins
n_rep <- 100 # number of simulation repeats
model_type <- "ternary" # either binary or ternary

# Set the path to error estimates (from Nathan)
err_fn <- 'data/shape_error_estimates_from_Nathan.csv'

# Finish setting the output folder
root_out <- paste(root_out,model_type,'/',sep='')
dir.create(root_out)
root_out_base <- root_out

# Load sources and data ---------------------------------------------------
# Get SHAPE and ct filename in the folder
name_list <- as.character(read.table(fn_name)$V1)
name_list <- name_list[1:23] # No Siegfield 16S RNA

src_fn <- list.files(source_dir) # List the files in the source directory
for (i in 1:length(src_fn)) source(paste0(source_dir,src_fn[i])) # Load the source files

# Loop across selected RNA -------------------------------------------------------
reac_range_store <- array(NA,dim=c(n_bin+1,length(rna_ix),n_rep))
for (i in 1:length(rna_ix)) {
  
  root_out <- paste(root_out_base,name_list[rna_ix[i]],'/',sep='')
  dir.create(root_out)
  root_out_rna_base <- root_out

  # loop across simulation repeats
  for (k in (1:n_rep)) {
    
    # Read the shape reactivity file
    if (k < 11) {
      shape_data <- read.table(paste(root_data,name_list[rna_ix[i]],'/',model_type,'_model/shape-00',toString(k-1),'.shape',sep=""))
    } else {
      shape_data <- read.table(paste(root_data,name_list[rna_ix[i]],'/',model_type,'_model/shape-0',toString(k-1),'.shape',sep=""))
    }

    #Set NaN reactivities
    nan_reac_filt <- shape_data$V2==-999
    Reactivities <- shape_data$V2
    Reactivities[nan_reac_filt] <- NA
    
    reac <- Reactivities
    
    # Compute the percentiles -------------------------------------------------
    reac_range <- quantile(reac,seq(0,1,1/n_bin),na.rm = T)
  
    reac_range_store[,i,k] <- reac_range
    
    ## With the percentile set to NA
    # Loop across ranges ------------------------------------------------------
    root_out_missing_base <- paste(root_out_rna_base,'missing_percentile/',sep='')
    dir.create(root_out_missing_base)
    root_out_missing_base2 <- paste(root_out_rna_base,'only_percentile/',sep='')
    dir.create(root_out_missing_base2)
    
    ## With the percentile set to NA
    for (l in 1:(length(reac_range)-1)) {
      
      root_out <- paste(root_out_missing_base,'percentile_',toString(l),'/',sep='')
      dir.create(root_out)
      
      reac_ix <- (reac>=reac_range[l]) & (reac<=reac_range[l+1]) & !is.na(reac)
      
      reac_noise <- reac
      reac_noise[reac_ix] <- NA
      
      # Write the shape file
      if (k < 11) {
        write_shape(fn = paste(root_out,'shape-00',toString(k-1),sep=""), x=reac_noise)
      } else {
        write_shape(fn = paste(root_out,'shape-0',toString(k-1),sep=""), x=reac_noise)
      }
    }
    
    ## Writing only the selected percentile
    for (l in 1:(length(reac_range)-1)) {
      
      root_out <- paste(root_out_missing_base2,'percentile_',toString(l),'/',sep='')
      dir.create(root_out)
      
      reac_ix <- (reac>=reac_range[l]) & (reac<=reac_range[l+1]) & !is.na(reac)
      
      reac_noise <- rep(NA,length(reac))
      reac_noise[reac_ix] <- reac[reac_ix]
      
      # Write the shape file
      if (k < 11) {
        write_shape(fn = paste(root_out,'shape-00',toString(k-1),sep=""), x=reac_noise)
      } else {
        write_shape(fn = paste(root_out,'shape-0',toString(k-1),sep=""), x=reac_noise)
      }
    }
  }
}


