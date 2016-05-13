## Generate data for the noise simulation and for the cross-validation analysis
# Author: Mirko Ledda | Date: 09/16/2015 | Version: v1.0

## Remove all active variables from the workspace
rm(list = ls())

# Set directory/file paths ------------------------------------------------
source_dir <- 'src/' # location of source scripts
root_ct <- 'data/ct_sequence_file/' # ct data directory
fn_name <- 'data/data_pool_name.txt' # RNA names
root_data <- 'data/shape_raw/' # shape data directory

# Output directory
root_out <- 'output/real_data_noise_crossval/'
dir.create(root_out)
root_out <- paste(root_out,'percentile_',toString(n_bin),'/',sep='')
dir.create(root_out)
root_out_base <- root_out

# Variable inputs --------------------------------------------------------
rna_ix <- 1:23 # Select RNAs based on name_list
err_multiplier <- c(0,1,2,3,4,5,6) # Multiplier for the variance of the normal distribution of random errors
n_bin <- 5 # number of bins
n_rep <- 10 # number of simulation repeats

# Set the path to error estimates (from Nathan)
err_fn <- 'data/shape_error_estimates_from_Nathan.csv'

# Load sources and data ---------------------------------------------------
# Get SHAPE and ct filename in the folder
name_list <- as.character(read.table(fn_name)$V1)
name_list <- name_list[1:23] # No Siegfield 16S RNA

src_fn <- list.files(source_dir) # List the files in the source directory
for (i in 1:length(src_fn)) source(paste0(source_dir,src_fn[i])) # Load the source files

err_estimate <- read.csv(err_fn)

# Linear regression of the reactivities vs the variances (log10-log10 scale) ------------------
# Take the absolute of the reactivities to estimate the variance
err_estimate <- abs(err_estimate)
filt <- !is.na(log10(err_estimate$boot1x.100.1100.react.var)) &
  !is.infinite(log10(err_estimate$boot1x.100.1100.react.var)) &
  !is.na(log10(err_estimate$boot1x.100.1100.reactivity)) &
  !is.infinite(log10(err_estimate$boot1x.100.1100.reactivity))
err_lm <- lm(log10(err_estimate$boot1x.100.1100.react.var[filt])~log10(err_estimate$boot1x.100.1100.reactivity[filt]))

# Loop across selected RNA -------------------------------------------------------
for (i in 1:length(rna_ix)) {
  
  root_out <- paste(root_out_base,name_list[rna_ix[i]],'/',sep='')
  dir.create(root_out)
  root_out_rna_base <- root_out
  
  # Read the shape reactivity file
  shape_data <- read.table(paste(root_data,name_list[rna_ix[i]],'.shape',sep=''))
  
  # Get the associated .ct amd .fa files
  ct_fn <- paste(root_ct,name_list[rna_ix[i]],'.ct',sep='')
  fa_fn <- paste(root_ct,name_list[rna_ix[i]],'.fa',sep='')
  
  #Set NaN reactivities
  nan_reac_filt <- shape_data$V2==-999
  Reactivities <- shape_data$V2
  Reactivities[nan_reac_filt] <- NA
  
  reac <- Reactivities

  # Compute the percentiles -------------------------------------------------
  reac_range <- quantile(reac,seq(0,1,1/n_bin),na.rm = T)
  
  # Loop across error sizes -------------------------------------------------
  for (j in 1:length(err_multiplier)) {
    
    root_out_err_base <- paste(root_out_rna_base,'err_multi_',toString(err_multiplier[j]),'/',sep='')
    dir.create(root_out_err_base)
    
    # Loop across repeats -----------------------------------------------------
    for (k in 1:n_rep) {
      
      # Generate a random white noise (Gaussian centered at 0 with unit variance) -------
      white_noise <- rnorm(length(reac),mean=0,sd = 1)
      
      # Compute the noise to add to the data as a function of the reactivity (using Nathan data for the fitting)
      # Take the absolute of the reactivity to estimate the variance
      # Assign zeros to the smallest value in the dataset
      est_var <- err_lm$coefficients[2]*log10(abs(reac)) + err_lm$coefficients[1]
      est_var[is.infinite(est_var)] <- min(est_var[!is.infinite(est_var)],na.rm=T)
      noise <- white_noise*sqrt(10^(est_var+err_multiplier[j]))
      
      # Write a shape with all reactivities perturbed ---------------------------
      root_out <- paste(root_out_err_base,'percentile_all/',sep='')
      dir.create(root_out)
      
      reac_noise <- reac + noise
      
      # Write the shape file
      if (k < 11) {
        write_shape(fn = paste(root_out,'shape-00',toString(k-1),sep=""), x=reac_noise)
      } else {
        write_shape(fn = paste(root_out,'shape-0',toString(k-1),sep=""), x=reac_noise)
      }
    }
  }
  
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
    k <- 1
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
    k <- 1
    if (k < 11) {
      write_shape(fn = paste(root_out,'shape-00',toString(k-1),sep=""), x=reac_noise)
    } else {
      write_shape(fn = paste(root_out,'shape-0',toString(k-1),sep=""), x=reac_noise)
    }
  }
  
}

save.image(paste(root_out_base,'results.RData',sep=''))


