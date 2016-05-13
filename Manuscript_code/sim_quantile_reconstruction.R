## Simulated SHAPE data reconstruction using a defined quantile ordering

## Remove all active variables from the workspace
rm(list = ls())

# Set directory/file paths ------------------------------------------------
source_dir <- 'src/' # location of source scripts
root_ct <- 'data/ct_sequence_file/' # ct data directory
fn_name <- 'data/data_pool_name.txt' # RNA names
root_data <- 'data/generated_data/sim_data_from_std_fit_distributions/' # shape data directory

# Output directory
root_out <- 'output/sim_data_quantile_reconstruction/'
dir.create(root_out)
root_out <- paste(root_out,'percentile_',toString(n_bin),'/',sep='')
dir.create(root_out)

# Variable inputs --------------------------------------------------------
rna_ix <- 1:23 # Select RNAs based on name_list
n_bin <- 5 # number of bins
model_type <- "ternary" # either binary or ternary
n_rep <- 20 # number of simulation repeats

if (n_bin==3) {
  bin_rank <- c(3,1,2)
} else if (n_bin==5) {
  bin_rank <- c(5,1,4,2,3)
} else if (n_bin==10) {
  bin_rank <- c(10,9,8,1,2,3,4,7,6,5)}

# Finish to initalize the output directory
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
    
    root_out_missing_base2 <- paste(root_out_rna_base,'only_percentile/',sep='')
    dir.create(root_out_missing_base2)
    
    reac_ix <- matrix(F,n_bin,length(reac))
    for (t in 1:n_bin) {
      reac_ix[t,(reac>=reac_range[t]) & (reac<=reac_range[t+1]) & !is.na(reac)] <- T
    }
    
    ## Writing only the selected percentile
    for (l in 1:(length(reac_range)-1)) {
      
      root_out <- paste(root_out_missing_base2,'percentile_',toString(l),'/',sep='')
      dir.create(root_out)
      
      reac_ixt <- NULL
      
      for (t in 1:l) {
        reac_ixt <- c(reac_ixt,which(reac_ix[bin_rank[t],]))
      }

      reac_noise <- rep(NA,length(reac))
      reac_noise[reac_ixt] <- reac[reac_ixt]
      
      # Write the shape file
      if (k < 11) {
        write_shape(fn = paste(root_out,'shape-00',toString(k-1),sep=""), x=reac_noise)
      } else {
        write_shape(fn = paste(root_out,'shape-0',toString(k-1),sep=""), x=reac_noise)
      }
    }
  }
}



