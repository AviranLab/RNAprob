# Fit probe simulation data and build binned files for RNAprob ----------------

## Remove all active variables from the workspace
rm(list = ls())

# Load libraries ----------------------------------------------------------
library(evd)

# Inputs ------------------------------------------------------------------
dir_out <- "output/" # Output folder
source_dir <- 'src/' # source directory
ctp <- "data/ct_sequence_file/" # Path to ct files
train_fn <- 'data/data_pool_list_train.txt' # RNA file names
root_data <- 'data/generated_data/mock_probe_simulation/'
nsim <- 10 # number of simulations
libr <- "evd" # Required libraries
bin_size <- 0.1 # Size of the bin for the empirical probability distributions

# Load sources and data ---------------------------------------------------
src_fn <- list.files(source_dir) # List the files in the source directory
for (i in 1:length(src_fn)) source(paste0(source_dir,src_fn[i])) # Load the source files

# Load training file
name_list <- as.character(read.table(train_fn)$V1)

# Scenario 1 --------------------------------------------------------------
sc_dir <- 'sim_probe_sc1/'
out_dir_sc <- paste0(dir_out,sc_dir)
dir.create(out_dir_sc)
root_data_sc <- paste0(root_data,sc_dir)
param_names <- c("lambda","mu","sigma","xi")

# Loop across simulations
for (k in 1:nsim) {
  ref_pair_all <- NULL
  dta_all <- NULL
  fitRes <- matrix(NA,3,length(param_names))
  colnames(fitRes) <- param_names
  rownames(fitRes) <- c("Unpaired","Helix-end","Stacked")
  
  # Loop across RNA
  for (rnai in 1:length(name_list)) {
    # read ct and sequence
    ct <- read.table(paste0(ctp,name_list[rnai],".ct"),skip=1)
    seq <- ct$V2
    
    # Get pairing
    ref_pair <- ct$V5
    ref_pair_all <- c(ref_pair_all,get_pair_3s(ref_pair))
    
    # Read the shape reactivity file
    if (k<11) {
      dta <- read.table(paste0(root_data_sc,name_list[rnai],'/shape-00',toString(k-1),'.shape',sep=''))$V2
    } else {
      dta <- read.table(paste(root_data_sc,name_list[rnai],'/shape-0',toString(k-1),'.shape',sep=''))$V2
    }
    
    dta_all <- c(dta_all,dta)
  }
  
  # Fit distributions
  upFit <- CFitExp(dta_all[ref_pair_all==0])
  phFit <- CFitExp(dta_all[ref_pair_all==1])
  psFit <- CFitGEV(dta_all[ref_pair_all==2])
  pFit <- CFitGEV(dta_all[ref_pair_all>0])
  
  # Store fitting results
  fitRes[1,1] <- upFit$lambda
  fitRes[2,1] <- phFit$lambda
  fitRes[3,2:4] <- c(psFit$mu,psFit$sigma,psFit$xi)
  
  # Write fitting and binned results to a file
  if (k<11) {
    write.table(fitRes,file = paste0(out_dir_sc,"fit-00",toString(k-1),'.txt'),quote=F)
    make_RNAprob_param_file_hist(x=dta_all,pairing = ref_pair_all,bin_size = bin_size,fn = paste0(out_dir_sc,"param-00",toString(k-1),'.txt'))
  } else {
    write.table(fitRes,file = paste0(out_dir_sc,"fit-0",toString(k-1),'.txt'),quote=F)
    make_RNAprob_param_file_hist(x=dta_all,pairing = ref_pair_all,bin_size = bin_size,fn = paste0(out_dir_sc,"param-0",toString(k-1),'.txt'))
  }
}


# Scenario 2 --------------------------------------------------------------
sc_dir <- 'sim_probe_sc2/'
out_dir_sc <- paste0(dir_out,sc_dir)
dir.create(out_dir_sc)
root_data_sc <- paste0(root_data,sc_dir)
param_names <- c("mean","sd")

# Loop across simulations
for (k in 1:nsim) {
  ref_pair_all <- NULL
  dta_all <- NULL
  fitRes <- matrix(NA,3,length(param_names))
  colnames(fitRes) <- param_names
  rownames(fitRes) <- c("Unpaired","Helix-end","Stacked")
  
  # Loop across RNA
  for (rnai in 1:length(name_list)) {
    # read ct and sequence
    ct <- read.table(paste0(ctp,name_list[rnai],".ct"),skip=1)
    seq <- ct$V2
    
    # Get pairing
    ref_pair <- ct$V5
    ref_pair_all <- c(ref_pair_all,get_pair_3s(ref_pair))
    
    # Read the shape reactivity file
    if (k<11) {
      dta <- read.table(paste0(root_data_sc,name_list[rnai],'/shape-00',toString(k-1),'.shape',sep=''))$V2
    } else {
      dta <- read.table(paste(root_data_sc,name_list[rnai],'/shape-0',toString(k-1),'.shape',sep=''))$V2
    }
    
    dta_all <- c(dta_all,dta)
  }
  
  # Fit distributions
  upFit <- CFitNormal(dta_all[ref_pair_all==0])
  phFit <- CFitNormal(dta_all[ref_pair_all==1])
  psFit <- CFitNormal(dta_all[ref_pair_all==2])
  pFit <- CFitNormal(dta_all[ref_pair_all>0])
  
  # Store fitting results
  fitRes[1,1:2] <- c(upFit$mean,upFit$sd)
  fitRes[2,1:2] <- c(phFit$mean,phFit$sd)
  fitRes[3,1:2] <- c(psFit$mean,psFit$sd)
  
  # Write fitting and binned results to a file
  if (k<11) {
    write.table(fitRes,file = paste0(out_dir_sc,"fit-00",toString(k-1),'.txt'),quote=F)
    make_RNAprob_param_file_hist(x=dta_all,pairing = ref_pair_all,bin_size = bin_size,fn = paste0(out_dir_sc,"param-00",toString(k-1),'.txt'))
  } else {
    write.table(fitRes,file = paste0(out_dir_sc,"fit-0",toString(k-1),'.txt'),quote=F)
    make_RNAprob_param_file_hist(x=dta_all,pairing = ref_pair_all,bin_size = bin_size,fn = paste0(out_dir_sc,"param-0",toString(k-1),'.txt'))
  }
}

