# Simulate data for a mock probe -------------------------------------------
# Scenario: Variance across nucleotides is not equal
# -> results in paired bases having a larger variance compared to unpaired
# Distributions: Normal distributions

## Remove all active variables from the workspace
rm(list = ls())

# Load libraries ----------------------------------------------------------
library(evd)

# Set directory/file paths ------------------------------------------------
shape_out <- "output/sim_probe_sc2/" # Output folder
dir.create(shape_out)
ctp <- "data/ct_sequence_file/" # Path to ct files
rnafn <- 'data/data_pool_name.txt' # RNA file names
source_dir <- 'src/' # location of source scripts

# Set density distributions parameters ----------------------------------------------------
nsim <- 10 # number of simulations

up <- list() # unpaired
ph <- list() # helix-end
ps <- list() # stacked
nuc.sd <- list()

up$mean <- 0.5
ph$mean <- 0.25
ps$mean <- 0

nuc.sd$A <- 0.1
nuc.sd$U <- 0.1
nuc.sd$G <- 1
nuc.sd$C <- 1

# Load sources ---------------------------------------------------
src_fn <- list.files(source_dir) # List the files in the source directory
for (i in 1:length(src_fn)) source(paste0(source_dir,src_fn[i])) # Load the source files

# Read RNA names ----------------------------------------------------------
name_list <- as.character(read.table(rnafn)$V1)
name_list <- name_list[1:23] # No Siegfield 16S RNA

# Create output directory
dir.create(shape_out)

# Initialize test result array 
fs <- array(NA,dim=c(length(iv_m),length(iv_b),length(name_list)))

# Loop across RNAs and simulations -------------------------------------
for (rnai in 1:length(name_list)) {
  # Read ct file ---------------------------------------------------------
  ct <- read.table(paste0(ctp,name_list[rnai],".ct"),skip=1)
  seq <- ct$V2
  
  # Get pairing
  ref_pair <- ct$V5
  ref_pair <- get_pair_3s(ref_pair)
  
  # Create output directory
  shape_out_rna <- paste0(shape_out,name_list[rnai],"/")
  dir.create(shape_out_rna)
  
  # Loop across simulations
  for (k in 1:nsim) {
    
    # Data simulation ---------------------------------------------------------
    sim_reac <- rep(NA,length(ref_pair))
    filtPair <- ref_pair!=0
    filtPairH <- ref_pair==1
    filtPairS <- ref_pair==2
    filtA <- seq=="A"
    filtU <- seq=="U"
    filtG <- seq=="G"
    filtC <- seq=="C"
    filtN <- seq=="N"
    
    # Unpaired
    mfilt <- !filtPair
    tfilt <- mfilt & filtA
    sim_reac[tfilt] <- rnorm(sum(tfilt),up$mean,nuc.sd$A)
    tfilt <- mfilt & filtU
    sim_reac[tfilt] <- rnorm(sum(tfilt),up$mean,nuc.sd$U)
    tfilt <- mfilt & filtG
    sim_reac[tfilt] <- rnorm(sum(tfilt),up$mean,nuc.sd$G)
    tfilt <- mfilt & filtC
    sim_reac[tfilt] <- rnorm(sum(tfilt),up$mean,nuc.sd$C)
    tfilt <- mfilt & filtN
    sim_reac[tfilt] <- rnorm(sum(tfilt),up$mean,mean(unlist(nuc.sd)))
    
    # Helix-end
    mfilt <- filtPairH
    tfilt <- mfilt & filtA
    sim_reac[tfilt] <- rnorm(sum(tfilt),ph$mean,nuc.sd$A)
    tfilt <- mfilt & filtU
    sim_reac[tfilt] <- rnorm(sum(tfilt),ph$mean,nuc.sd$U)
    tfilt <- mfilt & filtG
    sim_reac[tfilt] <- rnorm(sum(tfilt),ph$mean,nuc.sd$G)
    tfilt <- mfilt & filtC
    sim_reac[tfilt] <- rnorm(sum(tfilt),ph$mean,nuc.sd$C)
    tfilt <- mfilt & filtN
    sim_reac[tfilt] <- rnorm(sum(tfilt),ph$mean,mean(unlist(nuc.sd)))
    
    # Stacked
    mfilt <- filtPairS
    tfilt <- mfilt & filtA
    sim_reac[tfilt] <- rnorm(sum(tfilt),ps$mean,nuc.sd$A)
    tfilt <- mfilt & filtU
    sim_reac[tfilt] <- rnorm(sum(tfilt),ps$mean,nuc.sd$U)
    tfilt <- mfilt & filtG
    sim_reac[tfilt] <- rnorm(sum(tfilt),ps$mean,nuc.sd$G)
    tfilt <- mfilt & filtC
    sim_reac[tfilt] <- rnorm(sum(tfilt),ps$mean,nuc.sd$C)
    tfilt <- mfilt & filtN
    sim_reac[tfilt] <- rnorm(sum(tfilt),ps$mean,mean(unlist(nuc.sd)))
    
    # Write shape file
    x <- sim_reac
    n <- length(x)
    if (k<11) {
      fn_out <- paste0('shape-00',toString(k-1))
    } else {
      fn_out <- paste0('shape-0',toString(k-1))
    }
    fn_reac <- paste0(shape_out_rna,fn_out,'.shape')
    rdf <- data.frame("n"=seq(1,n),'SHAPE'=x)
    rdf[is.na.data.frame(rdf)] <- -999
    write.table(rdf,file=fn_reac,sep=" ",row.names = F, col.names = F)
    
  }
}

