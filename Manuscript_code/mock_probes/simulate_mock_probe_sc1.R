# Simulate data for a mock probe -------------------------------------------
# Scenario: Helix-end pairs react like Unpaired bases
# Distributions: SHAPE real distributions

## Remove all active variables from the workspace
rm(list = ls())

# Load libraries ----------------------------------------------------------
library(evd)

# Set directory/file paths ------------------------------------------------
shape_out <- "output/sim_probe_sc1/" # Output folder
dir.create(shape_out)
ctp <- "data/ct_sequence_file/" # Path to ct files
rnafn <- 'data/data_pool_name.txt' # RNA file names
source_dir <- 'src/' # location of source scripts

# Set parameters ----------------------------------------------------------
nsim <- 10 # number of simulations

# Set density distributions parameters
up <- list() # unpaired
ph <- list() # helix-end
ps <- list() # stacked

up$lambda <- 1.468            
ph$xi <- 0.821
ph$sigma <- 0.114
ph$mu <- 0.090
ps$xi <- 0.763
ps$sigma <- 0.049
ps$mu <- 0.040

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
    
    sim_reac[!filtPair] <- rexp(sum(!filtPair),up$lambda)
    sim_reac[filtPairH] <- rexp(sum(filtPairH),up$lambda)
    sim_reac[filtPairS] <- rgev(sum(filtPairS),loc=ps$mu,scale=ps$sigma,shape=ps$xi)
    
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

