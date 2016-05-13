## Simulate 3 states SHAPE using a Gaussian kernel density estimation (KDE)

## Remove all active variables from the workspace
rm(list = ls())

# Load libraries
lib <- c('ggplot2','reshape2')
lapply(lib, require, character.only=T)

# Set directory/file paths ------------------------------------------------
source_dir <- 'src/' # location of source scripts
fn_name <- 'data/data_pool_name.txt' # RNA names
root_out_base <- 'output/density_simulation/'
dir.create(root_out_base)
root_ct <- 'data/ct_sequence_file/' # ct data directory
root_data <- 'data/shape_raw/' # shape data directory

# User inputs -------------------------------------------------------------
sim_n <- 100 # number of simulation to generate

# Load sources and data ---------------------------------------------------
# Get SHAPE and ct filename in the folder
name_list <- as.character(read.table(fn_name)$V1)
name_list <- name_list[1:23] # No Siegfield 16S RNA

src_fn <- list.files(source_dir) # List the files in the source directory
for (i in 1:length(src_fn)) source(paste0(source_dir,src_fn[i])) # Load the source files

# Read data ---------------------------------------------------------------
all_pairing <- list()
all_pairing$pairing <- NULL
all_pairing$all <- NULL

for (i in 1:length(name_list)) {
  
  # Read the shape reactivity file
  shape_data <- read.table(paste0(root_data,name_list[i],'.shape'))
  nan_reac_filt <- shape_data$V2==-999
  Reactivities <- shape_data$V2
  Reactivities[nan_reac_filt] <- NA
  reac <- Reactivities
  
  x_pairing <- pairing_reactivities(x=reac,ref_fn = paste(root_ct,name_list[i],'.ct',sep=''),model = 'raw',states = 3)
  
  # concatenate the results
  all_pairing$pairing <- c(all_pairing$pairing,x_pairing$pairing)
  all_pairing$all <- c(all_pairing$all,reac)
  
}

# Separate data by state
dta <- list()
dta$un <- all_pairing$all[all_pairing$pairing==0 & !is.na(all_pairing$all)]
dta$ph <- all_pairing$all[all_pairing$pairing==1 & !is.na(all_pairing$all)]
dta$ps <- all_pairing$all[all_pairing$pairing==2 & !is.na(all_pairing$all)]

# Compute the KDE
dens <- list()
dens$un <- density(dta$un,bw="SJ")
dens$ph <- density(dta$ph,bw="SJ")
dens$ps <- density(dta$ps,bw="SJ")

# Get sample sizes
n <- list()
n$un <- sum(all_pairing$pairing==0) * sim_n
n$ph <- sum(all_pairing$pairing==1) * sim_n
n$ps <- sum(all_pairing$pairing==2) * sim_n

# Simulate new SHAPE values
reac.sim <- list()
reac.sim$un <- rnorm(n$un, sample(dta$un, size = n$un, replace = TRUE), dens$un$bw)
reac.sim$ph <- rnorm(n$ph, sample(dta$ph, size = n$ph, replace = TRUE), dens$ph$bw)
reac.sim$ps <- rnorm(n$ps, sample(dta$ps, size = n$ps, replace = TRUE), dens$ps$bw)

# cntrl plots  
plot(dens$un)
lines(density(reac.sim$un), col = "blue")
plot(dens$ph)
lines(density(reac.sim$ph), col = "blue")
plot(dens$ps)
lines(density(reac.sim$ps), col = "blue")

# Generate new SHAPE files ------------------------------------------------
cnt <- list()
cnt$un <- 0
cnt$ph <- 0
cnt$ps <- 0

for (i in 1:length(name_list)) {
  
  # create output directories
  root_out <- paste0(root_out_base,name_list[i],'/')
  dir.create(root_out)
  root_out <- paste0(root_out,'ternary_model/')
  dir.create(root_out)
  
  shape_data <- read.table(paste0(root_data,name_list[i],'.shape'))
  nan_reac_filt <- shape_data$V2==-999
  Reactivities <- shape_data$V2
  Reactivities[nan_reac_filt] <- NA
  reac <- Reactivities
  
  x_pairing <- pairing_reactivities(x=reac,ref_fn = paste(root_ct,name_list[i],'.ct',sep=''),model = 'raw',states = 3)
  
  for (k in 1:sim_n) {
    if (k<11) {
      fn_out <- paste0('shape-00',toString(k-1))
    } else {
      fn_out <- paste0('shape-0',toString(k-1))
    } 
    
    data.sim <- rep(NA,length(x_pairing$pairing))
    
    n <- sum(x_pairing$pairing==0)
    data.sim[x_pairing$pairing==0] <- reac.sim$un[(cnt$un+1):(cnt$un+n)]
    cnt$un <- cnt$un+n
    
    n <- sum(x_pairing$pairing==1)
    data.sim[x_pairing$pairing==1] <- reac.sim$ph[(cnt$ph+1):(cnt$ph+n)]
    cnt$ph <- cnt$ph+n
    
    n <- sum(x_pairing$pairing==2)
    data.sim[x_pairing$pairing==2] <- reac.sim$ps[(cnt$ps+1):(cnt$ps+n)]
    cnt$ps <- cnt$ps+n

    ## write .shape file
    write_shape(paste0(root_out,fn_out),data.sim)
  }
}
