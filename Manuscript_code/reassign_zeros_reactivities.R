# Random reassignement of negative values and zeros using negative value distribution --------
# Author: Mirko Ledda | Date: 06/08/2015 | v1.0

rm(list = ls())

# Load libraries ----------------------------------------------------------
library(PearsonDS)

# Set directory/file paths ------------------------------------------------
root_data <- 'data/shape_raw/' # shape data directory
root_out <- 'output/' # output directory
source_dir <- 'src/' # source directory
fn_name <- 'data/data_pool_name.txt' # RNA names

# Load sources and data ---------------------------------------------------
# Get SHAPE and ct filename in the folder
name_list <- as.character(read.table(fn_name)$V1)
name_list <- name_list[1:23] # No Siegfield 16S RNA

src_fn <- list.files(source_dir) # List the files in the source directory
for (i in 1:length(src_fn)) source(paste0(source_dir,src_fn[i])) # Load the source files

# create output directories -----------------------------------------------
root_out <- paste0(root_out,'reassign_zeros_reactivities/')
root_out_neg <- paste0(root_out,'random_negative/')
root_out_abs <- paste0(root_out,'random_absolute/')
root_out_abs_ln <- paste0(root_out,'random_absolute_ln/')
root_out_abs_boxcox <- paste0(root_out,'random_absolute_boxcox/')

dir.create(root_out)
dir.create(root_out_neg)
dir.create(root_out_abs)
dir.create(root_out_abs_ln)
dir.create(root_out_abs_boxcox)

# Get the negative values to construct the probability distribution -------
negative_reac <- NULL
for (i in 1:length(name_list)) {
  
  # Read the shape reactivity file
  shape_data <- read.table(paste0(root_data,name_list[i],'.shape'))
  nan_reac_filt <- shape_data$V2==-999
  Reactivities <- shape_data$V2
  Reactivities[nan_reac_filt] <- NA
  Reactivities <- Reactivities[!is.na(Reactivities)]
  if (sum(Reactivities<0,na.rm=T) > 0) negative_reac <- c(negative_reac,Reactivities[Reactivities<0])
}

# Compute distribution moments
negative_reac_moments <- empMoments(log(abs(negative_reac)))

# Reassign 0 and negative Reactivities to randomly generated values
for (i in 1:length(name_list)) {
  
  # Read the shape reactivity file
  shape_data <- read.table(paste0(root_data,name_list[i],'.shape'))
  nan_reac_filt <- shape_data$V2==-999
  Reactivities <- shape_data$V2
  Reactivities[nan_reac_filt] <- NA
  
  ix <- ((Reactivities<0) | (Reactivities==0)) & (!is.na(Reactivities))
  ndraw <- sum(ix)
  
  if (ndraw > 0) {
    # generate random number
    rand_gen <- -exp(rpearson(ndraw,moments=negative_reac_moments))
    # assign those random generated data
    Reactivities[ix] <- rand_gen
  }
  
  Reactivities[nan_reac_filt] <- NA
  
  ## write .shape files
  # with negative values
  write_shape(paste0(root_out_neg,name_list[i]),Reactivities)
  # absolute values
  write_shape(paste0(root_out_abs,name_list[i]),abs(Reactivities))
  # log-transformed absolute values
  write_shape(paste0(root_out_abs_ln,name_list[i]),log(abs(Reactivities)))
  # boxcox-transformed absolute values (lambda=0.1)
  write_shape(paste0(root_out_abs_boxcox,name_list[i]),boxcox.transform(abs(Reactivities),lambda=0.1))
}
