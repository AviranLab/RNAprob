# Bayesian posterior probability of base pairing --------------------------

# Author: Mirko Ledda | Date: 11/11/2015 | v1
#
# Uses the function `kde` from the package `ks` for kernel estimation

## Remove all active variables from the workspace
rm(list = ls())

# Set directory/file paths ------------------------------------------------
source_dir <- 'src/'
root_ct <- 'data/ct_sequence_file/' # root directory to the ct files
fn_name <- 'data/data_pool_name.txt' # RNA names

# User inputs -------------------------------------------------------------
dims <- c(1) # note that 0 correspond to the density estimation based on all reactivities
bayes_type = 'n5'
n_bins = 100
h=0.1 # either 'ls': least-square, 'lsc': least-square with covariance correction or a numerical value
data_transform <- 'raw' # either 'raw' or 'boxcox'

# Set other directories ---------------------------------------------------
if (data_transform == 'raw') {
  # raw data
  limits = c(-0.108895,2.220946)
} else if (data_transform == 'boxcox') {
  # transformed data
  limits <- c(-3.93168,0.8347476)
}

## set the shape data directory
if (data_transform == 'raw') {
  root_data <- 'data/shape_raw/'
} else if (data_transform == 'boxcox') {
  root_data <- 'data/generated_data/reassign_zeros_reactivities/random_absolute_boxcox/'
}

# output directory
root_out <- paste0('output/bayesian_prob_pairing_',data_transform,'_h_',toString(h),'_bin_',toString(n_bins),'/')

## Create needed directories
dir.create(root_out)
if (any(dims==0)) dir.create(paste0(root_out,'reactivity_density/'))
if (any(dims==1)) dir.create(paste0(root_out,'posterior_1D/'))
if (any(dims==2)) dir.create(paste0(root_out,'posterior_2D_',bayes_type,'/'))
if (any(dims==3)) dir.create(paste0(root_out,'posterior_3D/'))

# Automated Script --------------------------------------------------------
# Load source functions
src_fn <- list.files(source_dir)
for (i in 1:length(src_fn)) source(paste0(source_dir,src_fn[i]))

# Load sources and data ---------------------------------------------------
# Get SHAPE and ct filename in the folder
name_list <- as.character(read.table(fn_name)$V1)
name_list <- name_list[1:23] # exclude the ecoli_16S_siegfried data

# set list for concatenated data
all_pairing <- list()
all_pairing$un <- NULL; all_pairing$pe <- NULL; all_pairing$ps <- NULL; all_pairing$pairing <- NULL; all_pairing$all <- NULL

# Loop across datasets ----------------------------------------------------
for (i in 1:length(name_list)) {
  
  # Read the shape reactivity file
  shape_data <- read.table(paste0(root_data,name_list[i],'.shape'))
  nan_reac_filt <- shape_data$V2==-999
  Reactivities <- shape_data$V2
  Reactivities[nan_reac_filt] <- NA
  reac <- Reactivities
  
  # Pair the nucleotide according to 3 states, unpaired/stacked/helix-end
  x_pairing <- pairing_reactivities(x=reac,ref_fn = paste0(root_ct,name_list[i],'.ct'),model = 'raw',states = 3)
  
  # concatenate the results
  # exclude the ecoli_16S_siegfried data
  all_pairing$unpaired <- c(all_pairing$unpaired,x_pairing$unpaired)
  all_pairing$paired_end <- c(all_pairing$paired_end,x_pairing$paired_end)
  all_pairing$paired_stack <- c(all_pairing$paired_stack,x_pairing$paired_stack)
  all_pairing$pairing <- c(all_pairing$pairing,x_pairing$pairing)
  all_pairing$reac <- c(all_pairing$reac,reac)  
  all_pairing$flag <- c(all_pairing$flag,1,rep(0,length(reac)-2),1)
  
  # One-dimensional kernel estimation and posterior Bayesian probabilities---------------------------
  if (any(dims==1)) {
    ker <- kernel_3s(x=reac,xp = x_pairing,dims=1,nbin=n_bins,h=h,limits=limits)
    write.kernel.3s(ker=ker,dims = 1,fn = paste0(root_out,'posterior_1D/',name_list[i]))
    if (!any(c(is.na(ker$u$estimate),is.na(ker$e$estimate),is.na(ker$s$estimate)))) {
      post <- posterior.3s(ker,dims=1)
      if (!any(c(is.na(post$unpaired),is.na(post$paired_end),is.na(post$paired_stack)))) {
        plot1D_3s(post = post,fn=paste0(root_out,'posterior_1D/',name_list[i]))
        write.posterior.3s(post=post,dims = 1,fn = paste0(root_out,'posterior_1D/',name_list[i]))
      }
    }
  }
  # Two-dimensional kernel estimation and posterior Bayesian probabilities---------------------------
  if (any(dims==2)) {
    ker <- kernel_3s(x=reac,xp = x_pairing,dims=2,type = bayes_type,nbin=n_bins,h=h,limits=limits)
    write.kernel.3s(ker=ker,dims = 2,fn = paste0(root_out,'posterior_2D_',bayes_type,'/',name_list[i]))
    if (!any(c(is.na(ker$u$estimate),is.na(ker$e$estimate),is.na(ker$s$estimate)))) {
      post <- posterior.3s(ker,dims=2)
      if (!any(c(is.na(post$unpaired),is.na(post$paired_end),is.na(post$paired_stack)))) {
        plot2D_3s(post = post,fn=paste0(root_out,'posterior_2D_',bayes_type,'/',name_list[i]))
        write.posterior.3s(post=post,dims = 2,fn = paste0(root_out,'posterior_2D_',bayes_type,'/',name_list[i]))
      }
    }
  }
  # Three-dimensional kernel estimation and posterior Bayesian probabilities---------------------------
  if (any(dims==3)) {
    ker <- kernel_3s(x=reac,xp = x_pairing,dims=3,type = bayes_type,nbin=n_bins,h=h,limits=limits)
    write.kernel.3s(ker=ker,dims = 3,fn = paste0(root_out,'posterior_3D/',name_list[i]))
    if (!any(c(is.na(ker$u$estimate),is.na(ker$e$estimate),is.na(ker$s$estimate)))) {
      post <- posterior.3s(ker,dims=3)
      if (!any(c(is.na(post$unpaired),is.na(post$paired_end),is.na(post$paired_stack)))) {
        plot3D_3s(post = post,fn=paste0(root_out,'posterior_3D/',name_list[i]),data_transform=data_transform)
        write.posterior.3s(post=post,dims = 3,fn = paste0(root_out,'posterior_3D/',name_list[i]))
      }
    }
  }
  # Compute the density of the reactivities
  if (any(dims==0)) {
    dens <- density(reac[!is.na(reac)],from = limits[1],to = limits[2],n = n_bins,kernel="gaussian",adjust=1)
    result <- cbind(dens$x,dens$y)
    colnames(result) <- c('reactivities','density')
    write.table(x = result,row.names = F,col.names = T,quote = F,sep=' ',
                file = paste0(root_out,'reactivity_density/',name_list[i],'.reac.density'))
    png(filename = paste0(root_out,'reactivity_density/',name_list[i],'.png'),width = 800,height = 800)
    par(cex=1.3)
    p1 <- qplot(x=dens$x,y=dens$y) + geom_line() + xlab('Reactivities') + ylab('Density')
    grid.arrange(p1, nrow=1, ncol=1)
    dev.off()
  }
}

# One-dimensional kernel estimation and posterior Bayesian probabilities for all data
if (any(dims==1)) {
  ker <- kernel_3s(x=all_pairing$reac,xp = all_pairing,dims=1,nbin=n_bins,flag=1,h=h,limits=limits)
  write.kernel.3s(ker=ker,dims = 1,fn = paste0(root_out,'posterior_1D/all'))
  post <- posterior.3s(D=ker,dims=1)
  plot1D_3s(post = post,fn=paste0(root_out,'posterior_1D/all'))
  write.posterior.3s(post=post,dims = 1,fn = paste0(root_out,'posterior_1D/all'))
}
# Two-dimensional kernel estimation and posterior Bayesian probabilities for all data
if (any(dims==2)) {
  ker <- kernel_3s(x=all_pairing$reac,xp = all_pairing,dims=2,type = bayes_type,nbin=n_bins,flag=1,h=h,limits=limits)
  write.kernel.3s(ker=ker,dims = 2,fn = paste0(root_out,'posterior_2D_',bayes_type,'/all'))
  post <- posterior.3s(D=ker,dims=2)
  plot2D_3s(post = post,fn=paste0(root_out,'posterior_2D_',bayes_type,'/all'))
  write.posterior.3s(post=post,dims = 2,fn = paste0(root_out,'posterior_2D_',bayes_type,'/all'))
}
# Three-dimensional kernel estimation and posterior Bayesian probabilities for all data
if (any(dims==3)) {
  ker <- kernel_3s(x=all_pairing$reac,xp = all_pairing,dims=3,nbin=n_bins,flag=1,h=h,limits=limits)
  write.kernel.3s(ker=ker,dims = 3,fn = paste0(root_out,'posterior_3D/all'))
  post <- posterior.3s(D=ker,dims=3)
  plot3D_3s(post = post,fn=paste0(root_out,'posterior_3D/all'),data_transform=data_transform)
  write.posterior.3s(post=post,dims = 3,fn = paste0(root_out,'posterior_3D/all'))
}
# Compute the density of the reactivities
if (any(dims==0)) {
  dens <- density(all_pairing$reac[!is.na(all_pairing$reac)],from = limits[1],to = limits[2],n = n_bins,kernel="gaussian",adjust=1)
  result <- cbind(dens$x,dens$y)
  colnames(result) <- c('reactivities','density')
  write.table(x = result,row.names = F,col.names = T,quote = F,sep=' ',
              file = paste0(root_out,'reactivity_density/all.reac.density'))
  png(filename = paste0(root_out,'reactivity_density/all.png'),width = 800,height = 800)
  par(cex=1.3)
  qplot(x=dens$x,y=dens$y) + geom_line() + xlab('Reactivities') + ylab('Density')
  grid.arrange(p1, nrow=1, ncol=1)
  dev.off()
}

