## Analyse the results from the scheme robustness test when rejecting or including selected reactivities

# Author: Mirko Ledda
# Date: 09/15/2015
# Version: v1.0

## Remove all active variables from the workspace
rm(list = ls())

# Load libraries
lib <- c('ggplot2','reshape2','grid','gridExtra')
lapply(lib, require, character.only=T)

# Set directory/file paths ------------------------------------------------
source_dir <- 'src/' # location of source scripts
fn_name <- 'data/data_pool_name.txt' # RNA names
directory_out <- 'output/'
root_data <- 'data/shape_raw/'
root_ct <- 'data/ct_sequence_file/' # path to the ct files
fn_name_table <- 'data/data_pool_name_table.txt'

fd_ref <- 'data/folding_results/ref_folding.RData'
fd_reco_real <- 'data/folding_results/crossvalidation_real_folding.RData'
fd_reco_sim_out <- 'data/folding_results/crossvalidation_sim_leaveoneout_folding.RData'
fd_reco_sim_in <- 'data/folding_results/crossvalidation_sim_leaveonein_folding.RData'
fd_sim_ternary_fn <- 'data/folding_results/sim_shape_ternary_folding.RData'

# User inputs -------------------------------------------------------------
# Plot inputs
ati <- 12 # axis title font size
ate <- 12 # axis text font size

data_transform <- 'raw' # either 'raw' or 'ln' (Note: does not work with Box-Cox transformed data)
rna_ix <- 1:23 # Select RNAs based on name_list (1:23 would perform the analysis over ALL RNAs)
err_multiplier <- c(0,1,2,3,4,5,6) # Multiplier for the variance of the normal distribution of random errors
n_bin <- 5 # number of bins
n_rep <- 2 # number of simulation repeats

# Automated script starts here --------------------------------------------
# Load sources
src_fn <- list.files(source_dir)
for (i in 1:length(src_fn)) source(paste(source_dir,src_fn[i],sep=""))
root_out_base <- directory_out
load(fd_ref)

# Set matrices ------------------------------------------------------------
reac_range_store <- matrix(NA,n_bin+1,length(rna_ix))
noise_store <- array(NA,dim=c(10000,length(err_multiplier),n_rep,length(rna_ix)))
sim_data_store <- array(NA,dim=c(10000,length(err_multiplier),n_rep,length(rna_ix),n_bin))
rna_length <- rep(NA,length(rna_ix))
snr_store <- array(NA,dim=c(length(err_multiplier),n_rep,length(rna_ix),n_bin+1))
rank_store <- array(NA,dim=c(2,length(rna_ix),n_bin,5))
reac_store <- matrix(NA,10000,length(rna_ix))
pairw_ttest <- array(NA,dim=c(n_bin,5,2))

# Get SHAPE and ct filename in the RNA repository folder
name_list <- as.character(read.table(fn_name)$V1)
name_list <- name_list[1:23]

# Read result files -------------------------------------------------------
load(fd_reco_real)

# Loop across selected RNA -------------------------------------------------------
for (i in 1:length(rna_ix)) {
  
  # Read the shape reactivity file
  shape_data <- read.table(paste(root_data,name_list[rna_ix[i]],'.shape',sep=''))
  
  # Get the associated .ct amd .fa files
  ct_fn <- paste(root_ct,name_list[rna_ix[i]],'.ct',sep='')
  fa_fn <- paste(root_ct,name_list[rna_ix[i]],'.fa',sep='')
  
  #Set NaN reactivities
  nan_reac_filt <- shape_data$V2==-999
  Reactivities <- shape_data$V2
  Reactivities[nan_reac_filt] <- NA
  
  if (data_transform == 'raw') {
    reac <- Reactivities
  } else if (data_transform == 'ln') {
    reac <- exp(Reactivities)
  }
  
  # Compute the percentiles -------------------------------------------------
  reac_range <- quantile(reac,seq(0,1,1/n_bin),na.rm = T)
  reac_range_store[,i] <- reac_range
}

# Write the range file ----------------------------------------------------
name_list_table <- read.delim(fn_name_table,sep=';',header=F)
header_label <- "RNA"
for (l in 0:n_bin) {
  header_label <- c(header_label,toString(round((l) * 100/n_bin)))
}
reac_range_store <- t(reac_range_store)
reac_range_store <- as.data.frame(cbind(as.vector(name_list_table$V1[rna_ix]),round(reac_range_store,digits=2)))
colnames(reac_range_store) <- header_label
write.table(file = paste(root_out_base,'range.txt',sep = ''),x = reac_range_store,quote = F,sep = '\t',row.names = F)

# Plot the bars ------------------------------------------------
# Make a grey scale palette
glim <- c(80,240)/255
giv  <- (glim[2] - glim[1])/n_bin
gpalette <- NULL
for (j in 1:n_bin) {
  gpalette <- c(gpalette,rgb(glim[1]+giv*(j-1),glim[1]+giv*(j-1),glim[1]+giv*(j-1)))
}
ref_ix <- 25

# Averaged results --------------------------------------------------------
for (itype in 1:2) {
  data_label <- NULL
  for (l in 1:n_bin) {
    data_label <- c(data_label,paste(toString(round((l-1) * 100/n_bin)),'-',toString(round((l) * 100/n_bin)),'%   ',sep=''))
  }
  
  if (itype == 1) {
    data_sel <- 1:n_bin
  } else {
    data_sel <- (n_bin+1):(n_bin*2)
  }
  
  df <- matrix(NA,n_bin,5)
  df_sd <- matrix(NA,n_bin,5)
  
  df[,1] <- deigan[data_sel,ref_ix]
  df[,2] <- eddy_2c[data_sel,ref_ix]
  df[,4] <- eddy_2cSukosd[data_sel,ref_ix]
  df[,3] <- eddy_3c[data_sel,ref_ix]
  df[,5] <- eddy_3cSukosd[data_sel,ref_ix]
  
  dff <- data.frame(
    Scheme = factor(c(rep("Lin",n_bin),rep("P2",n_bin),rep("P3",n_bin),rep("P2s",n_bin),rep("P3s",n_bin))),
    resp = c(df),
    Quintiles = factor(rep(data_label,5))
  )
  
  dodge <- position_dodge(width=0.9)
  
  if (itype==2) {
    plot_t <- "Leave-one-in"
    plot1 <- ggplot(dff, aes(fill=Quintiles, y=resp, x=Scheme)) + 
      geom_bar(position="dodge", stat="identity",color="black") + 
      geom_hline(aes(yintercept=std_score[25,3]),size=0.9,linetype=2) +
      geom_segment(aes(x = 0.6, y = std_score[25,6], xend = 1.4, yend = std_score[25,6])) +
      geom_segment(aes(x = 1.6, y = std_score[25,9], xend = 2.4, yend = std_score[25,9])) +
      geom_segment(aes(x = 2.6, y = std_score[25,12], xend = 3.4, yend = std_score[25,12])) +
      geom_segment(aes(x = 3.6, y = std_score[25,15], xend = 4.4, yend = std_score[25,15])) +
      geom_segment(aes(x = 4.6, y = std_score[25,18], xend = 5.4, yend = std_score[25,18])) +
      ylab("") + xlab("") + theme_bw() + theme(legend.position="bottom") +
      ylim(0,100) + 
      scale_fill_manual(name="Quintiles   ",values=gpalette) +
      theme(axis.text=element_text(size=ate),
            axis.title=element_text(size=ati,face="bold"),
            legend.text=element_text(size=12),
            legend.title=element_text(size=14,face="bold"),
            plot.title=element_text(size=34,face="bold"))
  } else {
    plot_t <- "Leave-one-out"
    plot2 <- ggplot(dff, aes(fill=Quintiles, y=resp, x=Scheme)) + 
      geom_bar(position="dodge", stat="identity",color="black") + 
      geom_hline(aes(yintercept=std_score[25,3]),size=0.9,linetype=2) +
      geom_segment(aes(x = 0.6, y = std_score[25,6], xend = 1.4, yend = std_score[25,6])) +
      geom_segment(aes(x = 1.6, y = std_score[25,9], xend = 2.4, yend = std_score[25,9])) +
      geom_segment(aes(x = 2.6, y = std_score[25,12], xend = 3.4, yend = std_score[25,12])) +
      geom_segment(aes(x = 3.6, y = std_score[25,15], xend = 4.4, yend = std_score[25,15])) +
      geom_segment(aes(x = 4.6, y = std_score[25,18], xend = 5.4, yend = std_score[25,18])) +
      ylab("") + xlab("") + theme_bw() +
      ylim(0,100) + 
      scale_fill_manual(values=gpalette) +
      theme(axis.text=element_text(size=ati),
            axis.title=element_text(size=ate,face="bold"),
            legend.text=element_text(size=20,face="bold"),
            legend.title=element_text(size=26,face="bold"),
            plot.title=element_text(size=34,face="bold"))
  }
}


# SIMULATION DATA ---------------------------------------------------------
n_rep <- 20 # number of simulation repeats
itype <- c("missing","only")

load(fd_sim_ternary_fn)

# Loop across analysis type --------------------------------------------
for (ty in itype) {
  if (ty=="missing") {
    load(fd_reco_sim_out)
  } else {
    load(fd_reco_sim_in)
  }
  
  # Plot the averaged results -----------------------------------------------
  # Loop across bins
  for (j in 1:n_bin) {
    # Loop across repeats
    for (k in 1:n_rep) {
      deigan_mcc_avg[j,k] <- mcc.calc(TP = sum(deigan[,j,k,1],na.rm=T),FP = sum(deigan[,j,k,2],na.rm=T),TN = sum(deigan[,j,k,3],na.rm=T),FN = sum(deigan[,j,k,4],na.rm=T))
      eddy_2c_mcc_avg[j,k] <- mcc.calc(TP = sum(eddy_2c[,j,k,1],na.rm=T),FP = sum(eddy_2c[,j,k,2],na.rm=T),TN = sum(eddy_2c[,j,k,3],na.rm=T),FN = sum(eddy_2c[,j,k,4],na.rm=T))
      eddy_2cSukosd_mcc_avg[j,k] <- mcc.calc(TP = sum(eddy_2cSukosd[,j,k,1],na.rm=T),FP = sum(eddy_2cSukosd[,j,k,2],na.rm=T),TN = sum(eddy_2cSukosd[,j,k,3],na.rm=T),FN = sum(eddy_2cSukosd[,j,k,4],na.rm=T))
      eddy_3c_mcc_avg[j,k] <- mcc.calc(TP = sum(eddy_3c[,j,k,1],na.rm=T),FP = sum(eddy_3c[,j,k,2],na.rm=T),TN = sum(eddy_3c[,j,k,3],na.rm=T),FN = sum(eddy_3c[,j,k,4],na.rm=T))
      eddy_3cSukosd_mcc_avg[j,k] <- mcc.calc(TP = sum(eddy_3cSukosd[,j,k,1],na.rm=T),FP = sum(eddy_3cSukosd[,j,k,2],na.rm=T),TN = sum(eddy_3cSukosd[,j,k,3],na.rm=T),FN = sum(eddy_3cSukosd[,j,k,4],na.rm=T))
    }
  }
  
  df <- matrix(NA,n_bin,5)
  df_sd <- matrix(NA,n_bin,5)
  
  df[,1] <- rowMeans(deigan_mcc_avg,na.rm=T)
  df_sd[,1] <- RowSD(deigan_mcc_avg)
  df[,2] <- rowMeans(eddy_2c_mcc_avg,na.rm=T)
  df_sd[,2] <- RowSD(eddy_2c_mcc_avg)
  df[,4] <- rowMeans(eddy_2cSukosd_mcc_avg,na.rm=T)
  df_sd[,4] <- RowSD(eddy_2cSukosd_mcc_avg)
  df[,3] <- rowMeans(eddy_3c_mcc_avg,na.rm=T)
  df_sd[,3] <- RowSD(eddy_3c_mcc_avg)
  df[,5] <- rowMeans(eddy_3cSukosd_mcc_avg,na.rm=T)
  df_sd[,5] <- RowSD(eddy_3cSukosd_mcc_avg)
  
  dff <- data.frame(
    Scheme = factor(c(rep("Lin",n_bin),rep("P2",n_bin),rep("P3",n_bin),rep("P2s",n_bin),rep("P3s",n_bin))),
    resp = c(df*100),
    Bins = factor(rep(data_label,5)),
    se = c(df_sd*100)
  )
  
  se_limits <- aes(ymax = resp + se, ymin=resp - se)
  dodge <- position_dodge(width=0.9)
  
  
  if (ty=="only") {
    plot_t <- "Bin included"
    plot3 <- ggplot(dff, aes(fill=Bins, y=resp, x=Scheme)) + 
      geom_bar(position="dodge", stat="identity",color="black") + 
      geom_errorbar(se_limits,position=dodge, width=0.2) +
      geom_hline(aes(yintercept=std_score[25,3]),size=0.9,linetype=2) +
      geom_segment(aes(x = 0.6, y = mean(deigan_std_mcc_avg[1:n_rep],na.rm=T)*100, xend = 1.4, yend = mean(deigan_std_mcc_avg[1:n_rep],na.rm=T)*100)) +
      geom_segment(aes(x = 1.6, y = mean(eddy_2c_std_mcc_avg[1:n_rep],na.rm=T)*100, xend = 2.4, yend = mean(eddy_2c_std_mcc_avg[1:n_rep],na.rm=T)*100)) +
      geom_segment(aes(x = 2.6, y = mean(eddy_3c_std_mcc_avg[1:n_rep],na.rm=T)*100, xend = 3.4, yend = mean(eddy_3c_std_mcc_avg[1:n_rep],na.rm=T)*100)) +
      geom_segment(aes(x = 3.6, y = mean(eddy_2cSukosd_std_mcc_avg[1:n_rep],na.rm=T)*100, xend = 4.4, yend = mean(eddy_2cSukosd_std_mcc_avg[1:n_rep],na.rm=T)*100)) +
      geom_segment(aes(x = 4.6, y = mean(eddy_3cSukosd_std_mcc_avg[1:n_rep],na.rm=T)*100, xend = 5.4, yend = mean(eddy_3cSukosd_std_mcc_avg[1:n_rep],na.rm=T)*100)) +
      ylab("") + xlab("") + theme_bw() + 
      ylim(0,100) +
      scale_fill_manual(values=gpalette) +
      theme(axis.text=element_text(size=ate),
            axis.title=element_text(size=ati,face="bold"),
            legend.text=element_text(size=20,face="bold"),
            legend.title=element_text(size=26,face="bold"),
            plot.title=element_text(size=34,face="bold"))
  } else {
    plot_t <- "Bin excluded"
    plot4 <- ggplot(dff, aes(fill=Bins, y=resp, x=Scheme)) + 
      geom_bar(position="dodge", stat="identity",color="black") + 
      geom_errorbar(se_limits,position=dodge, width=0.2) +
      geom_hline(aes(yintercept=std_score[25,3]),size=0.9,linetype=2) +
      geom_segment(aes(x = 0.6, y = mean(deigan_std_mcc_avg[1:n_rep],na.rm=T)*100, xend = 1.4, yend = mean(deigan_std_mcc_avg[1:n_rep],na.rm=T)*100)) +
      geom_segment(aes(x = 1.6, y = mean(eddy_2c_std_mcc_avg[1:n_rep],na.rm=T)*100, xend = 2.4, yend = mean(eddy_2c_std_mcc_avg[1:n_rep],na.rm=T)*100)) +
      geom_segment(aes(x = 2.6, y = mean(eddy_3c_std_mcc_avg[1:n_rep],na.rm=T)*100, xend = 3.4, yend = mean(eddy_3c_std_mcc_avg[1:n_rep],na.rm=T)*100)) +
      geom_segment(aes(x = 3.6, y = mean(eddy_2cSukosd_std_mcc_avg[1:n_rep],na.rm=T)*100, xend = 4.4, yend = mean(eddy_2cSukosd_std_mcc_avg[1:n_rep],na.rm=T)*100)) +
      geom_segment(aes(x = 4.6, y = mean(eddy_3cSukosd_std_mcc_avg[1:n_rep],na.rm=T)*100, xend = 5.4, yend = mean(eddy_3cSukosd_std_mcc_avg[1:n_rep],na.rm=T)*100)) +
      ylab("") + xlab("") + theme_bw() + 
      ylim(0,100) +
      scale_fill_manual(values=gpalette) +
      theme(axis.text=element_text(size=ate),
            axis.title=element_text(size=ati,face="bold"),
            legend.text=element_text(size=20,face="bold"),
            legend.title=element_text(size=26,face="bold"),
            plot.title=element_text(size=34,face="bold"))
  }
  
}

# MAKE FINAL PLOT ---------------------------------------------------------
root_out_base <- directory_out
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
mylegend<-g_legend(plot1)

text.fs <- 8 # annotation text font size

setEPS()
postscript(file = paste0(root_out_base,"missing_summary.eps"))
par(cex=1)
grid.arrange(arrangeGrob(plot2 + theme(legend.position="none"),plot1 + theme(legend.position="none"),
                         plot4 + theme(legend.position="none"),plot3 + theme(legend.position="none"),
                         nrow=2, ncol=2,
                         left = textGrob("SLW-average MCC [%]", rot = 90, vjust = 1),
                         bottom = textGrob("Scheme"),
                         right=textGrob(""),
                         top=textGrob("")),
                         mylegend,nrow=2,heights=c(15, 1))

grid.text("Leave-one-out",x = unit(0.30, "npc"),y = unit(0.985, "npc"))
grid.text("Leave-one-in", x = unit(0.78, "npc"),y = unit(0.985, "npc"))

grid.text("Real", x = unit(0.985, "npc"),y = unit(0.78, "npc"),rot=270)
grid.text("Simulated", x = unit(0.985, "npc"),y = unit(0.35, "npc"),rot=270)

dev.off()



