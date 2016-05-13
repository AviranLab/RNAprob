## Analyse the results from the scheme robustness test when rejecting or including selected reactivities

# Author: Mirko Ledda
# Date: 09/15/2015
# Version: v1.0

## Remove all active variables from the workspace
rm(list = ls())

# Load libraries ----------------------------------------------------------
lib <- c('ggplot2','reshape2','gridExtra','grid')
lapply(lib,require,character.only = TRUE)

# Set directory/file paths ------------------------------------------------
# folding scores for simulations/real data
fd_ref <- 'data/folding_results/ref_folding.RData'
fd_knif_real <- 'data/folding_results/jackknife_real_folding.Rdata'
fd_knif_sim <- 'data/folding_results/jackknife_sim_folding.Rdata'

root_out <- 'output/' # output directory
source_dir <- 'src/' # location of source scripts
root_ct <- 'data/ct_sequence_file/' # path to the ct files
fn_name <- 'data/data_pool_name.txt' # RNA names
fn_name_clean <- 'data/data_pool_name_clean.txt' # RNA names

# Variables --------------------------------------------------------
n_bin <- 5 # number of bins
n_rep <- 20 # number of simulation repeats
atxt <- 14 # axis text font size
ltxt <- 16 # legend text font size
ltit <- 16 # legend title font size
rna_ix <- 1:23 # Select RNAs based on name_list (1:23 would perform the analysis over ALL RNAs)

# Set bin ranks
bin_rank <- c(5,1,4,2,3)
y_ax <- c(0,110)

# Automated script starts here --------------------------------------------
root_out_base <- root_out
# load source
src_fn <- list.files(source_dir)
for (i in 1:length(src_fn)) source(paste(source_dir,src_fn[i],sep=""))

# load data
load(fd_ref)

# Prep plots ----------------------------------------------------------
# Make a grey scale palette
glim <- c(80,240)/255
giv  <- (glim[2] - glim[1])/n_bin
gpalette <- NULL
for (j in 1:n_bin) {
  gpalette <- c(gpalette,rgb(glim[1]+giv*(j-1),glim[1]+giv*(j-1),glim[1]+giv*(j-1)))
}

ref_ix <- 25

data_label <- NULL
for (l in 1:n_bin) {
  data_label <- c(data_label,paste(toString(round((bin_rank[l]-1) * 100/n_bin)),'-',toString(round((bin_rank[l]) * 100/n_bin)),'%',sep=''))
}

dodge <- position_dodge(width=0.9)

# REAL ------------------------------------------------------------
load(fd_knif_real)

# Normalize the mcc data
deigan_mcc_norm <- array(NA,dim=c(1,n_bin))
eddy_2c_mcc_norm <- array(NA,dim=c(1,n_bin))
eddy_2cSukosd_mcc_norm <- array(NA,dim=c(1,n_bin))
eddy_3c_mcc_norm <- array(NA,dim=c(1,n_bin))
eddy_3cSukosd_mcc_norm <- array(NA,dim=c(1,n_bin))

for (i in 1:n_bin) {
  deigan_mcc_norm[,i] <- (deigan[i,25]-std_score[25,3])/(std_score[25,6]-std_score[25,3])
  eddy_2c_mcc_norm[,i] <- (eddy_2c[i,25]-std_score[25,3])/(std_score[25,9]-std_score[25,3])
  eddy_2cSukosd_mcc_norm[,i] <- (eddy_2cSukosd[i,25]-std_score[25,3])/(std_score[25,15]-std_score[25,3])
  eddy_3c_mcc_norm[,i] <- (eddy_3c[i,25]-std_score[25,3])/(std_score[25,12]-std_score[25,3])
  eddy_3cSukosd_mcc_norm[,i] <- (eddy_3cSukosd[i,25]-std_score[25,3])/(std_score[25,18]-std_score[25,3])
}

df <- matrix(NA,n_bin,5)
df_sd <- matrix(NA,n_bin,5)

df[,1] <- deigan_mcc_norm
df[,2] <- eddy_2c_mcc_norm
df[,4] <- eddy_2cSukosd_mcc_norm
df[,3] <- eddy_3c_mcc_norm
df[,5] <- eddy_3cSukosd_mcc_norm

dff <- data.frame(
  Scheme = factor(c(rep("RNALin   ",n_bin),rep("RNAProb-2   ",n_bin),rep("RNAProb-3   ",n_bin),rep("RNAProb-2s   ",n_bin),rep("RNAProb-3s   ",n_bin))),
  resp = c(df*100),
  Subset = factor(rep(data_label,5),levels=data_label)
)

plot1 <- ggplot(dff, aes(y=resp, x=Subset,colour=Scheme)) + 
  geom_line(aes(group=Scheme),size=1.2) + 
  geom_hline(aes(yintercept=100),size=1,linetype = 2) +
  geom_hline(aes(yintercept=0),size=1,linetype = 2) +
  ylab("") + xlab("") + theme_bw() + theme(legend.position="bottom") +
  ylim(y_ax) +
  scale_colour_discrete(name = "Scheme   ") +
  theme(axis.text=element_text(size=atxt),
        legend.text=element_text(size=ltxt),
        legend.title=element_text(size=ltit,face="bold"))


# SIMULATIONS ------------------------------------------------------------
load(fd_knif_sim)

df <- matrix(NA,n_bin,5)
df_sd <- matrix(NA,n_bin,5)

df[,1] <- rowMeans(deigan_mcc_avg_norm,na.rm=T)
df_sd[,1] <- RowSD(deigan_mcc_avg_norm)
df[,2] <- rowMeans(eddy_2c_mcc_avg_norm,na.rm=T)
df_sd[,2] <- RowSD(eddy_2c_mcc_avg_norm)
df[,4] <- rowMeans(eddy_2cSukosd_mcc_avg_norm,na.rm=T)
df_sd[,4] <- RowSD(eddy_2cSukosd_mcc_avg_norm)
df[,3] <- rowMeans(eddy_3c_mcc_avg_norm,na.rm=T)
df_sd[,3] <- RowSD(eddy_3c_mcc_avg_norm)
df[,5] <- rowMeans(eddy_3cSukosd_mcc_avg_norm,na.rm=T)
df_sd[,5] <- RowSD(eddy_3cSukosd_mcc_avg_norm)

dff <- data.frame(
  Scheme = factor(c(rep("RNALin",n_bin),rep("RNAProb-2",n_bin),rep("RNAProb-3",n_bin),rep("RNAProb-2s",n_bin),rep("RNAProb-3s",n_bin))),
  resp = c(df*100),
  Subset = factor(rep(data_label,5),levels=data_label),
  se = c(df_sd*100)
)

se_limits <- aes(ymax = resp + se, ymin=resp - se)

if (itype=="only") {
  plot_t <- "Signal reconstruction"
}

plot2 <- ggplot(dff, aes(y=resp, x=Subset,colour=Scheme)) + 
  geom_line(aes(group=Scheme),size=1.2) + 
  geom_errorbar(se_limits, width=0.2,size=1) +
  geom_hline(aes(yintercept=100),size=1,linetype = 2) +
  geom_hline(aes(yintercept=0),size=1,linetype = 2) +
  ylab("") + xlab("") + theme_bw() +
  ylim(y_ax) +
  theme(axis.text=element_text(size=atxt),
        legend.text=element_text(size=ltxt),
        legend.title=element_text(size=ltit,face="bold"))


# Combine plots -----------------------------------------------------------
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
mylegend<-g_legend(plot1)

text.fs <- 18 # annotation text font size


setEPS()
postscript(file = paste0(root_out_base,"reconstruction_summary.eps"),width = 10)
par(cex=1)
grid.arrange(arrangeGrob(plot1 + theme(legend.position="none"),plot2 + theme(legend.position="none"),
                         nrow=1, ncol=2,
                         left = textGrob("Normalized MCC [%]", rot = 90, vjust = 1,gp=gpar(fontsize=text.fs)),
                         bottom = textGrob("Quintile",gp=gpar(fontsize=text.fs)),
                         right=textGrob(""),
                         top=textGrob("")),
             mylegend,nrow=2,heights=c(15, 1))

grid.text("Real",x = unit(0.28, "npc"),y = unit(0.985, "npc"),gp=gpar(fontsize=text.fs))
grid.text("Simulated", x = unit(0.78, "npc"),y = unit(0.985, "npc"),gp=gpar(fontsize=text.fs))

dev.off()
