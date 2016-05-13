## Analyze Fei Naive simulation data

## Remove all active variables from the workspace
rm(list = ls())

# Load libraries
lib <- c('ggplot2','reshape2','grid')
lapply(lib, require, character.only=T)

# Set directory/file paths ------------------------------------------------
source_dir <- 'src/' # location of source scripts
fn_name <- 'data/data_pool_name.txt' # RNA names
root_out_base <- 'output/plot_all_folding/'

# folding scores for simulations/real data
fd_ref <- 'data/folding_results/ref_folding.RData'
fd_real <- 'data/folding_results/real_shape_folding.RData'
fd_sim <- 'data/folding_results/sim_shape_folding.RData'
fd_density <- 'data/folding_results/kde_density_folding.RData'

# User inputs -------------------------------------------------------------
sim_n <- 100 # number of simulations available
alpha <- 0.05 # Confidence interval

# Automated Script --------------------------------------------------------
# Load sources and data ---------------------------------------------------
dir.create(root_out_base)

# Get SHAPE and ct filename in the folder
name_list <- as.character(read.table(fn_name)$V1)
name_list <- name_list[1:23] # No Siegfield 16S RNA

# List the files in the source directory
src_fn <- list.files(source_dir)

# Load the source files
for (i in 1:length(src_fn)) source(paste0(source_dir,src_fn[i]))

# LOAD DATA ---------------------------------------------------------------
scheme_name <- c("deigan","eddy2c","eddy2cSukosd","eddy3c","eddy3cSukosd")

load(fd_density)
load(fd_sim)
load(fd_real)
load(fd_ref)

# process density sim data
data.sim.binary.mcc.all <- matrix(NA,5,sim_n)
data.sim.ternary.mcc.all <- matrix(NA,5,sim_n)
density.mcc.all <- matrix(NA,5,sim_n)
for (j in 1:5) {
  for (k in 1:sim_n) {
    data.sim.binary.mcc.all[j,k] <- mcc.calc(TP=sum(data.sim.binary[,j,k,1]),FP=sum(data.sim.binary[,j,k,3]),TN=sum(data.sim.binary[,j,k,2]),FN=sum(data.sim.binary[,j,k,4]))
    data.sim.ternary.mcc.all[j,k] <- mcc.calc(TP=sum(data.sim.ternary[,j,k,1]),FP=sum(data.sim.ternary[,j,k,3]),TN=sum(data.sim.ternary[,j,k,2]),FN=sum(data.sim.ternary[,j,k,4]))
    density.mcc.all[j,k] <- mcc.calc(TP=sum(data.density.sim[,j,k,1]),FP=sum(data.density.sim[,j,k,2]),TN=sum(data.density.sim[,j,k,3]),FN=sum(data.density.sim[,j,k,4]))
  }
}

# PLOT THE DATA -----------------------------------------------------------
ati <- 26 # axis title font size
ate <- 26 # axis text font size

glim <- c(80,240)/255
giv  <- (glim[2] - glim[1])/4
gpalette <- NULL
for (j in 1:4) {
  gpalette <- c(gpalette,rgb(glim[1]+giv*(j-1),glim[1]+giv*(j-1),glim[1]+giv*(j-1)))
}

data_label <- c("Real     ","Binary     ","Ternary     ","KDE     ")

df <- matrix(NA,4,5)
df_sd <- matrix(NA,4,5)

df[1,] <- data.real[25,]/100
df[2,] <- rowMeans(data.sim.binary.mcc.all)
df[3,] <- rowMeans(data.sim.ternary.mcc.all)
df[4,] <- rowMeans(density.mcc.all)

df_sd[1,] <- rep(0,5)
df_sd[2,] <- RowSD(data.sim.binary.mcc.all)
df_sd[3,] <- RowSD(data.sim.ternary.mcc.all)
df_sd[4,] <- RowSD(density.mcc.all)


dff <- data.frame(
  Scheme = factor(c(rep("RNAlin",4),rep("RNAProb-2",4),rep("RNAProb-2s",4),rep("RNAProb-3",4),rep("RNAProb-3s",4))),
  resp = c(df*100),
  Data = factor(rep(data_label,5),levels=data_label),
  se = c(df_sd*100)
)

se_limits <- aes(ymax = resp + se, ymin=resp - se)
dodge <- position_dodge(width=0.9)

ggplot(dff, aes(fill=Data, y=resp, x=Scheme)) + 
  geom_bar(position="dodge", stat="identity",color="black") + 
  geom_errorbar(se_limits,position=dodge, width=0.3) +
  geom_hline(aes(yintercept=std_score[25,3]),size=2,linetype=2) +
#  geom_segment(aes(x = 0.6, y = std_score[25,6], xend = 1.4, yend = std_score[25,6])) +
#  geom_segment(aes(x = 1.6, y = std_score[25,9], xend = 2.4, yend = std_score[25,9])) +
#  geom_segment(aes(x = 2.6, y = std_score[25,15], xend = 3.4, yend = std_score[25,15])) +
#  geom_segment(aes(x = 3.6, y = std_score[25,12], xend = 4.4, yend = std_score[25,12])) +
#  geom_segment(aes(x = 4.6, y = std_score[25,18], xend = 5.4, yend = std_score[25,18])) +
  ylab("SLW-average MCC [%]") + xlab("Scheme") + theme_bw() + theme(legend.position="bottom") +
  ylim(0,100) + 
  scale_fill_manual(name="Data     ",values=gpalette) +
  theme(axis.text=element_text(size=ate),
        axis.title=element_text(size=ati,face="bold"),
        legend.text=element_text(size=26),
        legend.title=element_text(size=30,face="bold"),
        plot.title=element_text(size=34,face="bold"),
        axis.title.x=element_text(vjust=-0.5)) +
  geom_segment(aes(x=3.5,y=95,xend=4,yend=95),size=2,linetype=2) +
  annotate("text", x = 4.3, y = 95, label = "No SHAPE",size=8) +
  geom_rect(aes(xmin = 3.4, xmax = 4.6, ymin = 92.5, ymax = 97.5),
            fill = "transparent", color = "black", size = 1)
  #annotate("rect", xmin = 3.4, xmax = 5, ymin = 90, ymax = 100,alpha = 0)

ggsave(filename = paste(root_out_base,'results_summary.eps',sep=''),width = 16,height = 11)




# STATS -------------------------------------------------------------------
p.binary <- matrix(NA,23,5)
p.ternary <- matrix(NA,23,5)
p.density <- matrix(NA,23,5)

for (i in 1:length(name_list)) {
  for (j in 1:length(scheme_name)) {
    p.binary[i,j] <- t.test(data.sim.binary.mcc[i,j,]*100,mu=data.real[i,j])$p.value
    p.ternary[i,j] <- t.test(data.sim.ternary.mcc[i,j,]*100,mu=data.real[i,j])$p.value
    p.density[i,j] <- t.test(density.mcc[i,j,]*100,mu=data.real[i,j])$p.value
  }
}

# Count the number of significant difference after Holm correction
n.binary <- colSums(matrix(p.adjust(p=c(p.binary),n = 23*5,method="fdr") < 0.05,23,5))
n.ternary <- colSums(matrix(p.adjust(p=c(p.ternary),n = 23*5,method="fdr") < 0.05,23,5))
n.density <- colSums(matrix(p.adjust(p=c(p.density),n = 23*5,method="fdr") < 0.05,23,5))
n.binary.05 <- colSums(p.binary < alpha)
n.ternary.05 <- colSums(p.ternary < alpha)
n.density.05 <- colSums(p.density < alpha)

# Test SLW-average MCC
p.binary.slw <- rep(NA,5)
p.ternary.slw <- rep(NA,5)
p.density.slw <- rep(NA,5)

for (j in 1:length(scheme_name)) {
    p.binary.slw[j] <- t.test(data.sim.binary.mcc.all[j,]*100,mu=data.real[25,j])$p.value
    p.ternary.slw[j] <- t.test(data.sim.ternary.mcc.all[j,]*100,mu=data.real[25,j])$p.value
    p.density.slw[j] <- t.test(density.mcc[i,j,]*100,mu=data.real[25,j])$p.value
}

p.binary.slw <- p.adjust(p.binary.slw,5,method = "fdr")
p.ternary.slw <- p.adjust(p.ternary.slw,5,method = "fdr")
p.density.slw <- p.adjust(p.density.slw,5,method = "fdr")


# Pairwise t-test between encoder/decoder results on simulations
scheme_name.binary <- c("deigan-2","eddy2c-2","eddy2cSukosd-2","eddy3c-2","eddy3cSukosd-2")
scheme_name.ternary <- c("deigan-3","eddy2c-3","eddy2cSukosd-3","eddy3c-3","eddy3cSukosd-3")

df <- c(data.sim.binary.mcc.all,data.sim.ternary.mcc.all)
df.grp <- c(rep(scheme_name.binary,sim_n),rep(scheme_name.ternary,sim_n))

ttest_p_value <- pairwise.t.test(df,df.grp,p.adjust.method = "fdr")$p.value
ttest_sig <- ttest_p_value<0.05

write.csv(ttest_p_value,file = paste0(root_out_base,"pairwise_ttest_sim_data.csv"))

# Comparison binary vs ternary
t.test(x=c(data.sim.binary.mcc),y=c(data.sim.ternary.mcc))
t.test(x=c(data.sim.binary.mcc.all),y=c(data.sim.ternary.mcc.all))

# Pairwise t-test between encoder/decoder results on real data
scheme_name.binary <- c("deigan-2","eddy2c-2","eddy2cSukosd-2","eddy3c-2","eddy3cSukosd-2")
scheme_name.ternary <- c("deigan-3","eddy2c-3","eddy2cSukosd-3","eddy3c-3","eddy3cSukosd-3")

df <- c(data.real[1:23,])
df.grp <- rep(scheme_name,each=23)

ttest_p_value <- pairwise.t.test(df,df.grp,p.adjust.method = "fdr")$p.value
ttest_sig <- ttest_p_value<0.05

write.csv(ttest_p_value,file = paste0(root_out_base,"pairwise_ttest_real_data.csv"))



