### Plot mock probe simulations

## Remove all active variables from the workspace
rm(list = ls())

# Load libraries ----------------------------------------------------------
lib <- c('ggplot2','gridExtra','grid','evd')
lapply(lib,require,character.only = TRUE)

# Set directory/file paths ------------------------------------------------
dir_out <- "output/" # Output folder
source_dir <- 'src/' # source directory
ctp <- "data/ct_sequence_file/" # Path to ct files
train_fn <- 'data/data_pool_list_train.txt' # RNA file names
nsim <- 10 # number of simulations

# Load data ---------------------------------------------------------------
# load sources
src_fn <- list.files(source_dir)
for (i in 1:length(src_fn)) source(paste(source_dir,src_fn[i],sep=""))

# Load training file
name_list <- as.character(read.table(train_fn)$V1)

# Scenario 1 --------------------------------------------------------------
sc_dir <- 'sim_probe_sc1/'
param_names <- c("mean","sd")
ref_pair_all <- NULL
seq_all <- NULL

# Loop across simulations
for (k in 1:nsim) {
  fitRes <- matrix(NA,3,length(param_names))
  colnames(fitRes) <- param_names
  rownames(fitRes) <- c("Unpaired","Helix-end","Stacked")
  
  # Loop across RNA
  for (rnai in 1:length(name_list)) {
    # read ct and sequence
    ct <- read.table(paste0(ctp,name_list[rnai],".ct"),skip=1)
    seq <- ct$V2
    seq_all <- c(seq_all,as.character(seq))
    
    # Get pairing
    ref_pair <- ct$V5
    ref_pair_all <- c(ref_pair_all,get_pair_3s(ref_pair))
  }
}

# get stats based on data
fAU <- seq_all == "A" | seq_all == "U"
fGC <- seq_all == "G" | seq_all == "C"
fup <- ref_pair_all == 0
fph <- ref_pair_all == 1
fps <- ref_pair_all == 2

# Plot data ---------------------------------------------------------------
# Plot options
trans_alpha <- 1 # transparency
ate <- 16
ati <- 16

# Build density functions
x0 <- 0
x1 <- 3
iv <- 0.01
x <- seq(x0,x1,iv)

Cdens1 <- function(x,iv,...) {
  y <- dgev(x,...)
  return(y)
}

Cdens2 <- function(x,iv,...) {
  y <- dexp(x,...)
  return(y)
}

df <- data.frame(
  x = x,
  up = Cdens2(x,iv,rate=1.468),
  ph = Cdens2(x,iv,rate=1.468),
  ph_real = Cdens1(x,iv,loc=0.090,scale=0.114,shape=0.821),
  ps = Cdens1(x,iv,loc=0.040,scale=0.049,shape=0.763)
)

# Plot the data
gg1 <- ggplot(data = df,aes(x=x)) +
  geom_ribbon(aes(ymax=up),ymin=0,fill="gray40",colour=NA,alpha=trans_alpha) +
  xlab("") + ylab("") +
  ylim(c(0,8)) +
  theme(axis.text=element_text(size=ate),
        axis.title=element_text(size=ati,face="bold"),
        legend.text=element_text(size=26),
        legend.title=element_text(size=30,face="bold"),
        plot.title=element_text(size=34,face="bold"),
        axis.title.x=element_text(vjust=-0.5)) +
  theme_bw()

gg2 <- ggplot(data = df,aes(x=x)) +
  geom_ribbon(aes(ymax=ph),ymin=0,fill="gray40",colour=NA,alpha=1) +
  geom_line(aes(y=ph_real),col="red",lwd=0.7,lty="dashed") +
  xlab("") + ylab("") +
  ylim(c(0,8)) +
  theme(axis.text=element_text(size=ate),
        axis.title=element_text(size=ati,face="bold"),
        legend.text=element_text(size=26),
        legend.title=element_text(size=30,face="bold"),
        plot.title=element_text(size=34,face="bold"),
        axis.title.x=element_text(vjust=-0.5)) +
  theme_bw()

gg3 <- ggplot(data = df,aes(x=x)) +
  geom_ribbon(aes(ymax=ps),ymin=0,fill="gray40",colour=NA,alpha=trans_alpha) +
  xlab("") + ylab("") +
  ylim(c(0,8)) +
  theme(axis.text=element_text(size=ate),
        axis.title=element_text(size=ati,face="bold"),
        legend.text=element_text(size=26),
        legend.title=element_text(size=30,face="bold"),
        plot.title=element_text(size=34,face="bold"),
        axis.title.x=element_text(vjust=-0.5)) +
  theme_bw()

# Scenario 2 --------------------------------------------------------------
sc_dir <- 'sim_probe_sc2/'
root_data_sc <- paste0(root_data,sc_dir)
param_names <- c("mean","sd")
ref_pair_all <- NULL
seq_all <- NULL

# Loop across simulations
for (k in 1:nsim) {
  fitRes <- matrix(NA,3,length(param_names))
  colnames(fitRes) <- param_names
  rownames(fitRes) <- c("Unpaired","Helix-end","Stacked")
  
  # Loop across RNA
  for (rnai in 1:length(name_list)) {
    # read ct and sequence
    ct <- read.table(paste0(ctp,name_list[rnai],".ct"),skip=1)
    seq <- ct$V2
    seq_all <- c(seq_all,as.character(seq))
    
    # Get pairing
    ref_pair <- ct$V5
    ref_pair_all <- c(ref_pair_all,get_pair_3s(ref_pair))
  }
}

# get stats based on data
fAU <- seq_all == "A" | seq_all == "U"
fGC <- seq_all == "G" | seq_all == "C"
fup <- ref_pair_all == 0
fph <- ref_pair_all == 1
fps <- ref_pair_all == 2

# Plot data ---------------------------------------------------------------
# Plot options
trans_alpha <- 1 # transparency
ate <- 16
ati <- 16
cAU <- "black"
cGC <- "gray"

# Build density functions
x0 <- 0
x1 <- 3
iv <- 0.01
x <- seq(x0,x1,iv)

Cdens1 <- function(x,iv,...) {
  y <- dnorm(x,...)
  return(y)
}

df <- data.frame(
  x = x,
  up_AU = Cdens1(x,iv,mean=0.5,sd=0.1)*sum(fAU & fup)/sum(fup),
  ph_AU = Cdens1(x,iv,mean=0.25,sd=0.1)*sum(fAU & fph)/sum(fph),
  ps_AU = Cdens1(x,iv,mean=0,sd=0.1)*sum(fAU & fps)/sum(fps),
  up_GC = Cdens1(x,iv,mean=0.5,sd=1)*sum(fGC & fup)/sum(fup),
  ph_GC = Cdens1(x,iv,mean=0.25,sd=1)*sum(fGC & fph)/sum(fph),
  ps_GC = Cdens1(x,iv,mean=0,sd=1)*sum(fGC & fps)/sum(fps)
)

# Plot the data
gg4 <- ggplot(data = df,aes(x=x)) +
  geom_ribbon(aes(ymax=up_GC),ymin=0,fill=cGC,colour=NA,alpha=trans_alpha) +
  geom_ribbon(aes(ymax=up_AU+up_GC,ymin=up_GC),fill=cAU,colour=NA,alpha=trans_alpha) +
  ylim(c(0,3)) +
  xlab("") + ylab("") +
  theme(axis.text=element_text(size=ate),
        axis.title=element_text(size=ati,face="bold"),
        legend.text=element_text(size=26),
        legend.title=element_text(size=30,face="bold"),
        plot.title=element_text(size=34,face="bold"),
        axis.title.x=element_text(vjust=-0.5)) +
  theme_bw()

gg5 <- ggplot(data = df,aes(x=x)) +
  geom_ribbon(aes(ymax=ph_GC),ymin=0,fill=cGC,colour=NA,alpha=trans_alpha) +
  geom_ribbon(aes(ymax=ph_AU+ph_GC,ymin=ph_GC),fill=cAU,colour=NA,alpha=trans_alpha) +
  ylim(c(0,3)) +
  xlab("") + ylab("") +
  theme(axis.text=element_text(size=ate),
        axis.title=element_text(size=ati,face="bold"),
        legend.text=element_text(size=26),
        legend.title=element_text(size=30,face="bold"),
        plot.title=element_text(size=34,face="bold"),
        axis.title.x=element_text(vjust=-0.5)) +
  theme_bw()

gg6 <- ggplot(data = df,aes(x=x)) +
  geom_ribbon(aes(ymax=ps_GC),ymin=0,fill=cGC,colour=NA,alpha=trans_alpha) +
  geom_ribbon(aes(ymax=ps_AU+ps_GC,ymin=ps_GC),fill=cAU,colour=NA,alpha=trans_alpha) +
  ylim(c(0,3)) +
  xlab("") + ylab("") +
  theme(axis.text=element_text(size=ate),
        axis.title=element_text(size=ati,face="bold"),
        legend.text=element_text(size=26),
        legend.title=element_text(size=30,face="bold"),
        plot.title=element_text(size=34,face="bold"),
        axis.title.x=element_text(vjust=-0.5)) +
  theme_bw()

# get probabilities external to the ploting window
Cdens_ext <- function(x,iv,...) {
  n <- length(x)
  y0 <- pnorm(0,...)
  yn <- 1-pnorm(x[n],...)
  return(c(y0,yn))
}

df_ext <- data.frame(
  up_AU = Cdens_ext(x,iv,mean=0.5,sd=0.1),
  ph_AU = Cdens_ext(x,iv,mean=0.25,sd=0.1),
  ps_AU = Cdens_ext(x,iv,mean=0,sd=0.1),
  up_GC = Cdens_ext(x,iv,mean=0.5,sd=1),
  ph_GC = Cdens_ext(x,iv,mean=0.25,sd=1),
  ps_GC = Cdens_ext(x,iv,mean=0,sd=1)
)


# Legend definition -------------------------------------------------------
# Make a junk plot to construct the legend
df_jk <- data.frame(
  x = c(0,1),
  y = c("A & U","G & C")
)
gg_jk <- ggplot(data = df_jk,aes(x=x,fill=y)) +
  geom_bar() + 
  scale_fill_manual(name="Bases",values=c("black","gray")) +
  theme(legend.position="bottom",
        legend.text = element_text(size=8),
        legend.key.height=unit(0.5,"cm"),
        legend.key.width=unit(0.5,"cm"))
  
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
mylegend<-g_legend(gg_jk)
lheight <- sum(mylegend$height)

# Plot all data -----------------------------------------------------------
setEPS()
postscript(file = paste0(dir_out,'mock_probe_dist.eps'),width=6,height=4)
par(cex=1)
grid.arrange(arrangeGrob(gg1,gg2,gg3,gg4,gg5,gg6,
                         nrow=2, ncol=3,
                         left = textGrob("Probability density", rot = 90, vjust = 1),
                         bottom = textGrob("Reactivity score"),
                         right=textGrob(""),
                         top=textGrob("")),
             mylegend,nrow=2,heights = unit.c(unit(1, "npc") - lheight, lheight))

grid.text("Unpaired",x = unit(0.22, "npc"),y = unit(0.970, "npc"))
grid.text("Helix-end", x = unit(0.54, "npc"),y = unit(0.970, "npc"))
grid.text("Stacked", x = unit(0.86, "npc"),y = unit(0.970, "npc"))

ptxt_size=6
grid.text(expression(paste("P(x<0|",x %in% "{A,U})" %~~% "0")),x = unit(0.25, "npc"),y = unit(0.45, "npc"),gp=gpar(fontsize=ptxt_size))
grid.text(expression(paste("P(x<0|",x %in% "{G,C}) = 0.31")),x = unit(0.25, "npc"),y = unit(0.40, "npc"),gp=gpar(fontsize=ptxt_size))
grid.text(expression(paste("P(x<0|",x %in% "{A,U})" %~~% "0")),x = unit(0.55, "npc"),y = unit(0.45, "npc"),gp=gpar(fontsize=ptxt_size))
grid.text(expression(paste("P(x<0|",x %in% "{G,C}) = 0.40")),x = unit(0.55, "npc"),y = unit(0.40, "npc"),gp=gpar(fontsize=ptxt_size))
grid.text(expression(paste("P(x<0|",x %in% "{A,U}) = 0.5")),x = unit(0.85, "npc"),y = unit(0.45, "npc"),gp=gpar(fontsize=ptxt_size))
grid.text(expression(paste("P(x<0|",x %in% "{G,C}) = 0.5")),x = unit(0.85, "npc"),y = unit(0.40, "npc"),gp=gpar(fontsize=ptxt_size))

grid.text("Scenario 1", x = unit(0.985, "npc"),y = unit(0.80 , "npc"),rot=270)
grid.text("Scenario 2", x = unit(0.985, "npc"),y = unit(0.42, "npc"),rot=270)

dev.off()

