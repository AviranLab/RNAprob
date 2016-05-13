
## Remove all active variables from the workspace
rm(list = ls())

# Load libraries ----------------------------------------------------------
lib <- c('ggplot2','reshape2','gridExtra','grid','nortest')
lapply(lib,require,character.only = TRUE)

# Set directory/file paths ------------------------------------------------
source_dir <- 'src/' # location of source scripts
root_ct <- 'data/ct_sequence_file/' # path to the ct files
rt_data <- 'data/generated_data/reassign_zeros_reactivities/' # root to generated data
root_raw_shape <- 'data/shape_raw/'
fn_name <- 'data/data_pool_name.txt' # RNA names
fn_name_clean <- 'data/data_pool_name_clean.txt'
root_out_base <- 'output/' # set output root directory

# Load sources and data ---------------------------------------------------
src_fn <- list.files(source_dir) # List the files in the source directory
for (i in 1:length(src_fn)) source(paste0(source_dir,src_fn[i])) # Load the source files

# Get SHAPE and ct filename in the RNA repository folder
name_list <- as.character(read.table(fn_name)$V1)
name_list <- name_list[1:23] # No Siegfield 16S RNA
name_list_clean <- as.character(read.table(fn_name_clean)$V1)

# Load data ---------------------------------------------------------------
for (k in 1:3) {
  
  # set list for concatenated data
  all_pairing <- list()
  all_pairing$pairing <- NULL
  all_pairing$all <- NULL
  
  if (k==1) {
    root_data <- root_raw_shape
  } else if (k==2) {
    root_data <- paste0(rt_data,'random_absolute_boxcox/')
  } else if (k==3) {
    root_data <- paste0(rt_data,'random_absolute_ln/')
  }
  
  # Loop across RNA types
  for (i in 1:length(name_list)) {
    
    # Read the shape reactivity file
    shape_data <- read.table(paste(root_data,name_list[i],'.shape',sep=''))
    
    # Get the associated .ct amd .fa files
    ct_fn <- paste(root_ct,name_list[i],'.ct',sep='')
    fa_fn <- paste(root_ct,name_list[i],'.fa',sep='')
    
    #Set NaN reactivities
    nan_reac_filt <- shape_data$V2==-999
    Reactivities <- shape_data$V2
    Reactivities[nan_reac_filt] <- NA
    reac <- Reactivities
    
    # Pair the nucleotide according to 3 states, unpaired/stacked/helix-end
    x_pairing <- pairing_reactivities(x=reac,ref_fn = paste(root_ct,name_list[i],'.ct',sep=''),model = 'raw',states = 3)
    
    # concatenate the results
    all_pairing$pairing <- c(all_pairing$pairing,x_pairing$pairing)
    all_pairing$all <- c(all_pairing$all,reac)
    
  }
  
  if (k==1) {
    raw <- list()
    raw$pairing <- all_pairing$pairing
    raw$reac <- all_pairing$all
  } else if (k==2) {
    boxcox <- list()
    boxcox$pairing <- all_pairing$pairing
    boxcox$reac <- all_pairing$all
  } else if (k==3) {
    ln <- list()
    ln$pairing <- all_pairing$pairing
    ln$reac <- all_pairing$all
  }
}

# Plot reactivities -------------------------------------------------------
tls <- 16 # title font size
ati <- 12 # axis title font size
ate <- 12 # axis text font size
lsi <- 1 # red line size

# Force digits in the axis
fmt <- function(){
  function(x) format(x,nsmall = 1,scientific = FALSE)
}

# Raw ---------------------------------------------------------------------
xl <- c(-1,5)
bw <- (xl[2]-xl[1])/100

raw$up <- raw$reac[raw$pairing == 0]
raw$ph <- raw$reac[raw$pairing == 1]
raw$ps <- raw$reac[raw$pairing == 2]

for1 <- function(x,lambda) {
  y = lambda*exp(-lambda*x)
  return(y)
}

for2 <- function(x,xi,sigma,mu) {
  y = (1/sigma)*((1+xi*(x-mu)/sigma)^-(1+1/xi))*exp(-(1+xi*(x-mu)/sigma)^(-1/xi))
  return(y)
}

pred <- data.frame(x=seq(0,xl[2],bw),y=for1(seq(0,xl[2],bw),lambda=1.468))
plot1 <- ggplot(data.frame(raw$up),aes(x=raw.up)) +
  geom_histogram(binwidth=bw,aes(y = ..density..)) +
  xlim(xl) + xlab("") + ylab("") +
  geom_line(data=pred,aes(x=x,y=y),colour="red",size=lsi) +
  theme_bw() +
  theme(axis.text=element_text(size=ate),
        axis.title=element_text(size=ati,face="bold"),
        plot.title=element_text(size=tls,face="bold"))


pred <- data.frame(x=seq(xl[1],xl[2],bw),y=for2(seq(xl[1],xl[2],bw),xi=0.821,sigma=0.114,mu=0.090))
plot2 <- ggplot(data.frame(raw$ph),aes(x=raw.ph)) +
  geom_histogram(binwidth=bw,aes(y = ..density..)) +
  xlim(xl) + xlab("") + ylab("") + scale_y_continuous(labels=fmt()) +
  geom_line(data=pred,aes(x=x,y=y),colour="red",size=lsi) +
  theme_bw() +
  theme(axis.text=element_text(size=ate),
        axis.title=element_text(size=ati,face="bold"),
        plot.title=element_text(size=tls,face="bold"))

pred <- data.frame(x=seq(xl[1],xl[2],bw),y=for2(seq(xl[1],xl[2],bw),xi=0.763,sigma=0.049,mu=0.040))
plot3 <- ggplot(data.frame(raw$ps),aes(x=raw.ps)) +
  geom_histogram(binwidth=bw,aes(y = ..density..)) +
  xlim(xl) + xlab("") + ylab("") +
  geom_line(data=pred,aes(x=x,y=y),colour="red",size=lsi) +
  theme_bw() +
  theme(axis.text=element_text(size=ate),
        axis.title=element_text(size=ati,face="bold"),
        plot.title=element_text(size=tls,face="bold"))

# Boxcox ---------------------------------------------------------------------
xl <- c(-8,4)
bw <- (xl[2]-xl[1])/100

boxcox$up <- boxcox$reac[boxcox$pairing == 0]
boxcox$ph <- boxcox$reac[boxcox$pairing == 1]
boxcox$ps <- boxcox$reac[boxcox$pairing == 2]

pred <- data.frame(x=seq(xl[1],xl[2],bw),y=dnorm(seq(xl[1],xl[2],bw),mean = mean(boxcox$up,na.rm=T),sd = sd(boxcox$up,na.rm=T)))
plot4 <- ggplot(data.frame(boxcox$up),aes(x=boxcox.up)) +
  geom_histogram(binwidth=bw,aes(y = ..density..)) +
  xlim(xl) + xlab("") + ylab("") +
  geom_line(data=pred,aes(x=x,y=y),colour="red",size=lsi) +
  theme_bw() +
  theme(axis.text=element_text(size=ate),
        axis.title=element_text(size=ati,face="bold"),
        plot.title=element_text(size=tls,face="bold"))

pred <- data.frame(x=seq(xl[1],xl[2],bw),y=dnorm(seq(xl[1],xl[2],bw),mean = mean(boxcox$ph,na.rm=T),sd = sd(boxcox$ph,na.rm=T)))
plot5 <- ggplot(data.frame(boxcox$ph),aes(x=boxcox.ph)) +
  geom_histogram(binwidth=bw,aes(y = ..density..)) +
  xlim(xl) + xlab("") + ylab("") +
  geom_line(data=pred,aes(x=x,y=y),colour="red",size=lsi) +
  theme_bw() +
  theme(axis.text=element_text(size=ate),
        axis.title=element_text(size=ati,face="bold"),
        plot.title=element_text(size=tls,face="bold"))

pred <- data.frame(x=seq(xl[1],xl[2],bw),y=dnorm(seq(xl[1],xl[2],bw),mean = mean(boxcox$ps,na.rm=T),sd = sd(boxcox$ps,na.rm=T)))
plot6 <- ggplot(data.frame(boxcox$ps),aes(x=boxcox.ps)) +
  geom_histogram(binwidth=bw,aes(y = ..density..)) +
  xlim(xl) + xlab("") + ylab("") +
  geom_line(data=pred,aes(x=x,y=y),colour="red",size=lsi) +
  theme_bw() +
  theme(axis.text=element_text(size=ate),
        axis.title=element_text(size=ati,face="bold"),
        plot.title=element_text(size=tls,face="bold"))

# Ln ---------------------------------------------------------------------
xl <- c(-10,3)
bw <- (xl[2]-xl[1])/100

ln$up <- ln$reac[ln$pairing == 0]
ln$ph <- ln$reac[ln$pairing == 1]
ln$ps <- ln$reac[ln$pairing == 2]

pred <- data.frame(x=seq(xl[1],xl[2],bw),y=dnorm(seq(xl[1],xl[2],bw),mean = mean(ln$up,na.rm=T),sd = sd(ln$up,na.rm=T)))
plot7 <- ggplot(data.frame(ln$up),aes(x=ln.up)) +
  geom_histogram(binwidth=bw,aes(y = ..density..)) +
  xlim(xl) + xlab("") + ylab("") +
  geom_line(data=pred,aes(x=x,y=y),colour="red",size=lsi) +
  theme_bw() +
  theme(axis.text=element_text(size=ate),
        axis.title=element_text(size=ati,face="bold"),
        plot.title=element_text(size=tls,face="bold"))

pred <- data.frame(x=seq(xl[1],xl[2],bw),y=dnorm(seq(xl[1],xl[2],bw),mean = mean(ln$ph,na.rm=T),sd = sd(ln$ph,na.rm=T)))
plot8 <- ggplot(data.frame(ln$ph),aes(x=ln.ph)) +
  geom_histogram(binwidth=bw,aes(y = ..density..)) +
  xlim(xl) + xlab("") + ylab("") +
  geom_line(data=pred,aes(x=x,y=y),colour="red",size=lsi) +
  theme_bw() +
  theme(axis.text=element_text(size=ate),
        axis.title=element_text(size=ati,face="bold"),
        plot.title=element_text(size=tls,face="bold"))

pred <- data.frame(x=seq(xl[1],xl[2],bw),y=dnorm(seq(xl[1],xl[2],bw),mean = mean(ln$ps,na.rm=T),sd = sd(ln$ps,na.rm=T)))
plot9 <- ggplot(data.frame(ln$ps),aes(x=ln.ps)) +
  geom_histogram(binwidth=bw,aes(y = ..density..)) +
  xlim(xl) + xlab("") + ylab("") +
  geom_line(data=pred,aes(x=x,y=y),colour="red",size=lsi) +
  theme_bw() +
  theme(axis.text=element_text(size=ate),
        axis.title=element_text(size=ati,face="bold"),
        plot.title=element_text(size=tls,face="bold"))


# Make final plot ---------------------------------------------------------
text.fs <- 8 # annotation text font size

setEPS()
postscript(file = paste0(root_out_base,'distributions.eps'))
par(cex=1)
grid.arrange(arrangeGrob(plot1,plot2, plot3, plot7, plot8, plot9, plot4, plot5, plot6,
                         nrow=3, ncol=3,
                         left = textGrob("Probability density", rot = 90, vjust = 1),
                         bottom = textGrob("Reactivity score"),
                         right=textGrob(""),
                         top=textGrob("")))

grid.text("Unpaired",x = unit(0.22, "npc"),y = unit(0.985, "npc"))
grid.text("Helix-end", x = unit(0.54, "npc"),y = unit(0.985, "npc"))
grid.text("Stacked", x = unit(0.86, "npc"),y = unit(0.985, "npc"))

# Raw parameters
grid.text(expression(paste("f(x; ",lambda,") = ",lambda,e^(-lambda*x))),x = unit(0.25, "npc"),y = unit(0.89, "npc"),gp=gpar(fontsize=text.fs))
grid.text(expression(paste(lambda," = 1.47")),x = unit(0.25, "npc"),y = unit(0.86, "npc"),gp=gpar(fontsize=text.fs))

grid.text(expression(paste("GEV(x; ",mu,",",sigma,",",xi,")")),x = unit(0.55, "npc"),y = unit(0.89, "npc"),gp=gpar(fontsize=text.fs))
grid.text(expression(paste(mu," = 0.09")),x = unit(0.57, "npc"),y = unit(0.86, "npc"),gp=gpar(fontsize=text.fs))
grid.text(expression(paste(sigma," = 0.11")),x = unit(0.57, "npc"),y = unit(0.84, "npc"),gp=gpar(fontsize=text.fs))
grid.text(expression(paste(xi," = 0.82")),x = unit(0.57, "npc"),y = unit(0.82, "npc"),gp=gpar(fontsize=text.fs))

grid.text(expression(paste("GEV(x; ",mu,",",sigma,",",xi,")")),x = unit(0.87, "npc"),y = unit(0.89, "npc"),gp=gpar(fontsize=text.fs))
grid.text(expression(paste(mu," = 0.04")),x = unit(0.89, "npc"),y = unit(0.86, "npc"),gp=gpar(fontsize=text.fs))
grid.text(expression(paste(sigma," = 0.05")),x = unit(0.89, "npc"),y = unit(0.84, "npc"),gp=gpar(fontsize=text.fs))
grid.text(expression(paste(xi," = 0.76")),x = unit(0.89, "npc"),y = unit(0.82, "npc"),gp=gpar(fontsize=text.fs))

# Ln parameters
grid.text(expression(paste(mu," = -0.89")),x = unit(0.18, "npc"),y = unit(0.59, "npc"),gp=gpar(fontsize=text.fs))
grid.text(expression(paste(sigma," = 1.31")),x = unit(0.18, "npc"),y = unit(0.57, "npc"),gp=gpar(fontsize=text.fs))
grid.text(expression(paste(mu," = -1.95")),x = unit(0.49, "npc"),y = unit(0.59, "npc"),gp=gpar(fontsize=text.fs))
grid.text(expression(paste(sigma," = 1.45")),x = unit(0.49, "npc"),y = unit(0.57, "npc"),gp=gpar(fontsize=text.fs))
grid.text(expression(paste(mu," = -2.55")),x = unit(0.80, "npc"),y = unit(0.59, "npc"),gp=gpar(fontsize=text.fs))
grid.text(expression(paste(sigma," = 1.38")),x = unit(0.80, "npc"),y = unit(0.57, "npc"),gp=gpar(fontsize=text.fs))

# BoxCox parameters
grid.text(expression(paste(mu," = -0.77")),x = unit(0.18, "npc"),y = unit(0.28, "npc"),gp=gpar(fontsize=text.fs))
grid.text(expression(paste(sigma," = 1.14")),x = unit(0.18, "npc"),y = unit(0.26, "npc"),gp=gpar(fontsize=text.fs))
grid.text(expression(paste(mu," = -1.69")),x = unit(0.49, "npc"),y = unit(0.28, "npc"),gp=gpar(fontsize=text.fs))
grid.text(expression(paste(sigma," = 1.14")),x = unit(0.49, "npc"),y = unit(0.26, "npc"),gp=gpar(fontsize=text.fs))
grid.text(expression(paste(mu," = -2.17")),x = unit(0.80, "npc"),y = unit(0.28, "npc"),gp=gpar(fontsize=text.fs))
grid.text(expression(paste(sigma," = 1.04")),x = unit(0.80, "npc"),y = unit(0.26, "npc"),gp=gpar(fontsize=text.fs))


grid.text("Helix-end", x = unit(0.54, "npc"),y = unit(0.985, "npc"))
grid.text("Stacked", x = unit(0.86, "npc"),y = unit(0.985, "npc"))

grid.text("Raw", x = unit(0.985, "npc"),y = unit(0.84, "npc"),rot=270)
grid.text("Ln", x = unit(0.985, "npc"),y = unit(0.53, "npc"),rot=270)
grid.text("Box-Cox", x = unit(0.985, "npc"),y = unit(0.21, "npc"),rot=270)

dev.off()

# Normality tests ---------------------------------------------------------
norm.test <- lillie.test
norm_p <- rep(NA,6)
norm_p[1] <- norm.test(ln$up)$p.value
norm_p[2] <- norm.test(ln$ph)$p.value
norm_p[3] <- norm.test(ln$ps)$p.value

norm_p[4] <- norm.test(boxcox$up)$p.value
norm_p[5] <- norm.test(boxcox$ph)$p.value
norm_p[6] <- norm.test(boxcox$ps)$p.value

norm_p <- cbind(c("ln_up","ln_ph","ln_ps","bc_up","bc_ph","bc_ps"),norm_p)
write.csv(norm_p,paste0(root_out_base,"normality_test.csv"),row.names=F,quote=F)

