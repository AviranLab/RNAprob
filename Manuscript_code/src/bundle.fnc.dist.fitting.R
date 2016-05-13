## Density function fitting

# Author: Mirko Ledda | Date: 03/03/2015 | v1
# WARNING: By default the number of bins in the histogram is defined by pretty(x,20)

#%%%%%%%%%%%%%%%%%%%%%%%%%%#
# List of fitted functions #
#%%%%%%%%%%%%%%%%%%%%%%%%%%#
# MODEL TYPE                          - ABRV    - FUNCTION    - FITTING TYPE
# 1. Normal (Gaussian)                - normal  - CFitNormal  - analogic
# 2. Gamma                            - gamma   - CFitGamma   - analogic
# 3. Exponential                      - exp     - CFitExp     - MLE
# 4. Generalized Extreme Value (GEV)  - GEV     - CFitGEV     - MLE (uses `evd`)
# 5. Pearson Distributions            - pearson - CFitPearson - Moments (uses `PearsonDS`)
#%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Set main functions
pkgTest <- function(x) {
  if (!require(x,character.only = TRUE)) {
    install.packages(x,dep=TRUE)
    require(x,character.only = TRUE) }}

JSD <- function(p, q) {
  # p: 1 by n - vector of data
  # q: 1 by n - vector of data
  
  # filter 0 out
  filt <- p>0 & q>0
  p <- p[filt]
  q <- q[filt]
  m <- 0.5 * (p + q)
  return(0.5 * (sum(p * log(p / m)) + sum(q * log(q / m))))
}

CFit_setX <- function(x,iv,n) {
  out <- list()
  
  lim <- range(x)
  
  if (is.null(iv)) {
    if (is.null(n)) {
      prty <- pretty(x,n=20)
      lim <- range(prty)
      n <- length(prty)-1
      iv <- diff(lim)/n
    } else {
      iv <- diff(lim)/n
    }
  }
  
  # Set x vectors
  out$hbreaks <- seq(lim[1],lim[2],iv)
  out$x0 <- seq(lim[1]+iv/2,lim[2],iv)
  
  return(out)
  
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# 1. Normal
CFitNormal <- function(x,iv=NULL,n=NULL,plot.flag=F) {
  # x:          1 by n  - Vector of data
  # iv:         Float   - Interval for the histogram bars (!!! prevales against n !!!)
  # n:          Integer - Number of histogram bars
  # plot.flag:  Boolean - Plot the fitting?
  
  out <- list()
  
  pdf <- dnorm
  cdf <- pnorm
  
  x[is.infinite(x)] <- NA
  x <- x[!is.na(x)]
  
  # Parameter estimation
  out$mean <- mean(x)
  out$sd <- sd(x)
 
  # Create the x-axis vector
  setX <- CFit_setX(x,iv,n)
  
  # Predicted values
  out$y <- pdf(setX$x0,mean = out$mean,sd = out$sd)
  
  # Histogram
  h <- hist(x,freq = F,breaks = setX$hbreaks,main=NULL,plot = plot.flag)
  if (plot.flag) lines(setX$x0,out$y,col="red",lwd=3)
  
  # Fitting statistics
  out$KS.p <- ks.test(x,cdf,mean = out$mean,sd = out$sd)$p.value # One-sample Kolmogorov-Smirnov test
  out$JSD <- JSD(out$y,h$density) # Jensen-Shannon divergence
  
  return(out)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# 2. Gamma
CFitGamma <- function(x,iv=NULL,n=NULL,plot.flag=F) {
  # x:          1 by n  - Vector of data
  # iv:         Float   - Interval for the histogram bars (!!! prevales against n !!!)
  # n:          Integer - Number of histogram bars
  # plot.flag:  Boolean - Plot the fitting?
  
  out <- list()
  
  pdf <- dgamma
  cdf <- pgamma
  
  x[is.infinite(x)] <- NA
  x <- x[!is.na(x)]
  
  # Parameter estimation
  mean.gam <- mean(x)
  var.gam <- var(x)
  
  out$rate <- mean.gam/var.gam # lambda
  out$shape <- (mean.gam^2)/var.gam # alpha
  
  # Create the x-axis vector
  setX <- CFit_setX(x,iv,n)
  
  # Predicted values
  out$y <- pdf(setX$x0,rate = out$rate,shape=out$shape)
  
  # Histogram
  h <- hist(x,freq = F,breaks = setX$hbreaks,main=NULL,plot = plot.flag)
  if (plot.flag) lines(setX$x0,out$y,col="red",lwd=3)
  
  # Fitting statistics
  out$KS.p <- ks.test(x,cdf,rate = out$rate,shape=out$shape)$p.value # One-sample Kolmogorov-Smirnov test
  out$JSD <- JSD(out$y,h$density) # Jensen-Shannon divergence
  
  return(out)
  
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
 
# 3. Exponential function
CFitExp <- function(x,iv=NULL,n=NULL,plot.flag=F) {
  # x:          1 by n  - Vector of data
  # iv:         Float   - Interval for the histogram bars (!!! prevales against n !!!)
  # n:          Integer - Number of histogram bars
  # plot.flag:  Boolean - Plot the fitting?
  
  out <- list()
  
  pdf <- dexp
  cdf <- pexp

  x[is.infinite(x)] <- NA
  x <- x[!is.na(x)]
  
  # Parameter estimation (by MLE)
  xfilt <- x>0
  x2 <- x[xfilt]
  out$lambda <- length(x2)/sum(x2)
  
  # Create the x-axis vector
  setX <- CFit_setX(x,iv,n)
  
  # Predicted values
  out$y <- pdf(setX$x0,rate = out$lambda)
  
  # Histogram
  h <- hist(x,freq = F,breaks = setX$hbreaks,main=NULL,plot = plot.flag)
  if (plot.flag) lines(setX$x0,out$y,col="red",lwd=3)
  
  # Fitting statistics
  out$KS.p <- ks.test(x,cdf,rate = out$lambda)$p.value # One-sample Kolmogorov-Smirnov test
  out$JSD <- JSD(out$y,h$density) # Jensen-Shannon divergence
  
  return(out)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# 4. Generalized extreme value
CFitGEV <- function(x,iv=NULL,n=NULL,plot.flag=F) {
  # x:          1 by n  - Vector of data
  # iv:         Float   - Interval for the histogram bars (!!! prevales against n !!!)
  # n:          Integer - Number of histogram bars
  # plot.flag:  Boolean - Plot the fitting?
  
  lapply("evd",pkgTest)
  out <- list()
  
  pdf <- dgev
  cdf <- pgev
  
  x[is.infinite(x)] <- NA
  x <- x[!is.na(x)]
  
  # Parameter estimation (by MLE)
  FitEst <- fgev(x,std.err = F)
  out$mu <- FitEst$estimate[1]
  out$sigma <- FitEst$estimate[2]
  out$xi <- FitEst$estimate[3]
  
  # Create the x-axis vector
  setX <- CFit_setX(x,iv,n)
  
  # Predicted values
  out$y <- pdf(setX$x0,loc = out$mu,scale = out$sigma,shape=out$xi)
  
  # Histogram
  h <- hist(x,freq = F,breaks = setX$hbreaks,main=NULL,plot = plot.flag)
  if (plot.flag) lines(setX$x0,out$y,col="red",lwd=3)
  
  # Fitting statistics
  out$KS.p <- ks.test(x,cdf,loc = out$mu,scale = out$sigma,shape=out$xi)$p.value # One-sample Kolmogorov-Smirnov test
  out$JSD <- JSD(out$y,h$density) # Jensen-Shannon divergence
  
  return(out)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# 5. Pearson distribution families
CFitPearson <- function(x,iv=NULL,n=NULL,plot.flag=F) {
  # x:          1 by n  - Vector of data
  # iv:         Float   - Interval for the histogram bars (!!! prevales against n !!!)
  # n:          Integer - Number of histogram bars
  # plot.flag:  Boolean - Plot the fitting?
  
  lapply("PearsonDS",pkgTest)
  out <- list()
  
  pdf <- dpearson
  cdf <- ppearson
  
  x[is.infinite(x)] <- NA
  x <- x[!is.na(x)]
  
  # Parameter estimation
  out$moments <- empMoments(x)
  
  # Create the x-axis vector
  setX <- CFit_setX(x,iv,n)
  
  # Predicted values
  out$y <- pdf(setX$x0,moments = out$moments)
  
  # Histogram
  h <- hist(x,freq = F,breaks = setX$hbreaks,main=NULL,plot = plot.flag)
  if (plot.flag) lines(setX$x0,out$y,col="red",lwd=3)
  
  # Fitting statistics
  out$KS.p <- ks.test(x,cdf,moments = out$moments)$p.value # One-sample Kolmogorov-Smirnov test
  out$JSD <- JSD(out$y,h$density) # Jensen-Shannon divergence
  
  return(out)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Fit all available distributions and plot fitting statistics to choose the best model

CFit <- function(x,iv=NULL,n=NULL,model=NULL) {
  # x:          1 by n  - Vector of data
  # iv:         Float   - Interval for the histogram bars (!!! prevales against n !!!)
  # n:          Integer - Number of histogram bars
  # model:      Vector  - Vector of string specifying tested models
  
  lapply("ggplot2",pkgTest)
  out <- list()
  
  fitKSp <- list()
  fitJSD <- list()
  
  # Available models
  if (is.null(model)) model <- c("normal","gamma","exp","GEV","pearson")
  
  # 1. Normal
  if (any(model=="normal")) {
    fitstat <- CFitNormal(x,iv,n)
    fitKSp$normal <- fitstat$KS.p
    fitJSD$normal <- fitstat$JSD  
  }
  
  # 2. Gamma
  if (any(model=="gamma")) {
    fitstat <- CFitGamma(x,iv,n)
    fitKSp$gamma <- fitstat$KS.p
    fitJSD$gamma <- fitstat$JSD  
  }
  
  # 3. Exponential
  if (any(model=="exp")) {
    fitstat <- CFitExp(x,iv,n)
    fitKSp$exp <- fitstat$KS.p
    fitJSD$exp <- fitstat$JSD  
  }
  
  # 4. Generalized Extreme Value (GEV)
  if (any(model=="GEV")) {
    fitstat <- CFitGEV(x,iv,n)
    fitKSp$GEV <- fitstat$KS.p
    fitJSD$GEV <- fitstat$JSD  
  }
  
  # 5. Pearson Distributions
  if (any(model=="pearson")) {
    fitstat <- CFitPearson(x,iv,n)
    fitKSp$pearson <- fitstat$KS.p
    fitJSD$pearson <- fitstat$JSD  
  }
  
  # Plot the fitting results
  df <- data.frame(
    x = unlist(fitKSp),
    y = unlist(fitJSD),
    model = model
  )
  
  ggplot(data=df,aes(x=x,y=y,label=model)) +
    xlim(0,1.2) + xlab("KS-test p-value") + ylab("JS distance") +
    geom_segment(x=0,y=0,xend=1,yend=0,linetype=2) +
    geom_segment(x=0,y=0,xend=0,yend=Inf,linetype=2) +
    geom_segment(x=1,y=0,xend=1,yend=Inf,linetype=2) +
    geom_point(x=1,y=0,col="red",size=2,pch=15) +
    expand_limits(y=0) +
    geom_label(size=5,hjust = 0, nudge_x = 0.02) +
    geom_point() +
    #scale_y_reverse() +
    theme_bw()
  
}
