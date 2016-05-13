# Compute kernel probability densities for 3 state reactivities ------------------------------

# Author: Mirko Ledda | Date: 06/16/2015 | v1
#
# Based on the `kde` function from the package `ks`

kernel_3s <- function(x,xp,dims=1,type='n5',nbin=50,limits=c(-0.5,1.5),flag=0,h='ls',cv=0,min.h=0.01) {
  # x:        1 by n  - Reactivities
  # xp:       list    - output of the pairing function
  # dims:     integer - Number of dimensions of the kernel
  # type:     char    - n3 = 3'-neighbor | n5 = 5'-neighbor | nmin = neighbor's min value | nmax = neighbor's max value
  # nbin:     integer - number of bins for x and y
  # limits:   1 by 2  - axis limits (same for x and y)
  # flag:     boolean - determine if the data where concatenated (1=yes, 0=no)
  # h:        float   - bandwith and covariance of the gaussian kernel
  # cv:       matrix  - covariance of the gaussian kernel
  
  ## Install or load required libraries
  lib <- c('ks')
  pkgTest <- function(x) {
    if (!require(x,character.only = TRUE)) {
      install.packages(x,dep=TRUE)
      require(x,character.only = TRUE) }}
  lapply(lib,pkgTest)
  
  ## Script starts
  ker <- list()
  
  iv <- seq(limits[1],limits[2],(limits[2]-limits[1])/nbin)
  
  # Call the flag if it exist
  if (flag==1) {
    xp$pairing[xp$flag] <- NA
  }
  
  ### Selection for the dimension
  if (dims==1) {
    ## 1D kernel
    H <- h
    
    ix <- list()
    ix$u <- which(xp$pairing==0)
    ix$e <- which(xp$pairing==1)
    ix$s <- which(xp$pairing==2)
    
    xk <- x[ix$u]
    filt <- !is.na(xk)
    ker$uc <- sum(filt)
    if (h=='ls' | h=='lsc') {
      H <- max(hlscv(unique(xk[filt])),min.h)
    }
    ker$u <- kde(xk[filt],xmin = limits[1],xmax=limits[2],gridsize = nbin,h=H)
    
    xk <- x[ix$e]
    filt <- !is.na(xk)
    ker$ec <- sum(filt)
    if (h=='ls' | h=='lsc') {
      H <- max(hlscv(unique(xk[filt])),min.h)
    }
    ker$e <- kde(xk[filt],xmin = limits[1],xmax=limits[2],gridsize = nbin,h=H)
    
    xk <- x[ix$s]
    filt <- !is.na(xk)
    ker$sc <- sum(filt)
    if (h=='ls' | h=='lsc') {
      H <- max(hlscv(unique(xk[filt])),min.h)
    }
    ker$s <- kde(xk[filt],xmin = limits[1],xmax=limits[2],gridsize = nbin,h=H)
    
  } else if (dims==2) {
    ## 2D kernel
    H <- matrix(c(h,cv,cv,h),nrow = 2)
    # select the kernel density type
    switch(type,
           # 5'-neighbor model -------------------------------------------------------
           n5={n<-1
               xp$pairing[n] <- NA
               # find pairing indexes
               ix <- list()
               ix$u <- which(xp$pairing==0)
               ix$e <- which(xp$pairing==1)
               ix$s <- which(xp$pairing==2)
               
               xk <- x[ix$u]
               yk <- x[ix$u-1]
               filt <- (!is.na(xk) & !is.na(yk))
               ker$uc <- sum(filt)
               if (h=='ls'){
                H <- Hlscv.diag(cbind(xk[filt],yk[filt]))
               } else if (h =='lsc') {
                 H <- Hlscv(cbind(xk[filt],yk[filt]))
               }
               ker$u <- kde(cbind(xk[filt],yk[filt]),xmin = rep(limits[1],2),xmax=rep(limits[2],2),gridsize = rep(nbin,2),H=H)
               
               xk <- x[ix$e]
               yk <- x[ix$e-1]
               filt <- (!is.na(xk) & !is.na(yk))
               ker$ec <- sum(filt)
               if (h=='ls'){
                 H <- Hlscv.diag(cbind(xk[filt],yk[filt]))
               }  else if (h =='lsc') {
                 H <- Hlscv(cbind(xk[filt],yk[filt]))
               }
               ker$e <- kde(cbind(xk[filt],yk[filt]),xmin = rep(limits[1],2),xmax=rep(limits[2],2),gridsize = rep(nbin,2),H=H)
               
               xk <- x[ix$s]
               yk <- x[ix$s-1]
               filt <- (!is.na(xk) & !is.na(yk))
               ker$sc <- sum(filt)
               if (h=='ls'){
                 H <- Hlscv.diag(cbind(xk[filt],yk[filt]))
               }  else if (h =='lsc') {
                 H <- Hlscv(cbind(xk[filt],yk[filt]))
               }
               ker$s <- kde(cbind(xk[filt],yk[filt]),xmin = rep(limits[1],2),xmax=rep(limits[2],2),gridsize = rep(nbin,2),H=H)},
           # 3'-neighbor model -------------------------------------------------------
           n3={n<-length(xp$pairing)
               xp$pairing[n] <- NA
               # find pairing indexes
               ix <- list()
               ix$u <- which(xp$pairing==0)
               ix$e <- which(xp$pairing==1)
               ix$s <- which(xp$pairing==2)
               
               xk <- x[ix$u]
               yk <- x[ix$u+1]
               filt <- (!is.na(xk) & !is.na(yk))
               ker$uc <- sum(filt)
               if (h=='ls'){
                 H <- Hlscv.diag(cbind(xk[filt],yk[filt]))
               } else if (h =='lsc') {
                 H <- Hlscv(cbind(xk[filt],yk[filt]))
               }
               ker$u <- kde(cbind(xk[filt],yk[filt]),xmin = rep(limits[1],2),xmax=rep(limits[2],2),gridsize = rep(nbin,2),H=H)
               
               xk <- x[ix$e]
               yk <- x[ix$e+1]
               filt <- (!is.na(xk) & !is.na(yk))
               ker$ec <- sum(filt)
               if (h=='ls'){
                 H <- Hlscv.diag(cbind(xk[filt],yk[filt]))
               } else if (h =='lsc') {
                 H <- Hlscv(cbind(xk[filt],yk[filt]))
               }
               ker$e <- kde(cbind(xk[filt],yk[filt]),xmin = rep(limits[1],2),xmax=rep(limits[2],2),gridsize = rep(nbin,2),H=H)
               
               xk <- x[ix$s]
               yk <- x[ix$s+1]
               filt <- (!is.na(xk) & !is.na(yk))
               ker$sc <- sum(filt)
               if (h=='ls'){
                 H <- Hlscv.diag(cbind(xk[filt],yk[filt]))
               } else if (h =='lsc') {
                 H <- Hlscv(cbind(xk[filt],yk[filt]))
               }
               ker$s <- kde(cbind(xk[filt],yk[filt]),xmin = rep(limits[1],2),xmax=rep(limits[2],2),gridsize = rep(nbin,2),H=H)},
           # Min neighbor model -------------------------------------------------------
           nmin={n<-c(1,length(xp$pairing))
                 xp$pairing[n] <- NA
                 # find pairing indexes
                 ix <- list()
                 ix$u <- which(xp$pairing==0)
                 ix$e <- which(xp$pairing==1)
                 ix$s <- which(xp$pairing==2)
                 
                 xk <- x[ix$u]
                 yk <- rowMin(cbind(x[ix$u-1],x[ix$u-1]))
                 filt <- (!is.na(xk) & !is.na(yk))
                 ker$uc <- sum(filt)
                 if (h=='ls'){
                   H <- Hlscv.diag(cbind(xk[filt],yk[filt]))
                 } else if (h =='lsc') {
                   H <- Hlscv(cbind(xk[filt],yk[filt]))
                 }
                 ker$u <- kde(cbind(xk[filt],yk[filt]),xmin = rep(limits[1],2),xmax=rep(limits[2],2),gridsize = rep(nbin,2),H=H)
                 
                 xk <- x[ix$e]
                 yk <- rowMin(cbind(x[ix$e-1],x[ix$e-1]))
                 filt <- (!is.na(xk) & !is.na(yk))
                 ker$ec <- sum(filt)
                 if (h=='ls'){
                   H <- Hlscv.diag(cbind(xk[filt],yk[filt]))
                 } else if (h =='lsc') {
                   H <- Hlscv(cbind(xk[filt],yk[filt]))
                 }
                 ker$e <- kde(cbind(xk[filt],yk[filt]),xmin = rep(limits[1],2),xmax=rep(limits[2],2),gridsize = rep(nbin,2),H=H)
                 
                 xk <- x[ix$s]
                 yk <- rowMin(cbind(x[ix$s-1],x[ix$s-1]))
                 filt <- (!is.na(xk) & !is.na(yk))
                 ker$sc <- sum(filt)
                 if (h=='ls'){
                   H <- Hlscv.diag(cbind(xk[filt],yk[filt]))
                 } else if (h =='lsc') {
                   H <- Hlscv(cbind(xk[filt],yk[filt]))
                 }
                 ker$s <- kde(cbind(xk[filt],yk[filt]),xmin = rep(limits[1],2),xmax=rep(limits[2],2),gridsize = rep(nbin,2),H=H)},
           # Max neighbor model -------------------------------------------------------
           nmax={n<-c(1,length(xp$pairing))
                 xp$pairing[n] <- NA
                 # find pairing indexes
                 ix <- list()
                 ix$u <- which(xp$pairing==0)
                 ix$e <- which(xp$pairing==1)
                 ix$s <- which(xp$pairing==2)
                 
                 xk <- x[ix$u]
                 yk <- rowMax(cbind(x[ix$u-1],x[ix$u-1]))
                 filt <- (!is.na(xk) & !is.na(yk))
                 ker$uc <- sum(filt)
                 if (h=='ls'){
                   H <- Hlscv.diag(cbind(xk[filt],yk[filt]))
                 } else if (h =='lsc') {
                   H <- Hlscv(cbind(xk[filt],yk[filt]))
                 }
                 ker$u <- kde(cbind(xk[filt],yk[filt]),xmin = rep(limits[1],2),xmax=rep(limits[2],2),gridsize = rep(nbin,2),H=H)
                 
                 xk <- x[ix$e]
                 yk <- rowMax(cbind(x[ix$e-1],x[ix$e-1]))
                 filt <- (!is.na(xk) & !is.na(yk))
                 ker$ec <- sum(filt)
                 if (h=='ls'){
                   H <- Hlscv.diag(cbind(xk[filt],yk[filt]))
                 } else if (h =='lsc') {
                   H <- Hlscv(cbind(xk[filt],yk[filt]))
                 }
                 ker$e <- kde(cbind(xk[filt],yk[filt]),xmin = rep(limits[1],2),xmax=rep(limits[2],2),gridsize = rep(nbin,2),H=H)
                 
                 xk <- x[ix$s]
                 yk <- rowMax(cbind(x[ix$s-1],x[ix$s-1]))
                 filt <- (!is.na(xk) & !is.na(yk))
                 ker$sc <- sum(filt)
                 if (h=='ls'){
                   H <- Hlscv.diag(cbind(xk[filt],yk[filt]))
                 } else if (h =='lsc') {
                   H <- Hlscv(cbind(xk[filt],yk[filt]))
                 }
                 ker$s <- kde(cbind(xk[filt],yk[filt]),xmin = rep(limits[1],2),xmax=rep(limits[2],2),gridsize = rep(nbin,2),H=H)}
    )
  } else if (dims==3) {
    H <- matrix(c(h,cv,cv,cv,h,cv,cv,cv,h),nrow = 3)
    
    n<-c(1,length(xp$pairing))
    xp$pairing[n] <- NA
    # find pairing indexes
    ix <- list()
    ix$u <- which(xp$pairing==0)
    ix$e <- which(xp$pairing==1)
    ix$s <- which(xp$pairing==2)
    
    xk <- x[ix$u]
    yk <- x[ix$u-1]
    zk <- x[ix$u+1]
    filt <- (!is.na(xk) & !is.na(yk) & !is.na(zk))
    ker$uc <- sum(filt)
    if (h=='ls'){
      H <- Hlscv.diag(cbind(xk[filt],yk[filt],zk[filt]))
    } else if (h =='lsc') {
      H <- Hlscv(cbind(xk[filt],yk[filt]))
    }
    ker$u <- kde(cbind(xk[filt],yk[filt],zk[filt]),xmin = rep(limits[1],3),xmax=rep(limits[2],3),gridsize = rep(nbin,3),H=H)
    
    xk <- x[ix$e]
    yk <- x[ix$e-1]
    zk <- x[ix$e+1]
    filt <- (!is.na(xk) & !is.na(yk) & !is.na(zk))
    ker$ec <- sum(filt)
    if (h=='ls'){
      H <- Hlscv.diag(cbind(xk[filt],yk[filt],zk[filt]))
    } else if (h =='lsc') {
      H <- Hlscv(cbind(xk[filt],yk[filt]))
    }
    ker$e <- kde(cbind(xk[filt],yk[filt],zk[filt]),xmin = rep(limits[1],3),xmax=rep(limits[2],3),gridsize = rep(nbin,3),H=H)
    
    xk <- x[ix$s]
    yk <- x[ix$s-1]
    zk <- x[ix$s+1]
    filt <- (!is.na(xk) & !is.na(yk) & !is.na(zk))
    ker$sc <- sum(filt)
    if (h=='ls'){
      H <- Hlscv.diag(cbind(xk[filt],yk[filt],zk[filt]))
    } else if (h =='lsc') {
      H <- Hlscv(cbind(xk[filt],yk[filt]))
    }
    ker$s <- kde(cbind(xk[filt],yk[filt],zk[filt]),xmin = rep(limits[1],3),xmax=rep(limits[2],3),gridsize = rep(nbin,3),H=H)
  }
  
  return(ker)
}


