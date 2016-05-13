## Plot 3 states Bayesian posterior probabilities --------------------------------------------

# Author: Mirko Ledda | Date: 06/17/2015 | v1

# Plot 3 states 1D-posterior probabilities ------------------------------
plot1D_3s <- function(post,fn) {
  # post: 1 by n  - posterior probabilities from `posterior_3s.R`
  # fn :  char    - filename to save the plot
  
  atxt <- 12 # axis text font size
  ptit <- 14 # plot title font size
  
  # Load required libraries
  lib <- c('ggplot2','reshape2','gridExtra','grid')
  pkgTest <- function(x) {
    if (!require(x,character.only = TRUE)) {
      install.packages(x,dep=TRUE)
      require(x,character.only = TRUE) }}
  lapply(lib,pkgTest)
  
  ## Script starts
  # Plot labels
  lab <- c('Unpaired','Helix-end','Stacked')
  
  # reshaspe the data into a data frame
  xy.df <- as.data.frame(cbind(post$x,post$unpaired,post$paired_end,post$paired_stack))
  xy.dfl <- melt(xy.df,id='V1')
  colnames(xy.df) <- c('x','u','e','s')
  
  # plot posterior probabilities 
  p1 <- qplot(x,u, data=xy.df) + ggtitle(lab[1]) + geom_line() +
    xlab('') + ylab('') +
    ylim(0,1) + theme_bw() +
    theme(axis.text=element_text(size=atxt),
          plot.title=element_text(size=ptit,face="bold"))
  p2 <- qplot(x,e, data=xy.df) + ggtitle(lab[2]) + geom_line() +
    xlab('') + ylab('') +
    ylim(0,1) + theme_bw() +
    theme(axis.text=element_text(size=atxt),
          plot.title=element_text(size=ptit,face="bold"))
  p3 <- qplot(x,s, data=xy.df) + ggtitle(lab[3]) + geom_line() +
    xlab('') + ylab('') +
    ylim(0,1) + theme_bw() +
    theme(axis.text=element_text(size=atxt),
          plot.title=element_text(size=ptit,face="bold"))
  
  xy.df <- as.data.frame(cbind(post$x,post$unpaired,post$paired_end+post$paired_stack))
  xy.dfl <- melt(xy.df,id='V1')
  colnames(xy.df) <- c('x','u','p')
  p4 <- qplot(x,p, data=xy.df) + ggtitle('Paired') + geom_line() +
    xlab('') + ylab('') +
    ylim(0,1) + theme_bw() +
    theme(axis.text=element_text(size=atxt),
          plot.title=element_text(size=ptit,face="bold"))
  
  setEPS()
  postscript(file = paste0(fn,'_1D_3s_post_prob.eps'))
  par(cex=1)
  grid.arrange(arrangeGrob(p1,p2, p4, p3,
                           nrow=2, ncol=2,
                           left = textGrob("Posterior probability", rot = 90, vjust = 1),
                           bottom = textGrob("Reactivity score"),
                           right=textGrob(""),
                           top=textGrob("")))
  dev.off()
  
}



# Plot 3 states 2D-posterior probabilities ------------------------------
plot2D_3s <- function(post,fn) {
  # post: n by n  - posterior probabilities from `posterior_3s.R`
  # fn :  char    - filename to save the plot
  
  # Load required libraries
  lib <- c('ggplot2','reshape2','gridExtra')
  pkgTest <- function(x) {
    if (!require(x,character.only = TRUE)) {
      install.packages(x,dep=TRUE)
      require(x,character.only = TRUE) }}
  lapply(lib,pkgTest)
  
  ## Script starts
  # Plot labels
  lab <- c('Unpaired','Helix-end','Stacked')
  
  # post: posterior probabilities
  maxp <- max(post$unpaired,post$paired_end,post$paired_stack)
  
  # compute normalized color levels
  normp <- list()
  normp$u <- (post$unpaired/maxp)
  normp$e <- (post$paired_end/maxp)
  normp$s <- (post$paired_stack/maxp)
  
  # Cluster the posterior probabilities into 3 categories for unpaired, helix-end and stacked
  ix_cluster <- matrix(NA,dim(post$unpaired)[1],dim(post$unpaired)[2])
  for (j in 1:dim(normp$u)[1]) {
    for (k in 1:dim(normp$u)[2]) {
      if (sum(is.na(c(post$unpaired[j,k],post$paired_end[j,k],post$paired_stack[j,k])))<3) {
        ix_cluster[j,k] <- which.max(c(post$unpaired[j,k],post$paired_end[j,k],post$paired_stack[j,k]))
      }
    }
  }
  
  # reshape the data into a data frame
  ix.melt <- melt(ix_cluster)
  image.df <- as.data.frame(cbind(post$x[ix.melt$Var1],post$y[ix.melt$Var2],ix.melt$value))
  colnames(image.df) <- c('x','y','z')
  
  u <- melt(post$unpaired)
  e <- melt(post$paired_end)
  s <- melt(post$paired_stack)
  post_p.df <- as.data.frame(cbind(post$x[u$Var1],post$y[u$Var2],u$value,e$value,s$value))
  colnames(post_p.df) <- c('x','y','u','e','s')
  
  # plot posterior probabilities
  p1 <- qplot(x,y, data=post_p.df, geom="tile", fill=u) + ggtitle(lab[1]) +
    scale_fill_gradient(limits=c(0, 1),low="white", high="red",name="p") +
    xlab('Reactivity') + ylab('Neighbor reactivity')
  
  p2 <- qplot(x,y, data=post_p.df, geom="tile", fill=e) + ggtitle(lab[2]) +
    scale_fill_gradient(limits=c(0, 1),low="white", high="red",name="p") +
    xlab('Reactivity') + ylab('Neighbor reactivity')
  p3 <- qplot(x,y, data=post_p.df, geom="tile", fill=s) + ggtitle(lab[3]) +
    scale_fill_gradient(limits=c(0, 1),low="white", high="red",name="p") +
    xlab('Reactivity') + ylab('Neighbor reactivity')
  
  image.df$Y1 <- cut(image.df$z,breaks = c(1,1.1:2.9,3.1),right = FALSE)
  
  p4 <- ggplot(data=image.df, aes(x=x, y=y)) + ggtitle("Merged") +
    geom_tile(aes(fill=Y1)) +
    xlab('Reactivity') +
    ylab('Neighbor reactivity') +
    scale_fill_discrete(name="Category",label=lab[sort(unique(image.df$z))])
  
  png(filename = paste(fn,'_2D_3s_post_prob.png',sep=''),width = 1200,height = 800)
  par(cex=1.3)
  grid.arrange(p1, p2, p3, p4, nrow=2, ncol=2)
  dev.off()
  
  # RGB continous plot for merged data
  u <- melt(normp$u)
  e <- melt(normp$e)
  s <- melt(normp$s)
  image.df <- as.data.frame(cbind(post$x[u$Var1],post$y[u$Var2],u$value,e$value,s$value))
  colnames(image.df) <- c('x','y','r','g','b')
  col_rgb <- rgb(image.df$r,image.df$g,image.df$b)
  
  p1 <- ggplot(data=image.df, aes(x=x, y=y)) + ggtitle("Merged") +
    geom_tile(fill=col_rgb) +
    xlab('Reactivity') +
    ylab('Neighbor reactivity')
  
  png(filename = paste0(fn,'_2D_3s_post_prob_merged.png'),width = 800,height = 800)
  par(cex=1.3)
  grid.arrange(p1, nrow=1, ncol=1)
  dev.off()
}

# Plot 3 states 3D-posterior probabilities ------------------------------
plot3D_3s <- function(post,fn,data_transform='raw') {
  # post:           array n   - posterior probabilities from `posterior_3s.R`
  # fn :            char      - filename to save the plot
  # data_transform: char      - Type of data transformation. raw | boxcox accepted
  
  # set axis and label parameters
  axis_lim <- c(min(post$x),mean(c(min(post$x),max(post$x))),max(post$x))
  lcolor <- "black"
  lab_alpha <- 0.6
  
  # Load required libraries
  lib <- c('ggplot2','reshape2','gridExtra')
  pkgTest <- function(x) {
    if (!require(x,character.only = TRUE)) {
      install.packages(x,dep=TRUE)
      require(x,character.only = TRUE) }}
  lapply(lib,pkgTest)
  
  ## Script starts
  # plot posterior probabilities --------------------------------------------
  # Compute the colorscale
  colRamp <- colorRamp(colors = c('white','orange','red'))(seq(0,1,0.01))/255
  
  colScale <- rep(NA,dim(colRamp)[1])
  for (j in 1:dim(colRamp)[1]) {
    colScale[j] <- rgb(red=colRamp[j,1],green=colRamp[j,2],blue=colRamp[j,3])
  }
  
  if (data_transform == 'raw') {
    iv <- seq(0,2,1)} else if (data_transform == 'boxcox') {
      iv <- seq(-3,1,1)}
  
  # Unpaired ---------------------------------------------------
  # Start 3D plotting
  rgl.open() 
  rgl.bg(color="white")
  
  rgl.bbox(
    xat = iv, xlab = iv, xunit = 0, xlen = length(iv),
    yat = iv, ylab = iv, yunit = 0, ylen = length(iv),
    zat = iv, zlab = iv, zunit = 0, zlen = length(iv),
    marklen.rel = TRUE,col=lcolor,expand = 1.1,alpha=lab_alpha)
  
  # plot the data
  image3d(post$unpaired,x=post$x,y=post$y,z=post$z,jitter=F,vlim=c(0,1),add=T,col=colScale,alpha=0.5)
  
  # Change the angle of the plot
  rgl.viewpoint(theta=310,phi=30,zoom = 1)
  
  # write axis label and main title
  rgl.texts(x=axis_lim[2],y=axis_lim[1],z=axis_lim[3],adj=c(-4,3),text='x',col=lcolor,alpha=lab_alpha)
  rgl.texts(x=axis_lim[1],y=axis_lim[2],z=axis_lim[1],adj=c(6,1),text='y',col=lcolor,alpha=lab_alpha)
  rgl.texts(x=axis_lim[1],y=axis_lim[1],z=axis_lim[2],adj=c(5,3),text='z',col=lcolor,alpha=lab_alpha)
  rgl.texts(x=axis_lim[3],y=axis_lim[3],z=axis_lim[1],adj=c(0.5,-2),text='Unpaired',col=lcolor,alpha=lab_alpha)
  
  # take a screenshot of the image and save it
  rgl.snapshot(filename=paste(fn,'_unpaired.png',sep=''), fmt = "png", top = TRUE )
  rgl.close()
  
  # Helix-end ---------------------------------------------------
  # Start 3D plotting
  rgl.open() 
  rgl.bg(color="white")
  
  rgl.bbox(
    xat = iv, xlab = iv, xunit = 0, xlen = length(iv),
    yat = iv, ylab = iv, yunit = 0, ylen = length(iv),
    zat = iv, zlab = iv, zunit = 0, zlen = length(iv),
    marklen.rel = TRUE,col=lcolor,expand = 1.1,alpha=lab_alpha)
  
  # plot the data
  image3d(post$paired_end,x=post$x,y=post$y,z=post$z,jitter=F,vlim=c(0,1),add=T,col=colScale,alpha=0.5)
  
  # Change the angle of the plot
  rgl.viewpoint(theta=310,phi=30,zoom = 1)
  
  # write axis label and main title
  rgl.texts(x=axis_lim[2],y=axis_lim[1],z=axis_lim[3],adj=c(-4,3),text='x',col=lcolor,alpha=lab_alpha)
  rgl.texts(x=axis_lim[1],y=axis_lim[2],z=axis_lim[1],adj=c(6,1),text='y',col=lcolor,alpha=lab_alpha)
  rgl.texts(x=axis_lim[1],y=axis_lim[1],z=axis_lim[2],adj=c(5,3),text='z',col=lcolor,alpha=lab_alpha)
  rgl.texts(x=axis_lim[3],y=axis_lim[3],z=axis_lim[1],adj=c(0.5,-2),text='Helix-end',col=lcolor,alpha=lab_alpha)
  
  # take a screenshot of the image and save it
  rgl.snapshot(filename=paste(fn,'_helix_end.png',sep=''), fmt = "png", top = TRUE )
  rgl.close()
  
  # Stacked ---------------------------------------------------
  # Start 3D plotting
  rgl.open() 
  rgl.bg(color="white")
  
  rgl.bbox(
    xat = iv, xlab = iv, xunit = 0, xlen = length(iv),
    yat = iv, ylab = iv, yunit = 0, ylen = length(iv),
    zat = iv, zlab = iv, zunit = 0, zlen = length(iv),
    marklen.rel = TRUE,col=lcolor,expand = 1.1,alpha=lab_alpha)
  
  # plot the data
  image3d(post$paired_stack,x=post$x,y=post$y,z=post$z,jitter=F,vlim=c(0,1),add=T,col=colScale,alpha=0.5)
  
  # Change the angle of the plot
  rgl.viewpoint(theta=310,phi=30,zoom = 1)
  
  # write axis label and main title
  rgl.texts(x=axis_lim[2],y=axis_lim[1],z=axis_lim[3],adj=c(-4,3),text='x',col=lcolor,alpha=lab_alpha)
  rgl.texts(x=axis_lim[1],y=axis_lim[2],z=axis_lim[1],adj=c(6,1),text='y',col=lcolor,alpha=lab_alpha)
  rgl.texts(x=axis_lim[1],y=axis_lim[1],z=axis_lim[2],adj=c(5,3),text='z',col=lcolor,alpha=lab_alpha)
  rgl.texts(x=axis_lim[3],y=axis_lim[3],z=axis_lim[1],adj=c(0.5,-2),text='Stacked',col=lcolor,alpha=lab_alpha)
  
  # take a screenshot of the image and save it
  rgl.snapshot(filename=paste0(fn,'_stacked.png'), fmt = "png", top = TRUE )
  rgl.close()
}


# Plot 3 states 2D-posterior probabilities ------------------------------
# Plot image mapping of the 3 state reactivity probabilities using a RGB-based scale
plot_diff_3s <- function(post,fn) {
  # post: n by n  - posterior probabilities from `posterior_3s.R`
  # fn :  char    - filename to save the image
  
  lib <- c('ggplot2','reshape2','gridExtra')
  pkgTest <- function(x) {
    if (!require(x,character.only = TRUE)) {
      install.packages(x,dep=TRUE)
      require(x,character.only = TRUE) }}
  lapply(lib,pkgTest)
  
  ## Script starts

    # Plot labels
  lab <- c('Unpaired','Helix-end','Stacked')
  
  u <- melt(post$unpaired)
  e <- melt(post$paired_end)
  s <- melt(post$paired_stack)
  post_p.df <- as.data.frame(cbind(post$x[u$Var1],post$y[u$Var2],u$value,e$value,s$value))
  colnames(post_p.df) <- c('x','y','u','e','s')
  
  # plot posterior probabilities
  p1 <- qplot(x,y, data=post_p.df, geom="tile", fill=u) + ggtitle(lab[1]) +
    scale_fill_gradient2(limits=c(-1, 1),name='diff')+
    xlab('Reactivity') + ylab('Neighbor reactivity')
  
  p2 <- qplot(x,y, data=post_p.df, geom="tile", fill=e) + ggtitle(lab[2]) +
    scale_fill_gradient2(limits=c(-1, 1),name='diff')+
    xlab('Reactivity') + ylab('Neighbor reactivity')
  p3 <- qplot(x,y, data=post_p.df, geom="tile", fill=s) + ggtitle(lab[3]) +
    scale_fill_gradient2(limits=c(-1, 1),name='diff')+
    xlab('Reactivity') + ylab('Neighbor reactivity')
  
  png(filename = paste0(fn,'_diff_3s_post_prob.png'),width = 1200,height = 800)
  par(cex=1.3)
  grid.arrange(p1, p2, p3, nrow=2, ncol=2)
  dev.off()
  
}



