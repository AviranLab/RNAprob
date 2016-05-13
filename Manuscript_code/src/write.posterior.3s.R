# Write 3 states posterior probabilities to a file ------------------------------

# Author: Mirko Ledda | Date: 06/17/2015 | v1
#
# Based on the custom build `posterior.3s` script

write.posterior.3s <- function(post,dims=1,fn) {
  # post: n by dims - posterior probabilities from `posterior.3s.R`
  # fn :  char      - filename to save the data
  # dims: integer   - Number of dimensions of the kernel
  
  ## Install or load required libraries
  lib <- c('ks')
  pkgTest <- function(x) {
    if (!require(x,character.only = TRUE)) {
      install.packages(x,dep=TRUE)
      require(x,character.only = TRUE) }}
  lapply(lib,pkgTest)
  
  ## Script starts
  if (dims == 1) {
    result <- cbind(post$x,post$unpaired,post$paired_end,post$paired_stack)
    colnames(result) <- c('x','unpaired','helix-end','stacked')
    write.table(x = result,row.names = F,col.names = T,quote = F,sep=' ',
                file = paste(fn,'_1D_posterior.density',sep=''))
  } else if (dims==2) {
    # reshape the data
    u <- melt(post$unpaired)
    e <- melt(post$paired_end)
    s <- melt(post$paired_stack)
    result <- as.data.frame(cbind(post$x[u$Var1],post$y[u$Var2],u$value,e$value,s$value))
    colnames(result) <- c('x','y','unpaired','helix-end','stacked')
    
    # write table
    write.table(x = result,row.names = F,col.names = T,quote = F,sep=' ',
                file = paste(fn,'_2D_posterior.density',sep=''))
  } else if (dims==3) {
    # reshape the data
    u <- melt(post$unpaired)
    e <- melt(post$paired_end)
    s <- melt(post$paired_stack)
    result <- as.data.frame(cbind(post$x[u$Var1],post$y[u$Var2],post$z[u$Var3],u$value,e$value,s$value))
    colnames(result) <- c('x','y','z','unpaired','helix-end','stacked')
    
    # write table
    write.table(x = result,row.names = F,col.names = T,quote = F,sep=' ',
                file = paste0(fn,'_3D_posterior.density'))
  }

}


