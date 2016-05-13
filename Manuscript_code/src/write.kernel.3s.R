# Write 3 states kernel densities to a file ------------------------------

# Author: Mirko Ledda | Date: 06/29/2015 | v1
#
# Based on the custom build `kernel.3s` script

write.kernel.3s <- function(ker,dims=1,fn) {
  # ker:  n by dims - kernel densities from `kernel.3s.R`
  # fn :  char      - filename to save the data
  # dims: integer   - Number of dimensions of the kernel
  
  ## Install or load required libraries
  lib <- c('reshape2')
  pkgTest <- function(x) {
    if (!require(x,character.only = TRUE)) {
      install.packages(x,dep=TRUE)
      require(x,character.only = TRUE) }}
  lapply(lib,pkgTest)
  
  ## Script starts
  if (dims == 1) {
    result <- cbind(ker$u$eval.points,ker$u$estimate,ker$e$estimate,ker$s$estimate)
    colnames(result) <- c('x','unpaired','helix-end','stacked')
    write.table(x = result,row.names = F,col.names = T,quote = F,sep=' ',
                file = paste0(fn,'_1D_ker.density'))
  } else if (dims==2) {
    # reshape the data
    u <- melt(ker$u$estimate)
    e <- melt(ker$e$estimate)
    s <- melt(ker$s$estimate)
    result <- as.data.frame(cbind(ker$u$eval.points[[1]][u$Var1],ker$u$eval.points[[2]][u$Var2],u$value,e$value,s$value))
    colnames(result) <- c('x','y','unpaired','helix-end','stacked')
    
    # write table
    write.table(x = result,row.names = F,col.names = T,quote = F,sep=' ',
                file = paste0(fn,'_2D_ker.density'))
  } else if (dims==3) {
    # reshape the data
    u <- melt(ker$u$estimate)
    e <- melt(ker$e$estimate)
    s <- melt(ker$s$estimate)
    result <- as.data.frame(cbind(ker$u$eval.points[[1]][u$Var1],ker$u$eval.points[[2]][u$Var2],ker$u$eval.points[[3]][u$Var3],u$value,e$value,s$value))
    colnames(result) <- c('x','y','z','unpaired','helix-end','stacked')
    
    # write table
    write.table(x = result,row.names = F,col.names = T,quote = F,sep=' ',
                file = paste0(fn,'_3D_ker.density'))
  }

}


