# Write SHAPE files -------------------------------------------------------

# Author: Mirko Ledda | Date: 11/10/2015 | v1

write_shape <- function(fn,x) {
  ## fn:  char    - filename. ".shape" will be used as extension
  ## x:   1 by n  - shape reactivities
  
  n <- length(x)
  
  fn_reac <- paste0(fn,'.shape')
  rdf <- data.frame("n"=seq(1,n),'SHAPE'=x)
  rdf[is.na.data.frame(rdf)] <- -999
  write.table(rdf,file=fn_reac,sep=" ",row.names = F, col.names = F)
}
