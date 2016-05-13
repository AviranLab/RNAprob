# Compute Bayesian posterior probabilities for 3 state reactivities ---------------------------------

# Author: Mirko Ledda | Date: 06/16/2015 | v1
#
# Based on usage in conjuction with the custom `kernel.3s` script

posterior.3s <- function(D,dims=1) {
  # D:    list    - kernel densities
  # dims: integer - Number of dimensions of the kernel
  
  # Note that u: unpaired, e:helix-end, s: stacked
  nu <- D$uc
  ne <- D$ec
  ns <- D$sc
  
  pu <- nu/(nu+ne+ns)
  pe <- ne/(nu+ne+ns)
  ps <- ns/(nu+ne+ns)
  
  den <- D$u$estimate*pu+D$e$estimate*pe+D$s$estimate*ps
  
  post <- list()
  post$unpaired <- (D$u$estimate*pu)/den
  post$paired_end <- (D$e$estimate*pe)/den
  post$paired_stack <- (D$s$estimate*ps)/den
    
  post$unpaired[is.infinite(post$unpaired) | is.na(post$unpaired)] <- 0
  post$paired_end[is.infinite(post$paired_end) | is.na(post$paired_end)] <- 0
  post$paired_stack[is.infinite(post$paired_stack) | is.na(post$paired_stack)] <- 0
  
  # extract the axis values
  if (dims==1) {
    post$x <- D$u$eval.points
  } else if (dims==2) {
    post$x <- D$u$eval.points[[1]]
    post$y <- D$u$eval.points[[2]]
  } else if (dims==3) {
    post$x <- D$u$eval.points[[1]]
    post$y <- D$u$eval.points[[2]]
    post$z <- D$u$eval.points[[3]]
  }
  
  return(post)
}
