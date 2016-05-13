# Disentangle stacked and helix-end paired bases
get_pair_3s <- function(rp) {
  n <- length(rp)
  new_pair <- rep(NA,n)
  for (i in 1:n) {
    if (rp[i] == 0) {
      new_pair[i] <- 0
    } else if (i==1 | i==n) {
      new_pair[i] <- 1
    } else if (rp[i-1]!=0 & rp[i+1]!=0 & rp[i-1]==rp[i]+1 & rp[i+1]==rp[i]-1) {
      new_pair[i] <- 2
    } else {
      new_pair[i] <- 1
    }
  }
  return(new_pair)
}