# Power transform using the Box-Cox formula -------------------------------

# Author: Mirko Ledda | Date: 06/02/2015 | v1
#
# If lambda is 0 then simply take the natural logarithm, box-cox otherwise

boxcox.transform <- function(x, lambda = 0,forw=T) {
  # x:      1 by n  - vector of data
  # lambda: float   - Box-Cox lambda parameter (default: 0)
  # forw:   boolean - forward or backward transformation?
  if (forw) {
    if (lambda == 0) {
      y <- log(x)
    } else {
      y <- ((x^lambda)-1)/lambda
    }
  } else {
    if (lambda == 0) {
      y <- exp(x)
    } else {
      y <- (lambda*x+1)^(1/lambda)
    }
  }
  
  return(y)
}