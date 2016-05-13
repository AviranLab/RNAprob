## Basic R functions

# Author: Mirko Ledda | Date: 11/10/2015 | v1

# is.even and is.odd implementation ---------------------------------------
is.even <- function(x) {x %% 2 == 0}
is.odd <- function(x) {x %% 2 != 0}

# row/col Min/Max implementations -----------------------------------------
rowMin <- function(x) apply(x,1,min,na.rm = TRUE)
rowMax <- function(x) apply(x,1,max,na.rm = TRUE)
colMin <- function(x) apply(x,2,min,na.rm = TRUE)
colMax <- function(x) apply(x,2,max,na.rm = TRUE)

# Row variance and sd calculation -----------------------------------------
RowVar <- function(x) rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
RowSD <- function(x) sqrt(RowVar(x))

# Compute MCC -------------------------------------------------------------
mcc.calc <- function(TP,FP,TN,FN) {
  TP <- as.numeric(TP); FP <- as.numeric(FP); TN <- as.numeric(TN); FN <- as.numeric(FN); # Needed in case of overflow
  (TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))}

# Concatenate factors -----------------------------------------------------
c.Factor <- function (x, y) {
  newlevels = union(levels(x), levels(y))
  m = match(levels(y), newlevels)
  ans = c(unclass(x), m[unclass(y)])
  levels(ans) = newlevels
  class(ans) = "factor"
  return(ans)
}
