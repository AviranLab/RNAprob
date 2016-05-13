## Classify reactivities according to a reference structure and under different models

# Author: Mirko Ledda | Date: 04/10/2015 | v1
#
# Handles 2 and 3 states scenarios
#
# Output a list of reactivities per states and a vector of pairing IDs

pairing_reactivities <- function(x,ref_fn,zerol=F,minl=F,model='raw',si=-0.6,sm=1.8,lambda=0.1,states=2) {
  # x:      1 by n  - vector of data
  # ref_fn: char    - reference .ct file
  # zerol:  float   - reactivity threshold for setting zeros
  # minl:   float   - reactivity minima
  # model:  char    - Transform model: raw = no transform | ln | log10 | boxcox | weeks = Weeks ln(a*x) + b formula
  # si:     float   - intercept in Weeks formula
  # sm:     float   - slope in Weeks formula
  # lambda: float   - Parameter for Box-Cox transform
  # states: integer - Number of status clustering. 2 = paired/unpaired | 3 = pair/helix-end/stacked
  
  # read ct
  ct <- read.table(ref_fn,skip=1)
  ref_pair <- ct$V5
  ref_pair[ref_pair>0] <- 1
  
  switch(model,
         raw={xt <- x},
         ln={x_zero <- x==0
         x[x_zero] <- NA
         xt <- log(x)},
         log10={x_zero <- x==0
         x[x_zero] <- NA
         xt <- log10(x)},
         weeks={xt <- log(x+1) * sm + si},
         boxcox={x_zero <- x==0
         x[x_zero] <- NA
         xt <- boxcox_transform(x,lambda)}
  )
  
  # Set zeros and set minimum if needed
  if (is.logical(zerol)) zerol <- min(xt[!is.infinite(xt)],na.rm=T)
  if (is.logical(minl)) minl <- min(xt[!is.infinite(xt)],na.rm=T)
  
  # Assign minima
  xt[xt<minl] <- minl
  
  # Treat zeros if required
  if (exists('x_zero')) xt[x_zero] <- zerol
  
  # Disentangle stacked and helix-end paired bases
  if (states==3) {
    for (i in 2:(length(ref_pair)-1)) {
      if (sum(!is.na(ref_pair[(i-1):(i+1)]))==3) {
        if (sum(ref_pair[(i-1):(i+1)]>0)==3) {
          ref_pair[i] <- 2
        }
      }
    }  
  }
  
  out <- list()
  out$unpaired <- xt[ref_pair==0]
  out$pairing <- ref_pair
  
  # nucleotides clustering
  if (states==2) {
    out$paired <- xt[ref_pair==1]
  } else {
    out$paired_stack <- xt[ref_pair==2]
    out$paired_end <- xt[ref_pair==1]
  }
  
  return(out)
}
