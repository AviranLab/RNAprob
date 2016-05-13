
make_RNAprob_param_file_hist <- function(x,pairing,bin_size,fn) {
  
  root_out_base <- "~/tmp/"
  
  match_range <- function(x,bin_size) {
    dx <- x/bin_size
    rm <- x-(floor(abs(dx))*sign(dx))*bin_size
    if (rm<0) {
      out <- x+abs(rm)
    } else {
      out <- x+(bin_size-rm)
    }
    return(out)
  }
  
  data_range <- range(x)
  data_range <- sapply(data_range,match_range,bin_size)
  data_range[1] <- data_range[1]-bin_size
  
  dta <- list()
  dta$un <- x[pairing==0 & !is.na(pairing)]
  dta$ph <- x[pairing==1 & !is.na(pairing)]
  dta$ps <- x[pairing==2 & !is.na(pairing)]
  
  # Compute the density of the data based on histograms
  xhist <- seq(data_range[1],data_range[2],bin_size)
  dens <- list()
  dens$un <- hist(dta$un,breaks=xhist,plot=F)$density
  dens$ph <- hist(dta$ph,breaks=xhist,plot=F)$density
  dens$ps <- hist(dta$ps,breaks=xhist,plot=F)$density
  
  # WRITE PARAM FILE --------------------------------------------------------
  fn.ph <- paste0(root_out_base,"ph.txt")
  fn.ps <- paste0(root_out_base,"ps.txt")
  fn.un <- paste0(root_out_base,"un.txt")
  
  # Write files
  ds <- as.data.frame(c(paste(min(xhist),bin_size),dens$ph))
  colnames(ds) <- ">SHAPE|helix_end|X"
  write.table(ds,file=fn.ph,quote=F,sep="",row.names = F)
  
  ds <- as.data.frame(c(paste(min(xhist),bin_size),dens$ps))
  colnames(ds) <- ">SHAPE|stacked|X"
  write.table(ds,file=fn.ps,quote=F,sep="",row.names = F)
  
  ds <- as.data.frame(c(paste(min(xhist),bin_size),dens$un))
  colnames(ds) <- ">SHAPE|unpaired|X"
  write.table(ds,file=fn.un,quote=F,sep="",row.names = F)
  
  system(paste0(paste("cat",fn.ph,fn.ps,fn.un,">"),fn,".txt"))
  system(paste("rm",fn.ph,fn.ps,fn.un))
  
}

