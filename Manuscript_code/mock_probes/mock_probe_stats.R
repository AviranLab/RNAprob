### Stats for mock probe simulation results

## Remove all active variables from the workspace
rm(list = ls())

# Inputs ------------------------------------------------------------------
data_dir <- "data/folding_results/mock_probes/" # master directory for individual RNA scores
root_out <- "output/" # outpout directory
n <- 10 # number of simulations

# Scenario 1 --------------------------------------------------------------
sc_name <- "scenario1"

# SWL-average MCC scores
dta <- read.table(paste0(data_dir,"simu-results-",sc_name,".csv"),header=T,sep=",")
dta <- dta[1:10,-1]
dta <- dta[1:10,2:3]
p1 <- t.test(dta[,1],dta[,2],paired = T)$p.value # paired t-test
mean1 <- colMeans(dta)
sd1 <- apply(dta,2,sd)

# Individual RNA scores
lin <- NULL
prob <- NULL

for (i in 1:n) {
  lin <- c(lin,read.table(paste0(data_dir,"RNAlin/",sc_name,"/simu",i-1,"/summary-mcc.score"))$V1[1:7])
  prob <- c(prob,read.table(paste0(data_dir,"RNAprob-3/",sc_name,"/simu",i-1,"/summary-mcc.score"))$V1[1:7])
}

cnt_best_lin1 <- sum(lin>prob)
cnt_best_prob1 <- sum(prob>lin)
cnt_same1 <- sum(prob==lin)

# Scenario 2 --------------------------------------------------------------
sc_name <- "scenario2"

# SWL-average MCC scores
dta <- read.table(paste0(data_dir,"simu-results-",sc_name,".csv"),header=T,sep=",")
dta <- dta[1:10,-1]
dta <- dta[1:10,2:3]
p2 <- t.test(dta[,1],dta[,2],paired = T)$p.value # paired t-test
mean2 <- colMeans(dta)
sd2 <- apply(dta,2,sd)

# Individual RNA scores
lin <- NULL
prob <- NULL

for (i in 1:n) {
  lin <- c(lin,read.table(paste0(data_dir,"RNAlin/",sc_name,"p/simu",i-1,"/summary-mcc.score"))$V1[1:7])
  prob <- c(prob,read.table(paste0(data_dir,"RNAprob-3/",sc_name,"p/simu",i-1,"/summary-mcc.score"))$V1[1:7])
}

cnt_best_lin2 <- sum(lin>prob)
cnt_best_prob2 <- sum(prob>lin)
cnt_same2 <- sum(prob==lin)

# Write output file
df <- cbind(mean1,sd1,p1,mean2,sd2,p2)
write.csv(df,paste0(root_out,"mock_probe_stats.csv"))

