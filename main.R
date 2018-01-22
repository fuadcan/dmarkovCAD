if(Sys.info()["sysname"]=="Linux"){setwd("~/Documents/CAD/")} else {setwd("~/CAD/")}

source("convDLV_dm.R")
source("dlvPath.R")
source("lnviDM2.R")
if(Sys.info()["sysname"]=="Linux")
  {library("parallel");mclapply.hack <- mclapply} else {source("mclapplyhack.R")}
dir.create("output",showWarnings = F)
dir.create("results",showWarnings = F)
dir.create("rejplots",showWarnings = F)
dir.create("pairplots",showWarnings = F)

# Estimating parameters for two state analysis
twostateress <- twostate("CAD_quarterly_BOP.rda")

# Estimating parameters for three state analysis
threestateress <- threestate("CAD_quarterly_BOP.rda")

# Reading outputs
ress_2state  <- get(load(paste0("output/","CAD_quarterly_BOP_2state_ress.rda")))

# Converting to tables
ress_2state <- totable(ress_2state)
# ress_dm <- lapply(ress_dm, function(r) r[,-8])

# Reformating results
ress_2state  <- correctRes(ress_2state[,1:7])
# ress_dm <- lapply(ress_dm, correctRes)

# data names and pair panels
# dnames     <- gsub("d_|_D_resALL.rda","",dir("output","_D_.*resALL.rda"))

dat <- get(load(paste0("CAD_quarterly_BOP",".rda")))
dat <- dat/10e10

# State switching series for dm and d
pathss_2state <- lapply(1:nrow(ress_2state), function(x) dlvPath_dm(ress_2state[x,],processdat(dat[,x])))
# pathss_D    <- lapply(1:length(ress_d), function(i) lapply(1:nrow(ress_d[[i]]), function(x) dlvPath_d(ress_d[[i]][x,-7],pdats[[i]][,x])))

# Correcting results (Reformatting results of pairs do not change state)
ischange  <- t(sapply(pathss_2state, function(p) {temp <- apply(p,1,sum) > 0; temp <- c(temp,temp,T,temp); return(temp)}))
parss_2state  <- ress_2state[,-8] * ischange
colnames(parss_2state) <- c("d_1","d_2","P_11","P_22","sigma","mu_1","mu_2")
write.csv(parss_2state, "CAD_quarterly_BOP_2state_results.csv")

# ischange  <- lapply(pathss_D, function(paths) t(sapply(paths, function(p) {temp <- apply(p,1,sum) > 0; temp <- c(temp,temp,T,T); return(temp)})))
# parss_d   <- lapply(1:length(ress_d),  function(i)  ress_d[[i]][,-7] * ischange[[i]])

# Reporting outputs
rep_dm <- do.call(rbind,lapply(parss_dm,report))
rep_dm <- data.frame(c(sapply(dnames, function(n) rep(n,5))),rownames(rep_dm),rep_dm)

rep_d  <- do.call(rbind,lapply(parss_d,report))
rep_d  <- data.frame(c(sapply(dnames, function(n) rep(n,5))),rownames(rep_d),rep_d)

colnames(rep_dm)  <- c("Data","Conv?",colnames(rep_dm)[-(1:2)])
colnames(rep_d)   <- c("Data","Conv?",colnames(rep_d)[-(1:2)])

# Writing reports
write.csv(rep_dm, "results/report_dm.csv")
write.csv(rep_d , "results/report_d.csv")

# Plotting all path graphs
source("plots_d.R")
source("plots_dm.R")
