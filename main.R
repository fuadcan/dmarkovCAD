if(Sys.info()["sysname"]=="Linux"){setwd("~/Documents/CAD/")} else {setwd("~/CAD/")}

source("lnviDM2.R")
source("lnviDS2.R")
source("lnviD2.R")
source("twostate_DM.R")
source("twostate_DS.R")
source("twostate_D.R")
source("dlvPath.R")

dir.create("output",showWarnings = F)
dir.create("results",showWarnings = F)
dir.create("rejplots",showWarnings = F)
dir.create("pairplots",showWarnings = F)

# Estimating parameters for two state analysis
twostateDMres <- twostate_DM("CAD_quarterly_BOP.rda")
twostateDSres <- twostate_DS("CAD_quarterly_BOP.rda")
twostateDres  <- twostate_D("CAD_quarterly_BOP.rda")

# Estimating parameters for three state analysis
threestateress <- threestate("CAD_quarterly_BOP.rda")

# Reading outputs
ress_2stateDM  <- get(load(paste0("output/","CAD_quarterly_BOP_2stateDM_ress.rda")))
ress_2stateDS  <- get(load(paste0("output/","CAD_quarterly_BOP_2stateDS_ress.rda")))
ress_2stateD   <- get(load(paste0("output/","CAD_quarterly_BOP_2stateD_ress.rda")))
ress_3state    <- get(load(paste0("output/","CAD_quarterly_BOP_ress.rda")))

# Converting to tables
ress_2stateDM <- totable(ress_2stateDM)
ress_2stateDS <- totable(ress_2stateDS)
ress_2stateD  <- totable(ress_2stateD)
ress_3state   <- totable(ress_3state)

# Reformating results
ress_2stateDM  <- correctRes(ress_2stateDM,"DM")
ress_2stateDS  <- correctRes(ress_2stateDS,"DS")
ress_2stateD   <- correctRes(ress_2stateD, "D")

# data
dat <- get(load(paste0("CAD_quarterly_BOP",".rda")))
dat <- dat/10e10

# State switching series for dm and d
pathss_2stateDM <- lapply(1:nrow(ress_2stateDM), function(x) dlvPath_DM(ress_2stateDM[x,],processdat(dat[,x])))
pathss_2stateDS <- lapply(1:nrow(ress_2stateDS), function(x) dlvPath_DS(ress_2stateDS[x,],processdat(dat[,x])))
pathss_2stateD  <- lapply(1:nrow(ress_2stateD),  function(x) dlvPath_D(ress_2stateD[x,],processdat(dat[,x])))
pathss_3state   <- lapply(1:nrow(ress_3state), function(x) dlvPath_dm3(ress_3state[x,1:13],processdat(dat[,x])))

# Correcting results (Reformatting results of pairs do not change state)
ischange  <- t(sapply(pathss_2stateDM, function(p) {temp <- apply(p,1,sum) > 0; temp <- c(temp,temp,T,temp,T,T); return(temp)}))
parss_2stateDM  <- ress_2stateDM * ischange
colnames(parss_2stateDM) <- c("d_1","d_2","P_11","P_22","sigma","mu_1","mu_2","Neg_lklh","Opt. Converged")
rownames(parss_2stateDM) <- colnames(dat)
write.csv(parss_2stateDM, "CAD_quarterly_BOP_2stateDM_results.csv")

ischange  <- t(sapply(pathss_2stateDS, function(p) {temp <- apply(p,1,sum) > 0; temp <- c(temp,temp,temp,T,T,T); return(temp)}))
parss_2stateDS  <- ress_2stateDS * ischange
colnames(parss_2stateDS) <- c("d_1","d_2","P_11","P_22","sigma_1","sigma_2","mu","Neg_lklh","Opt. Converged")
rownames(parss_2stateDS) <- colnames(dat)
write.csv(parss_2stateDS, "CAD_quarterly_BOP_2stateDS_results.csv")

ischange  <- t(sapply(pathss_2stateD, function(p) {temp <- apply(p,1,sum) > 0; temp <- c(temp,temp,T,T,T,T); return(temp)}))
parss_2stateD  <- ress_2stateD * ischange
colnames(parss_2stateD) <- c("d_1","d_2","P_11","P_22","sigma","mu","Neg_lklh","Opt. Converged")
rownames(parss_2stateD) <- colnames(dat)
write.csv(parss_2stateD, "CAD_quarterly_BOP_2stateD_results.csv")

ischange  <- t(sapply(pathss_3state, function(p) {temp <- apply(p,1,sum) > 0; temp <- c(temp,temp,temp,T,temp,T,T); return(temp)}))
parss_3state  <- ress_3state * ischange
colnames(parss_3state) <- c("d_1","d_2","d_3","P_11","P_22","P_33","P_12","P_21","P_31","sigma","mu_1","mu_2","mu_3","Neg_lklh","Opt. Converged")
rownames(parss_3state) <- colnames(dat)
write.csv(parss_3state, "CAD_quarterly_BOP_3state_results.csv")

# Plotting all path graphs
source("plots_d.R")
source("plots_dm.R")
