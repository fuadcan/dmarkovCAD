if(Sys.info()["sysname"]=="Linux"){setwd("~/Documents/dmarkovCAD/")} else {setwd("~/dMarkovCAD/")}

source("lnviDSM2.R")
source("lnviDM2.R")
source("lnviDS2.R")
source("lnviD2.R")
source("twostate_DM.R")
source("twostate_DS.R")
source("twostate_D.R")
source("dlvPath.R")
source("utils.R")

dir.create("output",showWarnings = F)
dir.create("results",showWarnings = F)
dir.create("rejplots",showWarnings = F)
dir.create("pairplots",showWarnings = F)

# Estimating parameters for two state analysis
twostateDSMres <- twostate_DSM("CAD_quarterly_BOP")
twostateDMres  <- twostate_DM("CAD_quarterly_BOP")
twostateDSres  <- twostate_DS("CAD_quarterly_BOP")
twostateDres   <- twostate_D("CAD_quarterly_BOP")

# Estimating parameters for three state analysis
threestateress <- threestate("CAD_quarterly_BOP.rda")

# Reading outputs
ress_2stateDSM <- get(load(paste0("output/","CAD_quarterly_BOP_2stateDSM_ress.rda")))
ress_2stateDM  <- get(load(paste0("output/","CAD_quarterly_BOP_2stateDM_ress.rda")))
ress_2stateDS  <- get(load(paste0("output/","CAD_quarterly_BOP_2stateDS_ress.rda")))
ress_2stateD   <- get(load(paste0("output/","CAD_quarterly_BOP_2stateD_ress.rda")))
# ress_3state    <- get(load(paste0("output/","CAD_quarterly_BOP_ress.rda")))

# Converting to tables
ress_2stateDSM <- totable(ress_2stateDSM)
ress_2stateDM  <- totable(ress_2stateDS)
ress_2stateDS  <- totable(ress_2stateDS)
ress_2stateD   <- totable(ress_2stateD)
# ress_3state    <- totable(ress_3state)

# Reformating results
ress_2stateDSM <- correctRes(ress_2stateDSM,"DSM")
ress_2stateDM  <- correctRes(ress_2stateDM,"DM")
ress_2stateDS  <- correctRes(ress_2stateDS,"DS")
ress_2stateD   <- correctRes(ress_2stateD, "D")

# data
dat <- get(load(paste0("CAD_quarterly_BOP",".rda")))
dat <- dat/10e10

# State switching series for dm and d
pathss_2stateDSM <- lapply(1:nrow(ress_2stateDSM), function(x) dlvPath_DSM(ress_2stateDSM[x,],processdat(dat[,x])))
pathss_2stateDM  <- lapply(1:nrow(ress_2stateDM),  function(x) dlvPath_DM(ress_2stateDM[x,],processdat(dat[,x])))
pathss_2stateDS  <- lapply(1:nrow(ress_2stateDS),  function(x) dlvPath_DS(ress_2stateDS[x,],processdat(dat[,x])))
pathss_2stateD   <- lapply(1:nrow(ress_2stateD),   function(x) dlvPath_D(ress_2stateD[x,],processdat(dat[,x])))
# pathss_3state    <- lapply(1:nrow(ress_3state),    function(x) dlvPath_dm3(ress_3state[x,1:13],processdat(dat[,x])))

# Correcting results (Reformatting results of pairs do not change state)
ischange  <- t(sapply(pathss_2stateDSM, function(p) {temp <- apply(p,1,sum) > 0; temp <- c(temp,T,T,temp,temp,T,T); return(temp)}))
parss_2stateDSM  <- ress_2stateDSM * ischange
colnames(parss_2stateDSM) <- c("d_1","d_2","P_11","P_22","sigma_1","sigma_2","mu_1","mu_2","Neg_lklh","Opt. Converged")
rownames(parss_2stateDSM) <- colnames(dat)
write.csv(parss_2stateDSM, "results/CAD_quarterly_BOP_2stateDSM_results.csv")
rep_2stateDSM    <- report(parss_2stateDSM,"DSM")

ischange  <- t(sapply(pathss_2stateDM, function(p) {temp <- apply(p,1,sum) > 0; temp <- c(temp,T,T,T,temp,T,T); return(temp)}))
parss_2stateDM  <- ress_2stateDM * ischange
colnames(parss_2stateDM) <- c("d_1","d_2","P_11","P_22","sigma","mu_1","mu_2","Neg_lklh","Opt. Converged")
rownames(parss_2stateDM) <- colnames(dat)
write.csv(parss_2stateDM, "results/CAD_quarterly_BOP_2stateDM_results.csv")
rep_2stateDM    <- report(parss_2stateDM,"DM")

ischange  <- t(sapply(pathss_2stateDS, function(p) {temp <- apply(p,1,sum) > 0; temp <- c(temp,T,T,temp,T,T,T); return(temp)}))
parss_2stateDS  <- ress_2stateDS * ischange
colnames(parss_2stateDS) <- c("d_1","d_2","P_11","P_22","sigma_1","sigma_2","mu","Neg_lklh","Opt. Converged")
rownames(parss_2stateDS) <- colnames(dat)
write.csv(parss_2stateDS, "results/CAD_quarterly_BOP_2stateDS_results.csv")
rep_2stateDS    <- report(parss_2stateDS,"DS")

ischange  <- t(sapply(pathss_2stateD, function(p) {temp <- apply(p,1,sum) > 0; temp <- c(temp,T,T,T,T,T,T); return(temp)}))
parss_2stateD  <- ress_2stateD * ischange
colnames(parss_2stateD) <- c("d_1","d_2","P_11","P_22","sigma","mu","Neg_lklh","Opt. Converged")
rownames(parss_2stateD) <- colnames(dat)
write.csv(parss_2stateD, "results/CAD_quarterly_BOP_2stateD_results.csv")
rep_2stateD    <- report(parss_2stateD,"D")

ischange  <- t(sapply(pathss_3state, function(p) {temp <- apply(p,1,sum) > 0; temp <- c(temp,T,T,temp,T,temp,T,T); return(temp)}))
parss_3state  <- ress_3state * ischange
colnames(parss_3state) <- c("d_1","d_2","d_3","P_11","P_22","P_33","P_12","P_21","P_31","sigma","mu_1","mu_2","mu_3","Neg_lklh","Opt. Converged")
rownames(parss_3state) <- colnames(dat)
write.csv(parss_3state, "results/CAD_quarterly_BOP_3state_results.csv")

# Postprocessing and binding reports

rep_2stateDM <- cbind(rep_2stateDM[,1:6],rep(NaN,5),rep_2stateDM[,7:8])
rep_2stateDS <- cbind(rep_2stateDS[,1:7],rep(NaN,5),rep_2stateDS[,8])
rep_2stateD  <- cbind(rep_2stateD[,1:6],rep(NaN,5),rep_2stateD[,7],rep(NaN,5))

rep_2state   <- do.call(rbind,list(rep_2stateDSM,rep_2stateDM,rep_2stateDS,rep_2stateD))
type_list    <- c("(d,mu,sigma)","(d,mu)","(d,sigma)","(d)")
rep_2state   <- data.frame(cbind(c(sapply(type_list, function(typ) rep(typ,5,)))),rownames(rep_2state),rep_2state)
colnames(rep_2state)[1:2] <- c("Model","Conv.")

write.csv(rep_2state,"results/report_2state.csv")

# Plotting all path graphs
source("plots_d.R")
source("plots_dm.R")
