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
ress_2stateD   <- correctRes(ress_2stateD,"D")

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

ischange  <- t(sapply(pathss_2stateDS, function(p) {temp <- apply(p,1,sum) > 0; temp <- c(temp,temp,temp,T,T,T); return(temp)}))
parss_2stateDS  <- ress_2stateDS * ischange

ischange  <- t(sapply(pathss_2stateD, function(p) {temp <- apply(p,1,sum) > 0; temp <- c(temp,temp,T,T,T,T); return(temp)}))
parss_2stateD  <- ress_2stateD * ischange

ischange  <- t(sapply(pathss_3state, function(p) {temp <- apply(p,1,sum) > 0; temp <- c(temp,temp,temp,T,temp); return(temp)}))
parss_3state  <- ress_3state[,1:13] * ischange

# d estimations vs year for each pair 
dsDM  <- lapply(1:nrow(ress_2stateDM), function(x) {matrix(ress_2stateDM[x,1:2],1,) %*% dlvPath_DM(ress_2stateDM[x,1:7],processdat(dat[,x]))})
dsDS  <- lapply(1:nrow(ress_2stateDS), function(x) {matrix(ress_2stateDS[x,1:2],1,) %*% dlvPath_DS(ress_2stateDS[x,1:7],processdat(dat[,x]))})
dsD   <- lapply(1:nrow(ress_2stateD),  function(x) {matrix(ress_2stateD[x,1:2],1,)  %*% dlvPath_D(ress_2stateD[x,1:6],  processdat(dat[,x]))})
ds3 <- lapply(1:nrow(ress_3state), function(x) {matrix(ress_3state[x,1:3],1,) %*% dlvPath_dm3(ress_3state[x,1:13],processdat(dat[,x]))})
# Standardizing the length of ds
dsDM  <- sapply(1:ncol(dat), function(i) deprocessds(dsDM[[i]],dat[,i]))
dsDS  <- sapply(1:ncol(dat), function(i) deprocessds(dsDS[[i]],dat[,i]))
dsD   <- sapply(1:ncol(dat), function(i) deprocessds(dsD[[i]],dat[,i]))
ds3   <- sapply(1:ncol(dat), function(i) deprocessds(ds3[[i]],dat[,i]))


tlabels <- unlist(lapply(1950:2018, function(y) paste0(y,c("Q1","Q2","Q3","Q4"))))
tlabels <- tlabels[1:nrow(dat)]

# Code for plotting d estimations vs year of given pair for given dataset
plot_specificDM <- function(serind){
  sername <- colnames(dat)[serind]
  ser     <- dsDM[,serind]; maintext <- paste0(sername, " Quarterly, CAD, 2 State, (d,mu) switching")
  temp    <- data.frame(Quarter=tlabels,d = ser)
  temp    <- temp[!is.na(temp$d),]
  
  p <- ggplot(temp, aes(x = Quarter, y = d,group=1))
  p <- p + geom_line() + scale_x_discrete("Quarter", breaks=temp$Quarter[seq(1,length(temp$Quarter),20)])
  p <- p + ggtitle(maintext) + if(nrow(temp)>200){theme(axis.text.x = element_text(angle=45))} else {NULL}
  ggsave(paste0("pairplots/","CAD_Q","_",sername,"_twostateDM.jpg"),p)
}

plot_specificDS <- function(serind){
  sername <- colnames(dat)[serind]
  ser     <- dsDS[,serind]; maintext <- paste0(sername, " Quarterly, CAD, 2 State, (d,sigma) switching")
  temp    <- data.frame(Quarter=tlabels,d = ser)
  temp    <- temp[!is.na(temp$d),]
  
  p <- ggplot(temp, aes(x = Quarter, y = d,group=1))
  p <- p + geom_line() + scale_x_discrete("Quarter", breaks=temp$Quarter[seq(1,length(temp$Quarter),20)])
  p <- p + ggtitle(maintext) + if(nrow(temp)>200){theme(axis.text.x = element_text(angle=45))} else {NULL}
  ggsave(paste0("pairplots/","CAD_Q","_",sername,"_twostateDS.jpg"),p)
}

plot_specificD <- function(serind){
  sername <- colnames(dat)[serind]
  ser     <- dsD[,serind]; maintext <- paste0(sername, " Quarterly, CAD, 2 State, d switching")
  temp    <- data.frame(Quarter=tlabels,d = ser)
  temp    <- temp[!is.na(temp$d),]
  
  p <- ggplot(temp, aes(x = Quarter, y = d,group=1))
  p <- p + geom_line() + scale_x_discrete("Quarter", breaks=temp$Quarter[seq(1,length(temp$Quarter),20)])
  p <- p + ggtitle(maintext) + if(nrow(temp)>200){theme(axis.text.x = element_text(angle=45))} else {NULL}
  ggsave(paste0("pairplots/","CAD_Q","_",sername,"_twostateD.jpg"),p)
}

lapply(1:74, plot_specificDM)
lapply(1:74, plot_specificDS)
lapply(1:74, plot_specificD)

plot_specific3 <- function(serind){
  sername <- colnames(dat)[serind]
  ser      <- ds3[,serind]; maintext <- paste0(sername, " Quarterly, CAD, 3 State")
  temp <- data.frame(Quarter=tlabels,d = ser)
  temp <- temp[!is.na(temp$d),]
  
  p <- ggplot(temp, aes(x = Quarter, y = d,group=1))
  p <- p + geom_line() + scale_x_discrete("Quarter", breaks=temp$Quarter[seq(1,length(temp$Quarter),20)])
  p <- p + ggtitle(maintext) + if(nrow(temp)>200){theme(axis.text.x = element_text(angle=45))} else {NULL}
  ggsave(paste0("pairplots/","CAD_Q","_",sername,"_threestate.jpg"),p)
}

lapply(1:74, plot_specific3)
