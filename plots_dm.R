# Reading outputs
ress_2state  <- get(load(paste0("output/","CAD_quarterly_BOP_2state_ress.rda")))
ress_3state  <- get(load(paste0("output/","CAD_quarterly_BOP_ress.rda")))

# Converting to tables
ress_2state <- totable(ress_2state)
ress_3state <- totable(ress_3state)

# Reformating results
ress_2state  <- correctRes(ress_2state[,1:7])

dat <- get(load(paste0("CAD_quarterly_BOP",".rda")))
dat <- dat/10e10

# State switching series for dm and d
pathss_2state <- lapply(1:nrow(ress_2state), function(x) dlvPath_dm(ress_2state[x,],processdat(dat[,x])))
pathss_3state <- lapply(1:nrow(ress_3state), function(x) dlvPath_dm3(ress_3state[x,1:13],processdat(dat[,x])))

# Correcting results (Reformatting results of pairs do not change state)
ischange  <- t(sapply(pathss_2state, function(p) {temp <- apply(p,1,sum) > 0; temp <- c(temp,temp,T,temp); return(temp)}))
parss_2state  <- ress_2state[,-8] * ischange

ischange  <- t(sapply(pathss_3state, function(p) {temp <- apply(p,1,sum) > 0; temp <- c(temp,temp,temp,T,temp); return(temp)}))
parss_3state  <- ress_3state[,1:13] * ischange

# d estimations vs year for each pair 
ds  <- lapply(1:nrow(ress_2state), function(x) {matrix(ress_2state[x,1:2],1,) %*% dlvPath_dm(ress_2state[x,1:7],processdat(dat[,x]))})
ds3 <- lapply(1:nrow(ress_3state), function(x) {matrix(ress_3state[x,1:3],1,) %*% dlvPath_dm3(ress_3state[x,1:13],processdat(dat[,x]))})
# Standardizing the length of ds
ds  <- sapply(1:ncol(dat), function(i) deprocessds(ds[[i]],dat[,i]))
ds3 <- sapply(1:ncol(dat), function(i) deprocessds(ds3[[i]],dat[,i]))


tlabels <- unlist(lapply(1950:2018, function(y) paste0(y,c("Q1","Q2","Q3","Q4"))))
tlabels <- tlabels[1:nrow(dat)]

# Code for plotting d estimations vs year of given pair for given dataset
plot_specific <- function(serind){
  sername <- colnames(dat)[serind]
  ser      <- ds[,serind]; maintext <- paste0(sername, " Quarterly, CAD")
  temp <- data.frame(Quarter=tlabels,d = ser)
  temp <- temp[!is.na(temp$d),]
  
  p <- ggplot(temp, aes(x = Quarter, y = d,group=1))
  p <- p + geom_line() + scale_x_discrete("Quarter", breaks=temp$Quarter[seq(1,length(temp$Quarter),20)])
  p <- p + ggtitle(maintext)
  ggsave(paste0("pairplots/","CAD_Q","_",sername,"_twostate.jpg"),p)
}

lapply(1:74, plot_specific)
