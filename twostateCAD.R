setwd("CAD/dMarkovCAD/")
source("lnviDM2.R")
# source("return_path.R")
# library("zoo")
# yearOrRegion <- "G7"

twostate <- function(dname){
  dat <- get(load(paste0(dname,".rda")))
  dat <- dat/10e10
  
  lowerV <- c(-2,-2,.8,.8,0.001,-5,-5)
  upperV <- c(2,2,.999,.999,5,5,5)
  const_mat <- matrix(0,length(lowerV),length(lowerV))
  diag(const_mat) <- 1
  const_mat <- rbind(const_mat,-const_mat)
  const_mat <- cbind(const_mat,c(lowerV,-upperV))

  inits  <- c(1.2,  1.7,  0.9,  0.9,  0.01, 0.1,0.1)

    processdat <- function(ser){
    datrange <- range(which(!is.na(ser)))
    ser <- ser[datrange[1]:datrange[2]]
    # x   <- zoo(ser)
    # ser <- na.approx(x,1:length(ser))
    ser
  }
  
  
  ress    <- lapply(1:ncol(dat), function(i) {cat(paste0(i,"\n")); tryCatch(constrOptim(inits, function(p) -lnviDM2(p,processdat(dat[,i])), NULL, ui = const_mat[,-8], const_mat[,8]),error=function(e) list(NA,NA))})

  save(ress, file = paste0("output/",dname,"_2state_ress.rda"))
  
  estms <- sapply(ress,function(res) c(res$par,res$value,res$convergence))
  colnames(estms) <- colnames(dat)
  estms <- t(estms)
  colnames(estms) <- c("d_1","d_2","d_3","P_11","P_22","P_33","P_12", "P_21", "P_31","sigma","mu_1","mu_2","mu_3","Neg_lklihood","Convergence")
  write.csv(estms, paste0("output/",dname,"2_state_ress.csv"))
  
}




