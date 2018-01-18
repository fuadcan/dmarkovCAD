setwd("~/../Dropbox/Ege - Fuat/dMarkovCADConvergence/")
save_file <- "~/../Dropbox/Ege - Fuat/dMarkovCADConvergence/results/"
data_file <- "~/../Dropbox/Ege - Fuat/dMarkovCADConvergence/data/"
source("lnviDM3.R")
source("return_path.R")
library("zoo")
# yearOrRegion <- "G7"
threestateCAD <- function(yearOrRegion){
  
  # dat <- read.csv(paste0("data/mds_",yearOrRegion,"-1950.csv"),sep=";")
  # dat <- readxl::read_excel("data/cad_gdp_oecd.xlsx")
  # dat <- as.data.frame(dat)
  
  dat <- get(load("CAD_quarterly_BOP.rda"))
  dat <- dat/10e10
  
  lowerV <- c( -2,-2,  -2, .6 , .6 , .6 , 0, 0, 0,0.001,-5,-5,-5)
  inits  <- c(1.2, 1, 0.7, .61, .61, .61,.1,.1,.1,   .4,.1,.1,.1) 
  upperV <- c(  2, 2,   2,.999,.999,.999,.3,.3,.3,    5, 5, 5, 5)
  
  const_mat <- matrix(0,13,13)
  diag(const_mat) <- 1
  const_mat <- rbind(const_mat,-const_mat)
  dmmvec <- -c(0,0,0,1,0,0,1,0,0,0,0,0,0)
  dmmvec <- c(dmmvec,0,dmmvec,0,dmmvec)[1:39]
  const_mat <- rbind(const_mat,matrix(dmmvec,3,,T))
  const_mat <- cbind(const_mat,c(lowerV,-upperV,rep(-1,3)))
  
  # ress    <- lapply(1:ncol(pdat), function(i) constrOptim(inits, function(p) -lnviD3(p,pdat[,i]), NULL, ui = const_mat[,-14], const_mat[,14]))
  
  processdat <- function(ser){
    datrange <- range(which(!is.na(ser)))
    ser <- ser[datrange[1]:datrange[2]]
    x   <- zoo(ser)
    ser <- na.approx(x,1:length(ser))
    ser
  }
  
  
  
  # lnviD3(inits,processdat(dat[,1]))
  # constrOptim(inits, function(p) -lnviD3(p,processdat(dat[,1])/10e10), NULL, ui = const_mat[,-14], const_mat[,14])
  # lnvid <- function(b,w){ res <- tryCatch(lnviD3(b,w),error= function(e) list(NA,NA)); return(res) }
  ress    <- lapply(1:ncol(dat), function(i) {cat(paste0(i,"\n")); constrOptim(inits, function(p) -lnviD3(p,processdat(dat[,i])), NULL, ui = const_mat[,-14], const_mat[,14])})
  save(ress, file = paste0(save_file,yearOrRegion,"_ress.rda"))
  
  return(ress) 
}

