setwd("../Desktop/CAD")
dat <- readxl::read_excel("CAD_quarterly_BOP_mined.xlsx")

dat <- data.frame(dat)
# dat <- reshape2::melt(dat, id.vars = c("Country.Name","Indicator.Name"))

holeexist <- function(ser){
  ser <- ser[-(1:2)]
  serbeg <- which(!is.na(ser))[1]
  serend <- tail(which(!is.na(ser)),1)
  ser    <- ser[serbeg:serend]
  sum(is.na(ser))
}

dat <- dat[apply(dat,1,holeexist)==0,]
dat[apply(!is.na(dat),1,sum) >= 60,1]
# ind1 <- apply(!is.na(dat),1,sum) >= 60
# ind2 <- apply(!is.na(dat),1,sum) >= 80
# ind3 <- apply(!is.na(dat),1,sum) >= 60 & apply(!is.na(dat),1,sum) < 80

dat <- dat[apply(!is.na(dat),1,sum) >= 80,]
CAD_dat <- dat[,-2]
rownames(CAD_dat) <- CAD_dat[,1]
CAD_dat <- CAD_dat[,-1]
CAD_dat <- t(CAD_dat)
