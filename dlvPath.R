
dlvPath_dm <- function(b,w){
  
  
  # b <- c(.4,.8,.99,.86,.4,.9,.1,-.3)
  
  size <- length(w)
  
  d1hat = b[1]      ; d2hat = b[2]
  p11hat = b[3]     ; p22hat = b[4];
  sigma1hat =  b[5] ; sigma2hat =  b[5];
  mu1hat = b[6]     ; mu2hat = b[7];
  
  rho1hat=0; theta1hat=0;
  
  tranmax = matrix(c(p11hat,1-p22hat,1-p11hat,p22hat),2) # CONTROL
  
  mu     <- c(mu1hat,mu2hat)
  sigma  <- c(sigma1hat,sigma2hat)
  coefd1 = matrix(cumprod((d1hat+(0:(size-2)))/(1:(size-1))),,1)
  coefd2 = matrix(cumprod((d2hat+(0:(size-2)))/(1:(size-1))),,1)
  coeff = cbind(coefd1,coefd2)
  
  # muhat is to store the selected path of state: from date 2 */
  eta1  = matrix(0,4,size); # row 1:2 for path1, 3:4 for path2 */
  eta   = matrix(0,4,size); # same as above */
  xihat = matrix(0,4,size); # same as above */
  ehat  = matrix(0,2,size); # row 1, for path1, row 2 for path2 */
  
  # to keep the long memory part, not iid */
  zsigmahat = matrix(0,2,size) # row 1, for path1, row 2 for path2 */
  zhat = matrix(0,2,size) # row 1, for path1, row 2 for path2 */
  lnlk <- matrix(0,2,size)
  
  path1 = matrix(0,2,size)
  path2 = matrix(0,2,size)
  path1[,1] <- c(1,0)
  path2[,1] <- c(0,1)
  #   path1[1,1:2] = matrix(c(1,0),1)
  #   path2[1,1:2] = matrix(c(0,1),1)
  zsigmahat[,1] = c(w[1] - mu[1], w[1] - mu[2])
  # zsigmahat is de-drifted series for each given path; or it is z*sigmahat
  
  zhat[,1] <- zsigmahat[,1] / sigma 
  
  ehat[,1] = zhat[,1]
  
  
  lnlk[,1] <- log(exp(-ehat[,1]^2/2 + 1e-50 ) / (sqrt(2*pi))*
                    c( (1-p22hat)/(2-p11hat-p22hat) +1e-50  , (1-p11hat)/(2-p11hat-p22hat)) +1e-50)
  
  j = 2
  
  
  ####################!!!!!!!!!!!!!!!!!!!!!
  while(j<=size){
    # tempsMUZ and temps are needed for likelihood possibility
    tempsMUZ <- c(t(matrix(c(rev(zsigmahat[1,1:(j-1)]),rev(zsigmahat[2,1:(j-1)])),,2))%*%coeff[1:(j-1),1:2])
    
    temps    <- rep(mu,2)[c(1,3,2,4)]+rho1hat*rep(sigma,2)[c(1,3,2,4)]*rep(zhat[1:2,j-1],2)+
      theta1hat*rep(sigma,2)[c(1,3,2,4)]*rep(ehat[1:2,j-1],2)
    
    # likelihood possibility
    eta1[,j]  <- exp(-(w[j] - temps - tempsMUZ)^2 / (2*rep(sigma,2)[c(1,3,2,4)]^2+1e-50))/
      (sqrt(2*pi*rep(sigma,2)[c(1,3,2,4)]^2))
    
    eta[,j] = eta1[,j] * c(tranmax)
    
    # Case is 1 for 1->1; 2 for 2->1; 3 for 1->2; 4 for 2->2
    case <- c(log(eta[1,j]+1e-50) + lnlk[1,j-1] > log(eta[2,j]+1e-50) + lnlk[2,j-1],
              log(eta[3,j]+1e-50) + lnlk[1,j-1] > log(eta[4,j]+1e-50) + lnlk[2,j-1])
    
    case <- c(sign(!case[1])+1,sign(!case[2])+3) # to convert cases to numbers
    
    # holding previous two path
    oldlnlk <- lnlk
    oldPath <- list(path1,path2)
    oldzsigmahat <- zsigmahat
    oldzhat <- zhat
    oldehat <- ehat
    
    # recalculating values
    # path ends up with state 1:
    lnlk[1,j]   <- log(eta[case[1],j]+1e-50)+oldlnlk[case[1],j-1]
    path1[,1:j] <- cbind(oldPath[[case[1]]][,1:(j-1)],c(1,0))
    zsigmahat[1,1:j] <- c(oldzsigmahat[case[1],1:(j-1)],  w[j] - mu[1] - tempsMUZ[case[1]])
    zhat[1,1:j] <- c(oldzhat[case[1],1:(j-1)], zsigmahat[1,j]/(sigma1hat))
    ehat[1,1:j] <- c(oldehat[case[1],1:(j-1)], (w[j] - tempsMUZ[case[1]] - temps[case[1]])/sigma1hat) # ONEMLI!!!
    
    # path ends up with state 2:
    lnlk[2,j]   = log(eta[case[2],j]+1e-50)+oldlnlk[case[2]-2,j-1];
    path2[,1:j] <- cbind(oldPath[[case[2]-2]][,1:(j-1)],c(0,1))
    zsigmahat[2,1:j] = c(oldzsigmahat[case[2]-2,1:(j-1)], w[j] - mu[2] - tempsMUZ[case[2]])
    zhat[2,1:j] = c(oldzhat[case[2]-2,1:(j-1)],zsigmahat[2,j]/(sigma2hat));
    ehat[2,1:j] = c(oldehat[case[2]-2,1:(j-1)],(w[j] - tempsMUZ[case[2]] -temps[case[2]])/sigma2hat)
    
    j=j+1
  }
  #   lnlipath1 = t(lnlk1);
  #   lnlipath2 = t(lnlk2);
  if(lnlk[1,size] > lnlk[2,size]){path <- path1} else {path <- path2}
  
  return(path)
}

dlvPath_dm3 <- function(b,w){
  
  # b      <- c(.4 ,.8,   1.1,    .3,    .3,.3,.3,.3,.3,.4,.1,.2,.3)
  # b <- inits
  # w <- rnorm(1000)
  size <- length(w)
  
  d1hat  = b[1] ; d2hat  = b[2] ; d3hat  = b[3]
  p11hat = b[4] ; p22hat = b[5] ; p33hat = b[6]
  p12hat = b[7] ; p21hat = b[8] ; p31hat = b[9] 
  sigma1hat <- sigma2hat <- sigma3hat <-  b[10];
  mu1hat <- b[11]; mu2hat <- b[12]; mu3hat <- b[13]
  
  rho1hat=0; theta1hat=0;
  tranmax <- matrix(,3,3)
  diag(tranmax) <- b[4:6]
  tranmax[c(4,2,3)] <- c(p12hat, p21hat, p31hat)
  tranmax[c(7,8,6)] <- apply(tranmax, 1, function(x) 1-sum(x,na.rm=T))
  
  mu     <- c(mu1hat,mu2hat,mu3hat)
  sigma  <- c(sigma1hat,sigma2hat,sigma3hat)
  coefd1 = matrix(cumprod((d1hat+(0:(size-2)))/(1:(size-1))),,1)
  coefd2 = matrix(cumprod((d2hat+(0:(size-2)))/(1:(size-1))),,1)
  coefd3 = matrix(cumprod((d3hat+(0:(size-2)))/(1:(size-1))),,1)
  coeff = cbind(coefd1,coefd2,coefd3)
  
  # muhat is to store the selected path of state: from date 2 */
  eta1  = matrix(0,9,size); # row 1:3 for path1, 4:6 for path2 */
  eta   = matrix(0,9,size); # same as above */
  xihat = matrix(0,9,size); # same as above */
  ehat  = matrix(0,3,size); # row 1, for path1, row 2 for path2 */
  
  # to keep the long memory part, not iid */
  zsigmahat = matrix(0,3,size) # row 1, for path1, row 2 for path2 */
  zhat      = matrix(0,3,size) # row 1, for path1, row 2 for path2 */
  lnlk <- matrix(0,3,size)
  
  path1 = matrix(0,3,size)
  path2 = matrix(0,3,size)
  path3 = matrix(0,3,size)
  path1[,1] <- c(1,0,0)
  path2[,1] <- c(0,1,0)
  path3[,1] <- c(0,0,1)
  
  zsigmahat[,1] <- sapply(mu, function(m) w[1] - m) 
  # zsigmahat is de-drifted series for each given path; or it is z*sigmahat
  
  zhat[,1] <- zsigmahat[,1] / sigma 
  
  ehat[,1] <- zhat[,1]
  
  limprobs <- (tranmax %*% tranmax %*% tranmax %*% tranmax %*% tranmax %*% tranmax )[1,]
  lnlk[,1] <- log(exp(-ehat[,1]^2/2 + 1e-50 ) / (sqrt(2*pi))*limprobs)
  
  j = 2
  
  
  ####################!!!!!!!!!!!!!!!!!!!!!
  while(j<=size){
    # tempsMUZ and temps are needed for likelihood possibility
    
    tempsMUZ <- c(t(matrix(c(rev(zsigmahat[1,1:(j-1)]),rev(zsigmahat[2,1:(j-1)]),rev(zsigmahat[3,1:(j-1)])),,3))%*%coeff[1:(j-1),])
    
    temps    <- c(matrix(rep(mu,3),,3,T)) + rho1hat*c(matrix(rep(sigma,3),,3,T))*rep(zhat[,j-1],3)+
      theta1hat*c(matrix(rep(sigma,3),,3,T))*rep(ehat[,j-1],3)
    
    # likelihood possibility
    eta1[,j]  <- exp(-(w[j] - temps - tempsMUZ)^2 / (2*c(matrix(rep(sigma,3),,3,T))^2+1e-50))/
      (sqrt(2*pi*c(matrix(rep(sigma,3),,3,T))^2))
    
    eta[,j] <- eta1[,j] * c(tranmax)
    
    # Case is 1 for 1->1; 2 for 2->1; 3 for 1->2; 4 for 2->2
    case  <- matrix(sapply(1:9, function(x) log(eta[x,j]+1e-50) + lnlk[(x-1) %% 3 + 1,j-1] ),,3)
    case  <- which(apply(case, 2, function(s) max(s) == s))
    case2 <- (case-1) %% 3 + 1
    
    # holding previous paths in temp memory
    oldlnlk <- lnlk
    oldPath <- list(path1,path2,path3)
    oldzsigmahat <- zsigmahat
    oldzhat <- zhat
    oldehat <- ehat
    
    # recalculating values
    # path ends up with state 1:
    lnlk[1,j]   <- log(eta[case[1],j]+1e-50) + oldlnlk[case2[1],j-1]
    path1[,1:j] <- cbind(oldPath[[case2[1]]][,1:(j-1)],c(1,0,0))
    zsigmahat[1,1:j] <- c(oldzsigmahat[case2[1],1:(j-1)], w[j] - mu[1] - tempsMUZ[case[1]])
    zhat[1,1:j] <- c(oldzhat[case2[1],1:(j-1)],zsigmahat[1,j]/(sigma1hat));
    ehat[1,1:j] <- c(oldehat[case2[1],1:(j-1)],(w[j] - tempsMUZ[case[1]] -temps[case[1]])/sigma1hat)
    
    # path ends up with state 2:
    lnlk[2,j]   <- log(eta[case[2],j]+1e-50) + oldlnlk[case2[2],j-1]
    path2[,1:j] <- cbind(oldPath[[case2[2]]][,1:(j-1)],c(0,1,0))
    zsigmahat[2,1:j] = c(oldzsigmahat[case2[2],1:(j-1)], w[j] - mu[2] - tempsMUZ[case[2]])
    zhat[2,1:j] <- c(oldzhat[case2[2],1:(j-1)],zsigmahat[2,j]/(sigma2hat));
    ehat[2,1:j] <- c(oldehat[case2[2],1:(j-1)],(w[j] - tempsMUZ[case[2]] -temps[case[2]])/sigma2hat)
    
    # path ends up with state 3:
    lnlk[3,j]   <- log(eta[case[3],j]+1e-50) + oldlnlk[case2[3],j-1]
    path3[,1:j] <- cbind(oldPath[[case2[3]]][,1:(j-1)],c(0,0,1))
    zsigmahat[3,1:j] = c(oldzsigmahat[case2[3],1:(j-1)], w[j] - mu[3] - tempsMUZ[case[3]])
    zhat[3,1:j] <- c(oldzhat[case2[3],1:(j-1)],zsigmahat[3,j]/(sigma3hat));
    ehat[3,1:j] <- c(oldehat[case2[3],1:(j-1)],(w[j] - tempsMUZ[case[3]] -temps[case[3]])/sigma3hat)
    
    # for(x in 0:2){
    #   lnlk[x+1,j]   <- log(eta[case[x+1],j]+1e-50) + oldlnlk[case2[x+1],j-1]
    #   paths[(3*x+1):(3*x+3),1:j] <- cbind(oldPath[(3*x+1):(3*x+3),1:(j-1)],(1:3 == (x+1)) * 1)
    #   zsigmahat[x+1,1:j] <- c(oldzsigmahat[case2[x+1],1:(j-1)],  w[j] - mu[case2[x+1]] - tempsMUZ[case[x+1]])
    #   zhat[x+1,1:j] <- c(oldzhat[case2[x+1],1:(j-1)], zsigmahat[x+1,j]/(sigma[x+1]))
    #   ehat[x+1,1:j] <- c(oldehat[case2[x+1],1:(j-1)], (w[j] - tempsMUZ[case[x+1]] - temps[case[x+1]])/sigma[x+1]) # ONEMLI!!!
    # }
    
    j=j+1
  }
  #   lnlipath1 = t(lnlk1);
  #   lnlipath2 = t(lnlk2);
  lnlkRes  <- lnlk[which(max(lnlk[,size]) == lnlk[,size])[1],size]
  bestpath <- get(paste0("path", which(max(lnlk[,size]) == lnlk[,size])[1]))
  
  # if(lnlk[1,size] > lnlk[2,size]){
  #   lnlkRes = lnlk[1,size]} else {lnlkRes = lnlk[2,size]}
  # 
  return(bestpath)
}
