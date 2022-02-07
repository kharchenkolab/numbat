#'@useDynLib numbat
#'@import Rcpp
NULL


## refer to https://github.com/evanbiederstedt/poilogcpp

#' @keywords internal  
rbipoilog <- function(S,mu1,mu2,sig1,sig2,rho,nu1=1,nu2=1,condS=FALSE,keep0=FALSE){
  
  sim <- function(nr){
    lamx <- rnorm(nr)
    lamy <- rho*lamx+sqrt(1-rho^2)*rnorm(nr)
    x <- rpois(nr,exp(sig1*lamx+mu1+log(nu1)))
    y <- rpois(nr,exp(sig2*lamy+mu2+log(nu2)))
    sel <- x+y 
    xy <- cbind(x,y)
    if (keep0) xy <- cbind(x,y) else xy <- cbind(x[sel>0],y[sel>0])
    return(xy)
  }
  
   
  if (S<1) stop('S is not positive')
  if (!is.finite(S)) stop('S is not finite')
  if ((S/trunc(S))!=1) stop('S is not an integer')
  if (sig1<0) stop('sig1 is not positive')
  if (sig2<0) stop('sig2 is not positive')
  if (nu1<0)    stop('nu1 is not positive')
  if (nu2<0)    stop('nu2 is not positive')
  if (-1>rho | rho>1) stop('rho must be in the range -1 to 1')
  if (!is.logical(keep0)) stop('keep0 should be TRUE or FALSE')
  if (length(c(mu1,mu2,sig1,sig2,rho,nu1,nu2))>7) stop('duplicate values of some parameters')
   
  if (condS){
    simMat <- matrix(NA,0,2)
    fac <- 2
    nr  <- S
    while (nrow(simMat)<S){
      simvals <- sim(nr*fac)
      simMat  <- rbind(simMat,simvals)
      fac     <- (1/(nrow(simvals)/(nr*fac)))*2
      fac     <- ifelse(is.finite(fac),fac,1000)
      nr      <- S-nrow(simvals)
    }
    simMat <- simMat[1:S,]
  }
   
  else simMat <- sim(S)
  return(simMat)
}

#' @keywords internal  
rpoilog <- function(S,mu,sig,nu=1,condS=FALSE,keep0=FALSE){
   
  sim <- function(nr){
    lamx <- rnorm(nr)
    x <- rpois(nr,exp(sig*lamx+mu+log(nu)))
    if (!keep0) x <- x[x>0]
    return(x)
  }
   
  if (S<1) stop('S is not positive')
  if (!is.finite(S)) stop('S is not finite')
  if ((S/trunc(S))!=1) stop('S is not an integer')
  if (sig<0) stop('sig is not positive')
  if (nu<0)    stop('nu is not positive')
   
  if (condS){
    simVec <- vector('numeric',0)
    fac <- 2
    nr  <- S
    while (length(simVec)<S){
      simvals <- sim(nr*fac)
      simVec  <- c(simVec,simvals)
      fac     <- (1/(length(simvals)/(nr*fac)))*2
      fac     <- ifelse(is.finite(fac),fac,1000)
      nr      <- S-length(simvals)
    }
    simVec <- simVec[1:S]
  }
   
  else simVec <- sim(S)
  return(simVec)
}

#' @keywords internal  
dbipoilog <- function(n1,n2,mu1,mu2,sig1,sig2,rho){
  
  if (length(n1)!=length(n2)) stop('n1 and n2 have unequal length')
  if (any((n1[n1!=0]/trunc(n1[n1!=0]))!=1)) stop('all values of n1 must be integers')
  if (any((n2[n2!=0]/trunc(n2[n2!=0]))!=1)) stop('all values of n2 must be integers')
  if (any(n1<0)) stop('one or several values of n1 are negative')
  if (any(n2<0)) stop('one or several values of n2 are negative')
  
  if (length(c(mu1,mu2,sig1,sig2,rho))>5) stop('vectorization of mu, sig and rho is currently not implemented')
  if (!all(is.finite(c(mu1,mu2,sig1,sig2,rho)))) stop('all parameters should be finite')
  if (-1>rho | rho>1) stop('rho must be in the range -1 to 1')
  if (sig1<0 | sig2<=0) stop('sig1 and/or sig2 is not larger than 0')
  
  poilog2(x=as.integer(n1), y=as.integer(n2), my1=as.double(mu1), my2=as.double(mu2),
    sig1=as.double(sig1^2), sig2=as.double(sig2^2), ro=as.double(rho), nrN=as.integer(length(n1)))
}




#' Returns the density for the Poisson lognormal distribution with parameters mu and sig
#' 
#' @param x vector of integers, the observations
#' @param mu mean of lognormal distribution
#' @param sig standard deviation of lognormal distribution
#' @param log boolean Return the log density if TRUE (default=FALSE)
#' @return NULL
#' @export
dpoilog <- function(x, mu, sig, log=FALSE){
  if (!(length(x) == length(mu) & length(x) == length(sig))) stop('All parameters must be same length') 
  if (any((x[x!=0]/trunc(x[x!=0]))!=1)) stop('all x must be integers')
  if (any(x<0)) stop('one or several values of x are negative')
  if (!all(is.finite(c(mu,sig)))) stop('all parameters should be finite')
  if (any(is.na(c(x,mu,sig)))) stop('Parameters cannot be NA')
  if (any(sig<=0)) {
      stop(c('sig is not larger than 0', unique(sig)))
  }

  p = poilog1(x=as.integer(x), my=as.double(mu), sig=as.double(sig^2), nrN=as.integer(length(x)))

  p[p == 0] = 1e-15

  if (log) {
    return(log(p))
  } else {
    return(p)
  }
}

#' @keywords internal  
poilogMLE <- function(n,startVals=c(mu=1,sig=2),nboot=0,zTrunc=TRUE,
                     method='BFGS',control=list(maxit=1000)){
  
  if (is.matrix(n) | (is.data.frame(n))) stop(paste('n has',ncol(n),'colums, supply a vector or use function bipoilogMLE',sep=' ')) 
  if (length(startVals)!=2) stop('length of startVals is not 2') 
  if (startVals[2]<=0) stop('start value of sig2 is not larger than 0')
  startVals[2] <- log(startVals[2])
  un <- unique(n)
  nr <- rep(NA,length(un))
  for (i in 1:length(un)){ nr[i] <- sum(n%in%un[i]) }
  
  lnL <- function(z,nr){
    if (z[2] < (-372)) z[2] <- -372
    if (z[2] >    354) z[2] <-  354
    b <- 0
    if (zTrunc) b <- log(1-dpoilog(0,z[1],exp(z[2])))
    logL <- -sum((log(dpoilog(un,z[1],exp(z[2])))-b)*nr)
    return(logL)
  }
  
  fit <- optim(startVals,lnL,nr=nr,control=control,method=method)
  
  if (fit$convergence!=0){
    if (fit$convergence==1) stop('the iteration limit has been reached!   try different startVals or increase maxit') 
    if (fit$convergence==10) stop('degeneracy of the Nelder Mead simplex ....')
    else stop(paste('unknown error in optimization', fit$message))
  } 
  
  fit$par <- c(as.numeric(fit$par),1-dpoilog(0,fit$par[1],exp(fit$par[2])))
  est <- list('par'=c('mu'=fit$par[1],'sig'=exp(fit$par[2])),'p'=fit$par[3],'logLval'=-fit$value,'gof'=NULL,boot=NULL)
  
  if (nboot>0){
    message(paste('estimates: mu: ',round(fit$par[1],3),' sig ',round(exp(fit$par[2]),3),sep=''),'\n')
    message('********     bootstrapping    ********\n')
    bMat <- matrix(NA,nboot,3)
    colnames(bMat) <- c('mu','sig2','logLval')
    count <- 0
    kat  <- seq(0,nboot,by=100)
    bStartVals <- fit$par
    while (count<nboot){
      bfit <- un <- nr <- NA
      # simulations are conditonal on the number of species in the observed data set :
      sim <- rpoilog(length(n),fit$par[1],exp(fit$par[2]),condS=TRUE,keep0=!zTrunc)
      un <- unique(sim)
      nr <- rep(NA,length(un))
      for (i in 1:length(un)){ nr[i] <- sum(sim%in%un[i]) }
      bfit <- try(optim(bStartVals,lnL,nr=nr,control=control,method=method),silent=TRUE)
      if (class(bfit)!='try-error'){
        count <- count+1
        bMat[count,] <- c(bfit$par[1],exp(bfit$par[2]),-bfit$value)
        if (count%in%kat) cat('   boot',count,'of',nboot,'\n')
      }
    }
    est$boot <- data.frame(bMat)
    est$gof  <- which(sort(c(est$logLval,bMat[,3]))==est$logLval)/nboot
  }  
  return(est)   
}


#' @keywords internal  
bipoilogMLE <- function(n1,n2=NULL,startVals=c(mu1=1,mu2=1,sig1=2,sig2=2,rho=0.5),
                       nboot=0,zTrunc=TRUE,file=NULL,
                       method='BFGS',control=list(maxit=1000)){
  if (is.null(n2)){
    if (!is.matrix(n1) & !is.data.frame(n1)) {stop('n2 is missing and n1 is not a matrix or data frame') }
    if (ncol(n1)!=2) stop('n1 should have 2 columns when n2 is not given')
    n2 <- n1[,2]
    n1 <- n1[,1]
  }
  if (!is.null(n2) & (is.matrix(n1) | is.data.frame(n1))) stop('n1 is a matrix and n2 is specified')
  if (length(n2)!=length(n2)) stop('n1 and n2 must be of equal length')
  
  z <- startVals
  if (length(z)!=5) stop('length of startVals is not 5') 
  if (z[3]<=0) stop('start value of sig1 is not larger than 0')
  if (z[4]<=0) stop('start value of sig2 is not larger than 0')
  if (-1>z[5] | z[5]>1) stop('start value of rho must be in the range -1 to 1')
  
  invRhoTrans <- function(w) {
     rho <- (1-exp(-1*w))/(1+exp(-1*w))
     if (rho>0.9999) rho <- 0.9999
     if (rho<(-0.9999)) rho <- -0.9999
     return(rho)
  } 
  rhoTrans <- function(rho) 1*log((1+rho)/(1-rho))
  
  z[3:4] <- log(z[3:4])
  z[5]   <- rhoTrans(z[5])

  xy <- cbind(n1,n2)
  un <- unique(xy)
  nr <- rep(NA,nrow(un))
  for (i in 1:nrow(un)){ nr[i] <- sum(apply(xy,1,function(x) x[1]==un[i,1] & x[2]==un[i,2])) }
  
  lnL <- function(z,nr){
    b <- 0
    if (z[3] < (-372)) z[3] <- -372
    if (z[4] < (-372)) z[4] <- -372
    if (z[3] >   354)  z[3] <-  354
    if (z[4] >   354)  z[4] <-  354
    if (zTrunc) b <- log(1-dbipoilog(0,0,z[1],z[2],exp(z[3]),exp(z[4]),invRhoTrans(z[5])))
    logL <- -sum((log(dbipoilog(un[,1],un[,2],z[1],z[2],exp(z[3]),exp(z[4]),invRhoTrans(z[5])))-b)*nr)
    return(logL)
  }
  
  fit <- optim(z,lnL,nr=nr,control=control,method=method)
  
  if (fit$convergence!=0){
    if (fit$convergence==1) stop('the iteration limit has been reached!   try different startVals or increase maxit') 
    if (fit$convergence==10) stop('degeneracy of the Nelder Mead simplex ....')
    else stop(paste('unknown error in optimization', fit$message))
  }
  
  fit$par <- c(as.numeric(fit$par),1-dpoilog(0,fit$par[1],exp(fit$par[3])),1-dpoilog(0,fit$par[2],exp(fit$par[4])))
  est <- list('par'=c('mu1'=fit$par[1],'mu2'=fit$par[2],'sig1'=exp(fit$par[3]),
         'sig2'=exp(fit$par[4]),'rho'=invRhoTrans(fit$par[5])),'p'=c(fit$par[6],fit$par[7]),'logLval'=-fit$value,'gof'=NULL,boot=NULL)
         
  if (nboot>0){
    cat(paste('estimates: mu1:',round(fit$par[1],3), ' mu2:',round(fit$par[2],3),
                ' sig1:',round(exp(fit$par[3]),3), ' sig2:',round(exp(fit$par[4]),3),
                ' rho:',  round(invRhoTrans(fit$par[5]),3), sep=''),'\n')
    cat('********     bootstrapping    ********\n')
    bMat <- matrix(NA,nboot,6)
    colnames(bMat) <- c('mu1','mu2','sig1','sig2','rho','logLval')
    count <- 0
    bStartVals <- fit$par
    while (count<nboot){
      bfit <- un <- nr <- NA
      # simulations are conditonal on the number of species in the observed data set :
      sim <- rbipoilog(length(n1),fit$par[1],fit$par[2],exp(fit$par[3]),exp(fit$par[4]),invRhoTrans(fit$par[5]),condS=TRUE,keep0=!zTrunc)
      un <- unique(sim)
      nr <- rep(NA,nrow(un))
      for (i in 1:nrow(un)){ nr[i] <- sum(apply(sim,1,function(x) x[1]==un[i,1] & x[2]==un[i,2])) }
      bfit <- try(optim(bStartVals,lnL,nr=nr,control=control,method=method),silent=TRUE)
      if (class(bfit)!='try-error'){
        count <- count+1
        bMat[count,] <- c(bfit$par[1],bfit$par[2],exp(bfit$par[3]),exp(bfit$par[4]),invRhoTrans(bfit$par[5]),-bfit$value)
        cat('   boot',count,'of',nboot,': ',c(bfit$par[1],bfit$par[2],exp(bfit$par[3]),
                                     exp(bfit$par[4]),invRhoTrans(bfit$par[5])),'\n')
      }
      if (!is.null(file)) write.table(bMat[1:count,],file,sep='\t')
    }
    est$boot <- data.frame(bMat)
    est$gof  <- which(sort(c(est$logLval,bMat[,6]))==est$logLval)/nboot
  }  
  
  return(est)   
}

 