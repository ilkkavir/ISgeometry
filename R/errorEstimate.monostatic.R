# 
# Error table calculationos for converting the ACF noise levels into
# plasma parameter error estimates. 
# 
# Reguires the ISspectrum pacakge
# 
# IV 2010-2023
# 
#  
#
# 
# 




divide <- function(x){return(x[2]/x[1])}

errorEstimate.monostatic <- function(p=c(1e11,300,1,0,0,.3),pm0=c(30.5,16),fradar=933.5e6,ind=seq(6),zeroLag=T,nLag=60,llag=933.5e6/fradar*1e-5,freq=seq(-100000,100000,by=10)*fradar/933.5e6,dp=1e-5,printLatex=T,noiseLevel=.01){
#
# Calculate the lookup table of plasma parameter errors for given plasma parameters, 
# radar carrier frequency, and noise level.
#
# INPUT:
#
#  p           plasma parameters: Ne [m^-3], Ti [K], Te/Ti, coll [s^-1], Vi [ms^-1], composition (O+/Ne if last element of pm0=16)
#  pm0         Ion masses in amu
#  fradar      Radar carrier frequency [Hz]
#  ind         Indices of plasma parameters to include in the parameter combinations.
#              For example, c(1,1,0,0,0,0) would select only Ne and Ti
#  zeroLag     Logical. If TRUE, also the zero-lag is included in the ACf
#  nLag        Number of lags
#  llag        Lag separation [s]
#  freq        Frequency axis in spectrum calculation [Hz]
#  printLatex  Logical. Should the error table be printed as a readily formatted LaTex table?
#  noiseLevel  The reference ACF noise level used in the calculations
#
# OUTPUT:
#   errtab (invisible) The final error table
#
#
#
#
# 
# 

    
  # scattering wave number
  kscatt <- 4*pi*fradar/299792458

  # frequency scale, THIS MUST BE THE SAME SCALE THAT IS USED IN ISspectrum!
  om0 <- kscatt*sqrt(2*1.3806503e-23/pm0[1]/1.66053886e-27)

  # parameter scales
  pscale              <- p
  pscale[4]           <- om0/2/pi
  pscale[5:length(p)] <- 1

  # length of frequency axis
  nf <- length(freq)

  # step size in frequency axis
  fstep <- mean(diff(freq))

  # number of plasma parameters
  np <- length(p)

  # finite differences to the parameters
  pdiff <- matrix(0,ncol=np,nrow=np)
  for(k in seq(1,np)){
    pdiff[k,]  <- p
    pdiff[k,k] <- pdiff[k,k] + dp*pscale[k]
  }

  # matrix for spectra
  s <- matrix(nrow=(np+1),ncol=nf)

    s[1,] <- ISspectrum::ISspectrum.guisdap(freq=freq,p=p,pm0=pm0,fradar=fradar)
    plot(s[1,])
  for(k in seq(2,(np+1))){
    s[k,] <- ISspectrum::ISspectrum.guisdap(freq=freq,p=pdiff[(k-1),],pm0=pm0,fradar=fradar)
  }

  # normalization to dimensionless units
  s <- s * om0 / p[1] / 9.978688e-29

  # autocorrelation function
  acf <- matrix(nrow=(np+1),ncol=nLag)
  tau <- (seq(nLag)-ifelse(zeroLag,1,0))*llag
  for(i in seq(np+1)){
    for(j in seq(nLag)){
      acf[i,j] <- sum(exp(1i*2*pi*freq*tau[j])*s[i,]) * fstep*2*pi/om0 
    }
  }

    plot(Re(acf[1,]))
    
  # ACF differences
  dacf <- t(acf[2:(np+1),] )
  for (k in seq(np)) dacf[,k] <- (dacf[,k] - acf[1,])/dp

  # linear theory matrix
  A <- matrix(0,ncol=np,nrow=2*nLag)
  A[1:nLag,] <- Re(dacf)
  A[(nLag+1):(2*nLag),] <- Im(dacf)

  # the frequency scale of Vallinkoski 1988
  om1 <- kscatt*sqrt(2*1.3806503e-23*p[2]/pm0[1]/1.66053886e-27)

  # time-lag scale
  tau0 <- pi/2/om1

  # number of lags within tau0
  Ntau <- round(tau0/llag)

  # variance of lagged products
  lpvar <- noiseLevel**2/4*Ntau

  Sigma <- diag(ncol=(2*nLag),nrow=(2*nLag),lpvar)

  # Fisher information 
  Q <- t(A)%*%solve(Sigma)%*%A

  # number of parameter indices
  nind <- length(ind)

  # number of rows in the final error table
  nr <- 2**nind-1

  errtab <- matrix(0,ncol=nind+1,nrow=nr)

  indtab <- matrix(F,ncol=nind,nrow=nr)

  for(i in seq(nind)){
    indtab[,nind-i+1] <- rev(rep(rep(c(T,F),each=2**(i-1)),length.out=nr))
  }


  # solve the standard deviations
  for(k in seq(nr)){
    errtab[k,which(indtab[k,])] <- sqrt(diag(solve(Q[ind[which(indtab[k,])],ind[which(indtab[k,])]])));
    errtab[k,nind+1] <- log10(
                          divide(
                              sort(
                               range(
                                 eigen(
                                   Q[ind[which(indtab[k,])],ind[which(indtab[k,])]],only.values=T,symmetric=T,EISPACK=T)$values
                                       )
                                     )
                                   )
                                 )

  if(printLatex){
    for(i in seq(nind)){
      if(indtab[k,i]){
        cat(sprintf(' %5.3f &',errtab[k,i]))
      }else{
        cat('       &')
      }
    }
    cat(sprintf(' %5.3f',errtab[k,nind+1]))
    cat('\\')
    cat('\\')
    cat('\n')
    }
  }

  errtab[!indtab] <- NA

  invisible(errtab)

} # errorEstimate.monostatic
