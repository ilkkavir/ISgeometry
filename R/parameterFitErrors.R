parameterFitErrors <- function(noiseLevel=.01,p=c(1e11,300,300,0,0,.3),pm0=c(30.5,16),fradar=233e6,ind=seq(6),zeroLag=FALSE,nLag=60,llag=233e6/fradar*1e-4,freq=seq(-100000,100000,by=100)*fradar/933.5e6,dp=1e-5){ 
#
# Calculate plasma parameter errors for the given ACF noise level plasma parameters, 
# radar carrier frequency.
#
# INPUT:
#
#  noiseLevel  The ACF noise level used in the calculations
#  p           plasma parameters: Ne [m^-3], Ti [K], Te [K], coll [s^-1], Vi [ms^-1], composition (O+/Ne if last element of pm0=16)
#  pm0         Ion masses in amu
#  fradar      Radar carrier frequency [Hz]
#  ind         Indices of plasma parameters to include in the parameter combinations.
#              For example, c(1,1,0,0,0,0) would select only Ne and Ti
#  zeroLag     Logical. If TRUE, also the zero-lag is included in the ACf
#  nLag        Number of lags
#  llag        Lag separation [s]
#  freq        Frequency axis in spectrum calculation [Hz]
#
# OUTPUT:
#   errrors    Standard deviations of the fited plasma parameters, NaN for those parameters that are not fited (ind=0)
#
#
#
# IV 2023
# 
# 

    # the simple spectrum calculation produces NaN at exact zero frequencies
    if(p[5]==0){p[5]<-1e-3}
    
    
    # scattering wave number
    kscatt <- 4*pi*fradar/299792458

    # frequency scale, THIS MUST BE THE SAME SCALE THAT IS USED IN ISspectrum!
    om0 <- kscatt*sqrt(2*1.3806503e-23/pm0[1]/1.66053886e-27)

    # parameter scales for finite differences
    pscale              <- p
    pscale[4]           <- om0/2/pi
    pscale[5:length(p)] <- 1

    # length of frequency axis
    nf <- length(freq)

    # step size in frequency axis
    fstep <- mean(diff(freq))

    # number of plasma parameters
    np <- length(p)

    # number of parameters to fit
    nfit <- length(ind)
    
    # finite differences to the parameters
    pdiff <- matrix(0,ncol=np,nrow=nfit)
    for(k in seq(1,nfit)){
        pdiff[k,]  <- p
        pdiff[k,ind[k]] <- pdiff[k,ind[k]] + dp*pscale[ind[k]]
    }
    
    # matrix for spectra
    s <- matrix(nrow=(nfit+1),ncol=nf)

#    pguisdap  <- p
#    pguisdap[3] <- p[3]/p[2]
#    s[1,] <- ISspectrum::ISspectrum.guisdap(freq=freq,p=pguisdap,pm0=pm0,fradar=fradar)
    s[1,] <- ISspectrumSimple(p=p,pm0=pm0,fradar=fradar,freq=freq)$s
    if(  (s[1,1]>.05*max(s[1,])) | (s[1,nf]>.05*max(s[1,])) | (sum(s[1,]>.1*max(s[1,]))<.05*nf) ){
        warning('The frequency axis is poorly selected, results may be inaccurate')
    }
    for(k in seq(2,(nfit+1))){
#        pguisdap  <- pdiff[(k-1),]
#        pguisdap[3] <- pdiff[(k-1),3]/pdiff[(k-1),2]
#        s[k,] <- ISspectrum::ISspectrum.guisdap(freq=freq,p=pguisdap,pm0=pm0,fradar=fradar)
        s[k,] <- ISspectrumSimple(p=pdiff[(k-1),],pm0=pm0,fradar=fradar,freq=freq)$s
    }

    # normalization to dimensionless units
    s <- s * om0 / p[1] / 9.978688e-29
#plot(freq,s[1,])
    # autocorrelation function
    acf <- matrix(nrow=(nfit+1),ncol=nLag)
    tau <- (seq(nLag)-ifelse(zeroLag,1,0))*llag
    for(i in seq(nfit+1)){
        for(j in seq(nLag)){
            acf[i,j] <- sum(exp(1i*2*pi*freq*tau[j])*s[i,]) * fstep*2*pi/om0
        }
    }
    # zero-lag power
    acf0  <- sum(s[1,]) * fstep*2*pi/om0
#plot(Re(acf[1,]))
#    print(max(abs(acf[1,])))
    
    # ACF differences
    dacf <- t(acf[2:(nfit+1),] )
    for (k in seq(nfit)) dacf[,k] <- (dacf[,k] - acf[1,])/dp

    # linear theory matrix
    A <- matrix(0,ncol=nfit,nrow=2*nLag)
    A[1:nLag,] <- Re(dacf)
    A[(nLag+1):(2*nLag),] <- Im(dacf)

    # the frequency scale of Vallinkoski 1988
    om1 <- kscatt*sqrt(2*1.3806503e-23*p[2]/pm0[1]/1.66053886e-27)

    # time-lag scale
    tau0 <- pi/2/om1

    # number of lags within tau0
    Ntau <- round(tau0/llag)

    # variance of lagged products
    lpvar <- noiseLevel**2*Ntau*acf0^2
#   lpvar <- noiseLevel**2/4*Ntau
#    lpvar <- noiseLevel**2*Ntau

    Sigma <- diag(ncol=(2*nLag),nrow=(2*nLag),lpvar)

    # Fisher information 
    Q <- t(A)%*%solve(Sigma)%*%A

    # parameter errors 
    errtab <- p*NA
    errtab[ind] <- sqrt(diag(solve(Q)))
    
    return(errtab*pscale)

} 
