parameterErrorEstimatesPrior <- function(lat,lon,alt,Ne,Ti,Te,Coll,Comp,fwhmRange,resR,intTime,NePr=1e12,TiPr=1e4,TePr=1e4,CollPr=0,ViPr=1e4,CompPr=0,pm0=c(30.5,16),hTeTi=110,Tnoise=300,Pt=3.5e6,locTrans=SKI,locRec=list(SKI,KAR,KAI),fwhmTrans=2.1,fwhmRec=c(1.2,1.7,1.7),RXduty=c(.75,1,1),mineleTrans=30,mineleRec=30,fradar=233e6,phArrTrans=TRUE,phArrRec=TRUE,fwhmIonSlab=50,dutyCycle=.25){
    #
    # Calculate plasma parameter error estimates in the given point with the given radar system
    # and plasma parameters
    #
    # INPUTS:
    #   lat          Geocentric latitude of the measurement volume [deg north]
    #   lon          Geocentric longitude of the measurement volume [deg east]
    #   alt          Altitude of the measurement volume [km]
    #   Ne           Electron density [m^-3]
    #   Ti           Ion temperature [K]
    #   Te           Electron temperature [K]
    #   Coll         Ion-neutral collision frequency [Hz]
    #   Comp         Ion composition (fraction of ion with mass pm0[2] out of the total ion number density)
    #   fwhmRange    Range resolution in ACF decoding [km]
    #   resR         Range resolution in the plasma parameter fit [km]
    #   intTime      Integration time [s], default 10
    #   pm0          Ion masses [amu], default c(30.5,16)
    #   hTeTi        Altitude, below which Te=Ti is assumed, [km], default 110
    #   Tnoise       Receiver noise temperature [K], default 300
    #   Pt           Transmitter power(s) [W], default 3.5e6
    #   locTrans     Transmitter location(s), list of lat, lon, [height] in degress [km], default list(c(69.34,20.21))
    #   locRec       Receiver locations, list of lat, lon, [height] in degrees [km], default list(c(69.34,20.21),c(68.48,22.52),c(68.27,19.45))
    #   fwhmTrans    Transmitter beam width(s) (at zenithg for phased-arrasy), [full width at half maximu, degrees], default c(2.1)
    #   fwhmRec      Receiver beam widhts (at zenith for phased-arrays) [full width at half maximum, degrees], default c(1.2,1.7,1.7)
    #   RXduty       "Receiver duty cycle" for each receiver, default c(.75,1,1)
    #   mineleTrans  Elevation limit of the transmitter [deg], default c(30)
    #   mineleRec    Elevation limit of the receivers [deg], default c(30,30,30)
    #   fradar       Radar system carrier frequency [Hz], default 233e6
    #   phArrTrans   Logical is (are) the transmitter(s) phased array(s)? Default c(T)
    #   phArrRec     Logical, are the receivers phased-arrays? A vector with a value for each receiver. Default c(T,T,T)
    #   fwhmIonSlab  Thickness of the ion slab that causes self-noise [km], default 50
    #   dutyCycle    Transmitter duty cycle, default 0.25
    #
    #
    # OUTPUT:
    #   errtabs      A list of error vectors
    #                  errtabs$los         is a list of error vectors for each individual site. 
    #                                      In the same order in which the receivers are listed in locRec, fwhmRec, etc. 
    #                  errtabs$multistatic is an error vector for "multistatic" analysis that merges data from all receivers.
    #                                      The velocity error estimate is the square root of the trace of the error covariance
    #                                      matrix of three orthogonal velocity vector components, and thus gives an upper limit 
    #                                      that the standard deviation of any projection of the velocity vector cannot exceed. If
    #                                      errors of all vector components are equal, the true error is sqrt(3) times the given error.
    #
    # IV 2023
    #

    # the ACF scale (we could do the calculations without defining acf0 and the noise levels, but this is easier to implement in ISgeometry)
    kscatt <- 4*pi*fradar/299792458
    om1 <- kscatt*sqrt(2*1.3806503e-23*Ti/pm0[1]/1.66053886e-27)
    tau0 <- pi/2/om1*1e6

    nlevRef <- 0.01

    
    if(!is.list(locTrans)){
        locTrans <- list(locTrans)
    }

    xyV <- sphericalToPlanar.geographic(lat=lat,lon=lon,zeroLatitude=locTrans[[1]][1],zeroMeridian=locTrans[[1]][2])

    
    noiseLevels <- multistaticNoiseLevels(refPoint=locTrans[[1]],locTrans=locTrans,locRec=locRec,locxy=FALSE,fwhmTrans=fwhmTrans,fwhmRec=fwhmRec,fwhmRange=fwhmRange,resR=resR,phArrTrans=phArrTrans,phArrRec=phArrRec,x=xyV$x,y=xyV$y,heights=alt,Tnoise=Tnoise,Pt=Pt,Ne=Ne,fwhmIonSlab=fwhmIonSlab,fradar=fradar,tau0=100,RXduty=RXduty,verbose=FALSE,mineleTrans=mineleTrans,mineleRec=mineleRec)

    dtau <- fwhmRange*1000/3e8*2


    errtabs <- list()
    errtabs$los <- list()
    # the plasma parameter vector. The velocity is practically independent from the other parameters, so not given as user input.
    # use 1e-3 to avoid problems at zero frequency with the simple spectrum calculation
    p <- c(Ne,Ti,Te,Coll,1e-3,Comp)
    pPr <- c(NePr,TiPr,TePr,CollPr,ViPr,CompPr)
    # errors with the first noise level
    # try to select a good frequency grid
    if(Coll<3e4){
        freq <- seq(-500,500,by=5)*1e4/tau0/log10(max(Coll,10))
    }else{ # with high collision frequencies the spectrum is very narrow (is this good enough?)
        freq <- seq(-100,100,by=1)*2e3*Ti/Coll*(fradar/233e6)**2
    }

    refErrs <- parameterFitErrorsPrior(noiseLevel=nlevRef,p=p,pPr=pPr,pm0=pm0,fradar=fradar,zeroLag=F,nLag=round(tau0/(dtau*1e6)*5),llag=dtau,freq=freq)
    # scale the errors with the noise levels at the other sites
    for(ii in seq(length(noiseLevels$noiseLevel.site))){
        # we have noise levels per one bit, scale to noise levels after integration and add contribution from Te/Ti
        nlev <- noiseLevels$noiseLevel.site[[ii]]$noiseLevel$noiseLevel[1,1,1]*sqrt(dtau/dutyCycle/intTime)*((1+Te/Ti)/2)
        errtabs$los[[ii]] <- nlev / nlevRef * refErrs
        names(errtabs$los[[ii]]) <- c('dNe','dTi','dTe','dColl','dVi','dComp')
    }
    # multistatic analysis, 
    nlevmultis <- noiseLevels$noiseLevel.isotropic$noiseLevel[1,1,1]*sqrt(dtau/dutyCycle/intTime)*((1+Te/Ti)/2)
    nlevVel <- noiseLevels$noiseLevel.velocity$noiseLevel[1,1,1]*sqrt(dtau/dutyCycle/intTime)*((1+Te/Ti)/2)
    errtabs$multistatic <- nlevmultis / nlevRef * refErrs
    errtabs$multistatic[5] <- nlevVel / nlevRef * refErrs[5]
    errtabs$multistatic <- errtabs$multistatic
    names(errtabs$multistatic) <- c('dNe','dTi','dTe','dColl','dVi','dComp')

    return(errtabs)
    
}

