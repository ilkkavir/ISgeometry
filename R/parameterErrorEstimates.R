parameterErrorEstimates <- function(lat,lon,alt,Ne,Ti,Te,Vi,Coll,Comp,fwhmRange,resR,intTime,pm0=c(30.5,16),hTeTi=110,Tnoise=300,Pt=3.5e6,locTrans=SKI,locRec=list(SKI,KAR,KAI),fwhmTrans=2.1,fwhmRec=c(1.2,1.7,1.7),RXduty=1,mineleTrans=30,mineleRec=30,fradar=233e6,phArrTrans=TRUE,phArrRec=TRUE,fwhmIonSlab=100,dutyCycle=.25){
    #
    # Calculate plasma parameter error estimates in the given point with the given radar system
    # and plasma parameters
    #
    # NOTICE: if we do not make the lookup tables, we can actually directly calculate the errors for
    #         a correct noise level!
    #         => do we even need to define the concept of ACF noise level?!?
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
    # the plasma parameter vector
    p <- c(Ne,Ti,Te,Coll,Vi,Comp)
    # errors with the first noise level
    # try to select a good frequency grid
    if(Coll<3e4){
        freq <- seq(-500,500,by=5)*1e4/tau0/log10(max(Coll,10))
    }else{ # with high collision frequencies the spectrum is very narrow (is this good enough?)
        freq <- seq(-100,100,by=1)*2e3*Ti/Coll
    }
    if(alt>hTeTi){
        parinds <- c(1,2,3,5)
    }else{
        parinds <- c(1,2,5)
    }
    refErrs <- parameterFitErrors(noiseLevel=nlevRef,p=p,pm0=pm0,fradar=fradar,ind=parinds,zeroLag=F,nLag=round(tau0/(dtau*1e6)*10),llag=dtau,freq=freq)
    # scale the errors with the noise levels at the other sites
    for(ii in seq(length(noiseLevels$noiseLevel.site))){
        # we have noise levels per one bit, scale to noise levels after integration
        nlev <- noiseLevels$noiseLevel.site[[ii]]$noiseLevel$noiseLevel[1,1,1]*sqrt(dtau/dutyCycle/intTime)
        errtabs$los[[ii]] <- nlev / nlevRef * refErrs
    }
    # multistatic analysis, 
    nlevmultis <- noiseLevels$noiseLevel.isotropic$noiseLevel[1,1,1]*sqrt(dtau/dutyCycle/intTime)
    nlevVel <- noiseLevels$noiseLevel.velocity$noiseLevel[1,1,1]*sqrt(dtau/dutyCycle/intTime)
    errtabs$multistatic <- nlevmultis / nlevRef * refErrs
    errtabs$multistatic[5] <- nlevVel / nlevRef * refErrs[5]

    return(errtabs)
    
}

