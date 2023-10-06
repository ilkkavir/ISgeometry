incoherentScatterSpectrum <- function(ele=c(1e11,300,0,0),ion=list(c(30.5,.7e11,300,0,0),c(16,.3e11,300,0,0),c(1,0,300,0,0)),fradar=933e6,scattAngle=180,freq=seq(-1000,1000)*10){
#
# a general incoherent scatter spectrum calculation, with max 7  singly-charged ions.
#
# This version does not use tabulated plasma dispersion function values
#
# INPUT:
#  ele         c(Ne,Te,nu_en,ve)
#  ion         list( c(m1,N1,T1,nu_1n,v1) , c(m2,N2,T2,nu_2n,v2) , ... )
#  fradar      radar frequency
#  scattAngle  scattering angle (the angle between incident and scattered wave vectors, 180 for backscattering)
#  freq        frequency axis
#
#  ion masses in amu
#  ion densities can be given as absolute values, but they are treated as relative abundances, so that sum(Ni) = Ne
#  only singly-charged ions at the moment
#  radar frequency in Hz
#  scattering angle in degrees
#  frequency axis points in Hz
#
#  OUTPUT:
#   a vector of power spectral densities at the given frequencies
#
# I. Virtanen 2015
#


    # transform the parameters
    if(is.list(ion)){
        nIon <- length(ion)
    }else{
        ion <- list(ion)
        nIon <- 1
    }

    # Boltzmann constant
    kb <- 1.3806503e-23

    # atomic mass unit
    amu <- 1.660538921e-27

    # electron mass in amu
    me <- 0.00054857990943

    # electron charge
    q <- 1.60217657e-19

    # permittivity of free space
    eps0 <- 8.85418782e-12

    # debye length
    D <- sqrt(eps0*kb*ele[2]/q**2/ele[1])

    # speed of light
    c <- 299792458.

    # scattering wave number for backscattering
    ks0 <- 4*pi*fradar/c

    # actual scattering wave number
    ks <- ks0 * sin(scattAngle/2*pi/180)

    # (k*D)^2
    kd2 <- ks**2 * D**2

    # sum of ion densities must match with the electron density
    nion <- sapply( ion , FUN=function(x){x[2]})
    nion <- nion/sum(nion) * ele[1]

    # number of frequency points
    nf <- length(freq)

    #######################
    # electron admittance #
    #######################
    # frequency scaling factor
    fscale <- 2*pi*sqrt(me*amu/2/kb/ele[2])/ks

    # normalized angular frequency, shifted by the electron velocity
    theta <- ( freq - 2 * fradar * ele[4] / c ) * fscale

    # normalized angular collision frequency (this is already angular frequency)
    psi <- ele[3] * fscale / 2 / pi

    # plasma dispersion function values for electrons
    pdfElectrons <- sapply( theta - 1i*psi , FUN=plasmaDispersionFunction )

    # finally the actual admittances (why theta instead of (theta-1iu*psi) ?)
    admittanceElectrons <- 1i * (1 - theta * pdfElectrons / ( 1 + 1i * psi * pdfElectrons ) )

    # electron admittance divided with the angular frequency
    scaledAdmittanceElectrons <- admittanceElectrons / ( 2 * pi * ( freq - 2 * fradar * ele[4] / c ))

    ###################
    # Ion admittances #
    ###################
    pdfIons <- list()
    admittanceIons <- list()
    scaledAdmittanceIons <- list()
    for(k in seq(nIon)){
        # frequency scaling factor
        fscale <- 2*pi*sqrt(ion[[k]][1]*amu/2/kb/ion[[k]][3])/ks

        # normalized angular frequency, shifted by the ion velocity
        theta <- ( freq - 2 * fradar * ion[[k]][5] / c ) * fscale

        # normalized angular collision frequency (these are angular frequencies from start)
        psi <- ion[[k]][4] * fscale / 2 /pi

        # plasma dispersion function values for ion k
        pdfIons[[k]] <- sapply( theta - 1i*psi, FUN=plasmaDispersionFunction )

        # finally the actual admittances (why theta instead of (theta-1i*psi) ??)
        admittanceIons[[k]] <- 1i * (1 - theta * pdfIons[[k]] / ( 1 + 1i * psi * pdfIons[[k]] ) )

        # the admittance divided with angular frequency
        scaledAdmittanceIons[[k]] <- admittanceIons[[k]] / ( 2 * pi * ( freq - 2 * fradar * ion[[k]][5] / c ))

    }

    ################
    # the spectrum #
    ################

    # ion - electron density ratios
    nine <- nion / ele[1]

    # ion - electron temperature ratio
    tite <- sapply( ion , FUN=function(x){x[[3]]}) / ele[2]

    # the final spectrum
    s <- (
        abs( admittanceElectrons )**2 *  colSums( apply( sapply( scaledAdmittanceIons , FUN=function(x){Re(x)}) , FUN=function(x,s){x*s},MARGIN=1,s=nine))
        +
        abs(
            colSums( apply( sapply( admittanceIons ,FUN=function(x){x}) , FUN=function(x,s){x*s},MARGIN=1,s=nine/tite))
            +
            1i * kd2
            )**2 * Re( scaledAdmittanceElectrons )
        ) / (
            abs(
                admittanceElectrons +
                colSums(apply( sapply( admittanceIons , FUN=function(x){x}) , FUN=function(x,s){x*s},MARGIN=1,s=nine/tite))
                + 1i*kd2
                )**2
            ) * ele[1] * 2.8179402894e-15**2 / pi * 4 * pi
term1 <-         abs( admittanceElectrons )**2 *  colSums( apply( sapply( scaledAdmittanceIons , FUN=function(x){Re(x)}) , FUN=function(x,s){x*s},MARGIN=1,s=nine))
term2 <-         abs(
            colSums( apply( sapply( admittanceIons ,FUN=function(x){x}) , FUN=function(x,s){x*s},MARGIN=1,s=nine/tite))
            +
            1i * kd2
            )**2 * Re( scaledAdmittanceElectrons )
term3 <-             abs(
                admittanceElectrons +
                colSums(apply( sapply( admittanceIons , FUN=function(x){x}) , FUN=function(x,s){x*s},MARGIN=1,s=nine/tite))
                + 1i*kd2
                )**2


    return(list(s=s,term1=term1,term2=term2,term3=term3,scaledAdmittanceIons=scaledAdmittanceIons,admittanceIons=admittanceIons,pdfIons=pdfIons))

} # incoherentScatterSpectrum
