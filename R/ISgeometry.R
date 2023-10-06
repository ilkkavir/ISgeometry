#
# R routines for measurement geometry evaluation
#
#
#
# I. Virtanen 2010, 2011
#
# M. Orispaa 2012
#   - Modifications to beam widening (May 2012)

# 20120605
# Modified levels. Now distinct levels for isotropic and velocity maps, iso.levels and vel.levels


# radius of Earth
EarthRadius <- function() return(6371)

# infinity when defining beam shapes
defaultInfinity <- function() return(1e3)

# some locations

# Tromso
TRO <-  c(69.58,19.23)

# Kiruna
KIR <- c(67.87,20.43)

# Sodankyla
SOD <- c(67.37,26.63)

# Andoya
#AND <- c(69.13,15.86)

# Andoya (Thomas)
AND <- c(69.31,16.12)

# Esrange
ESRA  <- c(67.88,21.12)

# Lokkavaara
LOK <- c(67.81,27.69)

# Muonio
MUO <- c(67.96,23.68)

# Kilpisjarvi
KIL <- c(69.02,20.86)

# Säytsjärvi
SAY <- c(69.35,27.20)

# Suolojavrras
SUO <- c(69.57,23.57)

# Miellejokk
MIE <- c(68.34,18.96)

# Parivierra
PAR <- c(67.06,19.59)

# Järämä
JAR <- c(68.38,21.08)

# Kautokeino
KAU <- c(69.01,23.04)

# Overbygd
OVE <- c(69.03,19.29)

# Abisko
ABI <- c(68.37,18.75)

## # Skibotn (by Cesar)
## SKI <- c(69 + 35/60,20 + 28/60)

# Coordinates by Thomas

# Longyearbyen
ESR <- c(78+ 9/60+11/3600,16+ 1/60+44/3600)

# Sekkujärvi
SKJ <- c(68.108444 ,20.996278)

# Järämä
JRM <- c(68.389160,21.014613)

# Karasjokk
KJK <- c(69.455051422,25.5342006683)

# Alta
ALT <- c(69+51.6/60 ,22+57.6/60)

# Masi
MAS <- c(69+26.51/60,23+38.884/60)

# Altaelva
AEV <- c(69+11.514/60,23+33.940/60)

# Hetta
HTA <- c(68+23.348/60,23+36.196/60)

# Palojoensuu
PLJ <- c(68+16.994/60,23+ 4.807/60 )

## # Karesuvanto
## KRS <- c(68+26.960/60 ,22+28.995/60)

# Peera
PEE <- c(68+53.081/60 ,21+ 3.452/60)

# updated E3D locations from Google maps 20231004
SKI <- c(69.34,20.31)
KAR <- c(68.48,22.52)
KAI <- c(68.27,19.45)



#############################################################
############# coordinate transformations ####################
#############################################################


sphericalToCartesian <- function(loc,degrees=T){
#
# conversion between spherical and Cartesian coordinate systems
#
# INPUT:
#
# loc      position in spherical coordinates c(elevation,azimuth,radius), e.g. c(latitude,longitude,earth_radius)
# degrees  logical, TRUE if the angles azimuth and elevation are in degrees, FALSE if radians are used
#
# OUPUT:
#
# a vector c(x,y,z) of the corresponding Cartesian coordinates
#

    if(degrees) loc[1:2] <- loc[1:2]/180*pi

    r <- ifelse(length(loc)<3,EarthRadius(),loc[3])

    x <- r*cos(loc[1])*cos(loc[2])
    y <- r*cos(loc[1])*sin(loc[2])
    z <- r*sin(loc[1])

    return(c(x,y,z))

} # sphericalToCartesian




#############################################################
####### some vector operations in Cartesian system ##########
#############################################################




vectorAngle.cartesian <- function(vec1,vec2,vecn=NULL,degrees=T){
#
# angle between vectors
#
# INPUT:
#
# vec1     c(x,y,z) vector 1 in Cartesian coordinates
# vec2     c(x,y,z) vector 2 in Cartesian coordinates
# vecn     c(x,y,z) normal of the plane of vec1 and vec2, determines the positive direction of rotation
# degrees  logical, if TRUE, the output is in degrees, otherwise in radians
#
# OUTPUT:
#
# angle between the vectors vec1 and vec2
#

    if(sum(abs(vec1))==0) return(0)
    if(sum(abs(vec2))==0) return(0)

    arg <- sum(vec1*vec2) / (sqrt(sum(vec1**2) * sum(vec2**2)))

    if(abs(arg)>1){
        if((abs(arg)-1)<1e-10) arg <- sign(arg)
    }

    s <- 1
    if(!is.null(vecn)){
        if(sum(vectorProduct.cartesian(vec1,vec2)*vecn)<0) s <- -1
    }

    if (degrees) return(s*180/pi*acos( arg ))

    return(s*acos( arg ))

} # vectorAngle.cartesian


vectorProduct.cartesian <- function(x,y){
#
# Vector (crossed) product of three-dimensional vectors x and y
#
# INPUT:
#
# x   c(x1,x2,x3)  a real vector
# y   c(y1,y2,y3)  a real vector
#
# OUTPUT:
#
# c(z1,z2,z3) a three dimensional vector x cross y
#

    if(length(x)!=3) error('length(x) must be 3')
    if(length(y)!=3) error('length(y) must be 3')

    return( c(x[2]*y[3]-x[3]*y[2],
              x[3]*y[1]-x[1]*y[3],
              x[1]*y[2]-x[2]*y[1]) )

} # vectorProduct.cartesian

normalUnitVector.cartesian <- function(vec1,vec2){
#
# A unit vector perpendicular to both vec1 and vec2
#
# INPUT:
#
# vec1 c(x1,y1,z1) a vector in cartesian coordinates
# vec2 c(x2,y2,z2)
#
# OUTPUT:
#
# c(xn,yn,zn) NaN:s if vec1 and vec2 are parallel
#

    nvec <- vectorProduct.cartesian(vec1,vec2)
    nvec <- nvec / sqrt(sum(nvec**2))
    return(nvec)

} # normalUnitVector.cartesian


#############################################################
##### gaussian pulse shapes and operations with them ########
#############################################################

pulseShape.gaussian.cartesian <- function(fwhmBeam,fwhmPulse,locAnt,locTarg,phArr){
#
# The 3D Gaussian distribution of a pulse with gaussian envelope.
# The antenna beam is assumed to be a gaussian as well, with the horizontal and
# "vertical" widths given as a function of beam width at zenith and zenith angle.
# (The beam is assumed to be symmetric at zenith)
#
#
# INPUT:
#
#  fwhmBeam  full width at half maximum of the beam cross section when the beam is pointed to zenith
#  fwhmPulse full width at half maximum of the Gaussian transmission envelope, in us
#  locAnt    c(x,y,z) position of the antenna in Cartesian system,  USE sphericalToCartesian() IF NECESSARY
#  locTarg   c(x,y,z) position of te target in Cartesian system
#  phArr     logical, if TRUE, phased-array antennas with beam widening as function of tilt angle are assumed, otherwise
#             beam width does not depent on pointing direction
#
#
# OUTPUT:
#  list with elements:
#    covar       covariance matrix defining the pulse shape in EFEC coordinates
#    zenithAngle zenith angle of the radar beam
#    vecAT       vector pointing from the antenna to the target
#

    # vector pointing from the antenna to the target
    vecAT <- locTarg - locAnt

    # distance from the antenna to the target
    distAT <- sqrt(sum(vecAT**2))

    # a unit vector to the target direction is used as the z-axis when defining pulse shapes
    zPulse <- vecAT/distAT

    # normal of the beam axis, horizontal at the antenna location (this is the x-axis direction when defining beam shape)
    # If the beam was pointed to zenith, select a random z-axis (this will be rotated to the geographic system anyway)
    if(abs(vectorAngle.cartesian(locTarg,locAnt))<.001){
        xPulse <- normalUnitVector.cartesian(locAnt,runif(3))
    }else{
        xPulse <- normalUnitVector.cartesian(locAnt,vecAT)
    }

    # the y-axis is perpendicular to x and z, and we have a right-handed system
    yPulse <- normalUnitVector.cartesian(zPulse,xPulse)

    # conversion of fwhm (in degrees) to standard deviations (in radians)
    stdBeam  <- fwhm2std(fwhmBeam)*pi/180
    # conversion of fwhm to standard deviation, both in km
    stdPulse <- fwhm2std(fwhmPulse)

    # zenith angle of the beam, in radians
    zenithAngle <- vectorAngle.cartesian(vecAT,locAnt,degrees=F)

    # The 3D covariance matrix of the pulse in the xyz-system calculated above
    covarPulse <- makePulseCovar(stdBeam,stdPulse,zenithAngle,distAT,phArr=phArr)

    # The same covariance matrix in the Cartesian GEOGRAPHIC system, where origin is in the centre of Earth
    # z-axis points to north pole, and x-axis to the zero-meridian at the equator
    covarPulseGeographic <- rotateGeographic(covarPulse,xPulse,yPulse,zPulse)

    # return the covariance in geographic coordinates
    return(list(covar=covarPulseGeographic,zenithAngle=zenithAngle,vecAT=vecAT))

} # pulseShape.gaussian.cartesian




# conversion of full width at half maximum to standard deviation
fwhm2std <- function(fwhm) return(fwhm/2.35482)

# conversion from  standard deviation to full width at half maximum
std2fwhm <- function(std) return(std*2.35482)


tiltedBeamWidth <- function(fwhmBeam,angle){
#
# Gives the fwhm of tilted beam
#
# INPUT:
#	fwhmBeam	beamwidth at Zenith
#	angle		tilting angle from zenith (in radians)
#
# OUTPUT:
#	fwhm of the tilted beam

    # If beamwidth is more than 180 degrees, return the original beamwidth
    if (fwhmBeam > pi)
	{
            return(fwhmBeam)
	}

    theta0 <- sin(fwhmBeam/2)
    sinAngle <- sin(angle)
    angle.left <- asin(sinAngle-theta0)
    if (sinAngle+theta0 < 1)
	{
            angle.right <- asin(sinAngle+theta0)
	}
    else
	{
            angle.right <- pi/2
	}
    tiltedBeamWidth <- angle.right - angle.left

    return(tiltedBeamWidth)

} #tiltedBeamWidth

makePulseCovar <- function(stdBeam,stdPulse,zenithAngle,distAT,degrees=F,phArr){
#
# covariance matrix defining the pulse shape in 3D.
# coorinates:
#        z-axis along the beam
#        x-axis horizontal near the antenna
#        y-axis perpendicular to x and z, such that the coordinate system is right-handed
#
# INPUT:
#    stdBeam     standard deviation of a symmetric zenith-pointed beam
#    stdPulse    standard deviation of the Gaussian envelope defining pulse shape in z-direction
#    zenithAngle zenith angle of the beam in radians
#    distAT      distance from the antenna to the target (scaling factor of beam widths
#     phArr      logical, if TRUE, phased-array antennas with beam widening as function of tilt angle are assumed, otherwise
#                beam width does not depent on pointing direction
#
# OUTPUT:
#    3 x 3 covariance matrix of the 3D pulse shape)
#
# Modification by M. Orispaa (May-4-2012)
#   - beam widening corrected
#   - For this, the fwhm (angle) is needed. In order to not mess with the original arguments, we have to do some extra work
#   - beam eccentrity in not taken into account!! The difference should be negligible.

    if (phArr)
	{
            # fwhm of the beam (in radians) (half of it)
            fwhmBeam <- std2fwhm(stdBeam)
            if (degrees) zenithAngle <- zenithAngle/180*pi

            Beam.widened <- tiltedBeamWidth(fwhmBeam,zenithAngle)
            stdBeam.widened <- fwhm2std(Beam.widened)

            return(diag( c( stdBeam * distAT, stdBeam.widened * distAT, stdPulse)**2))
	}
    else
	{
            return(diag( c( stdBeam * distAT,  stdBeam * distAT, stdPulse)**2))
	}

} # makePulseCovar

rotateGeographic <- function(covarPulse,xPulse,yPulse,zPulse){
#
# Coordinate system rotation to express the covariance matrix covarPulse,
# originally given in xPulse, yPulse, zPulse system, in geocentric coordinates,
# where x==c(1,0,0), y==c(0,1,0), and z==c(0,0,1)
#
# INPUT:
#   covarPulse  3 x 3 covariance matrix
#   xPulse, yPulse, zPulse  orthogonal unit basis vectors of the coordinate system, in which covarPulse is given
#
# OUTPUT:
#
#   3 x 3 covariance matrix in the geocentric coordinate system
#

    rotMat <- matrix(c(xPulse,yPulse,zPulse),byrow=F,nrow=3)

    return( (rotMat%*%covarPulse%*%t(rotMat)) )

} # rotateGeographic

beamIntersect.gaussian.cartesian <- function(locTrans,locRec,locScat,fwhmTrans,fwhmRec,infinity=defaultInfinity(),phArrTrans=TRUE,phArrRec=TRUE){
#
# An apporiximation of the intersection of two Gaussian beams.
#
# INPUT:
#
# locTrans   c(x,y,z) position of the transmitter antenna in Cartesian coordinates. Use the function sphericalToCartesian
#            to transform latitudes and longitudes to the Cartesian system
# locRec     c(x,y,z) position of the receiver antenna
# locScat    c(x,y,z) position of the target
# fwhmTrans  full width at half maximum of the transmitter beam when pointed to zenith, in degrees
# fwhmRec    full width at half maximum of the receiver beam when pointed to zenith, in degrees
# infinity   a large value used as standard deviation along the beam direction
# phArrTrans logical, if TRUE, phased-array antennas with beam widening as function of tilt angle are assumed, otherwise
#            beam width does not depent on pointing direction
# phArrRec   logical, if TRUE, phased-array antennas with beam widening as function of tilt angle are assumed, otherwise
#            beam width does not depent on pointing direction
#
# OUTPUT:
#   a list with elements:
#     covar 3 x 3     covarinace matrix of the approximated beam intersection, in cartesian geocentric coordinates (EFEC)
#     maxZenithAngle  larger of the two beam tilt angles
#     angleIntersect  angle between the two pointing directions
#     eletrans        transmitter beam elevation
#     eleRec          receiver beam elevation
#
#

    # The beams are approximated as the shapes of very long pulses transmitted with the antennas
    # i.e. 3D gaussians with very large variances in beam pointing directions
    covarTransBeam <<- pulseShape.gaussian.cartesian(fwhmTrans,infinity,locTrans,locScat,phArrTrans)
    covarRecBeam   <<- pulseShape.gaussian.cartesian(fwhmRec  ,infinity,locRec  ,locScat,phArrRec)

    # The above calls return the covariances in the same coordinate system, so we can simply multiply the distributions
    covarIntersect <- solve( solve(covarTransBeam$covar) + solve(covarRecBeam$covar) )

    # calculate also the angle between the beam pointing directions, this will be useful later
    angleIntersect <- vectorAngle.cartesian(covarTransBeam$vecAT,covarRecBeam$vecAT,degrees=F)

    return(list(covar=covarIntersect,maxZenithAngle=max(covarTransBeam$zenithAngle,covarRecBeam$zenithAngle),angleIntersect=angleIntersect,eleTrans=pi/2-covarTransBeam$zenithAngle,eleRec=pi/2-covarRecBeam$zenithAngle))

} # beamIntersect.gaussian.cartesian


scatteringVolume.gaussian.cartesian <- function(locTrans,locRec,locScat,fwhmTrans,fwhmRec,fwhmRange,fwhmIonSlab,infinity=defaultInfinity(),phArrTrans=TRUE,phArrRec=TRUE,mineleTrans=0,mineleRec=0){
#
# An apporiximation of the scattering volume
#
# INPUT:
#
# locTrans   c(x,y,z) position of the transmitter antenna in cartesian coordinates. Use the function sphericalToCartesian to transform
#            latitudes and longitudes
# locRec     c(x,y,z) position of the receiver
# locScat    c(x,y,z) position of the target
# fwhmTrans  full width at half maximum of the transmitter beam when pointed to zenith, in degrees
# fwhmRec    full width at half maximum of the receiver beam when pointed to zenith, in degrees
# fwhmRange  full width at half maximum of the gaussian range ambiguity function [km]
# fwhmIonSlab ionospheric slab thickness for self-noise calculations
# infinity   a large value used as standard deviation along the beam direction
# phArrTrans logical, if TRUE, phased-array antennas with beam widening as function of tilt angle are assumed, otherwise
#            beam width does not depent on pointing direction
# phArrRec   logical, if TRUE, phased-array antennas with beam widening as function of tilt angle are assumed, otherwise
#            beam width does not depent on pointing direction
# mineleTrans minimum allowed elevation of transmitter beam, in degrees
# mineleRec minimum allowed elevation of receiver beam, in degrees
#
# OUTPUT:
# 3 x 3 covarinace matrix of the approximated gaussian scattering volume, in cartesian geocentric coordinates
#

    # the beam intersection
    covarIntersect <- beamIntersect.gaussian.cartesian(locTrans,locRec,locScat,fwhmTrans,fwhmRec,infinity,phArrTrans,phArrRec)

    # Normal of the "scattering plane", range is measured along this direction
    zRange <- scatterPlaneNormal.cartesian(locTrans,locRec,locScat)

    # The "range covariance" is rotation symmetric with respect to the z-axis (zRange), we can thus choose arbitrary x- and y-axes,
    # as long as they form a right-handed right-angled coordinate system
    yRange <- normalUnitVector.cartesian(zRange,zRange+runif(3))
    xRange <- normalUnitVector.cartesian(yRange,zRange)

    # A 3D covariance matrix for defining the range resolution, in geographic coor
    # the cosine term originates from the geometry
    # phArr has no effect because zenithAngle==0
    covarRange <- rotateGeographic(
        makePulseCovar(stdBeam=infinity,stdPulse=fwhm2std(fwhmRange)*cos(covarIntersect$angleIntersect/2),
                       zenithAngle=0,distAT=1,phArr=TRUE) , xRange , yRange , zRange
        )

    # vertical direction at the scattering volume
    zVert <- scatterPlaneNormal.cartesian(locScat,locScat,2*locScat)
    yVert <- normalUnitVector.cartesian(zVert,zVert+runif(3))
    xVert <- normalUnitVector.cartesian(yVert,zVert)
    # covariance matrix for defining the thickness of the ionosphere, for self-noise calculations
    covarSlab <- rotateGeographic(makePulseCovar(stdBeam=infinity,stdPulse=fwhm2std(fwhmIonSlab),zenithAngle=0,distAT=1,phArr=TRUE),xVert,yVert,zVert)

    # covariance matrix defining the volume from which self-noise will originate
    if(fwhmIonSlab <= 0){
        covarNoise <<- matrix(0,ncol=3,nrow=3)
    }else{
        covarNoise <<- solve( solve(covarIntersect$covar) + solve(covarSlab) )
    }

    # multiply the beam intersection with the range distribution
    covarScattVol <<- solve( solve(covarIntersect$covar) + solve(covarRange) )

    # if either of the beams was pointed below the horizon, return zero-matrix
    if(covarIntersect$maxZenithAngle > (pi/2)) covarScattVol[,] <<- diag(c(0,0,0)) # NA

    # if elevation of transmitter or receiver beam is below minele we will also return a zero matrix
    if(covarIntersect$eleTrans < (mineleTrans*pi/180)) covarScattVol[,] <<- diag(c(0,0,0))
    if(covarIntersect$eleRec   < (mineleRec*pi/180)  ) covarScattVol[,] <<- diag(c(0,0,0))

    # return the scattering volume
    return(covarScattVol)

} # scatteringVolume.gaussian.cartesian



scatterPlaneNormal.cartesian <- function(locTrans,locRec,locScat){
#
# normal of the plane of constant range (which is an approximation of the true ellipsoidal shape)
#
# INPUT:
# locTrans   c(x,y,z) position of the transmitter in cartesian coordinates. Use the function sphericalToCartesian to transform from
#            latitudes and longitudes
# locRec     c(x,y,z) position of the receiver antenna
# locScat    c(x,y,z) position of the scatterer
#
# OUTPUT
# c(x,y,z)   a unit vector in cartesian geocentric coordinates
#


    # unit vector pointing from the transmitter to the target
    vecTS <- locScat - locTrans
    vecTS <- vecTS / sqrt(sum(vecTS**2))

    # unit vector pointing from the receiver to the target
    vecRS <- locScat - locRec
    vecRS <- vecRS / sqrt(sum(vecRS**2))

    # the normal of the scattering plane
    scatNorm <- vecTS + vecRS
    scatNorm <- scatNorm / sqrt(sum(scatNorm**2))

    return(scatNorm)

} # scatterPlaneNormal.cartesian


rotateHorizontal.covar.cartesian <- function(covar,pos){
#
# rotate the covariance matrix, given in cartesian geocentric coordinates, and with origin of the new coordinate system in pos
# to a coorinate system with x-axis directed to north, y-axis, to west, and z-axis upwards
#
# INPUT:
#  covar  3 x 3 covariance matrix
#  pos    c(x,y,z) coordinates of the target (origin of the coordinate system) in cartesian coordinates
#
# OUTPUT:
# 3 x 3 covariance matrix, with x-axis directed to north, y-axis to west, and z-axis upwards
#

    # z-axis direction in the new system
    z <- pos / sqrt(sum(pos**2))

    # y-axis is perpendicular to both the new z-axis and that of the original system
    if(pos[1]==90){
        y <- c(0,1,0)
    }else{
        y <- normalUnitVector.cartesian(z,c(0,0,1))
    }

    # x-axis is perpendicular to z and y, and the system is right-handed
    x <- normalUnitVector.cartesian(y,z)

    # rotation matrix
    rotMat <- matrix(c(x,y,z),ncol=3,byrow=T)

    # rotate the covariance
    covarRot <- rotMat%*%covar%*%t(rotMat)

    return(covarRot)

} # rotateHorizontal.covar.cartesian


rotateHorizontal.vector.cartesian <- function(vec,pos){
#
# rotate a vector matrix, given in cartesian geocentric coordinates, and with origin of the new coordinate system in pos
# to a coorinate system with x-axis directed to north, y-axis, to west, and z-axis upwards
#
# INPUT:
#  vec    c(x,y,z) the vector
#  pos    c(x,y,z) coordinates of the target (origin of the coordinate system) in cartesian coordinates
#
# OUTPUT:
# c(x',y',z') transformed input vector, with x-axis directed to north, y-axis to west, and z-axis upwards
#

    # z-axis direction in the new system
    z <- pos / sqrt(sum(pos**2))

    # y-axis is perpendicular to both the new z-axis and that of the original system
    if(pos[1]==90){
        y <- c(0,1,0)
    }else{
        y <- normalUnitVector.cartesian(z,c(0,0,1))
    }

    # x-axis is perpendicular to z and y, and the system is right-handed
    x <- normalUnitVector.cartesian(y,z)

    # rotation matrix
    rotMat <- matrix(c(x,y,z),ncol=3,byrow=T)

    # rotate the covariance
    vecRot <- c(matrix(vec,nrow=1)%*%t(rotMat))

    return(vecRot)

} # rotateHorizontal.vector.cartesian




get3Dresolution.cartesian <- function(locTrans,locRec,locScat,fwhmTrans,fwhmRec,fwhmRange,fwhmIonSlab,infinity=defaultInfinity(),phArrTrans=TRUE,phArrRec=TRUE,mineleTrans=0,mineleRec=0){
#
# spatial resolutions for the given measurement geometry
#
# INPUT:
#
# locTrans   c(x,y,z) position of the transmitter antenna in cartesian coordinates. Use sphericalToCartesian to transform from
#            latitudes and longitudes
# locRec     c(x,y,z) position of the receiver antenna in cartesian coordinates
# locScat    c(x,y,z) position of the scatterer (target) in cartesian coordinates
# fwhmTrans  full width at half maximum of the transmitter beam when pointed to zenith, in degrees
# fwhmRec    full width at half maximum of the receiver beam when pointed to zenith, in degrees
# fwhmRange  full width at half maximum of the gaussian range ambiguity function [km]
# fwhmIonSlab ionospheric slab thickness for self-noise calculations
# infinity   a large value used as standard deviation along the beam direction
# phArrTrans  logical, if TRUE, phased-array antennas with beam widening as function of tilt angle are assumed, otherwise
#             beam width does not depent on pointing direction
# phArrRec   logical, if TRUE, phased-array antennas with beam widening as function of tilt angle are assumed, otherwise
#             beam width does not depent on pointing direction
# mineleTrans minimum allowed elevation of transmitter beam, in degrees
# mineleRec minimum allowed elevation of receiver beam, in degrees
#
# OUTPUT:
#
# c(resNorth,resWest,resHeight) full width at half maximum in north-south, east-west, and height directions
#

    return(
        # conversion from standard deviation to full width at half maximum
        std2fwhm(
            # square root of variance to get the standard deviation
            sqrt(
                # diagonal contains variances of marginal distributions
                diag(
                    # coordinate system rotation to x pointing north, y to west and z upwards
                    rotateHorizontal.covar.cartesian(
                        # covariance matrix representing the true scattering volume with given
                        # antenna locations, beam widths, and range resolution
                        scatteringVolume.gaussian.cartesian(
                            locTrans=locTrans,
                            locRec=locRec,
                            locScat=locScat,
                            fwhmTrans=fwhmTrans,
                            fwhmRec=fwhmRec,
                            fwhmRange=fwhmRange,
                            fwhmIonSlab=fwhmIonSlab,
                            infinity=infinity,
                            phArrTrans=phArrTrans,
                            phArrRec=phArrRec,
                            mineleTrans=mineleTrans,
                            mineleRec=mineleRec
                            ),
                        pos = locScat
                        )
                    )
                )
            )
        )

} # get3Dresolution.cartesian



get3Dresolution.spherical <- function(locTrans,locRec,locScat,fwhmTrans,fwhmRec,fwhmRange,fwhmIonSlab,infinity=defaultInfinity(),phArrTrans=TRUE,phArrRec=TRUE,mineleTrans=0,mineleRec=0){
#
# spatial resolutions for the given measurement geometry. This is a wrapper which converts positions to cartesian system
# and calls get3Dresolution.cartesian
#
# INPUT:
#
# locTrans   c(lat,lon) or c(lat,lon,height) latitude, longitude (and height) of the transmitter antenna. Angles in degrees.
# locRec     c(lat,lon) or c(lat,lon,height) latitude, longitude (and height) of the receiver antenna. Angles in degrees.
# locScat    c(lat,lon,height) latitude, longitude and height of the scatterer (target), angles in degress. Height in km
# fwhmTrans  full width at half maximum of the transmitter beam when pointed to zenith, in degrees
# fwhmRec    full width at half maximum of the receiver beam when pointed to zenith, in degrees
# fwhmRange  full width at half maximum of the gaussian range ambiguity function [km]
# infinity   a large value used as standard deviation along the beam direction
# phArrTrans logical, if TRUE, phased-array antennas with beam widening as function of tilt angle are assumed, otherwise
#            beam width does not depent on pointing direction
# phArrRec   logical, if TRUE, phased-array antennas with beam widening as function of tilt angle are assumed, otherwise
#            beam width does not depent on pointing direction
# mineleTrans minimum allowed elevation of transmitter beam, in degrees
# mineleRec minimum allowed elevation of receiver beam, in degrees
#
# OUTPUT:
#
# c(resNorth,resWest,resHeight) full width at half maximum in north-south, east-west, and height directions
#
#
#

    return(get3Dresolution.cartesian(
        locTrans = sphericalToCartesian(locTrans),
        locRec = sphericalToCartesian(locRec),
        locScat = sphericalToCartesian(locScat),
        fwhmTrans = fwhmTrans,
        fwhmRec = fwhmRec,
        fwhmRange = fwhmRange,
        fwhmIonSlab = fwhmIonSlab,
        infinity = infinity,
        phArrTrans = phArrTrans,
        phArrRec = phArrRec,
        mineleTrans = mineleTrans,
        mineleRec = mineleRec
        )
           )

} # get3Dresolution.spherical




##############################################################
########### routines for actual evaluations ##################
##############################################################



bistaticResolutions.planar <- function(refPoint,locTrans,locRec,locxy=F,fwhmTrans=1,fwhmRec=1,fwhmRange=1,fwhmIonSlab=100,x=seq(-300,300,by=10),y=seq(-300,300,by=10),heights=c(100,200,500,1000),infinity=defaultInfinity(),phArrTrans=TRUE,phArrRec=TRUE,mineleTrans=0,mineleRec=0,verbose=FALSE){
#
# resolutions at points c(x,y,height), where x is measured along a circle of latitude, and y along a meridian. The distances are  measured from the point refPoint
#
# INPUT:
#
# refPoint   c(lat,lon) reference point, from which the distance x and y are measured
# locTrans   c(lat,lon) or c(lat,lon,height) latitude, longitude (and height) of the transmitter antenna. Angles in degrees.
#            OR c(x,y) position of the transmitter, IF locxy=T
# locRec     c(lat,lon) or c(lat,lon,height) latitude, longitude (and height) of the receiver antenna. Angles in degrees.
#            OR c(x,y)IF locxy=T
# locxy      logical, if T, the positions of the transmitter and receiver are given in cartesian coordinates
# fwhmTrans  full width at half maximum of the transmitter beam when pointed to zenith, in degrees
# fwhmRec    full width at half maximum of the receiver beam when pointed to zenith, in degrees
# fwhmRange  full width at half maximum of the gaussian range ambiguity function [km]
# fwhmIonSlab ionospheric slab thickness for self-noise calculations
# x          x coordinates
# y          y coordinates
# heights    heights of grid points in degrees
# infinity   a large value used as standard deviation along the beam direction
# phArrTrans  logical, if TRUE, phased-array antennas with beam widening as function of tilt angle are assumed, otherwise
#            beam width does not depent on pointing direction
# phArrRec   logical, if TRUE, phased-array antennas with beam widening as function of tilt angle are assumed, otherwise
#            beam width does not depent on pointing direction
# mineleTrans mimimum allowed elevation for transmitter(s) in degrees [0 90]
# mineleRec  minimum allowed elevation for receiver(s) in degrees [0 90]
#
# OUTPUT:
#  a list with elements:
#    res       a length(x) x length(y) x length(heights) x 3 array with resolutions in north-south, east-west, and vertical directions
#              resolutions are given as fwhm in km
#    x         same as input
#    y         same as input
#    heights   same as input
#    refPoint  same as input
#    fwhmTrans same as input
#    fwhmRec   same as input
#    fwhmRange same as input
#    locTrans  same as input
#    locRec    same as input
#    locxy     same as input
#    xyT       transmitter location in the planar xy-coordinates
#    xyR       receiver location in the planar xy-coordinates
#    lT        latitude and longitude of the transmitter
#    lR        latitude and longitude of the receiver
#    kscat     a length(x) x length(y) x length(heights) x 3 array of scattering wave vectors
#    distTS    a length(x) x length(y) x length(heights) array of distances from transmitter to target [km]
#    distRS    a length(x) x length(y) x length(heights) array of distances from receiver to target [km]
#    gainInt   a length(x) x length(y) x length(heights) array of integrals of the beam intersection over the whole measurement volume
#    gainIntNoise a length(x) x length(y) x length(heights) array of integrals of the beam intersection over the whole ion slab
#    beta      a length(x) x length(y) x length(heights) array of angles between incident and scattering (not scattered!) wave vectors
#    volume    a length(x) x length(y) x length(heights) array of volumes of the measurement volume [m^3]
#    vecTSh    a length(x) x length(y) x length(heights) x 3 array of vectors pointing from target to transmitter in
#              local horizontal coordinate systems
#

    if(locxy){
        latlonTrans <- planarToSpherical.geographic(x=locTrans[1],y=locTrans[2],refPoint=refPoint)
        latlonRec   <- planarToSpherical.geographic(x=locRec[1],y=locRec[2],refPoint=refPoint)
        lT <- c(latlonTrans$lat,latlonTrans$lon)
        lR <- c(latlonRec$lat,latlonRec$lon)
        xyT <- list(x=locTrans[1],y=locTrans[2])
        xyR <- list(x=locRec[1],  y=locRec[2])
    }else{
        latlonTrans <- locTrans
        latlonRec   <- locRec
        xyT <- sphericalToPlanar.geographic(lat=locTrans[1],lon=locTrans[2],zeroLatitude=refPoint[1],zeroMeridian=refPoint[2])
        xyR <- sphericalToPlanar.geographic(lat=locRec[1],lon=locRec[2],zeroLatitude=refPoint[1],zeroMeridian=refPoint[2])
        lT <- locTrans
        lR <- locRec
    }

    nx <- length(x)
    ny <- length(y)
    nh   <- length(heights)
    xyzT <- sphericalToCartesian(lT)
    xyzR <- sphericalToCartesian(lR)

    res          <- array(dim=c(nx,ny,nh,3))
    kscat        <- array(dim=c(nx,ny,nh,3))
    distTS       <- array(dim=c(nx,ny,nh))
    distRS       <- array(dim=c(nx,ny,nh))
    vecTSh       <- array(dim=c(nx,ny,nh,3))
    gainInt      <- array(dim=c(nx,ny,nh))
    gainIntNoise <- array(dim=c(nx,ny,nh))
    volume       <- array(dim=c(nx,ny,nh))
    volumeNoise  <- array(dim=c(nx,ny,nh))
    beta         <- array(dim=c(nx,ny,nh))

    for(i in seq(nx)){
        if(verbose){
            cat(sprintf('\r %5.0f /  %5.0f         ',i,nx))
        }
        for(j in seq(ny)){
            for(k in seq(nh)){
                scatPos        <- planarToSpherical.geographic(x=x[i],y=y[j],refPoint=refPoint)
                xyzScat        <<- sphericalToCartesian(c(scatPos$lat,scatPos$lon,EarthRadius()+heights[k]))
                res[i,j,k,]    <- get3Dresolution.cartesian(
                    locTrans=xyzT,locRec=xyzR,locScat=xyzScat,fwhmTrans=fwhmTrans,fwhmRec=fwhmRec,
                    fwhmRange=fwhmRange,fwhmIonSlab=fwhmIonSlab,infinity=infinity,phArrTrans=phArrTrans,phArrRec=phArrRec,
                    mineleTrans=mineleTrans,mineleRec=mineleRec)

                # the scattering wave vector
                kscatxyz       <- scatterPlaneNormal.cartesian(xyzT,xyzR,xyzScat)
                kscat[i,j,k,]  <- rotateHorizontal.vector.cartesian(kscatxyz,pos=xyzScat)

                # distances from the transmitter and the receiver to the target
                distTS[i,j,k]  <- sqrt(sum(abs(xyzT-xyzScat)**2))
                distRS[i,j,k]  <- sqrt(sum(abs(xyzR-xyzScat)**2))

                vecTSh[i,j,k,] <- rotateHorizontal.vector.cartesian( ((xyzT-xyzScat)/distTS[i,j,k]) , pos=xyzScat )

                # integral over the whole scattering distribution
                alphaT         <- vectorAngle.cartesian(xyzT,(xyzScat-xyzT),degrees=F)

                # Integral of transmission beam
                if (phArrTrans)
                    {
                        fwhmBeam <- fwhmTrans*pi/180
                        fwhmBeamWidened <- tiltedBeamWidth(fwhmBeam,alphaT)
                        intTbeam <- 2 * pi * distTS[i,j,k]**2 * fwhm2std(fwhmTrans*pi/180) * fwhm2std(fwhmBeamWidened)
                    }
                else
                    {
                        intTbeam <- 2 * pi * distTS[i,j,k]**2*fwhm2std(fwhmTrans*pi/180)**2
                    }


                alphaR         <- vectorAngle.cartesian(xyzR,(xyzScat-xyzR),degrees=F)

                # Integral of receiving beam
                if (phArrTrans)
                    {
                        fwhmBeam <- fwhmRec*pi/180
                        fwhmBeamWidened <- tiltedBeamWidth(fwhmBeam,alphaR)
                        intRbeam <- 2 * pi * distRS[i,j,k]**2 * fwhm2std(fwhmRec*pi/180) * fwhm2std(fwhmBeamWidened)
                    }
                else
                    {
                        intRbeam <- 2 * pi * distRS[i,j,k]**2*fwhm2std(fwhmRec*pi/180)**2
                    }


                volume[i,j,k]       <- (2*pi)**(3/2)*sqrt(det(covarScattVol))
                volumeNoise[i,j,k]  <- (2*pi)**(3/2)*sqrt(det(covarNoise))

                gainInt[i,j,k] <- 1e-3*volume[i,j,k]/intTbeam/intRbeam*ifelse(phArrTrans,cos(alphaT),1)

                gainIntNoise[i,j,k] <- 1e-3*volumeNoise[i,j,k]/intTbeam/intRbeam

                # angle between scattering wave vector and incident wave vector
                beta[i,j,k]    <- vectorAngle.cartesian((xyzScat-xyzT),kscatxyz,degrees=T)
            }
        }
    }

    if(verbose){
        cat('\r                                             \r')
    }

    return(list(res=res,x=x,y=y,heights=heights,refPoint=refPoint,fwhmTrans=fwhmTrans,fwhmRec=fwhmRec,fwhmRange=fwhmRange,locTrans=locTrans,
                locRec=locRec,locxy=locxy,xyT=c(xyT$x,xyT$y),xyR=c(xyR$x,xyR$y),lT=lT,lR=lR,kscat=kscat,distTS=distTS,distRS=distRS,
                gainInt=gainInt,gainIntNoise=gainIntNoise,beta=beta,volume=volume,vecTSh=vecTSh))


} # bistaticResolutions.planar


planarToSpherical.geographic <- function(refPoint=c(0,0),x=seq(-1000,1000,by=10),y=seq(-1000,1000,by=10),rEarth=EarthRadius()){
#
# conversion from distance measured along the Earth surface to latitudes and longitudes
#
# INPUT:
#   refPoint c(lat,lon) the reference point, from which the distance are measured, in degrees
#   x                   x values
#   y                   y values   length(x)==length(y) required!!
#   rEarth              Earth radius
#
# OUPUT:
#   a list with entries lat and lon, in degrees
#

    lat <- y/rEarth*180/pi + refPoint[1]

    lon <- x/(rEarth*cos(lat*pi/180))*180/pi + refPoint[2]

    return(list(lat=lat,lon=lon))

} # planarToSpherical.geographic


sphericalToPlanar.geographic <- function(lat,lon,zeroMeridian=19.23,zeroLatitude=69.58,rEarth=EarthRadius()){
#
# conversion from spherical coordinates (latitude,longitude) to system, where
#  x is distance to 0-meridian in km, measured along Earth surface
#  y is distance to equator in km, measured along Earth surface, negative in southern hemisphere
#
# INPUT:
#  lat                       latitude at Earht surface, in degrees
#  lon                       longitude at Earht surface, in degrees
#  zeroMeridian              the meridian we want to use as zero, in degrees
#  zeroLatitude              the latitude we want to use as zero, in degrees
#  rEarth                    Earth radius in km
#
# OUTPUT:
#   list, with components x and y

    x <- cos(lat*pi/180)*rEarth * (lon-zeroMeridian)*pi/180
    y <- (lat-zeroLatitude)*pi/180*rEarth
    names(x) <- NULL
    names(y) <- NULL

    return(list(x=x,y=y))

} # sphericalToPlanar.geographic




###########################################
############# graphics ####################
###########################################


savepdf <- function(devnum,fname){
##
## save the contents of graphics device number devnum into a pdf file name fname
##

    curdev <- dev.cur()
    dev.set(devnum)
    dev.print(device=pdf,file=fname)
    dev.set(curdev)

} ## savepdf

plotOnMap <- function(resData,f=function(x){return(x)},x,y,latlon=F,refPoint=c(69.58,19.23),xlim=range(x),ylim=range(y),zlim=range(f(resData),na.rm=T,finite=T),col=rev(gray(seq(1000)/1000)),xlab='[km]',ylab='[km]',cbarlab='Resolution [km]',main='',txloc=NA,rxloc=NA,levels=NA,gdev=x11,printInfo=F,printString=NULL,totalInfo=F,useRaster=useRaster){

    if(!latlon){
        maplims <- planarToSpherical.geographic(x=xlim,y=ylim,refPoint=refPoint)
        mapxlim <- maplims$x
        mapylim <- maplims$y
    }else{
        mapxlim <- xlim
        mapylim <- ylim
    }

    m <- map(database='worldHires',xlim=mapxlim,ylim=mapylim,plot=F,resolution=0,boundary=T)

    if(!latlon){
        mxy <- sphericalToPlanar.geographic(lat=m$y,lon=m$x,zeroMeridian=refPoint[2],zeroLatitude=refPoint[1])
    }else{
        mxy <- list(x=m$x,y=m$y)
    }

    speed <- sum(1/resData,na.rm=T)/length(x)/length(y)

    resd <- f(resData)
    resd[resd>zlim[2]] <- zlim[2]
    resd[is.nan(resd)] <- zlim[2]
    resd[is.na(resd)] <- zlim[2]

    if (printInfo)
        {
            gdev(width=6,height=10)
            layout(matrix(c(1,2,3),ncol=1),heights=c(.6,.05,.35))

        }
    else
        {
            gdev(width=6,height=6/.88)

            layout(matrix(c(1,2),ncol=1),heights=c(.88,.12))
        }
    par(mar=c(3,3,2,2),mgp=c(1.5,.5,0))

    if (totalInfo) main <-paste(main,', Speed: ',round(speed,4),' px/sec',sep='')

    # Draw filled contours
    image(x,y,resd,col=col,xlim=xlim,ylim=ylim,zlim=zlim,xlab=xlab,ylab=ylab,main=main,useRaster=useRaster)

    # Draw contour lines
    if(!all(is.na(levels))) contour(x,y,resData,add=T,levels=levels,lwd=2,col='red',labcex=1)

    par(mar=c(3,3,0,2))

    # Draw map
    lines(mxy$x,mxy$y,lwd=2)

    if(!is.na(txloc[1])){
        if(is.list(txloc)){
            for(k in seq(length(txloc))) points(txloc[[k]][1],txloc[[k]][2],col='blue',pch=24)
        }else{
            points(txloc[1],txloc[2],col='blue',pch=24)
        }
    }
    if(!is.na(rxloc[1])){
        if(is.list(rxloc)){
            for(k in seq(length(rxloc))) points(rxloc[[k]][1],rxloc[[k]][2],col='red',pch=19)
        }else{
            points(rxloc[1],rxloc[2],col='red',pch=19)
        }
    }

    image(seq(zlim[1],zlim[2],length.out=1000),c(1),matrix(seq(zlim[1],zlim[2],length.out=1000),ncol=1),col=col,yaxt='n',ylab='',xlab=cbarlab ,bty='L',useRaster=useRaster)

    if (printInfo)
        {
            plot(-1,xlim=c(0,1),ylim=c(0,1),ax=F,xlab='',ylab='')


            text(0.3,0.4,'TX sites:\nRX sites:\nfwhmTrans:\nfwhmRec:\nfwhmRange:\nResolution:\nPt:\nNe:\nTnoise:\nfradar:\ntau0:\ntargetNoiseLevel:\ndutyCycle:',pos=2,cex=1.5)


            text(0.31,0.4,printString,pos=4,cex=1.5)

        }


} # plotOnMap


plotkvectors <- function(resData,ix=ceiling(length(resData$x)/2),iy=ceiling(length(resData$y)/2),ih=1,xlim=range(resData$x),ylim=range(resData$y),zlim=c(0,max(resData$h)),theta=0,phi=15,r=sqrt(3),d=1,vlen=100){
#
# plot scattering wave vectors
#
# INPUT:
#   resData            an output list from the function multistaticIntegrationTimes or from multistaticNoiseLevels
#   ix, iy, ih         the point in the xyh-grid, whose k-vectors are plotted, 0,0,0 is at the lowest height in the southwest corner
#   xlim, ylim, zlim   axis limits in the plot
#
# OUTPUT:
#   none
#
#

    if(length(xlim)<2) xlim <- xlim + c(-10,10)
    if(length(ylim)<2) ylim <- ylim + c(-10,10)
    if(length(zlim)<2) zlim <- zlim + c(-10,10)

    if(xlim[1]==xlim[2]) xlim <- xlim + c(-10,10)
    if(ylim[1]==ylim[2]) ylim <- ylim + c(-10,10)
    if(zlim[1]==zlim[2]) zlim <- zlim + c(-10,10)

    # number of different k-vectors
    nk <- length(resData$noiseLevel.site)

    # a dummy call to create the viewing transformation matrix
    x <- seq(xlim[1],xlim[2])
    y <- seq(xlim[1],xlim[2])
    z <- matrix(NA,nrow=length(x),ncol=length(y))

    persp(x=x,y=y,z=z,xlim=xlim,ylim=ylim,zlim=zlim,xlab='East',ylab='North',zlab='Up',box=T,theta=theta,phi=phi,r=r,d=d,scale=T,asp=1) -> perspmat

    for(k in seq(nk)){
        x <- c(0,resData$noiseLevel.site[[k]]$kscat[ix,iy,ih,2]*vlen)+resData$x[ix]
        y <- c(0,-resData$noiseLevel.site[[k]]$kscat[ix,iy,ih,1]*vlen)+resData$y[iy]
        z <- c(0,-resData$noiseLevel.site[[k]]$kscat[ix,iy,ih,3]*vlen)+resData$h[ih]

        xy <- trans3d(x,y,z,perspmat)
        arrows(x0=xy$x[1],y0=xy$y[1],x1=xy$x[2],y1=xy$y[2],col=k,lwd=2)

        xy <- trans3d(x=resData$noiseLevel.site[[k]]$locRec[1],y=resData$noiseLevel.site[[k]]$locRec[[2]],0,perspmat)
        points(xy$x,xy$y,col=k,pch=20)

        xy <- trans3d(x=resData$noiseLevel.site[[k]]$locTrans[1],y=resData$noiseLevel.site[[k]]$locTrans[[2]],0,perspmat)
        points(xy$x,xy$y,col='black',pch='x')

    }

} # plotkvectors



#####################################################################
##### error analysis for multi-static ion velocity measurements #####
#####################################################################

vectorVelocityVariance <- function(kscats=list(),angle=list(),vars){
#
# given a set of velocity component directions  and optionally variances, estimate the variances of orthogonal
# vector velocity components
#
# INPUT:
#   kscats   a list of nx x ny x nh x 3 arrays, which contain the unit scattering wave vectors of each transmitter-receiver pair
#   angle    list of angles between   TX beam and scattering wave vector
#   vars     an optional nx x ny x nh array of variances
#
# OUPUT:
#   an nx x ny x nh x 3 array of variances of the vector velocity components
#

    # number of TX-RX-pairs
    if(is.list(kscats)){
        n <- length(kscats)
    }else{
        kscats <- list(kscats)
        n <- 1
    }
    kdims <- dim(kscats[[1]])
    nx <- kdims[1]
    ny <- kdims[2]
    nh <- kdims[3]


    # if the variances are given, then use them
    if(missing(vars)){
        vars <- list()
        for(k in seq(n)){
            vars[[k]] <- array(1,dim=c(nx,ny,nh))
        }
        nv <- n
    }else{
        if(is.list(vars)){
            nv <- length(vars)
        }else{
            vars <- list(vars)
            nv <- 1
        }
    }
  if(n!=nv) stop('Lengths of k and variance lists are different.')

    # an array for the vector element variances
    kvars <- kscats[[1]]*0

    A <- matrix(0,ncol=3,nrow=n)
    invsigma <- matrix(0,nrow=n,ncol=n)

    for(i in seq(nx)){
        for(j in seq(ny)){
            for(k in seq(nh)){
                A[,] <- 0
                invsigma[,] <- 0
                for(l in seq(n)){
                    A[l,] <- kscats[[l]][i,j,k,]
                    invsigma[l,l] <- 1/vars[[l]][i,j,k]
                }
                kvars[i,j,k,] <- tryCatch(diag(solve(t(A)%*%invsigma%*%A)),error=function(e){return(rep(NA,3))})
            }
        }
    }

    return(kvars)

} # vectorVelocityVariance


spatialIntegrationCoeffiecients <- function(resEW=10,resNS=10,resH=10,resolutions){
#
# Calculate the number of range-gate ACF samples that fit inside a resEW x resNS x resH spatial resolution cell
# INPUT:
#   resEW        resolution in east-west direction [km]
#   resNS        resolution in north-south direction [km]
#   resH         resolution in height [km]
#   resolutions  an output list from bistaticResolutions.planar or bistaticResolutions.spherical
#
# OUTPUT:
#   an nx x ny x nh array of integer factors telling how many range-gate estimates fit inside a voxel
#                   if  none fits, 0 is returned
#
# I. Virtanen 2011
#


    resdims <- dim(resolutions$res)
    nx <- resdims[1]
    ny <- resdims[2]
    nh <- resdims[3]
    intcoefs <- array(dim=c(nx,ny,nh))

    tres <- c(resNS,resEW,resH)

    for(i in seq(nx)){
        for(j in seq(ny)){
            for(k in seq(nh)){
                minint <- min( ( tres - resolutions$res[i,j,k,] ) / (resolutions$fwhmRange*abs(resolutions$vecTSh[i,j,k,])) )
                intcoefs[i,j,k] <- floor(ifelse( minint<0 , 0 , minint + 1 ))
            }
        }
    }

    return(intcoefs)


} # spatialIntegrationCoeffients




absoluteNoiseLevel <- function(resolutions,intCoefs,Pt=1e6,Ne=1e11,Tnoise=300,tau0,fradar=233e6,RXduty=1,verbose=F){
#
# absolute value of ACF noise level
#
#
# INPUT:
#   resolutions             an output list from bistaticResolutions.planar or bistaticResolutions.spherical
#   intCoefs    (optional)  an array of spatial integration coefficients (output array of spatialIntegrationCoefficients),
#                           unit values are used if missing
#   Pt                      transmitter power
#   Ne                      electron densities at each height
#   Tnoise                  receiver noise temperature
#   tau0                    acf time-scale (see vallinkoski 1988) [microseconds]
#   fradar                  radar carrier frequency
#   RXduty                  "Receiver duty cycle", Typically 1, but [0,1] for monostatic systems
#
# OUTPUT:
#   a list with elements:
#     noiseLevel   an nx x ny x nh array of noise level values
#     dtau         sample interval calculated from range resolution
#     lambda       radar wave length
#     ntau         number of lags between 0 and tau0
#     Pr           an array of scattering signal powers received from the measurement volumes
#     Pn           background noise power
#     Prn          an array of total received signal powers, including self-noise from outside the measurement volume
#     SNR          signal-to-noise ratio calculated as ratio Pr/Pn, i.e. neglecting the self-noise
#
#
#

    # Boltzmann constant
    kb <- 1.3806503e-23

    # speed of light
    c <- 299792458

    # single electron cross-section
    xi0 <- 1e-28

    # radar wave length
    lambda <- c/fradar

    # sampling interval from range resolution
    dtau <-resolutions$fwhmRange*1000/c*2

    # number of lags between 0 and tau0
    ntau <- round(tau0*1e-6/dtau/cos(resolutions$beta*pi/180))

    # if integration coefficients are missing, do not integrate
    if(missing(intCoefs)) intCoefs <- 1

    # noise power
    Pn <- kb*Tnoise/dtau

    # number of heights from dim(gainInt)
    nhei <- dim(resolutions$gainInt)[3]

    # received scattering signal power from the nominal range-gate
    Pr <- xi0/2 * 1/2*(1+cos(2*resolutions$beta*pi/180)**2) * Pt *  resolutions$gainInt * lambda**2/(4*pi)
    for( hh in seq(nhei)){
        Pr[,,hh] <- Pr[,,hh]*Ne[hh]
    }

    # total received scattering signal power for self-noise estimation
    Prn <- xi0/2 * 1/2*(1+cos(2*resolutions$beta*pi/180)**2) * Pt *  resolutions$gainIntNoise * lambda**2/(4*pi)
    for( hh in seq(nhei)){
        Prn[,,hh] <- Prn[,,hh]*Ne[hh]
    }

    # ACF noise level
    noiseLevel <- ( Pn  + Prn ) / Pr / sqrt(ntau) / sqrt(intCoefs) / sqrt(RXduty)

    if(verbose){
        cat(sprintf('Wave length [m] %19.9f      \n',lambda))
        cat(sprintf('Sample step [us]%19.9f      \n',dtau*1e6))
#        cat(sprintf('Number of lags  %19d        \n',ntau))
        cat('Mean values: \n')
        cat(sprintf('Signal power [W]%15.7fe-18  \n',mean(Pr*1e18)))
        cat(sprintf('Thermal noise power [W] %15.7fe-18  \n',mean(Pn*1e18)))
        cat(sprintf('Self-noise power [W] %15.7fe-18  \n',mean(Prn*1e18)))
        cat(sprintf('Signal-to-noise %19.9f      \n',mean(Pr/Pn)))
        cat(sprintf('ACF noise level %19.9f      \n',mean(noiseLevel)))
    }

    invisible( list(noiseLevel=noiseLevel,dtau=dtau,lambda=lambda,ntau=ntau,Pr=Pr,Pn=Pn,Prn=Prn,SNR=Pr/Pn) )

} # absoluteNoiseLevel


integrationTimes <- function(nLev, targetNoiseLevel=.01,dutyCycle=.25){

    return(nLev$noiseLevel**2/targetNoiseLevel**2 * nLev$dtau / dutyCycle)

} # integrationTimes



noiseLevels.planar <- function(refPoint,locTrans,locRec,locxy=F,fwhmTrans=1,fwhmRec=1,fwhmRange=1,fwhmIonSlab=100,x=seq(-300,300,by=10),y=seq(-300,300,by=10),heights=c(100,200,500,1000),infinity=defaultInfinity(),resNS,resEW,resH,Pt=1e6,Ne=1e11,Tnoise=300,fradar=233e6,tau0=100,RXduty=1,phArrTrans=TRUE,phArrRec=TRUE,verbose=FALSE){
#
# ACF noise levels at the grid points, with optional spatial integration, calculated by comparing to a reference measurement
#
# INPUT:
#   refPoint   c(lat,lon) reference point, from which the distance x and y are measured
#   locTrans   c(lat,lon) or c(lat,lon,height) latitude, longitude (and height) of the transmitter antenna. Angles in degrees.
#              OR c(x,y) position of the transmitter, IF locxy=T
#   locRec     c(lat,lon) or c(lat,lon,height) latitude, longitude (and height) of the receiver antenna. Angles in degrees.
#              OR c(x,y)IF locxy=T
#   locxy      logical, if T, the positions of the transmitter and receiver are given in cartesian coordinates
#   fwhmTrans  full width at half maximum of the transmitter beam when pointed to zenith, in degrees
#   fwhmRec    full width at half maximum of the receiver beam when pointed to zenith, in degrees
#   fwhmRange  full width at half maximum of the gaussian range ambiguity function [km]
#   fwhmIonSlab ionospheric slab thickness for self-noise calculations
#   latitudes  latitudes of grid points in degrees
#   longitudes longitudes of grid points in degrees
#   heights    heights of grid points in degrees
#   infinity   a large value used as standard deviation along the beam direction
#   resNS (optional) north-south resolution of the integrated resolution cells [km]
#   resEW (optional) east-west resolution of the integrated resolution cells [km]
#   resH  (optional) height resolution of the integrated resolution cells [km]
#   Pt         transmitter power [W]
#   Ne         electron density [m^-3]
#   Tnoise     system noise temperature [K]
#   fradar     radar carrier frequency [Hz]
#   tau0       ACF time-scale [us]
#   RXduty     "receiver duty cycle" [0,1]
#   phArrTrans logical, if TRUE, phased-array antennas with beam widening as function of tilt angle are assumed, otherwise
#              beam width does not depent on pointing direction
#   phArrRec   logical, if TRUE, phased-array antennas with beam widening as function of tilt angle are assumed, otherwise
#              beam width does not depent on pointing direction
#
# OUTPUT:
#   see absoluteNoiseLevel

    res <- bistaticResolutions.planar(refPoint=refPoint,locTrans=locTrans,locRec=locRec,locxy=locxy,fwhmTrans=fwhmTrans,fwhmRec=fwhmRec,fwhmRange=fwhmRange,fwhmIonSlab=fwhmIonSlab,x=x,y=y,heights=heights,infinity=infinity,phArrTrans=phArrTrans,phArrRec=phArrRec)

    if(any(missing(resNS),missing(resEW),missing(resH))){
        cat('missing spatial resolution, assuming analysis without spatial integration.\n')
        intCoefs <-1
    }else{
        intCoefs <- spatialIntegrationCoeffiecients(resEW=resEW,resNS=resNS,resH=resH,resolutions=res)
    }


    return(absoluteNoiseLevel(resolutions=res,intCoefs=intCoefs,Pt=Pt,Ne=Ne,Tnoise=Tnoise,fradar=fradar,tau0=tau0,RXduty=RXduty,verbose=verbose))


} # noiseLevels.planar



multistaticNoiseLevels <- function(refPoint,locTrans,locRec,locxy=F,fwhmTrans=c(1),fwhmRec=c(1),fwhmRange=c(1),x=seq(-300,300,by=10),y=seq(-300,300,by=10),heights=c(100,300,500),infinity=defaultInfinity(),resNS,resEW,resH,resR,Pt=c(1e6),Ne=1e11,fwhmIonSlab=100,Tnoise=c(300),fradar=233e6,tau0=100,phArrTrans=c(TRUE),phArrRec=c(TRUE),RXduty=c(1),verbose=FALSE,mineleTrans=0,mineleRec=0){
#
# ACF noise levels and velocity ACF noise levels at the grid points, with optional spatial integration,
# calculated by comparing to a reference measurement
#
# INPUT:
#   refPoint   c(lat,lon) reference point, from which the distance x and y are measured
#   locTrans   list(c(lat,lon)) or list(c(lat,lon,height)) latitudes, longitudes (and heights) of the transmitter antennas.
#              Angles in degrees.
#              OR list(c(x,y)) position of the transmitter, IF locxy=T
#   locRec     list(c(lat,lon)) or list(c(lat,lon,height)) latitudes, longitudes (and heights) of the receiver antennas.
#              Angles in degrees.
#              OR list(c(x,y)) IF locxy=T
#   locxy      logical, if T, the positions of the transmitter and receiver are given in cartesian coordinates
#   fwhmTrans  c(fwhmTrans1,fwhmTrans2,...) full width at half maximum of the transmitter beams when pointed to zenith, in degrees
#   fwhmRec    c(fwhmRec1,fwhmRec2,...) full width at half maximum of the receiver beams when pointed to zenith, in degrees
#   fwhmRange  c(fwhmRange1,fwhmRange2,...) full width at half maximum of the gaussian range ambiguity function [km]
#   latitudes  latitudes of grid points in degrees
#   longitudes longitudes of grid points in degrees
#   heights    heights of grid points in degrees
#   infinity   a large value used as standard deviation along the beam direction
#   resNS (optional) north-south resolution of the integrated resolution cells [km]
#   resEW (optional) east-west resolution of the integrated resolution cells [km]
#   resH  (optional) height resolution of the integrated resolution cells [km]
#   resR  (optional) integrated range-resolution [km], without checking the horizontal and vertical dimensions of the final resolution cells.
#                    used if none of the above three is given. Matched with fwhmRange if missing.
#   Pt         c(Pt1,Pt2,...) transmitter power [W]
#   Ne         electron density [m^-3]
#   fwhmIonSlab ionospheric slab thickness [km]
#   Tnoise     system noise temperatures of the receivers [K]
#   fradar     radar carrier frequency [Hz]
#   tau0       ACF time-scale [us]
#   phArrTrans logical, if TRUE, phased-array antennas with beam widening as function of tilt angle are assumed, otherwise
#               beam width does not depent on pointing direction
#   phArrRec    logical, if TRUE, phased-array antennas with beam widening as function of tilt angle are assumed, otherwise
#               beam width does not depent on pointing direction
#   RXduty      "Receiver duty cycles" [0,1]
#   verbose     optional printing of average snr etc.
#   mineleTrans mimimum allowed elevation for transmitter(s) in degrees [0 90]
#   mineleRec  minimum allowed elevation for receiver(s) in degrees [0 90]
#
#  OUTPUT:
#   a list with elements:
#     x                     same as input
#     y                     same as input
#     h                     the heights input
#     refPoint              same as input
#     xyT                   list of transmitter locations in xy-coordinates
#     xyT                   list of receiver locations in xy-coordinates
#     noiseLevel.isotropic  combined noise level for the isotropic parameters
#     noiseLevel.velocity   combined noise level for velocity
#     noiseLevel.site       outputs from bistaticResolutions.planar and absoluteNoise level for each T / R combination
#     velvars               velocity variances, see vectorVelocityVariance
#     kscats                scattering wave vectors
#     angles                angles between incident and scattering wave vectors
#
#




    if(!is.list(locTrans)) locTrans <- list(locTrans)
    if(!is.list(locRec)) stop('locRec should be a list of reeiver site locations.')

    # ACF noise levels for all transmitter-receiver pairs
    nTrans <- length(locTrans)
    nRec   <- length(locRec)
    nComb  <- nTrans*nRec

    Tnoise     <- rep(Tnoise,length.out=nRec)
    phArrTrans <- rep(phArrTrans,length.out=nTrans)
    phArrRec   <- rep(phArrTrans,length.out=nRec)
    fwhmTrans  <- rep(fwhmTrans,length.out=nTrans)
    fwhmRec    <- rep(fwhmRec,length.out=nRec)
    fwhmRange  <- rep(fwhmRange,length.out=nRec)
    Pt         <- rep(Pt,length.out=nTrans)
    mineleTrans <- rep(mineleTrans,length.out=nTrans)
    mineleRec <- rep(mineleRec,length.out=nRec)
    RXduty <- rep(RXduty,length.out=nRec)

    noiseLevels <- vector(mode='list',length=nComb)
    n <- 0
    for(i in seq(nTrans)){
        for(j in seq(nRec)){
            if(verbose){
                cat(sprintf('Transmitter %d / %d , receiver %d / %d\n',i,nTrans,j,nRec))
            }
            n <- (i-1)*nRec+j

            # gain integrals etc.
            noiseLevels[[n]] <- bistaticResolutions.planar(refPoint=refPoint,locTrans=locTrans[[i]],locRec=locRec[[j]],locxy=locxy,fwhmTrans=fwhmTrans[i],fwhmRec=fwhmRec[j],fwhmRange=fwhmRange[j],fwhmIonSlab=fwhmIonSlab,x=x,y=y,heights=heights,infinity=infinity,phArrTrans=phArrTrans[i],phArrRec=phArrRec[j],mineleTrans=mineleTrans[i],mineleRec=mineleRec[j])

            if(any(missing(resNS),missing(resEW),missing(resH))){
                if(missing(resR)){
                    intCoefs <- 1
                }else{
                    intCoefs <- resR / fwhmRange[j]
                }
            }else{
                intCoefs <- spatialIntegrationCoeffiecients(resEW=resEW,resNS=resNS,resH=resH,resolutions=noiseLevels[[n]])
            }

            # ACF noise level
            noiseLevels[[n]]$noiseLevel <- absoluteNoiseLevel(resolutions=noiseLevels[[n]],intCoefs=intCoefs,Pt=Pt[i],Ne=Ne,Tnoise=Tnoise[j],fradar=fradar,tau0=tau0,RXduty=RXduty[j],verbose=verbose)
        }
    }


    # combined noise level for uniform parameters
    noiseLevel.isotropic <- noiseLevels[[1]]$noiseLevel
    noiseLevel.isotropic$noiseLevel <- 1/noiseLevels[[1]]$noiseLevel$noiseLevel**2
    for(k in seq(2,nComb)){
        noiseLevel.isotropic$noiseLevel <- noiseLevel.isotropic$noiseLevel + 1/noiseLevels[[k]]$noiseLevel$noiseLevel**2
    }

    noiseLevel.isotropic$noiseLevel <- 1/sqrt(noiseLevel.isotropic$noiseLevel)


    # vector velocity noise level
    kscats <- vector(mode='list',length=nComb)
    angles <- vector(mode='list',length=nComb)
    vars   <- vector(mode='list',length=nComb)
    for(k in seq(nComb)){
        kscats[[k]] <- noiseLevels[[k]]$kscat
        angles[[k]] <- noiseLevels[[k]]$beta
        vars[[k]]   <- noiseLevels[[k]]$noiseLevel$noiseLevel**2
    }
    velvars <- vectorVelocityVariance(kscats,angles,vars)

    noiseLevel.velocity <- noiseLevels[[1]]$noiseLevel
    noiseLevel.velocity$noiseLevel<- sqrt(apply(velvars,MARGIN=c(1,2,3),FUN=sum))

    xyT <- xyR <- vector(mode='list',length=nComb)
    for(k in seq(nComb)){
        xyT[[k]] <- noiseLevels[[k]]$xyT
        xyR[[k]] <- noiseLevels[[k]]$xyR
    }
    xyT <- unique(xyT)
    xyR <- unique(xyR)




    return(list(x=x,y=y,h=heights,refPoint=refPoint,xyT=xyT,xyR=xyR,noiseLevel.isotropic=noiseLevel.isotropic,noiseLevel.velocity=noiseLevel.velocity,noiseLevel.site=noiseLevels,velvars=velvars,kscats=kscats,angles=angles))


} # multistaticNoiseLevels







multistaticIntegrationTimes <- function(refPoint=ISgeometry:::SKI,locTrans=ISgeometry:::SKI,locRec=list(ISgeometry:::SKI,ISgeometry:::KAR,ISgeometry:::KAI),locxy=FALSE,fwhmTrans=c(2.2),fwhmRec=c(1.3,1.8,1.8),fwhmRange=c(1),x=seq(-300,300,by=25),y=seq(-300,300,by=25),heights=c(300),infinity=defaultInfinity(),resNS,resEW,resH,resR,Pt=c(3.5e6),Ne=1e11,fwhmIonSlab=100,Tnoise=c(300),fradar=233e6,tau0=100,phArrTrans=c(TRUE),phArrRec=c(TRUE),targetNoiseLevel=.01,dutyCycle=c(.25),RXduty=1,gdev=NULL,zlim=c(0,5),zlimv=c(2,7),iso.levels=NULL,vel.levels=NULL,printInfo=F,plotResolution=F,NScut=0,heightProf=c(0,0),horCut=TRUE,verbose=FALSE,mineleTrans=30,mineleRec=30,useRaster=FALSE){
#
# integration times needed to reach target ACF noise level and velocity ACF noise level
# calculated by comparing to a reference measurement
#
# INPUT:
#   refPoint   c(lat,lon) reference point, from which the distance x and y are measured
#   locTrans   list(c(lat,lon)) or list(c(lat,lon,height)) latitudes, longitudes (and heights) of the transmitter antennas.
#              Angles in degrees.
#              OR list(c(x,y)) position of the transmitter, IF locxy=T
#   locRec     list(c(lat,lon)) or list(c(lat,lon,height)) latitudes, longitudes (and heights) of the receiver antennas.
#              Angles in degrees.
#              OR list(c(x,y)) IF locxy=T
#   locxy      logical, if T, the positions of the transmitter and receiver are given in cartesian coordinates
#   fwhmTrans  c(fwhmTrans1,fwhmTrans2,...) full width at half maximum of the transmitter beams when pointed to zenith, in degrees
#   fwhmRec    c(fwhmRec1,fwhmRec2,...) full width at half maximum of the receiver beams when pointed to zenith, in degrees
#   fwhmRange  c(fwhmRange1,fwhmRange2,...) full width at half maximum of the gaussian range ambiguity function [km]
#   x          grid point distance from refPoint in E-W direction [km]
#   y          grid point distance from refPoint in N-S direction [km]
#   heights    heights of grid points [km]
#   infinity   a large value used as standard deviation along the beam direction
#   resNS (optional) north-south resolution of the integrated resolution cells [km]
#   resEW (optional) east-west resolution of the integrated resolution cells [km]
#   resH  (optional) height resolution of the integrated resolution cells [km]
#   resR  (optional) integrated range-resolution [km], without checking the horizontal and vertical dimensions of the final resolution cells.
#                    used if some of the above three is not given. Matched with fwhmRange if missing.
#   Pt         c(Pt1,Pt2,...) transmitter power [W]
#   Ne         electron density [m^-3]. Either single value or a length(heights) vector. The last value will be repeated if length(Ne)<length(heights)
#   fwhmIonSlab ionospheric slab thickness [km]
#   Tnoise     system noise temperatures of the receivers [K]
#   fradar     radar carrier frequency [Hz]
#   tau0       ACF time-scale [us]
#   phArrTrans logical, if TRUE, phased-array antennas with beam widening as function of tilt angle are assumed, otherwise
#               beam width does not depent on pointing direction
#   phArrRec   logical, if TRUE, phased-array antennas with beam widening as function of tilt angle are assumed, otherwise
#               beam width does not depent on pointing direction
#   targetNoiseLevel
#              target noise level value, same target value used for both isotropic and velocity integration times
#   dutyCycle  transmitter duty cycle
#   RXduty     "Receiver duty cycles", 1 for remotes, [0,1] for monostatic sites
#   gdev       graphics device used for plotting
#   zlim       zlim for the plots of isotropic parameter integration times
#   zlimv      zlim for the plots of velocity integration times
#   iso.levels contour line levels for isometric parameter plot
#   vel.levels contour line levels for velocity plot
#   printInfo  plots used parameters below the map plot
#   plotResolution
#              plots resolution map in addition to the isotropic and velocity integration time maps.
#   NScut      a vector of indices of north-south directed vertical slices to plot. nothing will be plotted if NScut <= 0
#   heightProf indices (x,y) of a height profile of the speed estimates to plot. Plot is not produced if length(heightProf)<2,any(heightProf<=0), or is.null(heightProf).
#   horCut     if FALSE, the plots of horisontal cuts will not be produced. Useful e.g. when producing only a single vertical cut with x=c(0).
#   verbose    optional printing of average SNR / power etc.
#   mineleTrans mimimum allowed elevation for transmitter(s) in degrees [0 90]
#   mineleRec  minimum allowed elevation for receiver(s) in degrees [0 90]
#   useRaster  rasterize images?
#
#  OUTPUT:
#    nLevel (invisible) output from multistaticNoiseLevels padded with elements intTime.isotropic and intTime.velocity
#                       that are outputs from the function integrationTimes
#
#


    # If graphic device is not given choose one based on Sys.info
    if (is.null(gdev))
        {
            systemName <- tolower(Sys.info()['sysname'])
            if (systemName == "darwin")
                {
                    gdev <- quartz
                    cat('Graphic device: quartz\n')
                }
            else if (systemName == "linux" || systemName == "unix")
                {
                    gdev <- x11
                    cat('Graphic device: x11\n')
                }
            else if (systemName == "windows")
                {
                    gdev <- windows
                    cat('Graphic device: windows\n')
                }
            else
                {
                    stop('Please give preferred graphics device using argument keyword gdev.\n')
                }
        }

    # check that length(Ne) and length(heights) match, repeat lats element of Ne as necessary
    if ((lNe<-length(Ne))<(lhei<-length(heights))){
        Ne <- c(Ne,rep(Ne[lNe],lhei-lNe))
    }
    data.str <- NULL


    if(any(c(missing(resNS),missing(resEW),missing(resH)))){
        if(missing(resR)){
            nLevel <- multistaticNoiseLevels(refPoint=refPoint,locTrans=locTrans,locRec=locRec,
                                             locxy=locxy,fwhmTrans=fwhmTrans,fwhmRec=fwhmRec,
                                             fwhmRange=fwhmRange,phArrTrans=phArrTrans,phArrRec=phArrRec,
                                             x=x,y=y,heights=heights,Tnoise=Tnoise,infinity=infinity,
                                             Pt=Pt,Ne=Ne,fwhmIonSlab=fwhmIonSlab,
                                             fradar=fradar,tau0=tau0,RXduty=RXduty,verbose=verbose,
                                             mineleTrans=mineleTrans,mineleRec=mineleRec)
        }else{
            nLevel <- multistaticNoiseLevels(refPoint=refPoint,locTrans=locTrans,locRec=locRec,
                                             locxy=locxy,fwhmTrans=fwhmTrans,fwhmRec=fwhmRec,
                                             fwhmRange=fwhmRange,phArrTrans=phArrTrans,phArrRec=phArrRec,
                                             x=x,y=y,heights=heights,Tnoise=Tnoise,infinity=infinity,
                                             resR=resR,Pt=Pt,Ne=Ne,fwhmIonSlab=fwhmIonSlab,
                                             fradar=fradar,tau0=tau0,RXduty=RXduty,verbose=verbose,
                                             mineleTrans=mineleTrans,mineleRec=mineleRec)
        }
    }else{
        nLevel <- multistaticNoiseLevels(refPoint=refPoint,locTrans=locTrans,locRec=locRec,locxy=locxy,
                                         fwhmTrans=fwhmTrans,fwhmRec=fwhmRec,fwhmRange=fwhmRange,phArrTrans=phArrTrans,phArrRec=phArrRec,
                                         x=x,y=x,heights=heights,Tnoise=Tnoise,infinity=infinity,
                                         resNS=resNS,resEW=resEW,resH=resH,Pt=Pt,
                                         Ne=Ne,fwhmIonSlab=fwhmIonSlab,fradar=fradar,tau0=tau0,RXduty=RXduty,verbose=verbose,
                                         mineleTrans=mineleTrans,mineleRec=mineleRec)
    }

    nLevel$intTime.isotropic <- integrationTimes(nLevel$noiseLevel.isotropic,targetNoiseLevel=targetNoiseLevel,dutyCycle=dutyCycle)
    nLevel$intTime.velocity  <- integrationTimes(nLevel$noiseLevel.velocity, targetNoiseLevel=targetNoiseLevel,dutyCycle=dutyCycle)

    for(k in seq(length(heights))){
        # parse the information string
        if(any(c(missing(resNS),missing(resEW),missing(resH)))){
            if(missing(resR)){
                data.str <- paste(
                    paste(names(locTrans),collapse=', '),'\n',
                    paste(names(locRec),collapse=', '),'\n',
                    paste(fwhmTrans,collapse=', '),'\n',
                    paste(fwhmRec,collapse=', '),'\n',
                    paste(fwhmRange,collapse=', '),'\n',
                    'N/A\n',
                    paste(Pt/1e6,collapse=', '),' MW\n',
                    paste(Ne[k],collapse=', '),' m^-3\n',
                    paste(Tnoise,collapse=', '),' K\n',
                    paste(fradar/1e6,collapse=', '),' MHz\n',
                    paste(tau0,collapse=', '),' us\n',
                    paste(targetNoiseLevel,collapse=', '),'\n',
                    paste(dutyCycle,collapse=', ')
                    ,sep='')
            }else{
                data.str <- paste(
                    paste(names(locTrans),collapse=', '),'\n',
                    paste(names(locRec),collapse=', '),'\n',
                    paste(fwhmTrans,collapse=', '),'\n',
                    paste(fwhmRec,collapse=', '),'\n',
                    paste(fwhmRange,collapse=', '),'\n',
                    paste(resR,collapse=', '),'\n',
                    paste(Pt/1e6,collapse=', '),' MW\n',
                    paste(Ne[k],collapse=', '),' m^-3\n',
                    paste(Tnoise,collapse=', '),' K\n',
                    paste(fradar/1e6,collapse=', '),' MHz\n',
                    paste(tau0,collapse=', '),' us\n',
                    paste(targetNoiseLevel,collapse=', '),'\n',
                    paste(dutyCycle,collapse=', ')
                    ,sep='')
            }
        }else{
            data.str <- paste(
                paste(names(locTrans),collapse=', '),'\n',
                paste(names(locRec),collapse=', '),'\n',
                paste(fwhmTrans,collapse=', '),'\n',
                paste(fwhmRec,collapse=', '),'\n',
                paste(fwhmRange,collapse=', '),'\n',
                paste(resNS,resEW,resH,collapse=', '),'\n',
                paste(Pt/1e6,collapse=', '),' MW\n',
                paste(Ne[k],collapse=', '),' m^-3\n',
                paste(Tnoise,collapse=', '),' K\n',
                paste(fradar/1e6,collapse=', '),' MHz\n',
                paste(tau0,collapse=', '),' us\n',
                paste(targetNoiseLevel,collapse=', '),'\n',
                paste(dutyCycle,collapse=', ')
                ,sep='')

        }

  	if(is.null(iso.levels))
            {
                lev <- calc.levels(min(nLevel$intTime.isotropic[,,k],na.rm=T))
            }
  	else
            {
  		lev <- iso.levels
            }

        if(horCut)    plotOnMap(nLevel$intTime.isotropic[,,k],x=nLevel$x,y=nLevel$y,gdev=gdev,f=function(x){return(ceiling(log10(x)))},
                                col=rev(gray(seq(diff(zlim))/diff(zlim))),zlim=zlim,levels=lev,
                                cbarlab=expression(log[10]*"( integration time [s] )"),
                                refPoint=refPoint,main=paste('Isotropic parameters ',heights[k], ' km'),rxloc=nLevel$xyR,txloc=nLevel$xyT,
                                printInfo=printInfo,printString=data.str,totalInfo=F,useRaster=useRaster)

        if (is.null(vel.levels))
            {
                lev <- calc.levels(min(nLevel$intTime.velocity[,,k],na.rm=T))
            }
        else
            {
                lev <- vel.levels
            }

        if(horCut)    plotOnMap(nLevel$intTime.velocity[,,k],x=nLevel$x,y=nLevel$y,gdev=gdev,f=function(x){return(ceiling(log10(x)))},
                                col=rev(gray(seq(diff(zlimv))/diff(zlimv))),zlim=zlimv,levels=lev,
                                cbarlab=expression(log[10]*"( integration time [s] )"),
                                refPoint=refPoint,main=paste('Velocity ',heights[k],' km'),rxloc=nLevel$xyR,
                                txloc=nLevel$xyT,printInfo=printInfo,
                                printString=data.str,totalInfo=F,useRaster=useRaster)

        if (plotResolution)
            {

                max.res <- matrix(0,length(nLevel$x),length(nLevel$y))

                for (i in 1:length(nLevel$x))
                    {
                        for (j in 1:length(nLevel$y))
                            {
                                for (kk in 1:length(nLevel$noiseLevel.site))
                                    {
                                        max.res[i,j] <- max(max.res[i,j],nLevel$noiseLevel.site[[kk]]$res[,,k,3][i,j])
                                    }
                            }
                    }


                if(horCut)    plotOnMap(max.res,x=nLevel$x,y=nLevel$y,refPoint=refPoint,
                                        gdev=gdev,printInfo=T,main=paste('Height resolution ',heights[k],' km'),
                                        levels=c(seq(0,1,by=0.1),2:100),printString=data.str,useRaster=useRaster)

            }

    }


    # optional NS-cut
    if( any(NScut > 0 )){
        for( n in NScut){
            if( (n>0) & (n<=length(nLevel$x))){
                if(is.null(iso.levels))
                    {
                        lev <- calc.levels(min(nLevel$intTime.isotropic[n,,],na.rm=T))
                    }
                else
                    {
                        lev <- iso.levels
                    }
                plotNScut(nLevel$intTime.isotropic[n,,],f=function(x){return(ceiling(log10(x)))},
                          x=nLevel$x,y=nLevel$y,h=nLevel$h,refPoint=nLevel$refPoint,xlab='[km]',ylab='Height [km]',
                          cbarlab=expression(log[10]*"( integration time [s] )"),
                          main=paste('Isotropic parameters, x =',nLevel$x[n]),levels=lev,zlim=zlim,useRaster=useRaster,gdev=gdev)
                if(is.null(vel.levels))
                    {
                        lev <- calc.levels(min(nLevel$intTime.velocity[n,,],na.rm=T))
                    }
                else
                    {
                        lev <- iso.levels
                    }

                plotNScut(nLevel$intTime.velocity[n,,],f=function(x){return(ceiling(log10(x)))},
                          x=nLevel$x,y=nLevel$y,h=nLevel$h,refPoint=nLevel$refPoint,xlab='[km]',ylab='Height [km]',
                          cbarlab=expression(log[10]*"( integration time [s] )"),main=paste('Velocity, x =',nLevel$x[n]),
                          levels=lev,zlim=zlimv,useRaster=useRaster,gdev=gdev)
            }
        }
    }

    # optional height profile of the speed
    if( all(heightProf>0)){
        if(length(heightProf)==2){
            plotHeightProf(nLevel$intTime.isotropic[heightProf[1],heightProf[2],],f=function(x){return(log10(x))},h=nLevel$h,xlab=expression(log[10]*"( integration time [s] )"),ylab="Height [km]",xlim=zlim,main=paste('Isotropic parameters, x =',nLevel$x[heightProf[1]],' y =',nLevel$x[heightProf[2]]),gdev=gdev)
            plotHeightProf(nLevel$intTime.velocity[heightProf[1],heightProf[2],],f=function(x){return(log10(x))},h=nLevel$h,xlab=expression(log[10]*"( integration time [s] )"),ylab="Height [km]",xlim=zlimv,main=paste('Velocity, x =',nLevel$x[heightProf[1]],' y =',nLevel$x[heightProf[2]]),gdev=gdev)
        }
    }

    # warn if SNR > 1. The warning actually means that SNR is huge! Self-noise may dominate already at much lower values of SNR
    # that is calculated using only the thermal noise estimate
    maxsnr <- 0
    for( kk in seq(length(nLevel$noiseLevel.site))){
        maxsnr <- max( maxsnr , max(nLevel$noiseLevel.site[[kk]]$noiseLevel$SNR,na.rm=TRUE) )
    }
    if(maxsnr>1) cat("**********\n","   SNR > 1   ( ",maxsnr," )","\n**********\n",sep="")

    invisible(nLevel)

} # multistaticIntegrationTimes







calc.levels <- function(mm)
{
    ll <- NULL
    if (mm < 100000) ll <- c(10000,seq(20000,100000,by=20000),ll)
    if (mm < 10000) ll <- c(1000,seq(2000,10000,by=2000),ll)
    if (mm < 1000) ll <- c(100,seq(200,800,by=200),ll)
    if (mm < 100) ll <- c(10,seq(20,80,by=20),ll)
    if (mm < 10) ll <- c(1,seq(2,8,by=2),ll)
    if (mm < 1) ll <- c(0.1,seq(0.2,0.8,by=0.2),ll)
    if (mm < 0.1) ll <- c(0.01,seq(0.02,0.08,by=0.02),ll)
    if (mm < 0.01) ll <- c(0.001,seq(0.002,0.008,by=0.002),ll)

    return(ll)

}



integrationTimeVsArea <- function(A,levels=100,type='isotropic',plot=TRUE)
{
    # cell size in square kilometers
    x.step <- A$x[2] - A$x[1]
    y.step <- A$y[2] - A$y[1]
    cell.size <- x.step * y.step

    # Choose data
    if (type=='isotropic')
	{
            data <- A$intTime.isotropic
	}
    else if (type=='velocity')
	{
            data <- A$intTime.velocity
	}
    else
	{
            stop("Unknown type: ",type,sep="")
	}

    # Construct levels if given only the number of points
    if (length(levels)==1)
	{
            levels <- seq(min(data,na.rm=T),max(data,na.rm=T),len=levels)
	}

    areas <- rep(0,length(levels))

    for ( i in seq(levels))
	{
            areas[i] <- cell.size * sum(data<=levels[i],na.rm=T)
	}

    res <- list(levels=levels, areas=areas)

    if (!plot)
	{
            return(res)
	}
    else
	{
            plot(levels,areas,type='l',main='Integration time vs area',xlab='Integration time (seconds)',ylab='Area (square kilometers)')
            invisible(res)
	}
}




plotNScut <- function(resData,f=function(x){return(x)},x,y,h,refPoint,hlim=range(h),zlim=range(f(resData),na.rm=T,finite=T),col=rev(gray(seq(1000)/1000)),xlab='Latitude [degrees]', ylab='[km]',cbarlab='Resolution [km]',main='',levels=NA,gdev=x11,useRaster=FALSE){

    xlim <- range(y)

    resd <- f(resData)
    resd[resd>zlim[2]] <- zlim[2]
    resd[is.nan(resd)] <- zlim[2]
    resd[is.na(resd)] <- zlim[2]

    gdev(width=6,height=6/.88)
    layout(matrix(c(1,2),ncol=1),heights=c(.88,.12))
    par(mar=c(3,3,2,2),mgp=c(1.5,.5,0))

    # Draw filled contours
    image(y,h,resd,col=col,xlim=xlim,ylim=hlim,zlim=zlim,xlab=xlab,ylab=ylab,main=main,useRaster=useRaster)

    # Draw contour lines
    if(!is.null(levels)) contour(y,h,resData,add=T,levels=levels,lwd=2,col='red',labcex=1)

    # Draw the colorbar
    par(mar=c(3,3,0,2))

    image(seq(zlim[1],zlim[2],length.out=1000),c(1),matrix(seq(zlim[1],zlim[2],length.out=1000),ncol=1),col=col,yaxt='n',ylab='',xlab=cbarlab ,bty='L',useRaster=useRaster)

} # plotNScut

plotHeightProf <- function(resData,f=function(x){return(x)},h,xlab="( integration time [s] )",ylab="Height [km]",xlim=range(f(resData)),main='',gdev=x11){

    resd <- f(resData)

    gdev(width=6,height=6/.88)

    par(mar=c(3,3,2,2),mgp=c(1.5,.5,0))

    plot(resd,h,xlab=xlab,ylab=ylab,main=main,xlim=xlim)

} # plotHeightProf

