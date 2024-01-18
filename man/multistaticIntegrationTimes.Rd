\name{multistaticIntegrationTimes}
\alias{multistaticIntegrationTimes}
\title{Incoherent integration time estimation}
\description{
  Incoherent integration time required to reach a given ACF noise level
  in given points in space.

  The radar beams are modelled as 3D normal distributions to enable fast
  calculations.

  See references for definition of the ACF noise level.
}

\usage{
       multistaticIntegrationTimes(
        refPoint = SKI,
        locTrans = SKI,
          locRec = list(SKI,KAR,KAI),
           locxy = FALSE,
       fwhmTrans = c(2.2),
         fwhmRec = c(1.3,1.8,1.8),
       fwhmRange = c(1),
               x = seq(-300,300,by=25),
               y = seq(-300,300,by=25),
         heights = c(300),
        infinity = defaultInfinity(),
           resNS,
           resEW,
            resH,
            resR,
              Pt = c(3.5e6),
              Ne = 1e11,
     fwhmIonSlab = 100,
          Tnoise = c(300),
          fradar = 233e6,
            tau0 = 100,
      phArrTrans = TRUE,
        phArrRec = TRUE,
targetNoiseLevel = .01,
       dutyCycle = .25,
          RXduty = 1,
            gdev = NULL,
            zlim = c(0,5),
           zlimv = c(2,7),
      iso.levels = c(.1,.3,.7,1,3,7,10,33,67,100,333,667,1000,3333,6667,10000,33333,666667,1000000,3333333,6666667,1e7,3e7,7e7,1e8,3e8,7e8,1e9),
      vel.levels = c(.1,.3,.7,1,3,7,10,33,67,100,333,667,1000,3333,6667,10000,33333,666667,1000000,3333333,6666667,1e7,3e7,7e7,1e8,3e8,7e8,1e9),
       printInfo = FALSE,
  plotResolution = FALSE,
           NScut = 0,
          horCut = TRUE,
         verbose = FALSE,
     mineleTrans = 30,
       mineleRec = 30,
       useRaster = FALSE
       )
}

\arguments{
\item{refPoint}{c(lat,lon) reference point from which the distance x and y are measured}
\item{locTrans}{list(c(lat,lon)) or list(c(lat,lon,height)) latitudes, longitudes (and heights) of the transmitter antennas.
           Angles in degrees. OR list(c(x,y)) position of the transmitter, IF locxy=T}
\item{locRec}{     list(c(lat,lon)) or list(c(lat,lon,height)) latitudes, longitudes (and heights) of the receiver antennas.
           Angles in degrees. OR list(c(x,y)) IF locxy=T}
\item{locxy}{      logical, if T, the positions of the transmitter and
  receiver are given in cartesian coordinates, otherwise in latitudes
  and longitudes}
\item{fwhmTrans}{  c(fwhmTrans1,fwhmTrans2,...) full width at half maximum of the transmitter beams when pointed to zenith, in degrees. If only one value is given, it is used for all transmitters.}
\item{fwhmRec}{    c(fwhmRec1,fwhmRec2,...) full width at half maximum of the receiver beams when pointed to zenith, in degrees. If only one value is given, it is used for all receivers.}
\item{fwhmRange}{  c(fwhmRange1,fwhmRange2,...) full width at half maximum of the gaussian range ambiguity function [km]. If only one value is given, it is used for all transmitters.}
\item{x}{vector of grid points in West-East direction, relative to the reference point.}
\item{y}{vector of grid points in North-South direction, relative to the reference point.}
\item{heights}{vector of grid points in vertical direction, in kilometers. Separate maps are plotted for each height.}
\item{infinity}{   a large value used as standard deviation along the beam direction}
\item{resNS}{ (optional) north-south resolution of the integrated resolution cells [km]}
\item{resEW}{ (optional) east-west resolution of the integrated resolution cells [km]}
\item{resH}{  (optional) height resolution of the integrated resolution cells [km]}
\item{resR}{  (optional) integrated range resolution [km], used if some
  of the above three is not given. If missing, the resolution is matched
with fwhmRange}
\item{Pt}{         c(Pt1,Pt2,...) transmitter powers [W]. If only one value is given, it is used for all transmitters.}
\item{Ne}{         electron density [m^-3]}
\item{fwhmIonSlab}{Slab thickness of the ionosphere. Used in self-noise estimation}
\item{Tnoise}{     system noise temperatures of the receivers [K]. If only one value is given, it is used for all receivers.}
\item{fradar}{     radar carrier frequency [Hz]}
\item{tau0}{       ACF time-scale [us]}
\item{phArrTrans}{     logical, if TRUE, phased-array antennas with beam widening as function of tilt angle are assumed, otherwise beam width does not depend on pointing direction. If only one value is given, it is used for all transmitters.}
\item{phArrRec}{     logical, if TRUE, phased-array antennas with beam widening as function of tilt angle are assumed, otherwise beam width does not depend on pointing direction. If only one value is given, it is used for all receivers.}
\item{targetNoiseLevel}{  target noise level value, same target value used for both isotropic and velocity integration times}
\item{dutyCycle}{  transmitter duty cycle}
\item{RXduty}{"Receiver duty cycle"}
\item{gdev}{       graphics device used for plotting; default for Mac OS X is 'quartz', for Unix/Linux 'X11' and for Windows 'windows'.}
\item{zlim}{       zlim for the plots of isotropic parameter integration times}
\item{zlimv}{      zlim for the plots of velocity integration times}
\item{iso.levels}{     contour line levels for isometric parameter plot}
\item{vel.levels}{     contour line levels for velocity plot}
\item{printInfo}{ logical; if T, plots used parameters below the map plot}
\item{plotResolution}{logical; if T,plots resolution map in addition to the isotropic and velocity integration time maps.}
\item{plotResolution}{logical; if T, plots resolution map in addition to the isotropic and velocity integration time maps.}
\item{NScut}{a vector of indices of north-south directed vertical slices to plot. nothing will be plotted if NScut <= 0}
\item{horCut}{if FALSE, the plots of horisontal cuts will not be produced. Useful e.g. when producing only a single vertical cut with x=c(0).}
\item{verbose}{optional printing of average SNR / power etc.}
\item{mineleTrans}{mimimum allowed elevation for transmitter(s) in degrees [0 90]}
\item{mineleRec}{minimum allowed elevation for receiver(s) in degrees [0 90]}
\item{useRaster}{rasterize images? Works only with regular grids!}
	}

\value{Plots integration time maps for isotropic parameters and
  velocity. Optionally, if plotResolution is TRUE, also plots the resolution map. If called with a return variable, returns a list of noise levels, integration times, resolutions, etc. at the grid points.}

\details{}


%\note{}

\author{Ilkka Virtanen (University of Oulu, Finland)
  \email{ilkka.i.virtanen@oulu.fi}}

%\seealso{}

\references{

    1. Matti Vallinkosi, Statistics of incoherent scatter multiparameter
    fits, Journal of Atmospheric and Terrestrial Physics, 50 (9),
    839-851, 1988, https://doi.org/10.1016/0021-9169(88)90106-7
    
    2. Markku Lehtinen, Ilkka I. Virtanen, and Mikko R. Orispää,
    EISCAT_3D Measurement Methods Handbook, 2014,
    http://urn.fi/urn:isbn:9789526205854

}


\examples{

## EISCAT 3D system with the default values
multistaticIntegrationTimes()

## Four remotes plus a core transmit/receive site at Kiruna. Positions
## given in the xy-coordinate system:

multistaticIntegrationTimes(x=seq(-300,300,by=25),y=seq(-300,300,by=25),locTrans=c(0,0),locRec=list(c(0,0),c(-200,0),c(200,0),c(0,-200),c(0,200)),refPoint=KIR,locxy=T,Tnoise=c(100,40,40,40,40,40),heights=c(100,300,500),fwhmTrans=.6,fwhmRec=.6)


## save figure number 3 as pdf file "figure3.pdf"

savepdf(3,'figure3.pdf')


## evalution of the EISCAT UHF system assuming 10 x 10 x 10 km
## resolution cells
multistaticIntegrationTimes(x=seq(-300,300,by=25),y=seq(-300,300,by=25),locTrans=TRO,locRec=list(TRO,KIR,SOD),refPoint=KIR,locxy=F,Tnoise=c(100,40,40),heights=c(100,300,500),fwhmTrans=.6,fwhmRec=.6,phArr=F,resNS=10,resEW=10,resH=10)


## Using both steerable and phased arrays
multistaticIntegrationTimes(locTrans=list(KIR,SOD,KIL,JAR,KAU),locRec=list(MUO,OVE,KIR,SOD,KIL),phArrTrans=c(T,T,T,F,F),phArrRec=c(F,F,T,F,T),locxy=FALSE)

}


%\keyword{algebra}
