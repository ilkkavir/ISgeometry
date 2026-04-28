\name{parameterErrorEstimatesPrior}
\alias{parameterErrorEstimatesPrior}
\title{Plasma parameter error estimates for an ISR with prior information}
\description{
Calculates plasma parameter error estimates at a specific geocentric
location for a given radar system configuration, taking into account
prior information and realistic SNR.
}
\usage{
parameterErrorEstimatesPrior(lat, lon, alt, Ne, Ti, Te, Coll, Comp, 
    fwhmRange, resR, intTime, NePr = 1e+12, TiPr = 1e+04, TePr = 1e+04, 
    CollPr = 0, ViPr = 1e+04, CompPr = 0, pm0 = c(30.5, 16), Tnoise = 300, 
    Pt = 3.5e+06, locTrans = SKI, locRec = list(SKI, KAR, KAI), 
    fwhmTrans = 2.1, fwhmRec = c(1.2, 1.7, 1.7), RXduty = c(0.75, 1, 1), 
    mineleTrans = 30, mineleRec = 30, fradar = 233e+06, phArrTrans = TRUE, 
    phArrRec = TRUE, fwhmIonSlab = 50, dutyCycle = 0.25, maxLag)
}
\arguments{
  \item{lat}{Geocentric latitude of the measurement volume [deg north].}
  \item{lon}{Geocentric longitude of the measurement volume [deg east].}
  \item{alt}{Altitude of the measurement volume [km].}
  \item{Ne}{Electron density [m^-3].}
  \item{Ti}{Ion temperature [K].}
  \item{Te}{Electron temperature [K].}
  \item{Coll}{Ion-neutral collision frequency [Hz].}
  \item{Comp}{Ion composition (fraction of ion with mass pm0[2] out of the total ion number density).}
  \item{fwhmRange}{Range resolution in ACF decoding [km].}
  \item{resR}{Range resolution in the plasma parameter fit [km].}
  \item{intTime}{Integration time [s].}
  \item{NePr}{Prior standard deviation for electron density [m^-3].}
  \item{TiPr}{Prior standard deviation for ion temperature [K].}
  \item{TePr}{Prior standard deviavtion for electron temperature [K].}
  \item{CollPr}{Prior standard deviation for collision frequency [Hz].}
  \item{ViPr}{Prior standard deviation for ion velocity [m/s].}
  \item{CompPr}{Prior standard deviation for ion composition [dimensionless].}
  \item{pm0}{Ion masses [amu], default \code{c(30.5, 16)}.}
  \item{Tnoise}{Receiver noise temperature [K].}
  \item{Pt}{Transmitter power(s) [W].}
  \item{locTrans}{Transmitter location(s). A list containing lat, lon, and optionally height [deg, deg, km].}
  \item{locRec}{Receiver locations. A list of locations (lat, lon, height).}
  \item{fwhmTrans}{Transmitter beam width(s) at zenith [degrees FWHM].}
  \item{fwhmRec}{Receiver beam widths at zenith [degrees FWHM].}
  \item{RXduty}{Receiver duty cycle for each receiver.}
  \item{mineleTrans}{Minimum elevation limit of the transmitter [deg].}
  \item{mineleRec}{Minimum elevation limit of the receivers [deg].}
  \item{fradar}{Radar system carrier frequency [Hz].}
  \item{phArrTrans}{Logical; whether the transmitter is a phased array.}
  \item{phArrRec}{Logical; whether the receivers are phased arrays.}
  \item{fwhmIonSlab}{Full width at half maximum of the ion slab [km].}
  \item{dutyCycle}{Transmitter duty cycle.}
  \item{maxLag}{Maximum lag for the ACF [us]. 5 times the ACF lag
    corresponding to ion thermal velocity by default.}
}
\value{
Returns error estimates for the plasma parameters (Ne, Ti, Te, Coll, Vi,
Comp) based on the specified radar system geometry and plasma
parameters. The output list errtabs contains two lists:

errtabs$los is a list of error vectors for each individual site. The
sites are in the same order as in the input argument locRec.

errtabs$multistatic is an error vector for multistatic analysis. The
velocity error estimate is the square root of the trace of the error
covariance matrix of three orthogonal velocity vector components, and
thus gives an upper limit  that the standard deviation of any projection
of the velocity vector cannot exceed. If errors of all vector components
are equal, the true error is sqrt(3) times the given error. 
}
\author{
Ilkka Virtanen \email{ilkka.i.virtanen@oulu.fi}
}

\references{

    1. Matti Vallinkoski, Statistics of incoherent scatter multiparameter
    fits, Journal of Atmospheric and Terrestrial Physics, 50 (9),
    839-851, 1988, https://doi.org/10.1016/0021-9169(88)90106-7
    
    2. Markku Lehtinen, Ilkka I. Virtanen, and Mikko R. Orispää,
    EISCAT_3D Measurement Methods Handbook, 2014,
    http://urn.fi/urn:isbn:9789526205854

}

\examples{

## F region O+ ion fraction fit with electron density prior e.g. from
## plasma line data

FitErrors <- parameterErrorEstimatesPrior(SKI[1],SKI[2],250,1e11,2000,3000,1,.8,10,10,60,NePr=1e9,CompPr=.5,maxLag=500)

}