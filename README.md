# ISgeometry


[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6623186.svg)](https://doi.org/10.5281/zenodo.6623186)


Incoherent scatter radar performance evaluation


The package can be used for,

1. Estimating ACF noise levels in incoherent scatter radar measurements at given locations in space, when the radar system and electron density are known.

2. Estimating standard deviations of fitted Ne, Ti, Te, and Vi when the radar system and all plasma parameters are known.

See also the [e3doubt](https://github.com/Dartspacephysiker/e3doubt) python front-end, which implements also ionospheric models and beam scan pattern.

## Examples

### Incoherent integration time speed maps for the EISCAT3D system with default values

    multistaticIntegrationTimes()

### Incoherent integration time speed maps with a core site in Kiruna and four remote sites. Positions given in cartesian xy-coordinates in km.

    multistaticIntegrationTimes(x=seq(-300,300,by=25),y=seq(-300,300,by=25),locTrans=c(0,0),locRec=list(c(0,0),c(-200,0),c(200,0),c(0,-200),c(0,200)),refPoint=KIR,locxy=T,Tnoise=c(100,40,40,40,40,40),heights=c(100,300,500),fwhmTrans=.6,fwhmRec=.6)


### Standard deviations of plasma parameters
Errors of Ne, Ti, Te, and Vi at 300 km altitude above Skibotn, when Ne=1e11 m^-3, Ti=2000 K, Te=3000 K, ion-neutral collision frequency is negligibly small, only O+ ions, 5 km modulation bit length, post-integration to 20 km resolution and 60 s integration time.

    parameterErrorEstimates(lat=SKI[1],lon=SKI[2],alt=300,Ne=1e11,Ti=2000,Te=3000,Coll=0,Comp=1,fwhmRange=5,resR=20,intTime=60)


##References
Matti Vallinkosi, Statistics of incoherent scatter multiparameter fits, Journal of Atmospheric and Terrestrial Physics, 50 (9), 839-851, 1988, [https://doi.org/10.1016/0021-9169(88)90106-7](https://doi.org/10.1016/0021-9169(88)90106-7)

Markku Lehtinen, Ilkka I. Virtanen, and Mikko R. Orispää, EISCAT_3D Measurement Methods Handbook, 2014, [http://urn.fi/urn:isbn:9789526205854](http://urn.fi/urn:isbn:9789526205854)
