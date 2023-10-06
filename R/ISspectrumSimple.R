ISspectrumSimple <- function(p=c(1e11,300,300,0,0,0),pm0=c(30.5,16),fradar=223e6,scattAngle=180,freq=seq(-100000,100000,by=10)*fradar/933.5e6){


    ele <- c(p[1],p[3],p[4]*0.35714,p[5])
    ion <- list(
        c(pm0[1],(1-p[6]),p[2],p[4],p[5]),
        c(pm0[2],p[6],p[2],p[4],p[5])
    )

    return(ISgeometry:::incoherentScatterSpectrum(ele=ele,ion=ion,fradar=fradar,scattAngle=scattAngle,freq=freq))

}
