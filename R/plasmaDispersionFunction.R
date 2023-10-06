#
# General plasma dispersion function calculation for single scalar value z
#
plasmaDispersionFunction <- function(z){
    z1 <- -1i*z
    x  <- Re(z1)
    y  <- Im(z1)
#    if(abs(2*x*y)<30000){
        cs <- cos(2*x*y)
        sn <- sin(2*x*y)
#    }else{
#        cs <- 1
#        sn <- 0
#    }
    pisqr <- sqrt(pi)
    fn <- seq( max(1,floor(abs(2*y))-11) , (floor(abs(2*y))+11) )
    term1 <- exp(-fn*fn/4-y*y)
    if( max(abs(term1)) < 1e-100 ) term1 <- rep(0,length(fn))
    term2 <- exp( -(fn/2-y) * (fn/2-y) )/2
    term3 <- exp( -(fn/2+y) * (fn/2+y) )/2
    factor <- rep(1,length(fn)) / (fn*fn+4*x*x)
    sume <- sum(factor*(term1*cs-term2-term3))
    sumf <- sum(factor*(term1*2*x*sn+fn*(term2-term3)))
    if(abs(x*y)<1e-4){
        e1 <- (-exp(-y*y)*x*y*y/2+2*x*sume ) / pisqr
        f1 <- ( exp(-y*y)*y/2+sumf ) / pisqr
    }else{
        e1 <- ( exp(-y*y)*(cs-1)/4/x + 2*x*sume ) / pisqr
        f1 <- ( exp(-y*y)*sn/4/x+sumf ) / pisqr
    }
    if(abs(e1)<1e-100) e1 <- 0
    z1 <- 2i*(e1+1i*f1) + 1i*pisqr*exp(-y*y)*(cs+1i*sn)*cerfexp(-x)
    return(z1)

} # pldf
