# cerfexp.m: Approximation for the function  (1-erf(x))*exp(x**2), used in pldf
# GUISDAP v.1.60 96-05-27 Copyright Asko Huuskonen and Markku Lehtinen
#
# Conversion to R by Ilkka Virtanen 2011
#
#
# Approximation for the function  (1-erf(x))*exp(x**2) with
# accuracy 1e-8 for all x.  Method: power series for
# erf(x) for small x, asymptotic continued fractions expansion
# for large x (Abramovitz-Stegun formula 7.1.14)
#                                             Markku Lehtinen 2.7.1979
#
## res=cerfexp(x)
##
#  function res=cerfexp(x)
##
#  x1=abs(x);
#  sqx=x*x;
#     if (x1<1.) 
#        tail=0.;
#        term=sqx*sqx/10.;
#        erf=1.-sqx/3.+term;
#        for n=3:12 
#           term=-term*sqx*(2*n-1.)/(2*n+1.)/n;
#           tail=tail+term;
#        end 
#        erf=(erf+tail)*x*2./sqrt(pi);
#        res=(1-erf)*exp(sqx);
#      else
#        dum=0;
#        for n=(55:-1:1) 
#           dum=n/2/(dum+x1);
#         end 
#         res=1/(dum+x1)/sqrt(pi);
#         if (x<0), res=2*exp(sqx)-res; end 
#      end

cerfexp <- function(x){
  x1 <- abs(x)
  sqx <- x*x
  if(x1<1.){
    tail <- 0.
    term <- sqx*sqx/10
    erf <- 1.-sqx/3.+term
    for(n in seq(3,12)){
      term <- -term*sqx*(2*n-1.)/(2*n+1.)/n
      tail <- tail+term
    }
    erf <- (erf+tail)*x*2./sqrt(pi)
    res <- (1-erf)*exp(sqx)
  }else{
    dum <- 0
    for(n in seq(55,1,by=-1)){
      dum <- n/2/(dum+x1)
    }
    res <- 1/(dum+x1)/sqrt(pi)
    if(x<0) res <- 2*exp(sqx)-res
  }
  return(res)
}
