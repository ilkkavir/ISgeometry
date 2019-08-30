\name{integrationTimeVsArea}
\alias{integrationTimeVsArea}
\title{integrationTimeVsArea}


\description{Calculates and/or plots the covered area as function of integration time}

\usage{
integrationTimeVsArea(A,levels=100,type='isotropic',plot=TRUE)
}

\arguments{
\item{A}{A list given by \link{multistaticIntegrationTimes}}
\item{levels}{Either a single number giving the number of integration time levels calculated, or a vector
				giving the integration time levels explicitly.}
\item{type}{Either 'isotropic' or 'velocity'.}
\item{plot}{If true, a graph is drawn. If false, only data is returned.}

  }

  \value{A list containing the integration times (levels) in seconds and areas in square kilometers (areas).
  		If plot=T and no return variable is given, nothing is returned.}

  \details{}

  \author{Mikko Orispaa \email{mikko.orispaa@oulu.fi}}
  
  
  \examples{
  
  ## First calculate some integration times.
  ## Note the return variable A.
  A <- multistaticIntegrationTimes(fradar=235e6,Pt=2e6)
  
  ## Plot the integration times vs areas
  integrationTimeVsArea(A,type='isotropic')
  integrationTimeVsArea(A,levels = seq(1,100,by=5), type='velocity')
  
  ## Do not plot, only calculate
  B <- integrationTimeVsArea(A,type='isotropic',plot=F)
  
  }