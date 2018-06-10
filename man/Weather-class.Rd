\name{Weather-class}

\docType{class}

\alias{Weather-class}


\alias{show,Weather-method}
\alias{$,Weather-method}
\alias{$<-,Weather-method}
\alias{[[,Weather-method}
\alias{[[<-,Weather-method}


\title{Weather class}

\description{
Weather data
}


\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Weather", ...)}, or with the helper functions such as \code{weather}.
}

\section{Slots}{
Slots of Weather objects
	\describe{
    \item{\code{data}:}{data.frame with the weather data} 
    \item{\code{ID}:}{character} 
    \item{\code{name}:}{character} 
    \item{\code{country}:}{character} 
    \item{\code{longitude}:}{numeric} 
    \item{\code{latitude}:}{numeric} 
    \item{\code{elevation}:}{numeric} 
  }
}


\examples{
showClass("Weather")
}

\keyword{classes}
