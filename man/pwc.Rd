\name{pwc}

\alias{pwc}
\alias{pwc,numeric-method}
\alias{pwc,SpatRaster-method}


\title{pwc}

\description{
pwc
}

\usage{
\S4method{pwc}{numeric}(x, input="wbgt", adjust=TRUE)

\S4method{pwc}{SpatRaster}(x, input="wbgt", adjust=TRUE, filename="", overwrite=FALSE, ...)
}

\arguments{
  \item{x}{numeric or SpatRaster}
  \item{input}{character. One of "wbgt" or "utci"}
  \item{adjust}{logical. If \code{TRUE}, the Smallcombe et al. (2022) adjustment for a 7-hour workday is used}
  \item{filename}{character. Output filename}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}  
  \item{...}{additional arguments for writing files as in \code{\link[terra]{writeRaster}}}    
}

\value{
numeric or SpatRaster
}

\references{
Smallcombe et al.. 2022

Foster et al., 2022
}


\examples{
pwc(25)
}

