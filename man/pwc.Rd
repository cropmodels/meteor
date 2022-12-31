\name{pwc}

\alias{pwc}
\alias{pwc,data.frame-method}
\alias{pwc,SpatRaster-method}


\title{pwc}

\description{
pwc
}

\usage{
\S4method{pwc}{data.frame}(x, input="wbgt")

\S4method{pwc}{SpatRaster}(x, input="wbgt", filename="", overwrite=FALSE, ...)
}

\arguments{
  \item{x}{data.frame or SpatRaster}
  \item{input}{charcter. One of "wbgt" or "ucti"}
  \item{filename}{character. Output filename}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}  
  \item{...}{additional arguments for writing files as in \code{\link[terra]{writeRaster}}}    
}

\value{
numeric or SpatRaster
}

\examples{
pwc(25)
}

