\name{utci}

\alias{utci}
\alias{utci,data.frame-method}
\alias{utci,SpatRasterDataset-method}


\title{utci}

\description{
utci
}

\usage{
\S4method{utci}{data.frame}(x)
}

\arguments{
  \item{x}{data.frame}
 }

\value{
numeric or SpatRaster
}

\examples{
d <- data.frame(temp=20, rhum=80, tglb=22, wind=10)
utci(d)
}

