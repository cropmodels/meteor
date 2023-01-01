\name{globe}

\alias{Tg}
\alias{Tg,data.frame-method}
\alias{Tg,SpatRasterDataset-method}

\alias{Tnwb}
\alias{Tnwb,data.frame-method}
\alias{Tnwb,SpatRasterDataset-method}

\alias{WBGT}
\alias{WBGT,data.frame-method}
\alias{WBGT,SpatRasterDataset-method}


\title{globe temperature}

\description{
Globe temperature for either a \code{data.frame} or a \code{SpatRasterDataset}. These must have variables "temp" (C), "rhum" (\%), "wind" (m s-1), and "srad" (J s-1 m-2). The \code{data.frame} must also have a variable "date".
}

\usage{
\S4method{Tg}{data.frame}(x, latitude)
\S4method{Tg}{SpatRasterDataset}(x, filename="", overwrite=FALSE, ...)

\S4method{Tnwb}{data.frame}(x, latitude)
\S4method{Tnwb}{SpatRasterDataset}(x, filename="", overwrite=FALSE, ...)

\S4method{WBGT}{data.frame}(x, latitude)
\S4method{WBGT}{SpatRasterDataset}(x, filename="", overwrite=FALSE, ...)
}


\arguments{
  \item{x}{data.frame or SpatRasterDataset}
  \item{latitude}{numeric}
  \item{filename}{character. Output filename}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}  
  \item{...}{additional arguments for writing files as in \code{\link[terra]{writeRaster}}}    
}

\value{
numeric or SpatRaster
}

\examples{
wd <- data.frame(date=as.Date("2003-08-28") + 1:3, temp=c(19.1, 20.6, 19.4), rhum=c(66,71,73), wind=c(3.3, 1.9, 1.1), srad=c(168, 178, 125))
Tg(wd, 40.96) 

Tnwb(wd, 40.96) 

WBGT(wd, 40.96) 


library(terra)
r <- rast(ncol=2, nrow=2, nlyr=1)
temp <- setValues(r, 21:24)
time(temp) <- as.Date("2000-01-01")
rhum <- setValues(r, 81:84)
wind <- setValues(r, 9:12)
srad <- setValues(r, 100:103)
s <- sds(list(temp=temp, rhum=rhum, wind=wind, srad=srad))

x <- Tg(s)
y <- WBGT(s)
}
