% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/source.maps.r
\name{add.map.layer}
\alias{add.map.layer}
\title{Add Map Layer}
\usage{
add.map.layer(
  lon,
  lat,
  z,
  zlim = NULL,
  pal = "greyscale",
  col.na = NA,
  n = 255,
  refinement = NA,
  col.low = NA,
  col.high = NA,
  rev = FALSE,
  filled = FALSE,
  indicate = TRUE,
  verbose = FALSE
)
}
\arguments{
\item{lon}{longitude of the data}

\item{lat}{latitude of the data}

\item{z}{the data}

\item{zlim}{the limits of the z-axis. A value of NULL indicates that zlim should equal the range of z values.}

\item{pal}{the name of a palette generating function. Should be a string.}

\item{col.na}{the color with which any NAs should be drawn with. A value of NA skips this step.}

\item{n}{the number of distinct colors to request from the palette function.}

\item{refinement}{the level of bilinear refinement to apply to the image grid}
}
\description{
Add a image layer to the map!
}
\author{
Thomas Bryce Kelly
}
