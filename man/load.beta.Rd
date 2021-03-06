% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/source.beta.r
\name{load.beta}
\alias{load.beta}
\title{Load Beta Count Data}
\usage{
load.beta(filepath, tz = "US/Eastern", lookup = NULL)
}
\arguments{
\item{filepath}{The filepath to the raw count file.}

\item{tz}{The timezone of the computer that originally saved the datafile.}
}
\value{
It will return a list containing four objects: (1) _raw_ A list with the raw file contents, (2) meta A list containing the file metadata including instrument serial number and settings, (3) proc The processed and abridge dataset from the file, and (4) counts A dataframe containing the counting data from the file.
}
\description{
A function to read in and parse a raw beta coutner file from a Riso multicollector.
}
\author{
Thomas Bryce Kelly
}
\keyword{Beta}
\keyword{Thorium}
