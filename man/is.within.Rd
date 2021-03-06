% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/source.r
\name{is.within}
\alias{is.within}
\title{Is Within}
\usage{
is.within(x, bounds)
}
\arguments{
\item{x}{A vector of any type with strict ordering (i.e. where '>' or '<' are valid operators).}

\item{bounds}{A vector of two entries, where bounds[1] is the lower bound and bounds[2] is the upper.}
}
\description{
Filter a vector of entries with fall between a set of bounds (i.e. on a closed interval).
}
\author{
Thomas Bryce Kelly
}
