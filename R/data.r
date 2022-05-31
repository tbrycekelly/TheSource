#' @title Global Bathymetry
#'
#' A Matrix containing elevation data taken from NOAA via marmap.
#' @docType data
#' @format A matrix containing 2160 x 1080 values of elevation. Depth is given in negative. Resolution is 10 arc minutes (1/6th degree).
#' @source NOAA Bathymetric Database
"bathy.global"


#' @title Pacific Bathymetry
#'
#' A Matrix containing elevation data taken from NOAA via marmap.
#' @docType data
#' @format A list
#' @source NOAA Bathymetric Database
"bathy.pacific"


#' @title Arctic Bathymetry
#'
#' A Matrix containing elevation data taken from NOAA via marmap.
#' @docType data
#' @format A matrix containing 2160 x 1080 values of elevation. Depth is given in negative. Resolution is 10 arc minutes (1/6th degree).
#' @source NOAA Bathymetric Database
"bathy.arctic"


#' @title Antarctic Bathymetry
#'
#' A Matrix containing elevation data taken from NOAA via marmap.
#' @docType data
#' @format A matrix containing 2160 x 1080 values of elevation. Depth is given in negative. Resolution is 10 arc minutes (1/6th degree).
#' @source NOAA Bathymetric Database
"bathy.antarctic"


#' @title Gulf of Mexico Bathymetry
#'
#' A Matrix containing elevation data taken from NOAA via marmap.
#' @docType data
#' @format A matrix containing 2160 x 1080 values of elevation. Depth is given in negative. Resolution is 10 arc minutes (1/6th degree).
#' @source NOAA Bathymetric Database
"bathy.gom"

#' @title Geotraces' GA01 Salinity and Temperature Data
#' A dataframe containing T and S with various metadata.
#' @source See the Geotraces Intermediate Data Product and refences at http://www.geotraces.org/dp/idp2017
#' @format A Dataframe containing 5103 observations of 12 variables
"Geotraces.GA01.Sal"

#' @title Sample FRRF Datafile
#' A sample of an FRRF datafile for use in the demonstrations and tutorials.
#' @source Sven A Kranz, see Kranz et al. 2020 in JGR: Oceans.
#' @format a list containing the datastructure returned by TheSource::load.frrf()
"sample.frrf"


#' @title World Coastline Data
#' This is a coarse resolution coastline, copied from the OCE package, at scale 1:110M, with 10,696 points, suitable for world-scale plots plotted at a small size, e.g. inset diagrams.
#' @source Downloaded from https://www.naturalearthdata.com, in ne_110m_admin_0_countries.shp in July 2015, with an update on December 16, 2017. Copied from the OCE R package.
"coastlineWorld"


#' @title World Coastline Data (Medium Resolution)
#' This is a coarse resolution coastline, copied from the OCE package, at scale 1:110M, with 10,696 points, suitable for world-scale plots plotted at a small size, e.g. inset diagrams.
#' @source Downloaded from https://www.naturalearthdata.com, in ne_110m_admin_0_countries.shp in July 2015, with an update on December 16, 2017. Copied from the OCE R package.
"coastlineWorldMedium"


#' @title World Coastline Data (High Resolution)
#' This is a coarse resolution coastline, copied from the OCE package, at scale 1:110M, with 10,696 points, suitable for world-scale plots plotted at a small size, e.g. inset diagrams.
#' @source Downloaded from https://www.naturalearthdata.com, in ne_110m_admin_0_countries.shp in July 2015, with an update on December 16, 2017. Copied from the OCE R package.
"coastlineWorldFine"


#' @title World Coastline Data (L1)
#' This is a coarse resolution coastline, compiled from the global hierarchical coastline data product.
#' @source
"coastline1"


#' @title World Coastline Data (L2)
#' This is a coarse resolution coastline, compiled from the global hierarchical coastline data product.
#' @source
"coastline2"


#' @title World Coastline Data (L3)
#' This is a coarse resolution coastline, compiled from the global hierarchical coastline data product.
#' @source
"coastline3"


#' @title World Coastline Data (L4)
#' This is a coarse resolution coastline, compiled from the global hierarchical coastline data product.
#' @source
"coastline4"


#' @title World Coastline Data (L5)
#' This is a coarse resolution coastline, compiled from the global hierarchical coastline data product.
#' @source
"coastline5"


#' @title World EEZ Data
#' Global maritime eez boundaries. Creative Commons Attribution 4.0
#' @source Flanders Marine Institute (2019). Maritime Boundaries Geodatabase: Maritime Boundaries and Exclusive Economic Zones (200NM), version 11. Available online at https://www.marineregions.org/. https://doi.org/10.14284/386
"eez"
