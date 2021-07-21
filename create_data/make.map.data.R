library(TheSource)
library(marmap)


## Global
bathy = getNOAA.bathy(-180, 180, -90, 90, resolution = 15)
bathy.global = list(Lon = as.numeric(rownames(bathy)),
                    Lat = as.numeric(colnames(bathy)),
                    Z = bathy, res = 15)
map = make.map()
add.map.bathy(map, bathy.global, bathy.levels = -1000)

save(bathy.global, file = 'data/bathy.global.rdata')


## Save to netcdf
library(ncdf4)
lon1 <- ncdim_def("longitude", "degrees_east", bathy.global$Lon)
lat2 <- ncdim_def("latitude", "degrees_north", bathy.global$Lat)
mv <- -999999 #missing value to use
var_temp <- ncvar_def("depth", "meters", list(lon1, lat2), longname="NOAA bathymetry data pulled from gis.ngdc.noaa.gov", mv)

ncnew <- nc_create('./inst/extdata/bathy.global.nc', list(var_temp))
print(paste("The file has", ncnew$nvars,"variables"))
print(paste("The file has", ncnew$ndim,"dimensions"))

ncvar_put(ncnew, var_temp, bathy.global$Z, start = c(1,1),
          count = c(length(bathy.global$Lon), length(bathy.global$Lat)))
nc_close(ncnew)



## Arctic
bathy = getNOAA.bathy(-180, 180, 60, 90, resolution = 10)
bathy.arctic = list(Lon = as.numeric(rownames(bathy)),
                    Lat = as.numeric(colnames(bathy)),
                    Z = bathy, res = 10)

save(bathy.arctic, file = 'data/bathy.arctic.rdata')

map = make.map.arctic(dlon = 60)
add.map.bathy.shade(map, bathy.arctic, refinement = 0)


## Antarctic
bathy = getNOAA.bathy(-180, 180, -90, -45, resolution = 30)
bathy.antarctic = list(Lon = as.numeric(rownames(bathy)),
                       Lat = as.numeric(colnames(bathy)),
                       Z = bathy, res = 15)

save(bathy.antarctic, file = 'data/bathy.antarctic.rdata')

map = make.map(dlon = 60, lon.min = -180, lon.max = 180, lat.min = -90, lat.max = -60, p = make.proj(10, lat = -90,lon = 0))
add.map.bathy(map, bathy.antarctic, bathy.levels = c(-1000, -2500))


## Pacific
bathy = getNOAA.bathy(-70, 100, -60, 60, resolution = 15, antimeridian = TRUE)
bathy.pacific = list(Lon = as.numeric(rownames(bathy)),
                    Lat = as.numeric(colnames(bathy)),
                    Z = bathy, res = 15)

save(bathy.pacific, file = 'data/bathy.pacific.rdata')

map = make.map(coast = "coastlineWorld", p = "+proj=aea +lon_0=-170 +lat_1=0",
               lat.min = -10, lat.max = 54, lon.min = -206, lon.max = -119,
               dlat = 15, dlon = 15)
add.map.bathy(map, bathy.pacific)



## GOM
map = make.map(coast = "coastlineWorldFine", p = "+proj=merc",
               lat.min = 16, lat.max = 32,
               lon.min = -100, lon.max = -75,
               dlat = 5, dlon = 5)

bathy.gom = get.bathy(map, res = 5)
add.map.bathy(map, bathy.gom)
save(bathy.gom, file = 'data/bathy.gom.rdata')
