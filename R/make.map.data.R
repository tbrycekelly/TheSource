# #map = make.map(coast = "coastlineWorldFine", p = "+proj=merc",
# #               lat.min = 27, lat.max = 31,
# #               lon.min = -93, lon.max = -84,
# #               dlat = 5, dlon = 5)
# #add.map.bathy.shade(map, bathy.global)
#
# ## Global
# bathy = getNOAA.bathy(-180, 180, -90, 90, resolution = 15)
# bathy.global = list(Lon = as.numeric(rownames(bathy)),
#                     Lat = as.numeric(colnames(bathy)),
#                     Z = bathy, res = 15)
#
# save(bathy.global, file = 'data/bathy.global.rdata')
#
#
# ## Arctic
# bathy = getNOAA.bathy(-180, 180, 60, 90, resolution = 10)
# bathy.arctic = list(Lon = as.numeric(rownames(bathy)),
#                     Lat = as.numeric(colnames(bathy)),
#                     Z = bathy, res = 15)
#
# save(bathy.global, file = 'data/bathy.arctic.rdata')
#
# map = make.map(dlon = 60)
# add.map.bathy.shade(map, bathy.arctic)
#
# ## Antarctic
# bathy = getNOAA.bathy(-180, 180, -90, 60, resolution = 15)
# bathy.antarctic = list(Lon = as.numeric(rownames(bathy)),
#                     Lat = as.numeric(colnames(bathy)),
#                     Z = bathy, res = 15)
#
# save(bathy.antarctic, file = 'data/bathy.antarctic.rdata')
#
# map = make.map(dlon = 60, lat.min = -50, lat.max = -80, )
# add.map.bathy.shade(map, bathy.antarctic)
#
#
# ## Pacific
# bathy = getNOAA.bathy(-70, 100, -60, 60, resolution = 15, antimeridian = TRUE)
# bathy.pacific = list(Lon = as.numeric(rownames(bathy)),
#                     Lat = as.numeric(colnames(bathy)),
#                     Z = bathy, res = 15)
#
# save(bathy.pacific, file = 'data/bathy.pacific.rdata')
#
# map = make.map(coast = "coastlineWorld", p = "+proj=aea +lon_0=-170 +lat_1=0",
#                lat.min = -10, lat.max = 54, lon.min = -206, lon.max = -119,
#                dlat = 15, dlon = 15)
# add.map.bathy.shade(map, bathy.pacific, zlim = c(-6000, -100), filled = FALSE)
#
#
#
# ## GOM
# map = make.map(coast = "coastlineWorldFine", p = "+proj=merc",
#                lat.min = 16, lat.max = 32,
#                lon.min = -100, lon.max = -75,
#                dlat = 5, dlon = 5)
# bathy.gom = get.bathy(map, res = 5)
# add.map.bathy.shade(map, bathy.gom)
# save(bathy.gom, file = 'data/bathy.gom.rdata')
