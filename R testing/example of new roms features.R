
## Load data
grid.data = load.roms.grid('Z:/Data/Roms Southern Ocean/ocean_grd.nc')
time.data = load.roms.time('Z:/Data/Roms Southern Ocean/ocean1_000000.nc')
physics.data = load.roms.physics('Z:/Data/Roms Southern Ocean/ocean1_000000.nc')

## Do preliminary calculations (only do once per model)
grid.data$area = calc.roms.areas(grid.data)
grid.data$depth = calc.roms.depths(grid.data) ## Skip free surface estimates for now.
time.data$time = calc.roms.time(time.data)


## Calculate u and v from mass flux and area
physics.data$u = physics.data$Huon / grid.data$area$x
physics.data$v = physics.data$Hvom / grid.data$area$v




#### Example plots:
summary(as.numeric(grid.data$lat_psi))
summary(as.numeric(grid.data$lon_psi))

## Fast way to plot map:
plot.image(z = destagger.grid(physics.data$temp[,,20]), pal = 'ocean.thermal')


## Slow way to plot map (have to draw it as a bunch of points since lat and lon is not actually on a nice grid...)
map = make.map('coastlineWorldFine', lon.min = 235, lon.max = 260,
               lat.min = -75, lat.max = -69, p = make.proj('stere', lat = -70, lon = 250))
add.map.points(grid.data$lon_psi, grid.data$lat_psi,
               pch = 15, cex = 0.16, col = make.pal(pal = 'ocean.thermal', destagger.grid(physics.data$temp[,,20])))
redraw.map(map)




## Plot vertical profiles at particular lon,lat location:
lon = -115
lat = -72
index = calc.roms.location(lon, lat, grid.data)
plot(physics.data$salt[index$i, index$j,],
     grid.data$depth[index$i,index$j,],
     ylab = 'Depth', xlab = 'S')

plot(physics.data$salt[index$i, index$j,],
     exaggerate(-grid.data$depth[index$i,index$j,]),
     yaxt = 'n', ylab = 'Depth', xlab = 'S', ylim = exaggerate(c(500,0)), yaxs = 'i')
add.exaggerated.axis(2)


## Since diffusivity is on a staggered depth axis (21 depths vs the usual 20 for everything else), we can destagger it to match it up with the depths we're using:
plot(destagger.grid.z(physics.data$AKt[index$i, index$j,]),
     exaggerate(-grid.data$depth[index$i,index$j,]),
     yaxt = 'n', ylab = 'Depth', xlab = 'S', ylim = exaggerate(c(500,0)), yaxs = 'i')
add.exaggerated.axis(2)



