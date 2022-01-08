library(TheSource)

plot.new()
map = make.map.arctic(dlon = 20, dlat = 10)
map = make.map.cce()
map = make.map.nga()
map = make.map()

map = make.map.arctic(p = make.proj('stere', lat = 75, lon = 181, h = 1e8))

map = make.map('coastlineWorldMedium', lon.min = -190, lon.max = -140, lat.min = 50, p = make.proj('stere', lat = 60, lon = 180), dlat = 10, dlon = 20)
add.map.bathy(map, bathy.arctic, zlim = c(-1e3, 0), pal = 'ocean.haline', trim = T)
add.map.bathy(map, bathy.global, zlim = c(-1e3, 0), pal = 'ocean.haline', trim = T, subsample = 5)


map = test.map.arctic()
map = test.map.california()
map = test.map.mississippi()


map = make.map.nga(dlon = 1, dlat = 1)


map = make.map('coastlineWorldMedium',
               lon.min = 100,
               lon.max = 150,
               lat.min = -25,
               lat.max = 25,
               dlon = 10,
               dlat = 10,
               land.col = '#444464',
               p = make.proj(10, lat = 0, lon = 125))
add.map.bathy(map, bathy.global)

map = make.map(land.col = 'black', lat.max = 75, draw.grid = F, dlon = 60)
add.map.bathy(map, bathy.global, subsample = 4)



