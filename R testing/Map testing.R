library(TheSource)

{
  pdf('R testing/output/map testing - standard maps.pdf')
  map = make.map.arctic('coastline1')
  map = make.map.cce()
  map = make.map.nga()
  map = make.map()
  
  map = make.map.arctic(p = make.proj('nsper', lat = 90, lon = 181, h = 1e8))
  
  map = make.map('coastline4',
                 lon.min = -190,
                 lon.max = -140,
                 lat.min = 50,
                 p = make.proj('stere', lat = 60, lon = 180),
                 dlat = 10,
                 dlon = 20)
  add.map.bathy(map, bathy.arctic, zlim = c(-1e3, 0), pal = 'ocean.haline', trim = T, refine = -1)
  
  map = make.map('coastline1', lon.min = -190, lon.max = -140, lat.min = 50, p = make.proj('stere', lat = 60, lon = 180), dlat = 10, dlon = 20)
  add.map.bathy(map, bathy.global, zlim = c(-1e3, 0), pal = 'ocean.haline', trim = T, refine = -3)
  
  
  map = test.map.arctic()
  map = test.map.california()
  map = test.map.mississippi()
  
  
  map = make.map.nga(dlon = 1, dlat = 1)
  
  
  map = make.map('cosatline2',
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
  
  dev.off()
}


## Palau EEZ
map = make.map2('coastline4', scale = 1e2, lon = 131, lat = 2, dlon = 0.5, dlat = 0.5)
map = make.map2('coastline3', scale = 5e2, lon = 131, lat = 2, dlon = 5, dlat = 5)
map = make.map2('coastline3', scale = 7.5e2, lon = 134, lat = 5, dlon = 5, dlat = 5)

add.map.line(map, pos$Lon, pos$Lat)


for (i in 1:length(eez.data$Palau)) {
  add.map.line(map, eez.data$Palau[[i]]$lon, eez.data$Palau[[i]]$lat, col = 'red')
}


### Micro
map = make.map2('coastline3', scale = 5e2, lon = 134, lat = 5, dlon = 5, dlat = 5)
map = make.map2('coastline3', scale = 1e2, lon = 136.5, lat = 6.5, dlon = 1, dlat = 1)

add.map.line(map, pos$Lon, pos$Lat)

for (i in 1:length(eez.data$Micronesia)) {
  add.map.line(map, eez.data$Micronesia[[i]]$lon, eez.data$Micronesia[[i]]$lat, col = 'blue')
}







map = make.map2('coastline3', scale = 7.5e2, lon = 134, lat = 5, dlon = 5, dlat = 5)
add.map.eez(map, eez.data, 'Palau', col = 'red', lwd = 2)
add.map.line(map, pos$Lon, pos$Lat, lwd = 3)


map = make.map2('cosatline3', scale = 1e3, lon = 120, lat = -14, dlon = 5, dlat = 5)

add.map.eez(map, eez.data, 'Indonesia', col = 'red', lwd = 2)
add.map.eez(map, eez.data, 'Australia', col = 'orange', lwd = 2)
add.map.eez(map, eez.data, 'Highseas', col = 'blue', lwd = 2)




#### Testing various map projections

p = "+proj=aea +lon_0=0 +lat_0=0 +lat_1=30 +lat_2=40 +h=1e+08"
map = make.map2('coastline1', p = p, scale = 3e3)

p = "+proj=aeqd +lon_0=0 +lat_0=0 +lat_1=30 +lat_2=10 +h=1e+08"
map = make.map2('coastline1', p = p, scale = 1e4)

{
  res = data.frame(proj = rgdal::projInfo()$name, dt = NA)
  res$n = c(1:nrow(res))
  
  for (i in c(173, 160, 153, 152, 150, 144, 134, 116, 109, 107, 95, 94, 92, 83, 79, 71, 63, 62, 61, 51, 47, 30, 29, 28, 12)) {
    res = res[-i,]
  }
  
  
  for (i in 1:nrow(res)) {
    a = Sys.time()
    for (j in 1:50) {
      temp = get.map.line(list(p = paste0("+proj=", res$proj[i], " +lon_0=0 +lat_0=0 +lat_1=20 +lat_2=30 +h=1e+08 +sweep=y")), lon = c(-8, 5, -2), lat = c(-3, -3, 8))
    }
    res$dt[i] = (Sys.time() - a) / j
  }
}

res$f = res$dt / min(res$dt)
plot(res$n, res$dt / min(res$dt), pch = 16, ylim = c(0, 10)); grid();

res = res[order(res$dt),]


pdf('R testing/out/Map Projections.pdf')
for (proj in res$proj) {
  
  p = paste0("+proj=", proj, " +lon_0=0 +lat_0=0 +lat_1=30 +lat_2=10 +h=1e+08")
  map = make.map2('coastline1', p = p, scale = 2e3)
  mtext(paste0(proj, ' - ', round(res$f[res$proj == proj], digits = 2)))
  
}
dev.off()





