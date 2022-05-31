library(TheSource)

pdf('R testing/output/map testing - standard maps.pdf')
map = make.map.arctic('coastline1')
map = make.map.cce()
map = make.map.nga()
map = make.map()

map = make.map.arctic(p = make.proj('nsper', lat = 90, lon = 181, h = 1e8))

map = make.map('coastlineWorldFine',
               lon.min = -190,
               lon.max = -140,
               lat.min = 50,
               p = make.proj('stere', lat = 60, lon = 180),
               dlat = 10,
               dlon = 20)
add.map.bathy(map, bathy.arctic, zlim = c(-1e3, 0), pal = 'ocean.haline', trim = T, refine = -1)

map = make.map('coastlineWorld', lon.min = -190, lon.max = -140, lat.min = 50, p = make.proj('stere', lat = 60, lon = 180), dlat = 10, dlon = 20)
add.map.bathy(map, bathy.global, zlim = c(-1e3, 0), pal = 'ocean.haline', trim = T, refine = -3)


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


library(rgdal)
eez = readOGR('C:/Users/tom/Downloads/World_EEZ_v11_20191118/eez_boundaries_v11.shp')


## Palau EEZ
map = make.map2('coastlineWorldMedium', scale = 1e2, lon = 131, lat = 2, dlon = 0.5, dlat = 0.5)
map = make.map2('coastlineWorldMedium', scale = 5e2, lon = 131, lat = 2, dlon = 5, dlat = 5)
map = make.map2('coastlineWorldMedium', scale = 7.5e2, lon = 134, lat = 5, dlon = 5, dlat = 5)

add.map.line(map, pos$Lon, pos$Lat)



for (i in 1:length(eez.data$Palau)) {
  add.map.line(map, eez.data$Palau[[i]]$lon, eez.data$Palau[[i]]$lat, col = 'red')
}

k = which(sqrt((pos$Lon - 131.2)^2 + (pos$Lat - 2.04)^2) < 0.013)

pos[21775,]
pos[281130,]

### Micro
map = make.map2('coastlineWorldMedium', scale = 5e2, lon = 134, lat = 5, dlon = 5, dlat = 5)
map = make.map2('coastlineWorldMedium', scale = 1e2, lon = 136.5, lat = 6.5, dlon = 1, dlat = 1)

add.map.line(map, pos$Lon, pos$Lat)

for (i in 1:length(eez.data$Micronesia)) {
  add.map.line(map, eez.data$Micronesia[[i]]$lon, eez.data$Micronesia[[i]]$lat, col = 'blue')
}

136.45
6.53
k = which(sqrt((pos$Lon - 136.45)^2 + (pos$Lat - 6.53)^2) < 0.025)

pos$Datetime[12105]
pos$Datetime[289530]

### Micro
map = make.map2('coastlineWorldMedium', scale = 5e2, lon = 143, lat = 13, dlon = 5, dlat = 5)
map = make.map2('coastlineWorldMedium', scale = 1e2, lon = 142.3, lat = 11.8, dlon = 1, dlat = 1)

add.map.line(map, pos$Lon, pos$Lat)

for (i in 1:length(eez.data$`United States`)) {
  add.map.line(map, eez.data$`United States`[[i]]$lon, eez.data$`United States`[[i]]$lat, col = 'blue')
}

142.45
11.64
k = which(sqrt((pos$Lon - 142.45)^2 + (pos$Lat - 11.64)^2) < 0.012)

pos$Datetime[2756]



library(rgdal)
eez = readOGR('C:/Users/tom/Downloads/World_EEZ_v11_20191118/eez_boundaries_v11.shp')
eez@data$SOVEREIGN1[is.na(eez@data$SOVEREIGN1)] = 'Other'
eez@data$SOVEREIGN2[is.na(eez@data$SOVEREIGN2)] = 'Other'

highseas = readOGR('C:/Users/tom/Downloads/World_EEZ_v11_20191118/High_Seas_v1.shp')

eez.data = list()
for (name in unique(c(eez@data$SOVEREIGN1, eez@data$SOVEREIGN2))) {
  if (!is.na(name)) {
    l = which((name == eez@data$SOVEREIGN1 | name == eez@data$SOVEREIGN2) & eez@data$LINE_TYPE != 'Straight Baseline')

    eez.data[[name]] = list()
    for (i in 1:length(l)) {
      k = l[i]
      eez.data[[name]][[i]] = data.frame(lon = eez@lines[[k]]@Lines[[1]]@coords[,1],
                                         lat = eez@lines[[k]]@Lines[[1]]@coords[,2])
    }
  }
}

eez.data[['Highseas']] = list()
for (i in 1:length(highseas@polygons[[1]]@Polygons)) {
  for (j in 1:length(highseas@polygons[[1]]@Polygons[[i]])) {
    eez.data[['Highseas']][[length(eez.data[['Highseas']]) + 1]] = data.frame(lon = highseas@polygons[[1]]@Polygons[[i]]@coords[,1],
                                                                            lat = highseas@polygons[[1]]@Polygons[[i]]@coords[,2])
  }
}

save(eez.data, file = 'data/eez.rdata')

add.map.eez = function(map, eez, names = NULL, col = 'black', lwd = 1, lty = 1) {
  if (is.null(names)) {
    names = names(eez)
  }

  for (n in names) {
    l = which(names(eez) == n)
    if (length(l) > 0) {
      for (i in 1:length(eez[[n]])) {
        add.map.line(map, eez[[n]][[i]]$lon, eez[[n]][[i]]$lat, col = col, lwd = lwd, lty = lty)
      }
    }
  }
}

map = make.map2('coastlineWorldMedium', scale = 7.5e2, lon = 134, lat = 5, dlon = 5, dlat = 5)
add.map.eez(map, eez.data, 'Palau', col = 'red', lwd = 2)
add.map.line(map, pos$Lon, pos$Lat, lwd = 3)


map = make.map2('coastlineWorldMedium', scale = 1e3, lon = 120, lat = -14, dlon = 5, dlat = 5)

add.map.eez(map, eez.data, 'Indonesia', col = 'red', lwd = 2)
add.map.eez(map, eez.data, 'Australia', col = 'orange', lwd = 2)
add.map.eez(map, eez.data, 'Highseas', col = 'blue', lwd = 2)

for (i in 1:length(eez@lines)) {
  for (j in 1:length(eez@lines[[i]])) {
    add.map.line(map, eez@lines[[i]]@Lines[[j]]@coords[,1], eez@lines[[i]]@Lines[[j]]@coords[,2], greatCircle = F, col = 'red')
  }
}

eez.data[['Highseas']] = list()
for (i in 1:length(highseas@polygons[[1]]@Polygons)) {
  for (j in 1:length(highseas@polygons[[1]]@Polygons[[i]])) {
    eez.data[['Highseas']][length(eez.data[['Highseas']]) + 1] = data.frame(lon = highseas@polygons[[1]]@Polygons[[i]]@coords[,1],
                                                                            lat = highseas@polygons[[1]]@Polygons[[i]]@coords[,2])
    add.map.line(map,
                 highseas@polygons[[1]]@Polygons[[i]]@coords[,1],
                 highseas@polygons[[1]]@Polygons[[i]]@coords[,2],
                 greatCircle = F,
                 lwd = 2,
                 col = 'red')
  }
}



add.map.line(map, pos$Lon, pos$Lat, lwd = 1)
