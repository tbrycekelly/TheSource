

## Plot all eez lines (slow)
for (j in 1:length(eez)){
  for (i in 1:length(eez[[j]])) {
    add.map.line(map, eez[[j]][[i]]$lon, eez[[j]][[i]]$lat)
  }
}

## Plot all eez lines (fast)
for (j in 1:length(eez)){
  for (i in 1:length(eez[[j]])) {
    add.map.line(map, eez[[j]][[i]]$lon, eez[[j]][[i]]$lat, greatCircle = F)
  }
}


## Plot only Russia/other country EEZ
for (i in 1:length(eez$Russia)) {
  add.map.line(map, eez$Russia[[i]]$lon, eez$Russia[[i]]$lat, col = 'red')
}

## Plot only highseas/other lines
for (i in 1:length(eez$Highseas)) {
  add.map.line(map, eez$Highseas[[i]]$lon, eez$Highseas[[i]]$lat, col = 'blue')
}


get.transect.stations = function(start.lon, start.lat, end.lon, end.lat, n = 10) {
  
  points = data.frame(lon = c(start.lon, end.lon), lat = c(start.lat, end.lat))
  p = make.proj(projection = 'stere', lat = mean(points$lat), lon = mean(points$lon))
  
  points = rgdal::project(as.matrix(points), p)
  points = data.frame(lon = seq(points[1,1], points[2,1], length.out = n),
                      lon = seq(points[1,2], points[2,2], length.out = n))
  points = rgdal::project(as.matrix(points), proj = p, inv = T)
  
  data.frame(lon = points[,1], lat = points[,2])
}

