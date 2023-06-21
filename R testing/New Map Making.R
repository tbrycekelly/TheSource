library(microbenchmark)
library(TheSource)

R = 6.371e3


project.mercator = function(lon, lat, lon0 = 0, lat0 = 0, inv = F) {
  lon = (lon - lon0 + 180) %% 360 - 180
  lon = (lon / 180 * 3.14159)
  lat = lat / 180 * 3.14159
  lon0 = lon0 / 180 * 3.14159
  lat0 = lat0 / 180 * 3.14159
  
  x = lon
  y = log(tan(3.14159 * 0.25 + 0.5 * lat))
  #y = 0.5 * log( (1+sin(lat)) / (1 - sin(lat)) ) ## 0.5x
  
  
  list(x = x * R, y = y * R)
}


project.transmercator = function(lon, lat, lon0 = 0, lat0 = 0, inv = F) {
  lon = lon / 180 * 3.14159
  lat = lat / 180 * 3.14159
  lon0 = lon0 / 180 * 3.14159
  lat0 = lat0 / 180 * 3.14159
  
  
  B = cos(lat) * sin(lon - lon0)
  y = atan(tan(lat) / cos(lat - lat0)) - lat0
  x = atanh(B)
  
  
  list(x = x * R, y = y * R)
}


project.orthographic = function(lon, lat, lon0 = 0, lat0 = 0, inv = F) {
  lon = lon / 180 * 3.1415926535
  lat = lat / 180 * 3.1415926535
  lon0 = lon0 / 180 * 3.1415926535
  lat0 = lat0 / 180 * 3.1415926535
  
  x = cos(lat) * sin(lon - lon0)
  y = cos(lat0) * sin(lat) - sin(lat0) * cos(lat) * cos(lon - lon0)
  
  list(x = x * R, y = y * R)
}


project.equalarea = function(lon, lat, lon0 = 0, lat0 = 0, inv = F) {
  lon = lon / 180 * 3.1415926535
  lat = lat / 180 * 3.1415926535
  lon0 = lon0 / 180 * 3.1415926535
  lat0 = lat0 / 180 * 3.1415926535
  
  k = sqrt(2 / (1 + sin(lat0) * sin(lat) + cos(lat0) * cos(lat) * cos(lon - lon0)))
  x = k * cos(lat) * sin(lon - lon0)
  y = k * (cos(lat0) * sin(lat) - sin(lat0) * cos(lat) * cos(lon - lon0))
  
  list(x = x * R, y = y * R)
}


cut.coast = function(coast, lon0) {
  lon0 = lon0 %% 360 # 0 - 360
  lon1 = lon0 + 180
  
  for (i in 1:length(coast$data)) {
    coast$data[[i]]$longitude = (coast$data[[i]]$longitude + lon0) %% 360
    
    l = coast$data[[i]]$longitude >= 180
    l = which(diff(l) != 0)
    
    if (length(l) > 0) {
      l = c(l, nrow(coast$data[[i]]))

      for (k in seq(1, length(l)-1, by = 2)) {
        outside = (l[k]+1):l[k+1]
        
        coast$data[[length(coast$data) + 1]] = coast$data[[i]][outside,]
        coast$data[[i]] = coast$data[[i]][-outside,]
      }
    }
  }
  
  for (i in 1:length(coast$data)) {
    coast$data[[i]]$longitude = coast$data[[i]]$longitude - lon0
  }
  
  coast
}

plot(coast$data[[i]]$longitude[1:366], coast$data[[i]]$latitude[1:366], pch = 20)
points(coast$data[[i]]$longitude[l], coast$data[[i]]$latitude[l], pch = 16, col = 'red')
abline(v = lon0)
abline(v = lon1, col = 'red')

plot(NULL, NULL, xlim = c(-R, R)*pi, ylim = c(-R, R)*pi, axes = F, xlab = '', ylab = '')
#text.default(x = 0, y = 0, labels = '30', srt = 30)

coast = coastline1
coast$data = coast$data[1]
coast = cut.coast(coast, -50)
for (i in 1:length(coast$data)) {
  tmp = project.mercator(coast$data[[i]]$longitude, coast$data[[i]]$latitude, -40, 0)
  polygon(tmp$x, tmp$y, col = make.pal(i, min = 0, max = 10))
}

 coast = coastline2$data[[1]]


microbenchmark({
  tmp = project.mercator(coast$longitude, coast$latitude)
})

microbenchmark({
  tmp = project.orthographic(coast$longitude-120, coast$latitude-30)
})

microbenchmark({
  tmp = project.equalarea(coast$longitude-120, coast$latitude+20)
})

plot(tmp$x, tmp$y, type = 'l')
tmp = project.equalarea(c(-180:180), rep(75, 361))
lines(tmp)

microbenchmark({
  tmp = rgdal::project(cbind(coast$lon, coast$lat), proj = '+proj=merc')
})

microbenchmark({
  tmp = rgdal::project(cbind(coast$lon, coast$lat), proj = '+proj=ortho')
})

plot(tmp[,1], tmp[,2], type = 'l')
