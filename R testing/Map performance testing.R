library(TheSource)

map = make.map.arctic(coast = 'coastline2', dlon = 10, dlat = 10)
add.map.scale()
add.map.axis(map = map, N = 20, sides = 4)

map = make.map(coast = 'coastline2')
add.map.coastline(lakes2, p = map$p, land.col = 'blue')

map = make.map.nga()
add.map.coastline(lakes3, p = map$p, land.col = 'blue')



map = make.map2('coastline3', scale = 2000, lon = -100, lat = 50, p = make.proj(3))
temp = add.map.coastline(lakes2, p = map$p, land.col = 'blue')




a = c(1:4)


{
  a = Sys.time()
  v = 1 %in% a
  b = Sys.time() - a
}; b * 1000


{
  a = Sys.time()
  v = any(1 == a)
  b = Sys.time() - a
}; b * 1000


{
  a = Sys.time()
  v = sum(1 == a)
  b = Sys.time() - a
}; b * 1000




for (i in 1:100) {
  lon = runif(1, -180, 180)
  lat = runif(1, -90, 90)
  p = make.proj(round(runif(1, 1,9)), lon = lon, lat= lat)
  scale = runif(1, 1000, 5000)
  
  map = make.map2(coastline3, scale = scale, p = p, lon = lon, lat = lat)
  
  wait(10, F)
}

















