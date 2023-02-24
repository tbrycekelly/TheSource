library(TheSource)



{
  png('tests/Mapping Test 1.png', width = 2e3, height = 2e3)
  map = make.map2(scale = 3e3, lon = 0)
  add.map.layer(map, Aqua.Chl.8D$lon, Aqua.Chl.8D$lat, Aqua.Chl.8D$field, pal = 'inferno', zlim = c(0,5))
  
  dev.off()
}


{
  png('tests/Mapping Test 2.png', width = 2e3, height = 2e3)
  map = make.map2(scale = 1e3, lon = 180)
  add.map.layer(map, Aqua.Chl.8D$lon, Aqua.Chl.8D$lat, Aqua.Chl.8D$field, pal = 'inferno', zlim = c(0,5))
  
  dev.off()
}

{
  png('tests/Mapping Test 3.png', width = 2e3, height = 2e3)
  map = make.map.arctic()
  add.map.layer(map, Aqua.Zeu.8D$lon, Aqua.Zeu.8D$lat, Aqua.Zeu.8D$field, pal = 'cubicl', zlim = c(0,100), refine = -2)
  
  dev.off()
}


{
  png('tests/Mapping Test 4.png', width = 2e3, height = 2e3)
  map = make.map(lon.min = 90, lon.max = 360-90)
  add.map.layer(map, chl$lon, chl$lat, chl$field[,,1], pal = 'inferno', zlim = c(0,5), refine = -4)
  
  dev.off()
}


{
  png('tests/Mapping Test 5.png', width = 2e3, height = 2e3)
  map = make.map2(coast = 'antarctica3', lat = -90, p = make.proj('nsper', lat = -90), scale = 3e3)
  add.map.layer(map, Aqua.Zeu.8D$lon, Aqua.Zeu.8D$lat, Aqua.Zeu.8D$field, pal = 'cubicl', zlim = c(0,100), trim = F)
  
  dev.off()
}


{
  png('tests/Mapping Test 6.png', width = 2e3, height = 2e3)
  map = make.map2(coast = 'coastline3', lat = 35, lon = -119, p = make.proj('ortho', lat = 30, lon = -115), scale = 3e2)
  add.map.layer(map, Aqua.Zeu.8D$lon, Aqua.Zeu.8D$lat, Aqua.Zeu.8D$field, pal = 'cubicl', zlim = c(0,100))
  redraw.map(map)
  dev.off()
}
