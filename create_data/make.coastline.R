library(TheSource)

build.coastline = function(file) {

  coast = data.frame(longitude = NA, latitude = NA)
  empty = coast

  for (k in 1:length(file)) {
    data = rgdal::readOGR(file[k])
    message('Finished loading.')

    for (i in 1:length(data@polygons)) {
      if (i %% 1e3 == 0) {message(i)}
      for (j in 1:length(data@polygons[[i]]@Polygons)) {
        new = data.frame(longitude = data@polygons[[i]]@Polygons[[j]]@coords[,1],
                         latitude = data@polygons[[i]]@Polygons[[j]]@coords[,2])

        coast = rbind(coast, new) ##include NA row between each
      }
      coast = rbind(coast, empty)
    }
  }

  coastlineWorld@data = coast
  coastlineWorld
}


build.coastline = function(file) {

  coast = list()

  for (k in 1:length(file)) {
    data = rgdal::readOGR(file[k])
    message('Finished loading.')

    for (i in 1:length(data@polygons)) {
      if (i %% 1e3 == 0) {message(i)}
      for (j in 1:length(data@polygons[[i]]@Polygons)) {
        if (length(data@polygons[[i]]@Polygons[[j]]@coords[,1]) > 6) {
          coast[[length(coast) + 1]] = data.frame(longitude = data@polygons[[i]]@Polygons[[j]]@coords[,1],
                           latitude = data@polygons[[i]]@Polygons[[j]]@coords[,2])
        }

      }
    }
  }

  ## Return
  list(data = coast,
       meta = list(
         citation = NA,
         source.files = file,
         url = NA,
         time.creation = Sys.time(),
         size = NA,
         R.version = R.version.string,
         TheSource.version = TheSource.version()
       )
  )
}


coastline1 = build.coastline(list.files('C:/Users/Tom Kelly/Downloads/GSHHS_shp/c/', full.names = T, pattern = '.shp'))
coastline2 = build.coastline(list.files('C:/Users/Tom Kelly/Downloads/GSHHS_shp/l/', full.names = T, pattern = '.shp'))
coastline3 = build.coastline(list.files('C:/Users/Tom Kelly/Downloads/GSHHS_shp/i/', full.names = T, pattern = '.shp'))
coastline4 = build.coastline(list.files('C:/Users/Tom Kelly/Downloads/GSHHS_shp/h/', full.names = T, pattern = '.shp'))
coastline5 = build.coastline(list.files('C:/Users/Tom Kelly/Downloads/GSHHS_shp/f/', full.names = T, pattern = '.shp'))

save(coastline1, file = 'data/coastline1.rdata')
save(coastline2, file = 'data/coastline2.rdata')
save(coastline3, file = 'data/coastline3.rdata')
save(coastline4, file = 'data/coastline4.rdata')
save(coastline5, file = 'data/coastline5.rdata')


data('coastlineWorld', package = 'oce')
data('coastlineWorldMedium', package = 'ocedata')
data('coastlineWorldFine', package = 'ocedata')

temp = list(data = list(),
            meta = list(
              citation = NA,
              source.files = 'oce pacakge by Dan Kelley',
              url = NA,
              time.creation = Sys.time(),
              size = NA,
              R.version = R.version.string,
              TheSource.version = TheSource.version()
            )
          )

l = c(which(is.na(coastlineWorld@data$longitude)), length(coastlineWorld@data$longitude))
for (i in 1:(length(l)-1)) {
  temp$data[[i]] = data.frame(longitude = coastlineWorld@data$longitude[(l[i]+1):(l[i+1]-1)],
                              latitude = coastlineWorld@data$latitude[(l[i]+1):(l[i+1]-1)])
}
coastlineWorld = temp
save(coastlineWorld, file = 'data/coastlineWorld.rdata')

#### CoastlineWorldMedium
temp = list(data = list(),
            meta = list(
              citation = NA,
              source.files = 'oce pacakge by Dan Kelley',
              url = NA,
              time.creation = Sys.time(),
              size = NA,
              R.version = R.version.string,
              TheSource.version = TheSource.version()
            )
)

l = c(which(is.na(coastlineWorldMedium@data$longitude)), length(coastlineWorldMedium@data$longitude))
for (i in 1:(length(l)-1)) {
  temp$data[[i]] = data.frame(longitude = coastlineWorldMedium@data$longitude[(l[i]+1):(l[i+1]-1)],
                              latitude = coastlineWorldMedium@data$latitude[(l[i]+1):(l[i+1]-1)])
}
coastlineWorldMedium = temp
save(coastlineWorldMedium, file = 'data/coastlineWorldMedium.rdata')


#### CoastlineWorldFine
temp = list(data = list(),
            meta = list(
              citation = NA,
              source.files = 'oce pacakge by Dan Kelley',
              url = NA,
              time.creation = Sys.time(),
              size = NA,
              R.version = R.version.string,
              TheSource.version = TheSource.version()
            )
)

l = c(which(is.na(coastlineWorldFine@data$longitude)), length(coastlineWorldFine@data$longitude))
for (i in 1:(length(l)-1)) {
  temp$data[[i]] = data.frame(longitude = coastlineWorldFine@data$longitude[(l[i]+1):(l[i+1]-1)],
                              latitude = coastlineWorldFine@data$latitude[(l[i]+1):(l[i+1]-1)])
}
coastlineWorldFine = temp
save(coastlineWorldFine, file = 'data/coastlineWorldFine.rdata')






p = make.proj('stere', lon = -150, lat = 60)
map = make.map.nga('coastlineAlaska', lon.min = -156.4, lon.max = -157, lat.min = 71.1, lat.max = 71.4, dlon = 0.1, dlat = 0.1)
map = make.map.nga('coastlineWorldFine', lon.min = -151, lon.max = -150, lat.min = 59, lat.max = 60)
map = make.map2('coastlineL', lon = -150, lat = 60, scale = 2000, land.col = 'black', p = p)
map = make.map2('coastlineWorldFine', lon = -150, lat = 60, scale = 2000, land.col = 'black', p = p)




coastline = list(
  data = list(),
  meta = list(
    citation = NA,
    url = NA,
    time.creation = Sys.time(),
    size = NA,
    R.version = R.version.string,
    TheSource.version = package_version('TheSource')
  )
)








