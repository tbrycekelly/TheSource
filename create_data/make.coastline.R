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

#### Consider land areas:
area1 = data.frame(index = c(1:length(coastline1$data)), area = NA)
for (i in 1:nrow(area1)) {
  n = nrow(coastline1$data[[i]])
  area = abs(sum(sapply(c(2:n),
                function (x) {
                  coastline1$data[[i]]$longitude[x-1] * coastline1$data[[i]]$latitude[x] - coastline1$data[[i]]$longitude[x] * coastline1$data[[i]]$latitude[x-1]
                })) + coastline1$data[[i]]$longitude[n] * coastline1$data[[i]]$latitude[1] - coastline1$data[[i]]$longitude[1] * coastline1$data[[i]]$latitude[n])

  area1$area[i] = area
}

temp = coastline1
temp$data = temp$data[area1$area > quantile(area1$area, probs = 0.75)]
map = make.map(temp)

object.size(temp)
object.size(coastline1)

# Coastline2

get.coastline.area = function(coastline) {
  area = rep(NA, length(coastline$data))

  for (i in 1:length(area)) {
    n = nrow(coastline$data[[i]])
    area[i] = abs(sum(sapply(c(2:n),
                          function (x) {
                            coastline$data[[i]]$longitude[x-1] * coastline$data[[i]]$latitude[x] - coastline$data[[i]]$longitude[x] * coastline$data[[i]]$latitude[x-1]
                          })) + coastline$data[[i]]$longitude[n] * coastline$data[[i]]$latitude[1] - coastline$data[[i]]$longitude[1] * coastline$data[[i]]$latitude[n])
  }

  area
}

temp = coastline2
area = get.coastline.area(temp)
temp$data = temp$data[area > quantile(area, probs = 0.75)]
map = make.map2(temp, lat = 60, lon = -150)

object.size(temp) / object.size(coastline2)



## Coastline 3
temp = coastline3
area = get.coastline.area(temp)
temp$data = temp$data[area < quantile(area, probs = 0.75)]

map = make.map2(coastline3, lat = 60, lon = -150)
add.map.coastline(temp, p = map$p, land.col = 'red')

object.size(temp) / object.size(coastline3)


## Coastline 4
temp = coastline4
area = get.coastline.area(temp)
temp$data = temp$data[area < quantile(area, probs = 0.7)]

map = make.map2(coastline4, lat = 59, lon = -145, scale = 250)
add.map.coastline(temp, p = map$p, land.col = 'red')

object.size(temp) / object.size(coastline4)





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








