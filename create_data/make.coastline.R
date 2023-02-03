library(TheSource)

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


## Start the work!
file.url = 'https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/gshhg-shp-2.3.7.zip'
temp.file = tempfile()
temp.dir = tempdir()

## Download and extract coastline data
options(timeout = 600)
download.file(file.url, destfile = temp.file)
archive::archive_extract(temp.file, temp.dir)

coastline1 = build.coastline(list.files(paste0(temp.dir, '/GSHHS_shp/c/'), full.names = T, pattern = 'L1.shp'))
coastline2 = build.coastline(list.files(paste0(temp.dir, '/GSHHS_shp/l/'), full.names = T, pattern = 'L1.shp'))
coastline3 = build.coastline(list.files(paste0(temp.dir, '/GSHHS_shp/i/'), full.names = T, pattern = 'L1.shp'))
coastline4 = build.coastline(list.files(paste0(temp.dir, '/GSHHS_shp/h/'), full.names = T, pattern = 'L1.shp'))
coastline5 = build.coastline(list.files(paste0(temp.dir, '/GSHHS_shp/f/'), full.names = T, pattern = 'L1.shp'))

lakes1 = build.coastline(list.files(paste0(temp.dir, '/GSHHS_shp/c/'), full.names = T, pattern = 'L2.shp'))
lakes2 = build.coastline(list.files(paste0(temp.dir, '/GSHHS_shp/l/'), full.names = T, pattern = 'L2.shp'))
lakes3 = build.coastline(list.files(paste0(temp.dir, '/GSHHS_shp/i/'), full.names = T, pattern = 'L2.shp'))
lakes4 = build.coastline(list.files(paste0(temp.dir, '/GSHHS_shp/h/'), full.names = T, pattern = 'L2.shp'))
lakes5 = build.coastline(list.files(paste0(temp.dir, '/GSHHS_shp/f/'), full.names = T, pattern = 'L2.shp'))

antarctica1 = build.coastline(list.files(paste0(temp.dir, '/GSHHS_shp/c/'), full.names = T, pattern = 'L6.shp'))
antarctica2 = build.coastline(list.files(paste0(temp.dir, '/GSHHS_shp/l/'), full.names = T, pattern = 'L6.shp'))
antarctica3 = build.coastline(list.files(paste0(temp.dir, '/GSHHS_shp/i/'), full.names = T, pattern = 'L6.shp'))
antarctica4 = build.coastline(list.files(paste0(temp.dir, '/GSHHS_shp/h/'), full.names = T, pattern = 'L6.shp'))
antarctica5 = build.coastline(list.files(paste0(temp.dir, '/GSHHS_shp/f/'), full.names = T, pattern = 'L6.shp'))


## Coastline1
temp = coastline1
area = get.coastline.area(temp)
temp$data = temp$data[area > quantile(area, probs = 0.75)]
coastline1 = temp

# Coastline2

temp = coastline2
area = get.coastline.area(temp)
temp$data = temp$data[area > quantile(area, probs = 0.75)]
coastline2 = temp


## Coastline 3
temp = coastline3
area = get.coastline.area(temp)
temp$data = temp$data[area > quantile(area, probs = 0.75)]
coastline3 = temp


## Coastline 4
temp = coastline4
area = get.coastline.area(temp)
temp$data = temp$data[area > quantile(area, probs = 0.7)]
coastline4 = temp

## Coastline 5
temp = coastline5
area = get.coastline.area(temp)
temp$data = temp$data[area > quantile(area, probs = 0.7)]
coastline5 = temp


## Trim lakes
for (i in c(lakes1, lakes2, lakes3, lakes4, lakes5)) {
  area = get.coastline.area(i)
  i$data = i$data[area > quantile(area, probs = 0.5)]
}


## Add metadata
coastline1$meta$url = file.url
coastline2$meta$url = file.url
coastline3$meta$url = file.url
coastline4$meta$url = file.url
coastline5$meta$url = file.url
lakes1$meta$url = file.url
lakes2$meta$url = file.url
lakes3$meta$url = file.url
lakes4$meta$url = file.url
lakes5$meta$url = file.url
antarctica1$meta$url = file.url
antarctica2$meta$url = file.url
antarctica3$meta$url = file.url
antarctica4$meta$url = file.url
antarctica5$meta$url = file.url

## Add sizes
coastline1$meta$size = object.size(coastline1)
coastline2$meta$size = object.size(coastline2)
coastline3$meta$size = object.size(coastline3)
coastline4$meta$size = object.size(coastline4)
coastline5$meta$size = object.size(coastline5)
lakes1$meta$size = object.size(lakes1)
lakes2$meta$size = object.size(lakes2)
lakes3$meta$size = object.size(lakes3)
lakes4$meta$size = object.size(lakes4)
lakes5$meta$size = object.size(lakes5)
antarctica1$meta$size = object.size(antarctica1)
antarctica2$meta$size = object.size(antarctica2)
antarctica3$meta$size = object.size(antarctica3)
antarctica4$meta$size = object.size(antarctica4)
antarctica5$meta$size = object.size(antarctica5)

## Save
save(coastline1, file = 'data/coastline1.rdata')
save(coastline2, file = 'data/coastline2.rdata')
save(coastline3, file = 'data/coastline3.rdata')
save(coastline4, file = 'data/coastline4.rdata')
save(coastline5, file = 'data/coastline5.rdata')
save(lakes1, file = 'data/lakes1.rdata')
save(lakes2, file = 'data/lakes2.rdata')
save(lakes3, file = 'data/lakes3.rdata')
save(lakes4, file = 'data/lakes4.rdata')
save(lakes5, file = 'data/lakes5.rdata')
save(antarctica1, file = 'data/antarctica1.rdata')
save(antarctica2, file = 'data/antarctica2.rdata')
save(antarctica3, file = 'data/antarctica3.rdata')
save(antarctica4, file = 'data/antarctica4.rdata')
save(antarctica5, file = 'data/antarctica5.rdata')


map = make.map2(coastline1, scale = 200, lat = 5)
map = make.map2(coastline2, scale = 200, lat = 5)
map = make.map2(coastline3, scale = 200, lat = 5)
map = make.map2(coastline4, scale = 200, lat = 5)
map = make.map2(coastline5, scale = 200, lat = 5)


map = make.map2(coastline1, scale = 2000, lat = 60, lon = -150)
map = make.map2(coastline2, scale = 2000, lat = 60, lon = -150)
map = make.map2(coastline3, scale = 2000, lat = 60, lon = -150)
temp = add.map.coastline(lakes3, p = map$p, land.col = 'white')

map = make.map2(coastline1, scale = 200, lat = 60, lon = -150)
map = make.map2(coastline2, scale = 200, lat = 60, lon = -150)
map = make.map2(coastline3, scale = 200, lat = 60, lon = -150)
map = make.map2(coastline4, scale = 200, lat = 60, lon = -150)
map = make.map2(coastline5, scale = 200, lat = 60, lon = -150)


map = make.map2(coastline4, scale = 200, lat = 60, lon = -150)
temp = add.map.coastline(lakes2, p = map$p, land.col = 'white')


map = make.map2(lakes1, scale = 200, lat = 60, lon = -150)
map = make.map2(lakes2, scale = 200, lat = 60, lon = -150)
map = make.map2(lakes3, scale = 200, lat = 60, lon = -150)
map = make.map2(lakes4, scale = 200, lat = 60, lon = -150)
map = make.map2(lakes5, scale = 200, lat = 60, lon = -150)

