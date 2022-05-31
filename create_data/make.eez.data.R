library(rgdal)
eez = readOGR('inst/extdata/eez_boundaries_v11.shp')
eez@data$SOVEREIGN1[is.na(eez@data$SOVEREIGN1)] = 'Other'
eez@data$SOVEREIGN2[is.na(eez@data$SOVEREIGN2)] = 'Other'

highseas = readOGR('inst/extdata/High_Seas_v1.shp')

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
