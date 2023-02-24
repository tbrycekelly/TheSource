library(TheSource)

load('A:/Satellite/index.rdata')

index$File = gsub('Z:/Data/', 'A:/', index$File)

l = which(index$Param == 'CHL.chlor.a.4km' & index$Timeframe == '8D' & index$Year == 2020 & index$Julian.Day == 8 * 40 + 1)

Aqua.Chl.8D = read.satellite(index$File[l[1]])

temp = grid.subsample(z = Aqua.Chl.8D$field)
Aqua.Chl.8D$field = temp$z
Aqua.Chl.8D$lon = Aqua.Chl.8D$lon[seq(1, length(Aqua.Chl.8D$lon), by = 2)]
Aqua.Chl.8D$lat = Aqua.Chl.8D$lat[seq(1, length(Aqua.Chl.8D$lat), by = 2)]
Aqua.Chl.8D$grid = function(satellite) {expand.grid(lon = satellite$lon, lat = satellite$lat)}

save(Aqua.Chl.8D, file = 'data/Aqua.Chl.8D.rdata')




l = which(index$Param == 'ZLEE.Zeu.lee.4km' & index$Timeframe == 'MO' & index$Year == 2020)

Aqua.Zeu.8D = read.satellite(index$File[l[1]])

temp = grid.subsample(z = Aqua.Zeu.8D$field)
Aqua.Zeu.8D$field = temp$z
Aqua.Zeu.8D$lon = Aqua.Zeu.8D$lon[seq(1, length(Aqua.Zeu.8D$lon), by = 2)]
Aqua.Zeu.8D$lat = Aqua.Zeu.8D$lat[seq(1, length(Aqua.Zeu.8D$lat), by = 2)]
Aqua.Zeu.8D$grid = function(satellite) {expand.grid(lon = satellite$lon, lat = satellite$lat)}

save(Aqua.Zeu.8D, file = 'data/Aqua.Zeu.8D.rdata')
