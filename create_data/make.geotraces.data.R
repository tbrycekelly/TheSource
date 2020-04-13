library(TheSource)
library(openxlsx)


{## GA01.Sal
  Geotraces.GA01.Sal = read.xlsx('../TheSource.misc/GEOTRACES_IDP2017_v1/GA01.Salinity.xlsx')
  save(Geotraces.GA01.Sal, file = 'data/Geotraces.GA01.Sal.rdata')
}
