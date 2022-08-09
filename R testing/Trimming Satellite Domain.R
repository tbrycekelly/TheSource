library(TheSource)
library(ncdf4)

path = 'C:/Users/Tom Kelly/Desktop/T20160252016032.L3m_8D_CHL_chlor_a_4km.nc'

file = nc_open(path)

temp = file$var[[1]]
temp$size
temp$dimids
names(file$dim)

{
  a = Sys.time()
  ncdf4::file = nc_open(path)
  
  lats = c(-40, 0)
  lons = c(60, 100)
  
  if ('lat' %in% names(file$dim)) {
    file.lats = file$dim$lat$vals
    k = which(file.lats <= lats[2] & file.lats >= lats[1])
    #k = k[order(k)]
    
    if (max(diff(k)) == 1) {
      lat.count = k[length(k)] - k[1]
      lat.start = k[1]
    }
  }
  
  if ('lon' %in% names(file$dim)) {
    file.lons = file$dim$lon$vals
    k = which(file.lons <= lons[2] & file.lons >= lons[1])
    #k = k[order(kon
    
    if (max(diff(k)) == 1) {
      lon.count = k[length(k)] - k[1]
      lon.start = k[1]
    }
  }
  
  data = ncdf4::ncvar_get(file, names(file$var)[1], start = c(lon.start, lat.start), count = c(lon.count, lat.count))
  ncdf4::nc_close(file)
  
  message(Sys.time() - a)
}

object.size(data) / 1e6

{
  a = Sys.time()
  
  file = nc_open(path)
  data = ncvar_get(file, names(file$var)[1])
  nc_close(file)
  
  message(Sys.time() - a)
}
object.size(data) / 1e6


a = Sys.time()
test = ncdf4::ncvar_get(file, varid = 'npp', start = c(1,1), count = c(-1,1))
b = Sys.time()
test = ncdf4::ncvar_get(file, varid = 'npp', start = c(1,1), count = c(-1,-1))
c = Sys.time()

file$dim$fakeDim0$

close(file)

sate = read.satellite(path)
nsate = load.nc(path)
