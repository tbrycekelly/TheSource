library(TheSource)
p = make.proj()

p = "+proj=adams_hemi +lon_0=0 +lat_0=30 +lat_1=0 +lat_2=10 +h=1e+08"
map = make.map2('coastline1', p = p, scale = 1e4)

p = "+proj=aea +lon_0=0 +lat_0=0 +lat_1=30 +lat_2=10 +h=1e+08"
map = make.map2('coastline1', p = p, scale = 1e4)

p = "+proj=aeqd +lon_0=0 +lat_0=0 +lat_1=30 +lat_2=10 +h=1e+08"
map = make.map2('coastline1', p = p, scale = 1e4)

p = "+proj=vandg2 +lon_0=0 +lat_0=80 +lat_1=30 +lat_2=10 +h=1e+08"
map = make.map2('coastline2', lat = 80, p = p, scale = 1e4)

{
  res = data.frame(proj = rgdal::projInfo()$name, dt = NA)
  res$n = c(1:nrow(res))
  
  for (i in c(173, 160, 153, 152, 150, 144, 134, 116, 109, 107, 95, 94, 92, 83, 79, 71, 63, 62, 61, 51, 47, 30, 29, 28, 12)) {
    res = res[-i,]
  }
  
  
  for (i in 1:nrow(res)) {
    a = Sys.time()
    for (j in 1:50) {
      temp = get.map.line(list(p = paste0("+proj=", res$proj[i], " +lon_0=0 +lat_0=0 +lat_1=20 +lat_2=30 +h=1e+08 +sweep=y")), lon = c(-8, 5, -2), lat = c(-3, -3, 8))
    }
    res$dt[i] = (Sys.time() - a) / j
  }
}
res$f = res$dt / min(res$dt)
plot(res$n, res$dt / min(res$dt), pch = 16, ylim = c(0, 10)); grid();

head(res[order(res$dt),], 30)

