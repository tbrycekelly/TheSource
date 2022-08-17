library(tidync)


file = tidync::tidync(x = 'Z:/Data/Satellite/Raw/2012/V/V2012025.L4_8D_CBPM_4km.nc')

## Get dimension information
ndim = length(file$transforms)
dims = list()
for (i in 1:ndim) {
  dims[[names(file$transforms)[i]]] = file$transforms[[i]][[1]]
  
  file$transforms[[i]]$selected = F
  file$transforms[[i]]$selected[1:2] = T
  
  
}


tidync::hyper_filter
temp = as.data.frame(tidync::hyper_tibble(file))


temp = tidync::hyper_array(file, fakeDim0 = fakeDim0 > 80, fakeDim1 = fakeDim1 < 40)
tidync::