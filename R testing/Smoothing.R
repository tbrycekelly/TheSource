


grid = array(0, dim = c(15,16))

## Add Corners
grid[1,1] = 114
grid[1,16] = 121
grid[15,1] = 116
grid[15,16] = 121

delta.j1 = (grid[15,1] - grid[1,1]) / 14
delta.j6 = (grid[15,16] - grid[1,16]) / 14

## Fill in first and last column
for (i in c(1:15)) {
  grid[i,1] = grid[1,1] + delta.j1 * (i-1)
  grid[i,16] = grid[1,16] + delta.j6 * (i-1)
}


## Fill in all values
for (i in c(1:15)) {

  delta = (grid[i,16] - grid[i,1]) / 15

  for (j in 2:15) {
    grid[i,j] = grid[i,j-1] + delta
  }
}


smooth.gaussian = function(x, sd = 0.2, n = 5) {

  stencil = array(1, dim = rep(n, 2))

  for (i in 1:n) {
    stencil[i,] = stencil[i,] * dnorm(c(1:n), (n+1)/2, sd)
  }
  for (i in 1:n) {
    stencil[,i] = stencil[,i] * dnorm(c(1:n), (n+1)/2, sd)
  }
  stencil = stencil / sum(stencil, na.rm = T)


  new = x
  buffer = (n + 1) / 2
  for (i in buffer:(dim(x)[1] - buffer)) {
    for (j in buffer:(dim(x)[2] - buffer)) {
      new[i,j] = sum(stencil * x[(i-buffer+1):(i+buffer-1), (j-buffer+1):(j+buffer-1)])
    }
  }

  new
}


original = array(cumsum(runif(400, -1, 1)), dim = c(25,16))
image(original)
test = smooth.gaussian(original, sd = 5, n = 13)
image(test)


poc = read.satellite('C:/Data/Satellite/A20022492002256.L3m_8D_POC_poc_4km.nc', lon = c(-180, -140), lat = c(40, 80))
str(poc$field)
poc.smooth = smooth.gaussian(poc$field[[1]], 0.2, 3)

map = make.map.nga()
add.map.layer(map, poc$lon, poc$lat, poc$field[[1]], pal = 'parula', zlim = c(0, 1e3))
add.map.layer(map, poc$lon, poc$lat, poc.smooth, pal = 'parula', zlim = c(0, 1e3))
redraw.map(map)



frrf = load.frrf('C:/Users/Tom Kelly/Desktop/FRRF/CTG-Act2Run/Auto-saved FLC data files/20220308/')

plot(frrf[[1]]$A$E,
     frrf[[1]]$A$JPII,
     type = 'l',
     lwd = 3,
     yaxs = 'i',
     xaxs = 'i',
     xlim = c(0,700),
     xlab = 'E')

lines(frrf[[2]]$A$E,
      frrf[[2]]$A$JPII,
      lwd = 3,
      col = 'dark green')

lines(frrf[[3]]$A$E,
      frrf[[3]]$A$JPII,
      lwd = 3,
      col = 'dark green')
