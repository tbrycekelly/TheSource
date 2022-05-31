library(TheSource)
library(microbenchmark)

##### Gridder engine test

## gridding random noise
pdf('R testing/grid tests - runif.pdf', pointsize = 11)
par(mfrow = c(2,2))
for (n in c(10,100, 1e3)) {
  x = runif(n, 0, 10)
  y = runif(n, 0, 10)
  z = runif(n)
  z[sample(1:length(z), size = length(z) / 2)] = NA

  for (gridder in c('gridNN', 'gridNNI', 'gridODV', 'gridIDW')) {
    section = build.section(x, y, z, gridder = eval(parse(text = gridder))) # defaults
    plot.section(section, zlim = c(0,1), pal = 'inferno')
    mtext(adj = 0, gridder, cex = 0.8)
  }
}
dev.off()



pdf('R testing/grid tests - rnorm.pdf', pointsize = 11)
par(mfrow = c(2,2))
for (n in c(5, 25, 50)) {
  x = runif(n, 0, 10)
  y = runif(n, 0, 10)
  z = rnorm(n)
  z[sample(1:length(z), size = length(z) / 2)] = NA

  for (gridder in c('gridNN', 'gridNNI', 'gridODV', 'gridIDW')) {
    section = build.section(x, y, z, gridder = eval(parse(text = gridder))) # defaults
    plot.section(section, zlim = c(-1,1), pal = 'ocean.balance')
    mtext(adj = 0, gridder, cex = 0.8)
  }
}
dev.off()

pdf('R testing/grid tests - grid runif.pdf', pointsize = 11)
par(mfrow = c(2,2))
for (n in c(5, 25, 50)) {
  grid = expand.grid(x = c(1:n), y = c(1:n))
  z = runif(n^2)
  z[sample(1:length(z), size = length(z) / 2)] = NA

  for (gridder in c('gridNN', 'gridNNI', 'gridODV', 'gridIDW')) {
    section = build.section(grid$x, grid$y, z, gridder = eval(parse(text = gridder))) # defaults
    plot.section(section, zlim = c(0,1), pal = 'inferno')
    mtext(adj = 0, gridder, cex = 0.8)
  }
}
dev.off()


pdf('R testing/grid tests - grid rnorm.pdf', pointsize = 11)
par(mfrow = c(2,2))
for (n in c(5, 25, 50)) {
  grid = expand.grid(x = c(1:n), y = c(1:n))
  z = rnorm(n^2)
  z[sample(1:length(z), size = length(z) / 2)] = NA

  for (gridder in c('gridNN', 'gridNNI', 'gridODV', 'gridIDW')) {
    section = build.section(grid$x, grid$y, z, gridder = eval(parse(text = gridder))) # defaults
    plot.section(section, zlim = c(-1,1), pal = 'ocean.balance')
    mtext(adj = 0, gridder, cex = 0.8)
  }
}
dev.off()



pdf('R testing/grid tests - grid rnorm scaled with uncertainty.pdf', pointsize = 11)
par(mfrow = c(2,2))
for (n in c(5, 25, 50)) {
  grid = expand.grid(x = c(1:n), y = c(1:n))
  z = rnorm(n^2)
  z[sample(1:length(z), size = length(z) / 2)] = NA

  for (gridder in c('gridNN', 'gridNNI', 'gridODV', 'gridIDW')) {
    section = build.section(grid$x*100,
                            grid$y*100,
                            z,
                            gridder = eval(parse(text = gridder)),
                            nx = 100,
                            ny = 100,
                            uncertainty = 1) # defaults
    plot.section(section, zlim = c(-1,1), pal = 'ocean.balance')
    mtext(adj = 0, gridder, cex = 0.8)
  }
}
dev.off()

pdf('R testing/grid tests - grid rnorm scaled without uncertainty.pdf', pointsize = 11)
par(mfrow = c(2,2))
for (n in c(5, 25, 50)) {
  grid = expand.grid(x = c(1:n), y = c(1:n))
  z = rnorm(n^2)
  z[sample(1:length(z), size = length(z) / 2)] = NA

  for (gridder in c('gridNN', 'gridNNI', 'gridODV', 'gridIDW')) {
    section = build.section(grid$x*100,
                            grid$y*100,
                            z,
                            gridder = eval(parse(text = gridder)),
                            nx = 100,
                            ny = 100,
                            uncertainty = 1e-5) # defaults
    plot.section(section, zlim = c(-1,1), pal = 'ocean.balance')
    mtext(adj = 0, gridder, cex = 0.8)
  }
}
dev.off()



pdf('R testing/grid tests - x and y scale.pdf', pointsize = 11)
par(mfrow = c(2,2))
for (n in c(5, 25, 50)) {
  grid = expand.grid(x = c(1:n), y = c(1:n))
  z = rnorm(n^2)
  z[sample(1:length(z), size = length(z) / 2)] = NA

  for (gridder in c('gridNN', 'gridNNI', 'gridODV', 'gridIDW')) {
    section = build.section(grid$x,
                            grid$y*10,
                            z,
                            gridder = eval(parse(text = gridder)),
                            nx = 100,
                            ny = 100,
                            uncertainty = 1e-5,
                            x.factor = n)
    plot.section(section, zlim = c(-1,1), pal = 'ocean.balance')
    mtext(adj = 0, gridder, cex = 0.8)
  }
}
dev.off()



##############################################################
############## Neighborhood Tests ############################
##############################################################


pdf('R testing/grid tests gridODV - influence of neighbors on ODV.pdf', pointsize = 11)
par(mfrow = c(2,2))

grid = expand.grid(x = c(1:10), y = c(1:10))
z = rnorm(10^2)
z[sample(1:length(z), size = length(z) / 2)] = NA

for (n in c(2, 5, 10, 20)) {

  section = build.section(grid$x,
                            grid$y,
                            z,
                            gridder = gridODV,
                            nx = 100,
                            ny = 100,
                            uncertainty = 1,
                          neighborhood = n)
    plot.section(section, zlim = c(-1,1), pal = 'ocean.balance')
    mtext(adj = 0, paste0('gridODV w/ neighborhood = ', n), cex = 0.8)
}
dev.off()



pdf('R testing/grid tests gridIDW - influence of neighbors on IDW.pdf', pointsize = 11)
par(mfrow = c(2,2))

grid = expand.grid(x = c(1:10), y = c(1:10))
z = rnorm(10^2)
z[sample(1:length(z), size = length(z) / 2)] = NA

for (n in c(2, 5, 10, 20)) {

  section = build.section(grid$x,
                          grid$y,
                          z,
                          gridder = gridIDW,
                          nx = 100,
                          ny = 100,
                          uncertainty = 1,
                          neighborhood = n)
  plot.section(section, zlim = c(-1,1), pal = 'ocean.balance')
  mtext(adj = 0, paste0('gridIDW w/ neighborhood = ', n), cex = 0.8)
}
dev.off()




N = 5
x = runif(N)
x2 = x * 100
y = runif(N)
y2 = y * 100
z = runif(N)
z2 = x + y


s1 = build.section(x, y, z, gridder = gridODV)
s2 = build.section(x, y, z, gridder = gridIDW)
s3 = build.section(x, y, z2, gridder = gridODV)
s4 = build.section(x, y, z2, gridder = gridIDW)
s5 = build.section(x2, y, z2, gridder = gridODV)
s6 = build.section(x2, y, z2, gridder = gridIDW)
s7 = build.section(x2, y2, z, gridder = gridODV)
s8 = build.section(x2, y2, z, gridder = gridIDW)

s9 = build.section(x, y, z, gridder = gridNNI)
s10 = build.section(x2, y, z, gridder = gridNNI)
s11 = build.section(x2, y2, z, gridder = gridNNI)
s12 = build.section(x2, y2, z, gridder = gridNN)

{
  pdf('series of tests.pdf')
  par(mfrow = c(2,2))
  plot.section(s1, pal = 'inferno', zlim = c(0,1)); add.section.contour(s1)
  plot.section(s2, pal = 'inferno', zlim = c(0,1)); add.section.contour(s2)
  plot.section(s3, pal = 'inferno', zlim = c(0,2)); add.section.contour(s3)
  plot.section(s4, pal = 'inferno', zlim = c(0,2)); add.section.contour(s4)
  plot.section(s5, pal = 'inferno', zlim = c(0,2)); add.section.contour(s5)
  plot.section(s6, pal = 'inferno', zlim = c(0,2)); add.section.contour(s6)
  plot.section(s7, pal = 'inferno', zlim = c(0,1)); add.section.contour(s7)
  plot.section(s8, pal = 'inferno', zlim = c(0,1)); add.section.contour(s8)

  plot.section(s9, pal = 'inferno', zlim = c(0,1)); add.section.contour(s9)
  plot.section(s10, pal = 'inferno', zlim = c(0,1)); add.section.contour(s10)
  plot.section(s11, pal = 'inferno', zlim = c(0,1)); add.section.contour(s11)
  plot.section(s12, pal = 'inferno', zlim = c(0,1)); add.section.contour(s12)
  dev.off()
}

