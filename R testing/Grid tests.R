library(TheSource)
library(microbenchmark)

##### Gridder engine test

## gridding random noise
pdf('R testing/grid tests - runif.pdf')
par(mfrow = c(2,2))
for (n in c(10,100, 1e3)) {
  x = runif(n, 0, 10)
  y = runif(n, 0, 10)
  z = runif(n)

  for (gridder in c(gridNN, gridNNI, gridODV, gridIDW)) {
    section = build.section(x, y, z, gridder = gridder) # defaults
    plot.section(section, zlim = c(0,1), pal = 'inferno')
  }
}
dev.off()



pdf('R testing/grid tests - rnorm.pdf')
par(mfrow = c(2,2))
for (n in c(10,100, 1e3)) {
  x = runif(n, 0, 10)
  y = runif(n, 0, 10)
  z = rnorm(n)

  for (gridder in c(gridNN, gridNNI, gridODV, gridIDW)) {
    section = build.section(x, y, z, gridder = gridder) # defaults
    plot.section(section, zlim = c(-1,1), pal = 'ocean.balance')
  }
}
dev.off()

pdf('R testing/grid tests - grid runif.pdf')
par(mfrow = c(2,2))
for (n in c(5, 25, 125)) {
  grid = expand.grid(x = c(1:n), y = c(1:n))
  z = runif(n^2)

  for (gridder in c(gridNN, gridNNI, gridODV, gridIDW)) {
    section = build.section(grid$x, grid$y, z, gridder = gridder) # defaults
    plot.section(section, zlim = c(0,1), pal = 'inferno')
  }
}
dev.off()


pdf('R testing/grid tests - grid rnorm.pdf')
par(mfrow = c(2,2))
for (n in c(5, 25, 125)) {
  grid = expand.grid(x = c(1:n), y = c(1:n))
  z = rnorm(n^2)

  for (gridder in c(gridNN, gridNNI, gridODV, gridIDW)) {
    section = build.section(grid$x, grid$y, z, gridder = gridder) # defaults
    plot.section(section, zlim = c(-1,1), pal = 'ocean.balance')
  }
}
dev.off()



pdf('R testing/grid tests - grid rnorm scaled with uncertainty.pdf')
par(mfrow = c(2,2))
for (n in c(5, 25, 125)) {
  grid = expand.grid(x = c(1:n), y = c(1:n))
  z = rnorm(n^2)

  for (gridder in c(gridNN, gridNNI, gridODV, gridIDW)) {
    section = build.section(grid$x*100, grid$y*100, z, gridder = gridder, nx = 100, ny = 100, uncertainty = 1) # defaults
    plot.section(section, zlim = c(-1,1), pal = 'ocean.balance')
  }
}
dev.off()

pdf('R testing/grid tests - grid rnorm scaled without uncertainty.pdf')
par(mfrow = c(2,2))
for (n in c(5, 25, 125)) {
  grid = expand.grid(x = c(1:n), y = c(1:n))
  z = rnorm(n^2)

  for (gridder in c(gridNN, gridNNI, gridODV, gridIDW)) {
    section = build.section(grid$x*100, grid$y*100, z, gridder = gridder, nx = 100, ny = 100, uncertainty = 1e-5) # defaults
    plot.section(section, zlim = c(-1,1), pal = 'ocean.balance')
  }
}
dev.off()




#### Benchmarking

N = 25
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
}


times = data.frame(i = 0, t = 0)
for (i in 6) {
  message(i)
  N = 10^i

  x = runif(N)
  y = runif(N)
  z = runif(N) + x - y

  a = microbenchmark(build.section(x, y, z, gridder = gridODV, verbose = F), times = 20)

  temp = data.frame(i = i, t = as.numeric(a$time))
  times = rbind(times, temp)
}
times = times[-1,]

plot(NULL, NULL, xlim = c(1,6), ylim = c(0,7), xaxt = 'n', ylab = 'Calculations/s', xlab = 'N', yaxt = 'n')
add.log.axis(1)
add.log.axis(2)
add.violin(times$i, 10^times$i/times$t*1e12, col = 'grey', side = 2, scale = 0.03)
add.boxplot.box(times$i, log10(10^times$i/times$t*1e12), col = 'grey')
