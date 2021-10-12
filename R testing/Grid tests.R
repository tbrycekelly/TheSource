library(TheSource)
library(microbenchmark)

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


plot(NULL, NULL, xlim = c(0,6), ylim = c(0,100), xaxt = 'n')
add.log.axis(1)

for (i in 1:2) {
  N = 10^i

  x = runif(N)
  y = runif(N)
  z = runif(N) + x - y

  a = microbenchmark(build.section(x, y, z, gridder = gridODV, verbose = F))
  add.violin(i, a$time/1e6, col = 'grey', side = 2)

}


