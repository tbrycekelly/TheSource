## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
library(TheSource)
test.map.california()

## ------------------------------------------------------------------------
#par(mfrow=c(2,2))

for (i in 1:8) {
  p = make.proj(i)
  
  map = make.map(p = p)
  mtext(p, side = 3, line = 1, col = 'red', adj = 1, cex = 0.5)
  
}

