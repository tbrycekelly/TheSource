## ----setup, include = FALSE----------------------------------------------
library(TheSource)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
x = c(1:10)
y = x + rnorm(10)

plot(x, y)
add.error.bars(x, 1, y, rnorm(10))

