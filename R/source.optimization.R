t## Set of useful R functions for optimization problems
##
## Author: Thomas Bryce Kelly (tbk14 at fsu.edu)
## http://about.tkelly.org/
##
## Dept of Earth, Ocean & Atmospherical Sciences
## Florida State University
##
## Center for Ocean & Atmospheric Prediction Studies
## Florida State University
##
## National High Magnetic Field Laboratory
## Florida State University

#' @title Parameter Search
#' @author Thomas Bryce Kelly
#' @description Implements a recursive grid search routine to solve optimization problems in arbitrary dimensions.
#' @param n Number of recursions to perform
#' @param cost The cost function which must return a numeric value and accept parameter values as the first arguments and in the order they are provided.
#' @param ... Optional argument that is passed directly onto the cost function
#' @param bounds A dataframe containing the minimum and maximum values permitted of each parameter
#' @param splits The number of subdivisions to perform for each dimension (so grid size is n x splits ^ dimensionality)
#' @export
parameter.search = function(n, cost, ..., bounds, splits = 10) {

  if (splits <= 2) { stop('splits argument must be an integer greater than or equal to 3!')}

  splits = round(splits)

  ## How many dimensions
  dim = nrow(bounds)

  ## Setup search grid
  b = list()
  for (i in 1:nrow(bounds)) {
    b[[i]] = seq(bounds[i,1], bounds[i,2], length.out = splits)
  }
  grid = do.call('expand.grid', b)
  grid$cost = NA

  ## Calculate cost function at each grid location
  for (i in 1:nrow(grid)) {
    args = as.list(grid[i, 1:dim])
    args = c(args, list(...))
    names(args) = formalArgs(cost)[1:length(args)]
    grid$cost[i] = do.call(cost, args) # cost(grid[i, 1:dim])
  }

  ## Best grid location
  l = which.min(grid$cost)
  if (n == 1) {
    res = list(min = grid[l,], bounds = bounds, grid = grid, history = grid[l,])
  } else{

    ## Setup new bounding box
    loci = grid[l, 1:dim]
    bounds.new = data.frame(min = as.numeric(loci), max = as.numeric(loci))
    for (i in 1:dim) {
      bounds.new$min[i] = max(bounds$min[i], bounds.new$min[i] - (bounds$max[i] - bounds$min[i]) / splits)
      bounds.new$max[i] = min(bounds$max[i], bounds.new$max[i] + (bounds$max[i] - bounds$min[i]) / splits)
    }

    ## Call parameter.search recursively
    res = parameter.search(n-1, cost, ..., bounds = bounds.new, splits = splits)
    res$history = rbind(res$history, grid[l,])
  }
  ## Return
  res
}



#### Examples
# Find the local minimum of each function within a given bounded region
#
## A 1-D example function
#f1 = function(x) {
#  3*x^2 - 4*x + 5
#}
#
## A 2-D example
#f2 = function(x, y) {
#  3*x^2 - 4*x * y + y^2 +1
#}
#
## A 3-D example
#f3 = function(x, y, z) { 3*x * y + z^1.5 - 5}
#
## A 4-D example
#f4 = function(x, y, z, h) {
#  3*x * y + z^1.5 - 5 * sinh(h)
#}
#
#ans1 = parameter.search(4, f1, bounds = data.frame(min = c(-5), max = c(5)))
#ans2 = parameter.search(4, f3, z = 1, bounds = data.frame(min = c(-5, -5), max = c(5, 5)))
#ans3 = parameter.search(4, f3, bounds = data.frame(min = c(-5, -5, -5), max = c(5, 5, 5)))
#ans4 = parameter.search(4, f4, bounds = data.frame(min = c(-5, -5, -5, -1), max = c(5, 5, 5, 1)))

