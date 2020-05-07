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
#' @param progression The size of the new search-space to interrogate. A value between 1 and splits/2. Default value (NULL) will yeield a progression of max(1, splits/4), good for most problems.
#' @export
parameter.search = function(n, cost, grid = NULL, ..., bounds, splits = 10, progression = NULL) {

  if (splits <= 2) { stop('splits argument must be an integer greater than or equal to 3!')}
  if(is.null(progression)) { progression = max(1, ceiling(splits / 4))}

  splits = round(splits)

  ## How many dimensions
  dim = nrow(bounds)

  ## Setup search grid
  b = list()
  for (i in 1:nrow(bounds)) {
    b[[i]] = seq(bounds[i,1], bounds[i,2], length.out = splits)
  }
  grid = do.call('expand.grid', b)
  colnames(grid) = formalArgs(cost)[1:dim]
  grid$cost = NA

  ## Calculate cost function at each grid location
  for (i in 1:nrow(grid)) {
    args = as.list(grid[i, 1:dim])
    names(args) = formalArgs(cost)[1:length(args)]
    args = c(args, list(...))

    if (length(args) > length(formalArgs(cost))) { stop('Number of function arguments exceeds what function is expecting.') }
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
      bounds.new$min[i] = max(bounds$min[i], bounds.new$min[i] - progression * (bounds$max[i] - bounds$min[i]) / splits)
      bounds.new$max[i] = min(bounds$max[i], bounds.new$max[i] + progression * (bounds$max[i] - bounds$min[i]) / splits)
    }

    ## Call parameter.search recursively
    res = parameter.search(n-1, cost, ..., bounds = bounds.new, splits = splits, progression = progression)
    res$history = rbind(res$history, grid[l,])
  }
  ## Return
  res
}


#' @title Parameter Search using Random Walk
#' @author Thomas Bryce Kelly
#' @description Implements an MCMC algorithm (metropolis criterion) to determine a global optimum from within the stated bounds.
#' @param n Number of accepted solutions to seek
#' @param cost The cost function which must return a numeric value and accept parameter values as the first arguments and in the order they are provided.
#' @param ... Optional argument that is passed directly onto the cost function
#' @param bounds A dataframe containing the minimum and maximum values permitted of each parameter
#' @param progression The step length for each parameter (default is 1% of bounded range).
#' @param max.iter The maximum number of steps to take before returning (even if n is unsatisfied). **Goal for function to return via n rather than by max.iter.
#' @export
parameter.anneal = function (n, cost, ..., bounds, progression = NULL, max.iter = 1e6, k = 2) {

  if (is.null(progression)) {
    progression = as.numeric(bounds[,2] - bounds[,1]) / 100
    warning('Progression vector empty, defaulting to 1% of range.')
  }

  ## Setup History
  history = data.frame(Step = c(1:n))
  for (i in 1:nrow(bounds)) { history[formalArgs(cost)[i]] = NA }
  history$cost = NA

  ## Start in random location
  start = runif(nrow(bounds), min = bounds[,1], max = bounds[,2])
  args = as.list(start)
  names(args) = formalArgs(cost)[1:length(args)]
  args = c(args, list(...))
  if (length(args) > length(formalArgs(cost))) {
    stop("Number of function arguments exceeds what function is expecting.")
  }
  current.cost = do.call(cost, args)

  ## Update history
  history[1,] = c(1, start, current.cost)

  ## Setup loop
  i = 2
  for (j in 1:max.iter) {
    ## randomly perturbe the location
    param = start + rnorm(length(start), mean = 0, sd = progression)

    args = as.list(param)
    names(args) = formalArgs(cost)[1:length(args)]
    args = c(args, list(...))
    if (length(args) > length(formalArgs(cost))) {
      stop("Number of function arguments exceeds what function is expecting.")
    }
    new.cost = do.call(cost, args)

    if (new.cost < current.cost | exp( (new.cost - current.cost)/2 ) < runif(1)) {
      ## Update history
      history[i,] = c(i, param, new.cost)
      current.cost = new.cost
      start = param
      i = i + 1
      if (i > n) { break }
    }
  }

  list(min = history[i-1,], bounds = bounds, history = history, meta = list(j = j, max.iter = max.iter, progression = progression))
}



#' @title Parameter Search using Simulated Annealing
#' @author Thomas Bryce Kelly
#' @description Implements an MCMC algorithm (simulated annealing) to determine a global optimum from within the stated bounds.
#' @param n Number of accepted solutions to seek
#' @param cost The cost function which must return a numeric value and accept parameter values as the first arguments and in the order they are provided.
#' @param ... Optional argument that is passed directly onto the cost function
#' @param bounds A dataframe containing the minimum and maximum values permitted of each parameter
#' @param progression The step length for each parameter (default is 1% of bounded range).
#' @param max.iter The maximum number of steps to take before returning (even if n is unsatisfied). **Goal for function to return via n rather than by max.iter.
#' @param k The cooling coefficent. Larger values imply higher initial "temperature" and easier solution acceptance.
#' @export
parameter.anneal = function (n, cost, ..., bounds, progression = NULL, max.iter = 1e6, k = 2) {

  if (is.null(progression)) {
    progression = as.numeric(bounds[,2] - bounds[,1]) / 100
    warning('Progression vector empty, defaulting to 1% of range.')
  }

  ## Setup History
  history = data.frame(Step = c(1:n))
  for (i in 1:nrow(bounds)) { history[formalArgs(cost)[i]] = NA }
  history$cost = NA

  ## Start in random location
  start = runif(nrow(bounds), min = bounds[,1], max = bounds[,2])
  args = as.list(start)
  names(args) = formalArgs(cost)[1:length(args)]
  args = c(args, list(...))
  if (length(args) > length(formalArgs(cost))) {
    stop("Number of function arguments exceeds what function is expecting.")
  }
  current.cost = do.call(cost, args)

  ## Update history
  history[1,] = c(1, start, current.cost)

  ## Setup loop
  i = 2
  for (j in 1:max.iter) {
    ## randomly perturbe the location
    param = start + rnorm(length(start), mean = 0, sd = progression)

    args = as.list(param)
    names(args) = formalArgs(cost)[1:length(args)]
    args = c(args, list(...))
    if (length(args) > length(formalArgs(cost))) {
      stop("Number of function arguments exceeds what function is expecting.")
    }
    new.cost = do.call(cost, args)

    if (exp( (new.cost - current.cost) / ((n-i)/n * k)) < runif(1)) {
      ## Update history
      history[i,] = c(i, param, new.cost)
      current.cost = new.cost
      start = param
      i = i + 1
      if (i > n) { break }
    }
  }

  list(min = history[i-1,], bounds = bounds, history = history, meta = list(j = j, max.iter = max.iter, progression = progression))
}
