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
parameter.search = function(n, cost, grid = NULL, bounds, splits = 10, progression = NULL, verbose = T, ...) {

  if (verbose){
    message('Starting Parameter Search')
    message('Parameter search is an iterative grid search algorithm that seeks the parameters that minimize\nan objective cost function given a set of parameter bounds.')
  }

  ## Init Timer
  a = Sys.time()

  if (splits <= 2) { stop('splits argument must be an integer greater than or equal to 3!')}
  if(is.null(progression)) {
    progression = max(1, ceiling(splits / 4))
    if (verbose) { message(' No Progresssion given, default to:\t', progression) }
  }

  splits = round(splits)
  if (verbose) {message(' Number of splits set to:\t\t', splits)}
  ## How many dimensions
  dim = nrow(bounds)

  ## Setup search grid
  b = list()
  for (i in 1:nrow(bounds)) {
    if (bounds[i,1] == bounds[i,2]) {
      b[[i]] = bounds[i,1]
    } else {
        b[[i]] = seq(bounds[i,1], bounds[i,2], length.out = splits)
    }
  }
  argnames = formalArgs(cost)
  argnames = argnames[!argnames %in% names(list(...))] # don't include arguments passed through elipsis

  if (is.null(grid)) {
    grid = do.call('expand.grid', b)
    colnames(grid) = formalArgs(cost)[1:dim]
  }
  grid$cost = NA

  if (verbose) {
    message(' Mimimum bounds provided for ', paste(colnames(grid)[-ncol(grid)], collapse = ', '), ' are: ',
            paste(bounds[,1], collapse = ', '))
    message(' Maximum bounds provided for ', paste(colnames(grid)[-ncol(grid)], collapse = ', '), ' are: ',
            paste(bounds[,2], collapse = ', '))
  }

  ## Calculate cost function at each grid location
  best = Inf
  for (i in 1:nrow(grid)) {
    args = as.list(grid[i, 1:dim])
    names(args) = argnames[1:length(args)]
    args = c(args, list(...))

    #if (length(args) > length(formalArgs(cost))) { stop('Number of function arguments exceeds what function is expecting.') }
    grid$cost[i] = do.call(cost, args) # cost(grid[i, 1:dim])
    if (verbose & grid$cost[i] < best) {
      message(Sys.time(), ': New optimal parameter set found for n = ', n,': ', paste(grid[i,], collapse = ', '))
    }
  }

  ## Best grid location
  l = which.min(grid$cost)
  if (n == 1) {
    res = list(min = grid[l,], bounds = bounds, grid = grid, history = grid[l,], full.grid = grid)
  } else{

    ## Setup new bounding box
    loci = grid[l, 1:dim]
    bounds.new = data.frame(min = as.numeric(loci), max = as.numeric(loci))
    for (i in 1:dim) {
      bounds.new$min[i] = max(bounds$min[i], bounds.new$min[i] - progression * (bounds$max[i] - bounds$min[i]) / splits)
      bounds.new$max[i] = min(bounds$max[i], bounds.new$max[i] + progression * (bounds$max[i] - bounds$min[i]) / splits)
    }

    ## Call parameter.search recursively
    res = parameter.search(n-1, cost, ..., bounds = bounds.new, splits = splits, progression = progression, verbose = F)
    res$history = rbind(res$history, grid[l,])
    res$full.grid = rbind(res$full.grid, grid)
  }
  ## Return
  res
}

#' @title Parameter Search for Parallel Processing
#' @author Thomas Bryce Kelly
#' @description Implements a recursive grid search routine to solve optimization problems in arbitrary dimensions.
#' @param n Number of recursions to perform
#' @param cost The cost function which must return a numeric value and accept parameter values as the first arguments and in the order they are provided.
#' @param ... Optional argument that is passed directly onto the cost function
#' @param bounds A dataframe containing the minimum and maximum values permitted of each parameter
#' @param splits The number of subdivisions to perform for each dimension (so grid size is n x splits ^ dimensionality)
#' @param progression The size of the new search-space to interrogate. A value between 1 and splits/2. Default value (NULL) will yeield a progression of max(1, splits/4), good for most problems.
#' @export
parameter.search.parallel = function(n, cost, grid = NULL, ..., bounds, splits = 10, progression = NULL) {

  if (splits <= 2) { stop('splits argument must be an integer greater than or equal to 3!')}
  if(is.null(progression)) { progression = max(1, ceiling(splits / 4))}

  splits = round(splits)

  ## How many dimensions
  dim = nrow(bounds)

  ## Setup search grid
  b = list()
  for (i in 1:nrow(bounds)) {
    if (bounds[i,1] == bounds[i,2]) {
      b[[i]] = bounds[i,1]
    } else {
      b[[i]] = seq(bounds[i,1], bounds[i,2], length.out = splits)
    }
  }
  grid = do.call('expand.grid', b)
  colnames(grid) = formalArgs(cost)[1:dim]

  cn = parallel::detectCores() - 1
  cl = parallel::makeCluster(cn)
  ### PASS THE OBJECT FROM MASTER PROCESS TO EACH NODE
  parallel::clusterExport(cl, varlist = c('grid', 'cost', names(optional.args), ls(envir = globalenv())), envir = environment())

  ### DIVIDE THE DATAFRAME BASED ON # OF CORES
  sp = parallel::parLapply(cl, parallel::clusterSplit(cl = cl, seq = seq(nrow(grid))), function(c) {grid[c,]})

  grid$cost = Reduce(c, parallel::parLapply(cl, sp, fun = function(s) {
    s = data.frame(s)
    res = rep(NA, nrow(s))
    if (length(s) > 0) {
      for (i in 1:nrow(s)) {
        args = as.list(s[i,])
        names(args) = formalArgs(cost)[1:length(args)]
        args = c(args, list(...))

        res[i] = do.call(cost, args) # cost(grid[i, 1:dim])
      }
    } else {
      return()
    }
    return(res)
  }, chunk.size = 1))
  parallel::stopCluster(cl)

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
parameter.walk = function (n, cost, ..., bounds, progression = NULL, max.iter = 1e6, k = 2) {

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
    l = param < bounds[,1]
    ll = param > bounds[,2]
    param[l] = bounds[l,1]
    param[ll] = bounds[ll,2]

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

  list(min = history[i-1,], bounds = bounds, history = history[!is.na(history$cost)], meta = list(j = j, max.iter = max.iter, progression = progression))
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

  list(min = history[i-1,], bounds = bounds, history = history[!is.na(history$cost),], meta = list(j = j, max.iter = max.iter, progression = progression))
}


#' @title A Gradient Descent Algorithm
#' @author Thomas Bryce Kelly
#' @description Implements a recursive grid search routine to solve optimization problems in arbitrary dimensions.
#' @param n Number of recursions to perform
#' @param cost The cost function which must return a numeric value and accept parameter values as the first arguments and in the order they are provided.
#' @param ... Optional argument that is passed directly onto the cost function
#' @param bounds A dataframe containing the minimum and maximum values permitted of each parameter
#' @param splits The number of subdivisions to perform for each dimension (so grid size is n x splits ^ dimensionality)
#' @param progression The size of the new search-space to interrogate. A value between 1 and splits/2. Default value (NULL) will yeield a progression of max(1, splits/4), good for most problems.
#' @export
parameter.descent.old = function(cost, guess = NULL, ..., bounds = NULL, max.iter = 1e6, progression = NULL, step = 1e-6, tol = 1e-18) {

  ## Determine cost function arguemnts
  argnames = formalArgs(cost)
  argnames = argnames[!argnames %in% names(list(...))] # don't include arguments passed through elipsis


  if (is.null(guess)) {
    message('No initial guess supplied, starting at the origin: guess = c(0,...,0).')
    guess = rep(0, length(argnames))
  }
  if(is.null(progression)) {
    message('No progression length provided, using 1.')
    progression = 0.5
  }

  iter = 0
  ans = guess
  hist = as.data.frame(matrix(guess, ncol = length(guess)))
  colnames(hist) = argnames[1:length(guess)]

  while (T) {
    iter = iter + 1

    if (iter > max.iter) {
      break
    }

    grad = rep(0, length(guess))

    ## Calculate cost at current guess point
    arg0 = as.list(guess)
    names(arg0) = argnames[1:length(guess)]
    arg0 = c(arg0, list(...))
    cost0 = do.call(cost, arg0)

    ## Calculate gradient at current guess
    for (j in 1:length(guess)) {
      arg = arg0
      arg[[j]] = arg0[[j]] + step
      grad[j] = (do.call(cost, arg) - cost0) / step
    }

    ## Update guess
    guess = guess - progression * grad

    arg1 = as.list(guess)
    names(arg1) = argnames[1:length(guess)]
    arg1 = c(arg1, list(...))
    cost1 = do.call(cost, arg1)

    ## Add to History
    hist = rbind(hist, guess)

    if (abs(cost0 - cost1) <= tol) {
      break
    }
  }

  ## Return
  list(min = hist[nrow(hist),], bounds = bounds, guess = hist[1,], history = hist)
}

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
parameter.descent = function(cost, guess = NULL, ..., bounds = NULL, max.iter = 100, progression = NULL, step = 1, tol = 1e-8, verbose = T) {

  if (verbose){
    message('Starting Parameter Descent')
    message('Parameter search is an iterative grid search algorithm that seeks the parameters that minimize\nan objective cost function given a set of parameter bounds.')
  }

  ## Init Timer
  a = Sys.time()

  ## How many dimensions
  argnames = formalArgs(cost)
  argnames = argnames[!argnames %in% names(list(...))] # don't include arguments passed through elipsis
  dim = length(argnames)

  ## Setup simplex
  b = array(NA, dim = c(dim+1, dim))
  if (!is.null(guess)) {
    if (length(guess) == dim) {
      b[1,] = guess
    } else {
      message('Improper guess provided, ignoring.')
      b[1,] = 0 ## default is to start at origin.
    }
  } else {
    b[1,] = 0 ## default is to start at origin.
  }

  for (i in 1:dim) {
    b[i+1,] = b[1] + step / dim
    b[i+1,i] = b[1,i] + step ##
  }

  ## Setup simplex and history object
  simplex = as.data.frame(b)
  history = data.frame(cost = rep(NA, max.iter))
  for (n in argnames) {
    history[[n]] = NA
  }
  simplex = data.frame(cost = NA, simplex)

  ## Calculate cost function at each node
  for (j in which(is.na(simplex$cost))) {
    args = as.list(simplex[j,-1])
    names(args) = argnames[1:length(args)]
    args = c(args, list(...))

    simplex$cost[j] = do.call(cost, args) # cost(grid[i, 1:dim])
  }

  for (i in 1:max.iter) {
    ## Order points
    simplex = simplex[order(simplex$cost, decreasing = T),]
    history[i,] = simplex[nrow(simplex),]

    ## Add reflected point
    centroid = colMeans(simplex[-1,])
    test.point = centroid + 1 * (centroid - simplex[1,]) ## Reflected point

    ## Calculate cost
    args = as.list(test.point[-1])
    names(args) = argnames[1:length(args)]
    args = c(args, list(...))
    test.point[1] = do.call(cost, args) # cost(grid[i, 1:dim])

    ## Expansion
    if (max(simplex$cost) > test.point[1]) {
      simplex[1,] = centroid + 2 * (centroid - simplex[1,]) ## Replace worse case

      args = as.list(simplex[1,-1])
      names(args) = argnames[1:length(args)]
      args = c(args, list(...))
      simplex[1,1] = do.call(cost, args) # cost(grid[i, 1:dim])

      if (simplex[1,1] > test.point[1]) {
        simplex[1, ] = test.point # revert to test point if expansion fails (optional step)
        if (verbose) { message(' Simplex reflextion.') }
      } else {
        if (verbose) { message(' Simplex expansion.') }
      }

    } else {
      ## Contraction
      simplex[1,] = centroid - 0.5 * (centroid - simplex[1,]) ## replace reflected point
      args = as.list(simplex[1,-1])
      names(args) = argnames[1:length(args)]
      args = c(args, list(...))
      simplex[1,1] = do.call(cost, args)
      if (verbose) { message(' Simplex contraction') }
    }

    if (i > 2 * dim + 1) {
      if (abs(history$cost[i] - history$cost[i-2*dim-1]) < tol) {
        if (verbose) { message('Cost tolerance met, ending.')}
        break
      }
    }
  }
  res = list(min = simplex[which.min(simplex$cost),], simplex = simplex, history = history[!is.na(history$cost),])

  ## Return
  res
}




#' @title Get Optimization Test Function
#' @author Thomas Bryce Kelly
#' @references Wikipedia: https://en.wikipedia.org/wiki/Test_functions_for_optimization
#' @param n the number of the test function desired.
#' @export
test.function = function(n, verbose = T) {

  Rastrigin = function(...) {
    params = list(...)
    ans = 10 * length(params)
    for (i in 1:length(params)) { ans = ans + params[[i]]^2 - 10 * cos(2*pi * params[[i]]) }
    ans
  }

  Ackley = function(x, y) {
    -10 * exp(-0.2 * sqrt(0.5 * (x^2 + y^2))) - exp(-0.5 * (cos(2 * pi * x) + cos(2 * pi * y))) + 2.7182818 + 20
  }

  Sphere = function(...) {
    params = list(...)
    ans = 0
    for (i in 1:length(params)) {ans = ans + params[[i]]^2}
    return(ans)
  }

  Rosenbrock = function(...) {
    params = list(...)
    ans = 0
    for (i in 1:(length(params) - 1)) {ans = ans + 100 * (params[[i+1]] - params[[i]]^2)^2 + (1 - params[[i]])^2}
    return(ans)
  }

  Beale = function(x, y) {
    (1.5 - x + x*y)^2 + (2.25 - x + x * y^2)^2 + (2.625 - x + x * y^3)^2
  }

  Goldstein.Price  = function(x, y) {
    a = 1 + (x + y + 1)^2 * (19 - 14*x + 3 * x^2 - 14 * y + 6 * x * y + 3 * y^2)
    b = (30 + (2 * x - 3 * y)) * (18 - 32 * x + 12 * x^2 + 48 * y - 36 * x * y + 27 * y^2)
    a * b
  }

  Booth = function(x, y) {
    (x + 2 * y - 7)^2 + (2 * x + y - 5)^2
  }

  Bukin = function(x, y) {
    100 * sqrt(abs(y - 0.01 * x^2)) + 0.01 * abs(x + 10)
  }

  Matyas = function(x, y) {
    0.26 * (x^3 + y^3) - 0.48 * x * y
  }

  ## Start return series
  if (n == 1) {
    if (verbose) { message('Returning Rastrigin test function. Global minimum at f(0,...,0) = 0') }
    return(Rastrigin)
  }
  if (n == 2) {
    if (verbose) { message('Returning Ackley test function. Global minimum at f(0,0) = 0') }
    return(Ackley)
  }
  if (n == 3) {
    if (verbose) { message('Returning Sphere test function. Global minimum at f(0,...,0) = 0')}
    return(Sphere)
  }
  if (n == 4) {
    if (verbose) { message('Returning Rosenbrock test function. Global minimum at f(1,...,1) = 0') }
    return(Rosenbrock)
  }
  if (n == 5) {
    if (verbose) { message('Returning Beale test function. Global minimum at f(3, 0.5) = 0') }
    return(Beale)
  }

  if (n == 6) {
    if (verbose) { message('Returning Goldstein.Price test function. Global minimum at f(0, -1) = 3') }
    return(Goldstein.Price)
  }
  if (n == 7) {
    if (verbose) { message('Returning Booth test function. Global minimum at f(1, 3) = 0') }
    return(Booth)
  }
  if (n == 8) {
    if (verbose) { message('Returning Bukin test function. Global minimum at f(-10,1) = 0') }
    return(Bukin)
  }
  if (n == 9) {
    if (verbose) { message('Returning Matyas test function. Global minimum at f(0, 0) = 0') }
    return(Matyas)
  }

  message('Argument n should be between 1 and 9. Returning trivial test function')
  return(function(...) {0})
}
