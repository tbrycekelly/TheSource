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
