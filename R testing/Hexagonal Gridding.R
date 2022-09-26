make.grid = function(x, y, z = NULL, fun = length, type = 1, res = NULL) {
  
  if (is.null(res)) {
    res = min(diff(range(x)), diff(range(y)))/10
  }
  
  if (type = 1) {
    xx = seq(min(x), max(x), by = res)
    yy = seq(min(y), max(y), by = res)
  } else if (type == 2) {
    xx = seq(min(x), max(x), by = res)
    yy = seq(min(y), max(y), by = res)
    nx = length(xx)
    ny = length(yy)
    
    xx = matrix(xx, nrow = nx, ncol = ny)
    yy = t(matrix(rep(c(yy, yy + res/2), nx/2), nrow = ny, ncol = nx))
    
  }
}

