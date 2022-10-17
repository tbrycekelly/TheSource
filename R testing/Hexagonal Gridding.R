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



calc.vertex.poly = function(x, y) {
  dx = 
}




calc.vertex = function(x, y) {
  
  ## Diffs
  dx.dx = t(diff(t(x))) / 2
  dx.dy = diff(x) / 2
  dy.dx = t(diff(t(y))) / 2
  dy.dy = diff(y) / 2
  
  poly = vector("list", length(x) * length(y))
  
  
  
  ## Vertex
  vertex.x = matrix(NA, nrow = dim(x)[1]+1, ncol = dim(x)[2]+1)
  vertex.y = matrix(NA, nrow = dim(y)[1]+1, ncol = dim(y)[2]+1)
  
  
  ## Field
  for (i in 2:(dim(vertex.x)[1] - 1)) {
    for (j in 2:(dim(vertex.x)[2] - 1)) {
      ii = max(i-1, 1)
      jj = max(j-1, 1)
      
      vertex.x[i,j] = x[ii,jj] + dx.dx[ii,jj] + dx.dy[ii,jj]
      vertex.y[i,j] = y[ii,jj] + dy.dx[ii,jj] + dy.dy[ii,jj]
      
      poly[i * j] = data.frame(x = c(x[ii,jj] + dx.dx[ii,jj] + dx.dy[ii,jj], x[ii,jj] + dx.dx[ii,jj] + dx.dy[ii,jj], x[ii,jj] + dx.dx[ii,jj] + dx.dy[ii,jj], x[ii,jj] + dx.dx[ii,jj] + dx.dy[ii,jj]),
                               y = y[ii,jj] + dy.dx[ii,jj] + dy.dy[ii,jj])
      
    }
  }
  
  
  ## Fill in perimeter
  # i = 1
  for (j in 2:(dim(vertex.x)[2] - 1)) {
    jj = max(j-1, 1)
    
    vertex.x[1,j] = x[1,jj] + dx.dx[1,jj] - dx.dy[1,jj]
    vertex.y[1,j] = y[1,jj] + dy.dx[1,jj] - dy.dy[1,jj]
  }
  
  # j = 1
  for (i in 2:(dim(vertex.x)[1] - 1)) {
    ii = max(i-1, 1)
    
    vertex.x[i,1] = x[ii,1] - dx.dx[ii,1] + dx.dy[ii,1]
    vertex.y[i,1] = y[ii,1] - dy.dx[ii,1] + dy.dy[ii,1]
  }
  
  # j = dim(vertex.x)[2]
  for (i in 1:(dim(vertex.x)[1] - 1)) {
    ii = max(i-1, 1)
    j = dim(vertex.x)[2]
    
    vertex.x[i,j] = x[ii,j-1] + dx.dx[ii,j-2] + dx.dy[ii,j-1]
    vertex.y[i,j] = y[ii,j-1] + dy.dx[ii,j-2] + dy.dy[ii,j-1]
  }
  
  # i = dim(vertex.x)[2]
  for (j in 1:(dim(vertex.x)[2] - 1)) {
    jj = max(j-1, 1)
    i = dim(vertex.x)[1]
    
    vertex.x[i,j] = x[i-1,jj] + dx.dx[i-1,jj] + dx.dy[i-2,jj]
    vertex.y[i,j] = y[i-1,jj] + dy.dx[i-1,jj] + dy.dy[i-2,jj]
  }
  
  ## Fill in corners
  ## both = 1
  vertex.x[1,1] = x[1,1] - dx.dx[1,1] - dx.dy[1,1]
  vertex.y[1,1] = y[1,1] - dy.dx[1,1] - dy.dy[1,1]
  
  i = dim(vertex.x)[1]
  vertex.x[i,1] = x[i-1,1] - dx.dx[i-1,1] + dx.dy[i-2,1]
  vertex.y[i,1] = y[i-1,1] - dy.dx[i-1,1] + dy.dy[i-2,1]
  
  i = dim(vertex.x)[1]
  j = dim(vertex.x)[2]
  vertex.x[i,j] = x[i-1,j-1] + dx.dx[i-1,j-2] + dx.dy[i-2,j-1]
  vertex.y[i,j] = y[i-1,j-1] + dy.dx[i-1,j-2] + dy.dy[i-2,j-1]
  
  j = dim(vertex.x)[2]
  vertex.x[1,j] = x[1,j-1] + dx.dx[1,j-2] - dx.dy[1,j-1]
  vertex.y[1,j] = y[1,j-1] + dy.dx[1,j-2] - dy.dy[1,j-1]
  
  list(x = vertex.x, y = vertex.y)
}
