#' @title Load NetCDF
#' @import ncdf4
#' @param file The location of a .nc file.
#' @export
load.nc = function(file, var = NULL, test = F, verbose = T) {
  
  if (verbose) { message('Attempting to load data from ', file)}
  file = ncdf4::nc_open(file)
  if (verbose) { message(' File openned.')}
  
  all.var = names(file$var)
  for (name in names(file$dim)) {
    if (file$dim[[name]]$dimvarid$id > 0) {
      all.var = c(all.var, name)
    }
  }
  all.var = all.var[all.var != 'nv']
  
  
  if (is.null(var)) {
    var = all.var
    if (verbose) { message(' No variables specified, loading all ', length(var), ' entires.') }
  }
  
  ## Allow number of variable to be used as well
  if (is.numeric(var)) {
    var = names(file$var)[var]
  }
  
  ## Load data
  data = list()
  
  for (v in var) {
    
    if (v %in% all.var) {
      a = Sys.time()
      data[[v]] = ncdf4::ncvar_get(nc = file, varid = v)
      if (test) { data[[v]] = data[[v]][1] }
      if (verbose) {
        d = dim(data[[v]])
        if (is.null(d)) { d = length(data[[v]])}
        message(' Variable ', v, ' (', paste(d, collapse = 'x'),') loaded from file. \t', round(difftime(Sys.time(), a, units = 'secs')), 's')
      }
    } else {
      a = Sys.time()
      data[[v]] = NA
      if (verbose) {
        message(' Variable ', v, ' (X) does not exist in file. \t', round(difftime(Sys.time(), a, units = 'secs')), 's')
      }
    }
  }
  
  if (verbose) {message(' Closing file. Finished.')}
  ncdf4::nc_close(file)
  
  ##return
  data
}


#' @title Load NetCDF
#' @import ncdf4
#' @param file The location of a .nc file.
#' @export
load.nc.subset = function(file, var = NULL, x = NULL, y = NULL, test = F, verbose = T) {
  
  if (verbose) { message('Attempting to load data from ', file)}
  file = ncdf4::nc_open(file)
  if (verbose) { message(' File openned.')}
  
  all.var = names(file$var)
  for (name in names(file$dim)) {
    if (file$dim[[name]]$dimvarid$id > 0) {
      all.var = c(all.var, name)
    }
  }
  all.var = all.var[all.var != 'nv']
  
  
  if (is.null(var)) {
    var = all.var
    if (verbose) { message(' No variables specified, loading all ', length(var), ' entires.') }
  }
  
  ## Allow number of variable to be used as well
  if (is.numeric(var)) {
    var = names(file$var)[var]
  }
  
  ## Load data
  data = list()
  
  for (v in var) {
    
    if (v %in% all.var) {
      a = Sys.time()
      data[[v]] = ncdf4::ncvar_get(nc = file, varid = v, )
      if (test) { data[[v]] = data[[v]][1] }
      if (verbose) {
        d = dim(data[[v]])
        if (is.null(d)) { d = length(data[[v]])}
        message(' Variable ', v, ' (', paste(d, collapse = 'x'),') loaded from file. \t', round(difftime(Sys.time(), a, units = 'secs')), 's')
      }
    } else {
      a = Sys.time()
      data[[v]] = NA
      if (verbose) {
        message(' Variable ', v, ' (X) does not exist in file. \t', round(difftime(Sys.time(), a, units = 'secs')), 's')
      }
    }
  }
  
  if (verbose) {message(' Closing file. Finished.')}
  ncdf4::nc_close(file)
  
  ##return
  data
}


#' @title Read NetCDF Variables
#' @import ncdf4
#' @param file The location of a .nc file.
#' @export
load.nc.vars = function(file, verbose = F) {
  
  if (verbose) { message('Attempting to load variables from ', file)}
  file = ncdf4::nc_open(file)
  if (verbose) { message(' File openned.')}
  
  var = c(names(file$var), names(file$dim))
  
  if (verbose) {message(' Closing file. Finished.')}
  ncdf4::nc_close(file)
  
  ##return
  var
}


#' @title Read NetCDF Variables
#' @import ncdf4
#' @param file The location of a .nc file.
#' @export
load.nc.dims = function(file, verbose = F) {
  
  if (verbose) { message('Attempting to load dims from ', file)}
  file.con = ncdf4::nc_open(file)
  if (verbose) { message(' File openned.')}
  
  dim = list()
  var = names(file.con$dim)
  for (v in var) {
    dim[[v]] = file.con$dim[[v]]$vals
  }
  
  vars = list()
  var = names(file.con$var)
  
  for (v in var) {
    vars[[v]] = file.con$var[[v]]$size
    names(vars[[v]]) = file.con$var[[v]]$dimids + 1
  }
  
  if (verbose) {message(' Closing file. Finished.')}
  ncdf4::nc_close(file.con)
  
  ##return
  list(dim = dim, var = vars)
}