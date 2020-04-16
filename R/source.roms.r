## Set of useful R functions for working with Regional Ocean
## Model System (ROMS) data and netcdf files.
##
## Including LTRANS files.
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


#' @title Print NetCDF
#' @export
#' @import ncdf4
print.nc = function(file) {
    file = ncdf4::nc_open(file)
    print(file)
    ncdf4::nc_close(file)
}


#' @title Convert ROMS Time
#' @author Thomas Bryce Kelly
#' @param x The datetime numeric to be converted
#' @param tz The timezone for the conversion, almost always UTC
#' @export
conv.time.roms = function(x, tz='UTC') {
    as.POSIXct(x, origin="1900-01-01")
}


#' @title Get ROMS Time
#' @export
#' @author Thomas Bryce Kelly
#' @import ncdf4
get.roms.time = function(path, convert = TRUE) {

    ## Open the nc file
    roms = ncdf4::nc_open(path, write=FALSE)

    nhist = ncdf4::ncvar_get(roms, varid = 'nHIS') # dt * nHist = time between records
    dstart = ncdf4::ncvar_get(roms, varid = 'dstart') # days since 1900-01-01 00:00:00
    dt = ncdf4::ncvar_get(roms, varid = 'dt') # Sec
    ntimes = ncdf4::ncvar_get(roms, varid = 'ntimes') # number of time frames

    ## Close the nc file
    ncdf4::nc_close(roms)

    ## Calculate the times (seconds from 1900-01-01)
    roms.times = seq(from = dstart, by = dt * nhist / 86400, length.out = ntimes/nhist) * 86400

    if (convert) {
        roms.times = conv.time.roms(roms.times, tz='GMT')
    }
    roms.times
}


#' @title Load in ROMS Data
#' @export
#' @author Thomas Bryce Kelly
#' @import ncdf4
load.roms = function(path) {
    file = ncdf4::nc_open(path, write=FALSE)

    #### Load variables
    ## Vertical Grid
    h = ncdf4::ncvar_get(file, 'h')
    hc = ncdf4::ncvar_get(file, 'hc')
    theta = ncdf4::ncvar_get(file, 'theta_s')
    N = 42 # Hard coded

    h = (h[2:dim(h)[1],] + h[1:(dim(h)[1] - 1),]) / 2 ## interpolate to center of each grid cell
    h = (h[,2:dim(h)[2]] + h[,1:(dim(h)[2] - 1)]) / 2

    ## Horizontal
    lat = ncdf4::ncvar_get(file, 'lat_psi')
    lon = ncdf4::ncvar_get(file, 'lon_psi')
    grid = list(lat = c(lat), lon = c(lon))

    ## Advect
    u = ncdf4::ncvar_get(file, 'u')
    u = (u[,2:dim(u)[2],,] + u[,1:(dim(u)[2] - 1),,]) / 2

    v = ncdf4::ncvar_get(file, 'v')
    v = (v[2:dim(v)[1],,,] + v[1:(dim(v)[1] - 1),,,]) / 2

    w = ncdf4::ncvar_get(file, 'w')
    w = (w[,,2:dim(w)[3],] + w[,,1:(dim(w)[3] - 1),]) / 2 # depth averaged w
    w = (w[2:dim(w)[1],,,] + w[1:(dim(w)[1] - 1),,,]) / 2
    w = (w[,2:dim(w)[2],,] + w[,1:(dim(w)[2] - 1),,]) / 2

    dims = c(dim(u)[1] * dim(u)[2], dim(u)[3], dim(u)[4])
    dim(u) = dims
    dim(v) = dims
    dim(w) = dims

    ## Diffusivity
    AKt = ncdf4::ncvar_get(file, 'AKt')

    ## Fields
    temp = ncdf4::ncvar_get(file, 'temp')
    if ('rho' %in% names(file$var)) {rho = ncvar_get(file, 'rho')} else {rho = NULL}
    salt = ncvar_get(file, 'salt')

    temp = (temp[2:dim(temp)[1],,,] + temp[1:(dim(temp)[1] - 1),,,]) / 2
    temp = (temp[,2:dim(temp)[2],,] + temp[,1:(dim(temp)[2] - 1),,]) / 2

    rho = (rho[2:dim(rho)[1],,,] + rho[1:(dim(rho)[1] - 1),,,]) / 2
    rho = (rho[,2:dim(rho)[2],,] + rho[,1:(dim(rho)[2] - 1),,]) / 2

    salt = (salt[2:dim(salt)[1],,,] + salt[1:(dim(salt)[1] - 1),,,]) / 2
    salt = (salt[,2:dim(salt)[2],,] + salt[,1:(dim(salt)[2] - 1),,]) / 2

    AKt = (AKt[2:dim(AKt)[1],,,] + AKt[1:(dim(AKt)[1] - 1),,,]) / 2
    AKt = (AKt[,2:dim(AKt)[2],,] + AKt[,1:(dim(AKt)[2] - 1),,]) / 2
    AKt = (AKt[,,2:dim(AKt)[3],] + AKt[,,1:(dim(AKt)[3] - 1),]) / 2

    dim(temp) = dims
    dim(rho) = dims
    dim(salt) = dims
    dim(AKt) = dims

    dim(lon) = dims[1]
    dim(lat) = dims[1]
    dim(h) = dims[1]

    #### Calculated
    z = get.zlevel.roms(h, hc, theta, N=N)
    time = get.roms.time(path, TRUE)

    #### Return object
    list(path=path, lat=lat, lon=lon, grid=grid, depth= function(x) {-get.zlevel.roms(x, hc, theta, N=N)},
         h=h, u=u, v=v, w=w, AKt=AKt, T=temp, rho=rho, S=salt, time=time)
}


#' @title Load LTRANS Data
#' @author Thomas Bryce Kelly
#' @export
#' @import ncdf4
load.ltrans = function(path, roms.path, backward = TRUE) {
    cycle = ncdf4::nc_open(path, write=FALSE)
    lat = ncdf4::ncvar_get(cycle, 'lat')
    lon = ncdf4::ncvar_get(cycle, 'lon')
    S = ncdf4::ncvar_get(cycle, 'salinity')
    T = ncdf4::ncvar_get(cycle, 'temperature')
    age = ncdf4::ncvar_get(cycle, 'age')
    dob = ncdf4::ncvar_get(cycle, 'dob')
    depth = ncdf4::ncvar_get(cycle, 'depth')
    model.time = ncdf4::ncvar_get(cycle, 'model_time')
    ncdf4::nc_close(cycle)

    rt = get.roms.time(roms.path, FALSE)

    if (backward) {
        model.time = max(rt) - model.time
        dob = max(rt) - dob
    }

    ## Convert to POSIX
    model.time = conv.time.roms(model.time, tz='GMT')
    dob = conv.time.roms(dob, tz='GMT')
    rt = conv.time.roms(rt, tz='GMT')

    ## Return
    model = list(file = path, lat = lat, lon = lon, sal = S, temp = T, n.times = length(model.time),
                 age = age, dob = dob, depth = depth, time = model.time, n.particles = length(lon[,1]),
                 roms.times = rt)
}


#' @title Get ROMS z-levels
#' @author Thomas Bryce Kelly
get.zlevel.roms = function(h, hc, theta, b=0.6, N) {
    ds = 1 / N
    s = seq(-1 + ds/2, -ds/2, by = ds)
    A = sinh(theta * s) / sinh(theta)
    B = (tanh(theta * (s + 0.5)) - tanh(theta * 0.5)) / (2 * tanh(theta * 0.5))
    C = (1 - b) * A + b * B

    hc * s + (h - hc) * C
}
