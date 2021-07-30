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



##########################
##### LOADING FUNCTIONS
##########################

#' @title Load in ROMS physics
#' @export
#' @author Thomas Bryce Kelly
#' @import ncdf4
load.roms.physics = function(path, verbose = T) {

    model = load.nc(file = path,
                    var = c('u', 'Huon', 'v', 'Hvom', 'w', 'AKt', 'Akt_bak',
                            'temp', 'rho', 'salt', 'hice', 'f', 'nl_tnu2', 'snow_thick',
                            'uice', 'vice', 'zeta'),
                    verbose = verbose)

    ## Return
    model
}

#' @title Load in ROMS grid
#' @export
#' @author Thomas Bryce Kelly
#' @import ncdf4
load.roms.grid = function(path, verbose = T) {

    model = load.nc(file = path,
                    var = c('h', 'hc', 'theta_s', 'theta_b',
                            'lat_psi', 'x_psi', 'lon_psi', 'y_psi',
                            'zice', 'Cs_r', 'Cs_w',
                            'mask_psi', 'mask_rho', 'mask_u', 'mask_v',
                            'pm', 'pn', 's_rho', 's_w',
                            'x_psi', 'x_rho', 'x_u', 'x_v',
                            'y_psi', 'y_rho', 'y_u', 'y_v',
                            'Vtransform', 'Vstretching'),
                    verbose =  verbose)
    ## Return
    model
}



#' @title Load in ROMS temporal info
#' @export
#' @author Thomas Bryce Kelly
#' @import ncdf4
load.roms.time = function(path, verbose = T) {

    model = load.nc(file = path,
                    var = c('dstart', 'dt', 'dtfast',
                            'nAVG', 'ndefAVG', 'ndefHIS', 'ndtfast',
                            'nHIS', 'nRST', 'ntimes', 'ntsAVG', 'ocean_time'),
                    verbose =  verbose)
    ## Return
    model
}



###########################
#### CONVERSION FUNCTIONS
###########################

#' @title Convert ROMS Time
#' @author Thomas Bryce Kelly
#' @param x The datetime numeric to be converted
#' @param tz The timezone for the conversion, almost always UTC
#' @export
conv.time.roms = function(x, origin = "1900-01-01", tz = 'UTC') {
    as.POSIXct(x, origin = origin, tz = tz)
}


#' @title Get ROMS Time
#' @export
#' @author Thomas Bryce Kelly
calc.roms.time = function(time.data, tz = 'UTC') {

    ## If we have ocean_time, then use it, if not try and calculate it based on time steps.
    if (!is.na(time.data$ocean_time[1])) {
        time = conv.time.roms(time.data$ocean_time)
    } else {
        nhist = time.data$nHIS # dt * nHist = time between records
        dstart = time.data$dstart # days since 1900-01-01 00:00:00
        dt = time.data$dt # Sec
        ntimes = time.data$ntimes # number of time frames


        ## Calculate the times (seconds from 1900-01-01)
        time = seq(from = dstart, by = dt * nhist / 86400, length.out = ntimes/nhist) * 86400
    }

    time
}


#' @export
destagger.grid = function(x) {
    nx = dim(x)[1]
    ny = dim(x)[2]

    if (length(dim(x)) == 2) {
        x= (x[1:(nx-1), 1:(ny-1)] + x[2:nx, 1:(ny-1)] + x[1:(nx-1), 2:ny] + x[2:nx, 2:ny]) / 4
    } else if (length(dim(x)) == 3) {
        x = (x[1:(nx-1), 1:(ny-1),] + x[2:nx, 1:(ny-1),] + x[1:(nx-1), 2:ny,] + x[2:nx, 2:ny,]) / 4
    } else if (length(dim(x)) == 4) {
        x = (x[1:(nx-1), 1:(ny-1),,] + x[2:nx, 1:(ny-1),,] + x[1:(nx-1), 2:ny,,] + x[2:nx, 2:ny,,]) / 4
    } else if (length(dim(x)) == 4) {
        x = (x[1:(nx-1), 1:(ny-1),,,] + x[2:nx, 1:(ny-1),,,] + x[1:(nx-1), 2:ny,,,] + x[2:nx, 2:ny,,,]) / 4
    }

    ## Return x
    x
}

#' @export
calc.roms.depths = function(grid.data, zeta = NULL, verbose = T) {

    ## Setup
    nx = dim(grid.data$h)[1]
    ny = dim(grid.data$h)[2]
    h = destagger.grid(grid.data$h)
    if (!is.na(grid.data$zice[1])) {
        zice = destagger.grid(grid.data$zice)
    } else {
        zice = matrix(0, nrow = nx - 1, ncol = ny - 1)
    }


    ## Free surface
    if (is.null(zeta)) {
        zeta = matrix(0, nrow = nx - 1, ncol = ny - 1)
    } else {
        if (length(dim(zeta)) > 2) {
            message(' Zeta has dimensions greater than 2. Ignoring zeta...')
            zeta = matrix(0, nrow = nx - 1, ncol = ny - 1)
        } else {
            zeta = destagger.grid(zeta)
        }
    }

    ## Setup depths and calculate
    z = array(0, dim = c(nx-1, ny-1, length(grid.data$Cs_r)))

    if (!is.na(grid.data$Vtransform) & grid.data$Vtransform == 1) { ## Equation 1
        for (i in 1:length(grid.data$Cs_r)) {
            z0 = (grid.data$s_rho[i] - grid.data$Cs_r[i]) * grid.data$hc + grid.data$Cs_r[i] * (h + zice)
            z[,,i] = z0 + zeta * (1 + z0 / (h + zice)) + zice
        }
    } else if (!is.na(grid.data$Vtransform) & grid.data$Vtransform == 2) { ## Equation 2
        for (i in 1:length(grid.data$Cs_r)) {
            z0 = (grid.data$hc * grid.data$s_rho[i] + grid.data$Cs_r[i]* (h + zice)) / (grid.data$hc + h + zice)
            z[,,i] = zeta + (zeta + h + zice) * z0
        }
    } else {
        message(' No valid Vtransform given. Returning zeros!')
    }

    z
}


#' @title Get ROMS z-levels
#' @author Thomas Bryce Kelly
calc.roms.depths.old = function(h, hc, theta, b = 0.6, N) {
    z = array(NA, dim = c(dim(h), N))
    ds = 1 / N
    s = seq(-1 + ds/2, -ds/2, by = ds)
    A = sinh(theta * s) / sinh(theta)
    B = (tanh(theta * (s + 0.5)) - tanh(theta * 0.5)) / (2 * tanh(theta * 0.5))
    C = (1 - b) * A + b * B

    for (i in 1:dim(h)[1]) {
        for (j in 1:dim(h)[2]) {
            z[i,j,] = hc * s + (h[i,j] - hc) * C
        }
    }
    z
}



#####################
#### OLD STUFF
#####################

#' @title Get ROMS Time old
#' @export
#' @author Thomas Bryce Kelly
#' @import ncdf4
get.roms.time.old = function(path, convert = TRUE) {

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


#' @export
load.roms.old = function(path, verbose = T) {
    if (verbose) { message(Sys.time(), ': Attempting to load ROMS file ', path, appendLF = F)}
    file = ncdf4::nc_open(path, write=FALSE)
    if (verbose) { message(' success.')}

    vars = names(file$var)

    #### Load variables
    ## Vertical Grid
    if ('h' %in% vars) {
        if (verbose) { message(' Loading h...')}
        h = ncdf4::ncvar_get(file, 'h')
    } else {
        stop('Error, no h variable found.')
    }

    if ('hc' %in% vars) {
        if (verbose) { message(' Loading hc...')}
        hc = ncdf4::ncvar_get(file, 'hc')
    } else {
        stop('Error, no hc variable found.')
    }

    if ('theta_s' %in% vars) {
        if (verbose) { message(' Loading theta.s...')}
        theta.s = ncdf4::ncvar_get(file, 'theta_s')
    } else {
        stop('Error, no theta_s variable found.')
    }

    if ('theta_b' %in% vars) {
        if (verbose) { message(' Loading theta.b...')}
        theta.b = ncdf4::ncvar_get(file, 'theta_b')
    } else {
        message('Warning, no theta_b variable found, setting to zero.')
        theta.b = 0
    }

    h = (h[2:dim(h)[1],] + h[1:(dim(h)[1] - 1),]) / 2 ## interpolate to center of each grid cell
    h = (h[,2:dim(h)[2]] + h[,1:(dim(h)[2] - 1)]) / 2

    ## Horizontal
    if ('lat_psi' %in% vars) {
        if (verbose) { message(' Loading lat_psi...')}
        lat = ncdf4::ncvar_get(file, 'lat_psi')

        if (verbose) { message(' Loading lon_psi...')}
        lon = ncdf4::ncvar_get(file, 'lon_psi')
    } else if ('x_psi' %in% vars) {
        if (verbose) { message(' Loading x_psi...')}
        lon = ncdf4::ncvar_get(file, 'x_psi')
        if (verbose) { message(' Loading y_psi...')}
        lat = ncdf4::ncvar_get(file, 'y_psi')
    } else {
        message('No positional coordiantes found!')
    }

    grid = list(lon = lon, lat = lat)
    if (verbose) { message(' Formed grid ', length(grid$lon), 'x', length(grid$lat))}

    ## Advect
    if ('u' %in% vars) {
        if (verbose) { message(' Loading u... ', appendLF = F); a = Sys.time()}
        u = ncdf4::ncvar_get(file, 'u')
        if (verbose) { message('(',Sys.time() - a, ')')}
        u = (u[,2:dim(u)[2],,] + u[,1:(dim(u)[2] - 1),,]) / 2
    } else if ('Huon' %in% vars) {
        if (verbose) { message(' Loading Huon... ', appendLF = F); a = Sys.time()}
        u = ncdf4::ncvar_get(file, 'Huon')
        if (verbose) { message('(',Sys.time() - a, ')')}
        u = (u[,2:dim(u)[2],,] + u[,1:(dim(u)[2] - 1),,]) / 2
    } else {
        stop('Error, no u velocity found.')
    }

    if ('v' %in% vars) {
        if (verbose) { message(' Loading v... ', appendLF = F); a = Sys.time()}
        v = ncdf4::ncvar_get(file, 'v')
        if (verbose) { message('(',Sys.time() - a, ')')}
        v = (v[2:dim(v)[1],,,] + v[1:(dim(v)[1] - 1),,,]) / 2
    } else if ('Hvom' %in% vars) {
        if (verbose) { message(' Loading Hvom... ', appendLF = F); a = Sys.time()}
        v = ncdf4::ncvar_get(file, 'Hvom')
        if (verbose) { message('(',Sys.time() - a, ')')}
        v = (v[2:dim(v)[1],,,] + v[1:(dim(v)[1] - 1),,,]) / 2
    } else {
        stop('Error, no v velocity found.')
    }

    if ('w' %in% vars) {
        if (verbose) { message(' Loading w... ', appendLF = F); a = Sys.time()}
        w = ncdf4::ncvar_get(file, 'w')
        if (verbose) { message('(',Sys.time() - a, ')')}
        w = (w[,,2:dim(w)[3],] + w[,,1:(dim(w)[3] - 1),]) / 2 # depth averaged w
        w = (w[2:dim(w)[1],,,] + w[1:(dim(w)[1] - 1),,,]) / 2
        w = (w[,2:dim(w)[2],,] + w[,1:(dim(w)[2] - 1),,]) / 2
    } else {
        stop('Error, no w velocity found.')
    }

    if (verbose) { message(' Unifying grid layout...')}
    dims = c(dim(u)[1] * dim(u)[2], dim(u)[3], dim(u)[4])
    dim(u) = dims
    dim(v) = dims
    dim(w) = dims
    dim(lon) = dims[1]
    dim(lat) = dims[1]
    dim(h) = dims[1]
    N = dims[2]

    ## Diffusivity
    if ('AKt' %in% vars) {
        if (verbose) { message(' Loading AKt...')}
        AKt = ncdf4::ncvar_get(file, 'AKt')
    } else {
        AKt = NULL
    }

    ## Fields
    if ('temp' %in% vars) {
        if (verbose) { message(' Loading pot temp...')}
        temp = ncdf4::ncvar_get(file, 'temp')
    } else {
        temp = NULL
    }

    if ('rho' %in% vars) {
        if (verbose) { message(' Loading rho...')}
        rho = ncvar_get(file, 'rho')
    } else {
        rho = NULL
    }
    if ('salt' %in% vars) {
        if (verbose) { message(' Loading salinity...')}
        salt = ncvar_get(file, 'salt')
    } else {
        salt = NULL
    }

    if (verbose) { message(' Adjusting tracer grid...')}
    if (!is.null(temp)) {
        temp = (temp[2:dim(temp)[1],,,] + temp[1:(dim(temp)[1] - 1),,,]) / 2
        temp = (temp[,2:dim(temp)[2],,] + temp[,1:(dim(temp)[2] - 1),,]) / 2
        dim(temp) = dims
    }

    if (!is.null(rho)) {
        rho = (rho[2:dim(rho)[1],,,] + rho[1:(dim(rho)[1] - 1),,,]) / 2
        rho = (rho[,2:dim(rho)[2],,] + rho[,1:(dim(rho)[2] - 1),,]) / 2
        dim(rho) = dims
    }

    if (!is.null(salt)) {
        salt = (salt[2:dim(salt)[1],,,] + salt[1:(dim(salt)[1] - 1),,,]) / 2
        salt = (salt[,2:dim(salt)[2],,] + salt[,1:(dim(salt)[2] - 1),,]) / 2
        dim(salt) = dims
    }

    if (!is.null(AKt)) {
        AKt = (AKt[2:dim(AKt)[1],,,] + AKt[1:(dim(AKt)[1] - 1),,,]) / 2
        AKt = (AKt[,2:dim(AKt)[2],,] + AKt[,1:(dim(AKt)[2] - 1),,]) / 2
        AKt = (AKt[,,2:dim(AKt)[3],] + AKt[,,1:(dim(AKt)[3] - 1),]) / 2
        dim(AKt) = dims
    }

    #### Calculated
    z = get.zlevel.roms(h, hc, theta.s, N=N)
    time = get.roms.time(path, TRUE)
    meta = list(path = path, theta.b = theta.b, theta.s = theta.s, hc = hc, N = N, dims = dims)

    #### Return object
    list(lat = lat,
         lon = lon,
         grid = grid,
         depth = function(x) {-get.zlevel.roms(x, hc, theta.s, N = N)},
         z = z,
         h = h,
         u = u,
         v = v,
         w = w,
         AKt = AKt,
         T = temp,
         rho = rho,
         S = salt,
         time = time,
         meta = meta)
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
    list(file = path, lat = lat, lon = lon, sal = S, temp = T, n.times = length(model.time),
                 age = age, dob = dob, depth = depth, time = model.time, n.particles = length(lon[,1]),
                 roms.times = rt)
}



