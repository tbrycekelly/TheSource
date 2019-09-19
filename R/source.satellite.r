###############################
#### Conversion functions #####
###############################

#' @title Convert Data to SST (Kahru)
#' @author Thomas Bryce Kelly
#' @description
#' @param
#' @export
conv.kahru.sst = function(x) {
    ## Apply fix for signed vs unsigned int values
    if (any(x < 0)) {
        x[which(x < 0)] = x[which(x < 0)] + 256
    }
    ## Remove invalid data
    l = which(x <= 0 | x >= 255)
    x[l] = NA

    0.15 * x - 3
}

#' @title Convert Data to Chl (Kahru)
#' @author Thomas Bryce Kelly
#' @description
#' @param
#' @export
conv.kahru.chl = function(x) {
    ## Apply fix for signed vs unsigned int values
    if (any(x < 0)) {
        x[which(x < 0)] = x[which(x < 0)] + 256
    }
    ## Remove invalid data
    l = which(x <= 0 | x >= 255)
    x[l] = NA

    10 ^ (0.015 * x - 2)
}

#' @title Convert Data to NPP (Kahru)
#' @author Thomas Bryce Kelly
#' @description
#' @param
#' @export
conv.kahru.npp = function(x) {
    x[x < 0] = NA
    x
}

#' @title Convert Data to Chl (NASA)
#' @author Thomas Bryce Kelly
#' @description
#' @param
#' @export
conv.nasa.chl = function(x) {
    10 ^ x
}



###########################
###### File SECTION #######
###########################

#' @title Read Satellite Data
#' @author Thomas Bryce Kelly
#' @description
#' @param
#' @export
## Read nc file and do preliminary parsing/conversion
read.satellite = function(input.dir, file, conv = function(x){x}, entry = 1) {
    nc.file = nc_open(paste0(input.dir, file))

    x = ncvar_get(nc.file, varid = names(nc.file$var)[entry])

    ## Lat and Lon
    l.lat = which(names(nc.file$dim) == 'lat')
    l.lon = which(names(nc.file$dim) == 'lon')

    if (length(l.lat) > 0 & length(l.lon) > 0) {
        lat = ncvar_get(nc.file, varid = 'lat')
        lon = ncvar_get(nc.file, varid = 'lon')
    } else {
        lat = NULL
        lon = NULL
    }

    attr = ncatt_get(nc = nc.file, varid = names(nc.file$var)[entry])

    nc_close(nc.file)

    x = conv(x)

    list(field = x, file = file, dir = input.dir, grid = NULL, lon = lon, lat = lat,
         times = get.satellite.times(file), conv = conv, attr = attr)
}

#' @title Get Satellite Grid (Kahru)
#' @author Thomas Bryce Kelly
#' @description
#' @param
#' @export
get.kahru.grid = function(x) {
    kahru.lat = 0
    kahru.lon = 0

    if (nrow(x$field) == 540) {  # Standard 4k
        x$lat = seq(45, 30.03597, length.out = ncol(x$field))
        x$lon = seq(-140, -115.5454, length.out = nrow(x$field))
        x$grid = expand.grid(lon = x$lon, lat = x$lat)
    }
    if (nrow(x$field) == 588) {  ## Calcofi Fit
        x$lat = seq(37, 29.51349, length.out = ncol(x$field))
        x$lon = seq(-126.125, -116.6828, length.out = nrow(x$field))
        x$grid = expand.grid(lon = x$lon, lat = x$lat)
    }
    if (nrow(x$field) == 3840) {
        load('../Data/Satellite/Kahru Grid 1k.rdata')
        x$lat = grid$Lat
        x$lon = grid$Lon
        x$grid = data.frame(lon = as.numeric(x$lon), lat = as.numeric(x$lat))
    }

    return(x)
}

#' @title Get Satellite Grid (NASA)
#' @author Thomas Bryce Kelly
#' @description
#' @param x
#' @export
get.nasa.grid = function(x, rm.cce = FALSE) {

    kd.lat = seq(90, -90, length.out = 4320)
    kd.lon = seq(-180, 180, length.out = 8640)


    if (rm.cce) {
        l = which(kd.lat >= cce.lat[1] & kd.lat <= cce.lat[2])
        kd.lat = kd.lat[l]
        x$field = x$field[,l]

        l = which(kd.lon >= cce.lon[1] & kd.lon <= cce.lon[2])
        kd.lon = kd.lon[l]
        x$field = x$field[l,]

    }

    x$lons = kd.lon
    x$lats = kd.lat
    x$grid = expand.grid(lon = kd.lon, lat = kd.lat)

    return(x)
}



###########################
###### TIME SECTION #######
###########################


#' @title Get Satellite Times
#' @author Thomas Bryce Kelly
#' @description
#' @param
#' @export
get.satellite.times = function(x) {
    start = rep(conv.time.unix(1), length(x))
    mid = rep(conv.time.unix(1), length(x))
    end = rep(conv.time.unix(1), length(x))

    for (i in 1:length(x)) {

        # Kahru NPP
        if(grepl('_VGPM', x[i])) {
            start[i] = as.POSIXct(as.numeric(substr(x[1], 6, 8))*86400, origin = paste0(substr(x[i], 2, 5), '-01-01'), tz='GMT')
            mid[i] = start[i]
            end[i] = start[i]
        } else{
            if (substr(x[i], 1, 1) == 1 | substr(x[i], 1, 1) == 2) {
                start[i] = as.POSIXct(as.numeric(substr(x[i], 5, 7))*86400,
                                   origin = paste0(substr(x[i], 1, 4), '-01-01'), tz='GMT')
                #end[i] = as.POSIXct(as.numeric(substr(x[i], 12, 14))*86400,
                #                 origin = paste0(substr(x[i], 1, 4), '-01-01'), tz='GMT')
                #mid[i] = as.POSIXct((as.numeric(substr(x[i], 12, 14)) + as.numeric(substr(x[i], 5, 7)))/2*86400,
                #                 origin = paste0(substr(x[i], 1, 4), '-01-01'), tz='GMT')
            }
            else {
                start[i] = as.POSIXct(as.numeric(substr(x[i], 6, 8))*86400,
                                   origin = paste0(substr(x[i], 2, 5), '-01-01'), tz='GMT')
                #end[i] = as.POSIXct(as.numeric(substr(x[i], 13, 15))*86400,
                #                 origin = paste0(substr(x[i], 2, 5), '-01-01'), tz='GMT')
                #mid[i] = as.POSIXct((as.numeric(substr(x[i], 13, 15)) + as.numeric(substr(x[i], 6, 8)))/2*86400,
                #                 origin = paste0(substr(x[i], 2, 5), '-01-01'), tz='GMT')
            }
        }
    }

    #mid[is.na(mid)] = start[is.na(mid)]
    #end[is.na(end)] = start[is.na(end)]

    list(start = start, mid = mid, end = end)
}

#' @title Print NetCDF
#' @import ncdf4
#' @export
print.nc = function(file) {
  file = nc_open(file)
  print(file)
  nc_close(file)
}

#' @title Load NetCDF
#' @import ncdf4
#' @export
load.nc = function(file) {
  file = nc_open(file)
  data = file
  nc_close(file)

  ##return
  data
}



###########################
###### Index SECTION #######
###########################

#' @title
#' @author Thomas Bryce Kelly
#' @description
#' @param
#' @export
rebuild.satellite.index = function(dir = '../../Data_Satellite/', verbose = FALSE) {
    files = list.files(pattern = '.nc', path = paste0(dir, 'Raw/'))

    summary = data.frame(File = files, Sensor = NA, Datetime = NA, Year = NA,
                         Julian.Day = NA, Month = NA, Day = NA, Level = NA,
                         Timeframe = NA, Param = NA, Resolution = NA, stringsAsFactors = FALSE)

    for (i in 1:length(files)) {
        if (verbose) {cat(i)}

        if (any(grepl('UKMO', files[i]))) {
            prime = strsplit(files[i], split = '-')[[1]]

            datetime = substr(prime[1], 1,14)
            summary$Year[i] = as.numeric(substr(datetime, 1, 4))
            summary$Month[i] = as.numeric(substr(datetime, 5, 6))
            summary$Day[i] = as.numeric(substr(datetime, 7, 8))
            summary$Datetime[i] = paste0(as.POSIXct(
                paste0(substr(datetime, 1, 4), '-', substr(datetime, 5, 6), '-', substr(datetime, 7, 8), ' 12:00:00'), tz = 'GMT'
            ))

            summary$Level[i] = 'L4'
            summary$Timeframe[i] = 'DAY'
            summary$Param[i] = 'GHRSST'
            summary$Resolution[i] = '1/20'

            summary$Sensor[i] = 'Composite'

        } else {
            prime = strsplit(files[i], split = '\\.')[[1]]
            second = strsplit(prime[2], split = '_')[[1]]
            datetime = substr(prime[1], 2, 15)

            summary$Sensor[i] = substr(prime[1], 1, 1)

            summary$Year[i] = as.numeric(substr(datetime, 1, 4))
            summary$Julian.Day[i] = as.numeric(substr(datetime, 5, 7))

            summary$Datetime[i] = paste0(as.POSIXct(summary$Julian.Day[i]*86400,
                                                    origin = paste0(summary$Year[i], '-01-01'), tz='GMT'))
            summary$Month[i] = as.numeric(substr(summary$Datetime[i], 6, 7))
            summary$Day[i] = as.numeric(substr(summary$Datetime[i], 9, 10))

            summary$Level[i] = second[1]
            summary$Timeframe[i] = second[2]
            summary$Param[i] = paste(second[3], second[4], sep = '.')
            #summary$Param[i] = second[4]
            summary$Resolution[i] = second[length(second)]
        }

    }
    summary$Datetime = as.POSIXct(summary$Datetime)

    ## SST
    summary$Param[summary$Param == 'SST.sst'] = 'SST'
    summary$Param[summary$Param == 'SNPP.SST'] = 'SST'
    summary$Param[summary$Param == 'SNPP.NSST'] = 'NSST'
    summary$Param[summary$Param == 'NSST.sst'] = 'SST'
    summary$Param[summary$Param == 'SST4.sst4'] = 'SST4'

    ## CHL
    summary$Param[summary$Param == 'SNPP.CHL'] = 'CHL'
    summary$Param[summary$Param == 'CHL.chlor'] = 'CHL'

    ## Other
    summary$Param[summary$Param == 'KD490.Kd'] = 'KD490'
    summary$Param[summary$Param == 'SNPP.KD490'] = 'KD490'
    summary$Param[summary$Param == 'JPSS1.KD490'] = 'KD490'
    summary$Param[summary$Param == 'PAR.par'] = 'PAR'

    if (verbose) {
        print(unique(summary$Year))
        print(unique(summary$Sensor))
        print(unique(summary$Level))
        print(unique(summary$Timeframe))
        print(unique(summary$Param))
        print(unique(summary$Resolution))
    }

    summary
}
