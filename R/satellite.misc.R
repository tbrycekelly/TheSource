
#' @title Rebuild Satellite Index
#' @author Thomas Bryce Kelly
#' @description rebuild
#' @param verbose Flag to print diagnostics
#' @export
rebuild.satellite.index = function(dir = 'Z:/Data/Data_Satellite/Raw', verbose = FALSE) {

  #### Find files
  if (verbose) { message(Sys.time(), ': Looking for available satellite files...')}

  # List files
  a = Sys.time()
  files = list.files(pattern = '.nc', path = dir, recursive = F)
  b = Sys.time()

  if (verbose) { message(Sys.time(), ': Found ', length(files), ' files! Search took ', round(as.numeric(b) - as.numeric(a)), ' seconds.')}


  #### Populate data fields
  summary = data.frame(File = c(files, NA), Sensor = NA, Datetime = NA, Year = NA,
                       Julian.Day = NA, Month = NA, Day = NA, Level = NA,
                       Timeframe = NA, Param = NA, Resolution = NA, stringsAsFactors = FALSE)

  if (length(files) > 0) {
    for (i in 1:length(files)) {
      if (file.size(paste0(dir, files[i])) > 0) {
        if (verbose) { message(Sys.time(), ': Attempting to scrub metadata from ', files[i], ' (', i, '/', length(files), ')... ', appendLF = F) }

        summary$File[i] = paste0(dir, files[i])

        prime = strsplit(files[i], split = '\\.')[[1]]
        second = strsplit(prime[2], split = '_')[[1]]
        datetime = substr(prime[1], 2, 8)

        summary$Sensor[i] = substr(prime[1], 1, 1)

        summary$Year[i] = as.numeric(substr(datetime, 1, 4))
        summary$Julian.Day[i] = as.numeric(substr(datetime, 5, 7))

        summary$Datetime[i] = paste0(as.POSIXct(summary$Julian.Day[i]*86400,
                                                origin = paste0(summary$Year[i], '-01-01'), tz = 'GMT'))
        summary$Month[i] = as.numeric(substr(summary$Datetime[i], 6, 7))
        summary$Day[i] = as.numeric(substr(summary$Datetime[i], 9, 10))

        summary$Level[i] = second[1]
        summary$Timeframe[i] = second[2]
        summary$Resolution[i] = second[length(second)]
        if (length(second) > 4) {
          summary$Param[i] = paste(second[-c(1:2)], collapse = '.')
        } else {
          summary$Param[i] = second[3]
        }

        if (verbose) { message('Success.') }
      }
    }
  }

  ## Remove bad files or folders
  summary = summary[!is.na(summary$Param),]

  c = Sys.time()

  #### Finish up
  summary$Datetime = as.POSIXct(summary$Datetime, tz = 'GMT')

  ## Recursively delve into subfolders:
  dirs = list.dirs(path = dir, recursive = F, full.names = T)
  if (length(dirs) > 0) {
    for (i in 1:length(dirs)) {
      if (verbose) { message(' Delving into ', dirs[i])}
      summary = rbind(summary, rebuild.satellite.index(dir = paste0(dirs[i], '/'), verbose = verbose))
    }
  }
  ## Return
  summary
}


###############################
#### Conversion functions #####
###############################

## Depreciated

#' @title Convert Data to SST (Kahru)
#' @author Thomas Bryce Kelly
#' @description Convert
#' @param x The pixel values (or a vector or matrix)
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
#' @description Convert Kahru's data to mg Chl m-2
#' @param x Raw value
#' @export
conv.kahru.chl = function(x) {
  ## Apply fix for signed vs unsigned int values
  if (any(x < 0)) {
    x[which(x < 0)] = x[which(x < 0)] + 256
  }
  ## Remove invalid data
  l = which(x <= 0 | x >= 255)
  x[l] = NA

  10^(0.015 * x - 2)
}


#' @title Convert Data to NPP (Kahru)
#' @author Thomas Bryce Kelly
#' @description Converts kahru's data to npp
#' @param x raw value
#' @export
conv.kahru.npp = function(x) {
  x[x < 0] = NA
  x
}


#' @title Convert Data to Chl (NASA)
#' @author Thomas Bryce Kelly
#' @description convert
#' @param x raw value
#' @export
conv.nasa.chl = function(x) { 10 ^ x }

