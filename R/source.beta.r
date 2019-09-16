#' @title Load Beta Count Data
#' @author Thomas Bryce Kelly
#' @description A function to read in and parse a raw beta coutner file from a Riso multicollector.
#' @keywords Beta Thorium
#' @return It will return a list containing four objects: (1) _raw_ A list with the raw file contents, (2) meta A list containing the file metadata including instrument serial number and settings, (3) proc The processed and abridge dataset from the file, and (4) counts A dataframe containing the counting data from the file.
#' @params filepath The filepath to the raw count file.
#' @params tz The timezone of the computer that originally saved the datafile.
#' @export
load.beta = function(filepath, tz = 'US/Atlantic') {
  con = file(filepath, "r")
  raw = readLines(con) ## Read in raw file
  close(con)

  interface = substr(strsplit(raw[1], 'Interface ')[[1]][2], 1, 7)
  temp = strsplit(raw[2], '\\s+')[[1]]
  file.time = as.POSIXct(strptime(paste(temp[1], temp[2]), format = '%d/%m/%y %H:%M'), tz = tz)
  preset.time = as.numeric(temp[5])
  cycles = as.numeric(temp[10])


  data = list(raw = raw,
              meta = list(file = filepath,
                          cycles = cycles,
                          preset.time = preset.time,
                          interface = interface,
                          cummulative = NA,
                          file.time = file.time,
                          analysis.time = Sys.time()),
              proc = data.frame(Sample.ID = rep(NA, 6),
                                Time.Start = NA,
                                Time.Mid = NA,
                                Time.End = NA,
                                Counts = NA,
                                Durration = NA,
                                DPM = NA,
                                s.DPM = NA),
              counts = data.frame(Time = NA, Cycle = NA)
  )

  temp = strsplit(strsplit(raw[3], ':')[[1]], '\\s+')

  if (length(temp) == 6) { names = c(temp[[2]][2], temp[[3]][2], temp[[4]][2], temp[[5]][2], temp[[6]][2], 'Guard')
  } else { names = c('D1', 'D2', 'D3', 'D4', 'D5', 'Guard') }

  for (n in names) { data$counts[[n]] = NA }

  for (i in 7:length(raw)) {
    line = strsplit(raw[i], ',')[[1]]
    if (grepl('stopped', line[1])) { break }
    if (length(line) == 9) { ## Valid line
      time = as.POSIXct(strptime(paste(line[1], line[2]), format = '%d/%m/%y %H:%M'), tz = tz)
      cycle = line[3]

      temp = data.frame(Time = time, Cycle = as.numeric(cycle),
                        D1 = as.numeric(line[4]),
                        D2 = as.numeric(line[5]),
                        D3 = as.numeric(line[6]),
                        D4 = as.numeric(line[7]),
                        D5 = as.numeric(line[8]),
                        Guard = as.numeric(line[9]))
      colnames(temp) = colnames(data$counts)
      data$counts = rbind(data$counts, temp)
    }
  }
  data$counts = data$counts[-1,]
  data$counts$Time = conv.time.unix(data$counts$Time)

  if (data$counts[nrow(data$counts),3] > mean(data$counts[1,3]) * 3) { data$meta$cummulative = TRUE }

  for (i in 3:ncol(data$counts)) {
    data$proc$Sample.ID[i-2] = colnames(data$counts)[i]
    ## Time
    data$proc$Time.Start[i-2] = data$counts$Time[1] - 60 * data$meta$preset.time
    data$proc$Time.End[i-2] = data$counts$Time[nrow(data$counts)]
    data$proc$Time.Mid[i-2] = mean(data$counts$Time)
    ## Add decays
    if (data$meta$cummulative) { data$proc$Counts[i-2] = data$counts[nrow(data$counts),i] }
    else { data$proc$Counts[i-2] = sum(data$counts[,i]) }
    ## Durration in minutes
    data$proc$Durration[i-2] = data$meta$preset.time * (sum(diff(data$counts$Cycle)) + 1)
    ## DPM
    data$proc$DPM[i-2] = data$proc$Counts[i-2] / data$proc$Durration[i-2]
    data$proc$s.DPM[i-2] = 1 / sqrt(data$proc$Counts[i-2])
  }

  data$proc$Time.Start = conv.time.unix(data$proc$Time.Start)
  data$proc$Time.Mid = conv.time.unix(data$proc$Time.Mid)
  data$proc$Time.End = conv.time.unix(data$proc$Time.End)

  data
}

#' @title Plot Beta Counts
#' @author Thomas Bryce Kelly
#' @description A helper function to quickly plot up raw data from a beta count data file.
#' @keywords Thorium
#' @params data A beta counter object as returned by _load.beta()_.
#' @params ylim The ylim to use. The default uses the limits of the data.
#' @params cols A set of 5 colors to use in the plot.
#' @params show.legend A boolean flag to turn the legend on and off.
#' @export
plot.beta.counts = function(data, ylim = NULL, cols = pals::glasbey(5), show.legend = TRUE) {

  if (is.null(ylim)) { ylim = c(0, max(data$counts[,3:7])) }

  plot(NULL, NULL, ylim = ylim, xlim = c(1, max(data$counts$Cycle)), xaxs='i', yaxs='i',
       main=data$meta$file, xlab = "Cycle", ylab = "Counts")

  lines(data$counts$Cycle, data$counts[,3], col=cols[1], lwd=2)
  lines(data$counts$Cycle, data$counts[,4], col=cols[2], lwd=2)
  lines(data$counts$Cycle, data$counts[,5], col=cols[3], lwd=2)
  lines(data$counts$Cycle, data$counts[,6], col=cols[4], lwd=2)
  lines(data$counts$Cycle, data$counts[,7], col=cols[5], lwd=2)

  if (show.legend) {
    legend(1.2, max(data$counts[,3:7], na.rm = TRUE) * 0.9, colnames(data$counts)[3:7], col = cols, lwd = 2, cex = 0.7)
  }
  mtext(paste('Cummulative:', data$meta$cummulative), adj = 1)
}

