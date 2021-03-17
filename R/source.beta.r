## Author: Thomas Bryce Kelly (tkelly@alaska.edu)
## http://about.tkelly.org/
##
## College of Fisheries and Ocean Sciences
## University of Alaska Fairbanks
##
## Dept of Earth, Ocean & Atmospheric Sciences
## Florida State University


#' @title Load Beta Count Data
#' @author Thomas Bryce Kelly
#' @description A function to read in and parse a raw beta coutner file from a Riso multicollector.
#' @keywords Beta Thorium
#' @return It will return a list containing four objects: (1) _raw_ A list with the raw file contents, (2) meta A list containing the file metadata including instrument serial number and settings, (3) proc The processed and abridge dataset from the file, and (4) counts A dataframe containing the counting data from the file.
#' @param filepath The filepath to the raw count file.
#' @param tz The timezone of the computer that originally saved the datafile.
#' @export
load.beta = function(filepath, tz = 'US/Eastern', lookup = NULL) {
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
                          cummulative = FALSE,
                          file.time = file.time,
                          analysis.time = Sys.time(),
                          file.timezone = tz),
              proc = data.frame(Cruise = rep(NA, 6),
                                Sample.ID = NA,
                                Time.Start = NA,
                                Time.Mid = NA,
                                Time.End = NA,
                                Counts = NA,
                                Guard = NA,
                                Durration = NA,
                                DPM = NA,
                                s.DPM = NA),
              counts = data.frame(Time = NA, Cycle = NA)
  )

  ## Detector line to pull sample names from:
  temp = strsplit(strsplit(raw[3], ':')[[1]], '\\s+')

  if (length(temp) == 6) { names = c(temp[[2]][2], temp[[3]][2], temp[[4]][2], temp[[5]][2], temp[[6]][2], 'Guard')
  } else { names = c('D1', 'D2', 'D3', 'D4', 'D5', 'Guard') }

  ## Fix names
  if (!is.null(lookup)) {
    look = read.xlsx(lookup)
    for (i in 1:length(names)) {
      l = which(names[i] == look[,1])
      if (length(l) > 0) { names[i] = look[l[1],2]}
    }
  }
  names = make.unique(names) ## Keep duplicate names from collapsing number of columns

  for (n in names) { data$counts[[n]] = NA }

  ## Determine start line in case there are notes in the header
  start.line = which(grepl('Date,Time,Cycle', raw))[1] + 1
  end.line = which(grepl('Counter stopped', raw))[1]
  if (is.na(end.line)) {end.line = length(raw)}

  for (i in start.line:end.line) {
    line = strsplit(raw[i], ',')[[1]]
    #if (grepl('stopped', line[1])) { break }
    if (length(line) == 9) { ## Valid line
      time = as.POSIXct(strptime(paste(line[1], line[2]), format = '%d/%m/%y %H:%M'), tz = tz)
      cycle = line[3]

      temp = data.frame(Time = time,
                        Cycle = as.numeric(cycle),
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

  if (data$counts$Guard[nrow(data$counts)] > data$counts$Guard[1] * 2) { data$meta$cummulative = TRUE }

  for (i in 3:ncol(data$counts)) {
    data$proc$Sample.ID[i-2] = colnames(data$counts)[i]
    data$proc$Cruise[i-2] = strsplit(data$proc$Sample.ID[i-2], split = '-|_')[[1]][1]
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
    data$proc$DPM[i - 2] = data$proc$Counts[i-2] / data$proc$Durration[i - 2]
    data$proc$s.DPM[i - 2] = 1 / sqrt(data$proc$Counts[i - 2])
  }

  data$proc$Time.Start = conv.time.unix(data$proc$Time.Start)
  data$proc$Time.Mid = conv.time.unix(data$proc$Time.Mid)
  data$proc$Time.End = conv.time.unix(data$proc$Time.End)

  data.frame(file = data$meta$file,
             cycles = data$meta$cycles,
             preset.time = data$meta$preset.time,
             cummulative = data$meta$cummulative,
             file.time = data$meta$file.time,
             analysis.time = data$meta$analysis.time,
             file.timezone = data$meta$file.timezone,
             interface = data$meta$interface,
             position = c(1:6),
             cruise = data$proc$Cruise,
             sample.id = data$proc$Sample.ID,
             time.start = data$proc$Time.Start,
             time.mid = data$proc$Time.Mid,
             time.end = data$proc$Time.End,
             counts = data$proc$Counts,
             durration = data$proc$Durration,
             DPM = data$proc$DPM,
             sDPM = data$proc$s.DPM,
             stringsAsFactors = FALSE)
}


#' @title Correct Beta Sample IDs
#' @author Thomas Bryce Kelly
#' @export
correct.beta.names = function(summary, lookup.file) {
  lookup = read.xlsx(lookup.file)
  for (i in 1:nrow(summary)) {
    l = which(summary$Sample.ID[i] == lookup[,1])
    if (length(l) == 1) {summary$Sample.ID[i] = lookup[l,2]}
    if (length(l) > 1) {stop('Name occurs more than once in lookup list: ', summary$Sample.ID[i])}
  }
}

#' @title Plot Beta Counts
#' @author Thomas Bryce Kelly
#' @description A helper function to quickly plot up raw data from a beta count data file.
#' @keywords Thorium
#' @param data A beta counter object as returned by _load.beta()_.
#' @param ylim The ylim to use. The default uses the limits of the data.
#' @param cols A set of 5 colors to use in the plot.
#' @param show.legend A boolean flag to turn the legend on and off.
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

#' @title Load Yield Analysis ICPMS Data
#' @author Thomas Bryce Kelly
#' @export
load.beta.icpms = function(file) {
  ## Load Data
  raw = readLines(file)
  raw = sub('yes', '', raw)
  raw = raw[1:(length(raw)-3)]
  head = strsplit(raw[1:13], '\t+')
  rp = as.numeric(strsplit(head[[9]][2],'\\s+')[[1]][c(1,3,5,7,9,11)])


  ## Parse Data
  l = which(grepl('IS before BS', raw))
  body = strsplit(raw[l:length(raw)], '\t+|,+')
  if (length(body[[1]]) > 4) {
    names = body[[1]][c(2,4:length(body[[1]]))]
    names = sub('\\)', '', sub('\\(', '.', names))

    dat = data.frame(Time = rep(NA, length(body) - 10))
    for (n in names) {
      dat[[n]] = NA
    }

    for (i in 11:length(body)) {
      dat$Time[i-10] = as.numeric(body[[i]][1])
      for (k in 2:ncol(dat)) {
        dat[i-10,k] = as.numeric(body[[i]][2*k-1])
      }
    }
    dat = dat[!is.na(dat$Time),]

    proc = data.frame(Element = colnames(dat)[2:ncol(dat)],
                      Mean = apply(dat[,2:ncol(dat)], 2, mean),
                      SEM = apply(dat[,2:ncol(dat)], 2, sd) / sqrt(nrow(dat))
    )
    rownames(proc) = c(1:nrow(proc))
    rownames(dat) = c(1:nrow(dat))
  } else {
    dat = NA
    proc = NA
    warning('File does not contain individual run information: ', file)
  }

  ## initialize data structure
  data = list(file = list(file = file,
                          load.time = Sys.time()
  ),
  meta = list(Data.File = head[[3]][2],
              Error = head[[4]][2],
              Analysis.Date = substr(head[[5]][2], 6, 100),
              Sample.Name = head[[6]][2],
              Tune.Parameters = head[[7]][2],
              Method = head[[8]][2],
              Runs = list(Low = rp[1], Med = rp[3], High = rp[5]),
              Passes = list(Low = rp[2], Med = rp[4], High = rp[6]),
              Switch.Delay = as.numeric(head[[10]][2]),
              Washtime = as.numeric(head[[11]][2]),
              Takeup = as.numeric(head[[12]][2]),
              Deadtime = as.numeric(head[[13]][2])
  ),
  data = dat,
  proc = proc
  )
  data
}
