## Functions for analyzing FRRF Output
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


#' @title Load FRRF Datafiles
#' @author Thomas Bryce Kelly
#' @param input.dir the input directory
#' @param files.names the name of the file in the input directory
#' @export
load.frrf = function(input.dir, file.names) {
    result = list()

    for (file.name in file.names) {

    con = file(paste0(input.dir, file.name), "r")

    col.names = c('Saq', 'E', 'Start', 's', 'Chl', 'ADC', 'rP.obs', 'rP.fit', 'JPII', 'JVPII', 'Fo', 'Fm', 'Fv.Fm',
                  'C', 'p', 'RSigma', 'Sigma', 'CSQ', 'TauES', 'NPQ', 'NSV', 'QR', 'Qo', 'Qm', 'Qo.SE', 'Qm.SE',
                  'Q.SE', 'Q.SE.ratio', 'Qo.points', 'Qo.slope', 'Qo.int', 'Qm.points', 'Qm.slope', 'Qm.int')

    data.a = data.frame(Saq = NA, E = NA, Start = NA, s = NA, Chl = NA, ADC = NA,
                        rP.obs = NA, rP.fit = NA, JPII = NA, JVPII = NA, Fo = NA,
                        Fm = NA, Fv.Fm = NA, C = NA, p = NA, RSigma = NA, Sigma = NA, CSQ = NA,
                        TauES = NA, NPQ = NA, NSV = NA, QR = NA, Qo = NA, Qm = NA,
                        Qo.SE = NA, Qm.SE = NA, Q.SE = NA, Q.SE.ratio = NA, Qo.points = NA,
                        Qo.slope = NA, Qo.int = NA, Qm.points = NA, Qm.slope = NA, QM.int = NA)

    data.b = data.a
    data.c = data.a
    data.d = data.a

    data.s = data.frame(LED = c('A','B','C','D'), Alpha = NA, Ek = NA, Pm = NA, Em = NA, AlphaW = NA, SErP = NA,
                        nPSII = NA, RCII = NA, stringsAsFactors = FALSE)

    parms = data.frame(Ka = NA, Chl.multiplier = NA, QR.threshold = NA, C.derivation = NA, Tau.PQ = NA, Tau.SQ = NA)
    suppressWarnings({
    while (TRUE) {
        ## Load file line by line
        line = readLines(con, n = 1)
        if (length(line) == 0) {
          break
        }
        line = strsplit(line, split = ',')[[1]]

        if (length(line) > 1) {
            ## Check for PARAM values
            if (line[1] == 'FRRf3 Ka:') {
                parms$Ka = as.numeric(line[2])
            }
            else if (line[1] == '[Chl] multiplier:') {
                parms$Chl.multiplier = as.numeric(line[2])
            }
            else if (line[1] == 'QR threshold:') {
                parms$QR.threshold = line[2]
            }
            else if (line[1] == 'C derivation:') {
                parms$C.derivation = line[2]
            }
            else if (line[1] == 'TauPQ:') {
                parms$Tau.PQ = as.numeric(line[2])
            }
            else if (line[1] == 'TauSQ:') {
                parms$Tau.SQ = as.numeric(line[2])
            }


            else if (line[1] == 'Alpha:') {
                data.s$Alpha = as.numeric(line[2:5])
            }
            else if (line[1] == 'Ek:') {
                data.s$Ek = as.numeric(line[2:5])
            }
            else if (line[1] == 'Pm:') {
                data.s$Pm = as.numeric(line[2:5])
            }
            else if (line[1] == 'Em:') {
                data.s$Em = as.numeric(line[2:5])
            }
            else if (line[1] == 'AlphaW:') {
                data.s$AlphaW = as.numeric(line[2:5])
            }
            else if (line[1] == 'SErP:') {
                data.s$SErP = as.numeric(line[2:5])
            }

            ## Data A Values
            else if (length(line) > 1 & line[2] == 'LED combination A (450 nm alone)') {
                line = readLines(con, n = 1) ## Next line
                while (TRUE) {
                    line = readLines(con, n = 1) ## Next line
                    if (length(line) == 0) {
                        break
                    }
                    line = as.numeric(strsplit(line, split = ',')[[1]])

                    if (is.na(line[2])) {
                        break
                    }
                    if (line[2] > 0 & line[2] < 25) {
                        data.a[line[2],] = line[c(2:22,24:30,32:37)]
                    }

                }
            }

            ## Data B Values
            else if (length(line) > 1 & line[2] == 'LED combination B (450 + 530 nm)') {
                line = readLines(con, n = 1) ## Next line
                while (TRUE) {
                    line = readLines(con, n = 1) ## Next line
                    if (length(line) == 0) {
                        break
                    }
                    line = as.numeric(strsplit(line, split = ',')[[1]])

                    if (is.na(line[2])) {
                        break
                    }
                    if (length(line) > 1 & line[2] > 0 & line[2] < 25) {
                        data.b[line[2],] = line[c(2:22,24:30,32:37)]
                    }

                }
            }

            ## Data C Values
            else if (length(line) > 1 & line[2] == 'LED combination C (450 + 624 nm)') {
                line = readLines(con, n = 1) ## Next line
                while (TRUE) {
                    line = readLines(con, n = 1) ## Next line
                    if (length(line) == 0) {
                        break
                    }
                    line = as.numeric(strsplit(line, split = ',')[[1]])

                    if (is.na(line[2])) {
                        break
                    }
                    if (line[2] > 0 & line[2] < 25) {
                        data.c[line[2],] = line[c(2:22,24:30,32:37)]
                    }

                }
            }

            ## Data D Values
            else if (length(line) > 1 & line[2] == 'LED combination D (450 + 530 + 624 nm)') {
                line = readLines(con, n = 1) ## Next line
                while (TRUE) {
                    line = readLines(con, n = 1) ## Next line
                    if (length(line) == 0) {
                        break
                    }
                    line = as.numeric(strsplit(line, split = ',')[[1]])

                    if (length(line) > 1 & is.na(line[2])) {
                        break
                    }
                    if (length(line) > 1 & line[2] > 0 & line[2] < 25) {
                        data.d[line[2],] = line[c(2:22,24:30,32:37)]
                    }

                }
            }

        }
    }
    })

    close(con)


    ## Generate Datetime stamp:
    year = substr(file.name, 1, 4)
    month = substr(file.name, 5, 6)
    day = substr(file.name, 7, 8)
    hour = substr(file.name, 10, 11)
    minute = substr(file.name, 12, 13)
    second = substr(file.name, 14, 15)

    datetime = as.POSIXct(paste0(year, '-', month, '-', day, ' ', hour, ':', minute, ':', second))

    result[[paste0('Y', year, month, day, hour, minute, second)]] = list(Params = parms, S = data.s,
         A = data.a, B = data.b, C = data.c, D = data.d,
         Datetime = datetime, File = paste0(input.dir, file.name))
    }
    ## Return
    result
}
