## Set of useful R functions for the calculation of some physical ocean properties
## or other physically oriented features
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


#####################################
#####################################
###### Equations of State  ##########
#####################################
#####################################


#' @title Calculate Argon Solubility
#' @author Thomas Bryce Kelly
#' @description Calculate Argon solubility in seawater. Should return a value of 13.46 in 10 degree water in 35 PSU water.
#' @param S Salinity in PSU
#' @param T Temperature
#' @export
calc.ArSol = function(S = 35, T = 10) {
  ## Check value at 35PSU, 10C
  #### 13.46

  # convert T to scaled temperature
  Ts = log((298.15 - T) / (273.15 + T))

  A.0 = 2.79150
  A.1 = 3.17609
  A.2 = 4.13116
  A.3 = 4.90379
  B.0 = -6.96233e-3
  B.1 = -7.66670e-3
  B.2 = -1.16888e-2

  ## Calcualte ln(C)
  lnC = A.0 + A.1 * Ts + A.2 * Ts^2 + A.3 * Ts^3 + S * (B.0 + B.1 * Ts + B.2 * Ts^2)

  exp(lnC) # umol/kg
}

#' @title Calculate N2 Solubility
#' @author Thomas Bryce Kelly
#' @description Calculate the N2 solubility in seawater.
#' @param S Salinity (numeric or a vector) in PSU
#' @param T Temperature (numeric or a vector) in centigrade
#' @export
calc.N2Sol = function(S = 35, T = 10) {
  ## Check value at 35PSU, 10C
  #### 500.885

  # convert T to scaled temperature
  Ts = log((298.15 - T) / (273.15 + T))

  A.0 = 6.42931
  A.1 = 2.92704
  A.2 = 4.32531
  A.3 = 4.69149
  B.0 = -7.44129e-3
  B.1 = -8.02566e-3
  B.2 = -1.46775e-2

  ## Calcualte ln(C)
  lnC = A.0 + A.1 * Ts + A.2 * Ts^2 + A.3 * Ts^3 + S * (B.0 + B.1 * Ts + B.2 * Ts^2)

  exp(lnC) # umol/kg
}


#' @title Calculate Neon Solubility
#' @author Thomas Bryce Kelly
#' @description Calculates the Neon solubility in seawater based on salinity and temperature
#' @references Hamme & Emerseon 2004 (Table 4)
#' @param S Salinity in PSU
#' @param T Temperature in centigrade
#' @export
calc.NeSol = function(S = 35, T = 10) {
  ## Check value at 35PSU, 10C
  #### 7.34121e-3

  # convert T to scaled temperature
  Ts = log((298.15 - T) / (273.15 + T))

  ## Coefficients for Ar from Hamme & Emerseon 2004 (Table 4)
  A.0 = 2.18156
  A.1 = 1.29108
  A.2 = 2.12504
  A.3 = 0
  B.0 = -5.94737e-3
  B.1 = -5.13896e-3
  B.2 = 0

  ## Calcualte ln(C)
  lnC = A.0 + A.1 * Ts + A.2 * Ts^2 + A.3 * Ts^3 + S * (B.0 + B.1 * Ts + B.2 * Ts^2) # nmol/kg

  exp(lnC) * 1e-3 # umol/kg
}


#' @title Calcualte O2 Solubility
#' @author Thomas Bryce Kelly
#' @param S Salinity in PSU
#' @param T Temperature in centigrade
#' @export
calc.O2sol = function(S = 35, T = 10, verbose = F) { ## umol Kg-1
    ## Check value at 35PSU, 10C
    #### 274.61

    # convert T to scaled temperature
    Ts = log((298.15 - T) / (273.15 + T))

    # constants from Table 1 of Garcia & Gordon for the fit to Benson and Krause (1984)
    A.0 = 5.80871
    A.1 = 3.20291
    A.2 = 4.17887
    A.3 = 5.10006
    A.4 = -9.86643e-2
    A.5 = 3.80369

    B.0 = -7.01577e-3
    B.1 = -7.70028e-3
    B.2 = -1.13864e-2
    B.3 = -9.51519e-3

    C.0 = -2.75915e-7

    A.calc = A.0 + A.1*Ts + A.2*Ts^2 + A.3*Ts^3 + A.4*Ts^4 + A.5*Ts^5
    B.calc = B.0 + B.1*Ts + B.2*Ts^2 + B.3*Ts^3

    if (verbose) { message(Sys.time(), ': Calculating O2 Solubility (umol O2 kg-1) with check value (S=35, T=10) of 274.61') }

    ## (umol / kg)
    exp(A.calc + S * B.calc + C.0 * S^2)
}


#' @title Calculate Density
#' @author Thomas Bryce Kelly
#' @param S Salinity in PSU
#' @param T Temperature in centigrade, provide potential temperature for potential density.
#' @param p Pressure in db, set p=0 for potential density
#' @references Massel 2015
#' @export
calc.rho = function(S = 30, T = 15, P = 0, verbose = F) {
    a.0 = 999.842594
    a.1 = 6.793953e-2
    a.2 = -9.095290e-3
    a.3 = 1.001685e-4
    a.4 = -1.120083e-6
    a.5 = 6.536332e-9
    rho.smow = a.0 + a.1*T + a.2*T^2 + a.3*T^3 + a.4*T^4 + a.5*T^5

    b.0 = 8.2449e-1
    b.1 = -4.0899e-3
    b.2 = 7.6438e-5
    b.3 = -8.2467e-7
    b.4 = 5.3875e-9
    B.1 = b.0 + b.1*T + b.2*T^2 + b.3*T^3 + b.4*T^4

    c.0 = -5.7246e-3
    c.1 = 1.0227e-4
    c.2 = -1.6546e-6
    d.0 = 4.8314e-4
    C.1 = c.0 + c.1*T + c.2*T^2

    if (verbose) { message(Sys.time(), ': Calculating seawater density based on in situ T, S and pressure. For potential density set p = 0 regardless of in situ pressure and provide potential temperature. Value returned is in kg m-3 with paramterization given in Massel 2015.') }

    rho = rho.smow + B.1 * S + C.1 * S^1.5 + d.0 * S^2
    rho / (1 - 0.1 * P / calc.seawater.compressibility(S, T, P/10))
}


#' @titel Calculate Seawater Compressibility (module)
#' @author Thomas Bryce Kelly
#' @param P Pressure in bar
#' @references https://link.springer.com/content/pdf/bbm%3A978-3-319-18908-6%2F1.pdf
calc.seawater.compressibility = function(S = 8, T = 10, P = 10, verbose = F) {
  e0 = 19652.210000
  e1 = 148.420600
  e2 = -2.327105
  e3 = 1.360477e-2
  e4 = -5.155288e-5

  Kw = e0 + e1 * T + e2 * T^2 + e3 * T^3 + e4 * T^4

  f0 = 54.674600
  f1 = -0.603459
  f2 = 1.099870e-2
  f3 = -6.167000e-5

  F1 = f0 + f1 * T + f2 * T^2 + f3 * T^3

  g0 = 7.9440e-2
  g1 = 1.6483e-2
  g2 = -5.3009e-4

  G1 = g0 + g1 * T + g2 * T^2

  ## Calculate module at surface
  K0 = Kw + F1 * S + G1 * S^1.5

  #### Pressure Correction
  h0 = 3.23990
  h1 = 1.43713e-3
  h2 = 1.16092e-4
  h3 = -5.77905e-7
  i0 = 2.28380e-3
  i1 = -1.09810e-5
  i2 = -1.6078e-6
  j0 = 1.91075e-4

  Aw = h0 + h1 * T + h2 * T^2 + h3 * T^3
  A1 = Aw + (i0 + i1 * T + i2 * T^2) * S + j0 * S^1.5

  k0 = 8.50935e-5
  k1 = -6.12293e-6
  k2 = 5.27870e-8

  m0 = -9.9348e-7
  m1 = 2.0816e-8
  m2 = 9.1697e-10

  Bw = k0 + k1*T + k2*T^2
  B2 = Bw + (m0 + m1*T + m2*T^2) * S

  ## Module at pressure P
  K = K0 + A1 * P + B2 * P^2

  if (verbose) {
    message('Module of compressibility for Seawater calcualted from UNESCO 1981. Check value for 8PSU, 10C, 0bar = 21351.408.')
    message('\tKw = ', Kw)
    message('\tA1 = ', A1)
    message('\tB2 = ', B2)
    message('\tK = ', K)
  }

  K
}

#' @title Calculate Adiabatic Temperature Gradient
#' @author Thomas Bryce Kelly
#' @param P Pressure in db
#' @param S Salinity in PSU
#' @param T Temperature in centigrade
calc.adiabatic.temp.grad = function(S, T, P) {
  T68 = 1.00024 * T

  a0 =  3.5803e-5
  a1 = 8.5258e-6
  a2 = -6.836e-8
  a3 =  6.6228e-10

  b0 = 1.8932e-6
  b1 = -4.2393e-8

  c0 = 1.8741e-8
  c1 = -6.7795e-10
  c2 = 8.733e-12
  c3 = -5.4481e-14

  d0 = -1.1351e-10
  d1 =  2.7759e-12

  e0 = -4.6206e-13
  e1 = +1.8676e-14
  e2 = -2.1687e-16

  ADTG = a0 + (a1 + (a2 + a3 * T68) * T68) * T68 + (b0 + b1 * T68) * (S-35)
  + ((c0 + (c1 + (c2 + c3 * T68) * T68) * T68) + (d0 + d1 * T68) * (S-35)) * P
  + (e0 + (e1 + e2 * T68) * T68 ) * P^2

  ADTG
}


#' @title Calculate Potential Temperature
#' @author Thomas Bryce Kelly
#' @param S Salinity in PSU
#' @param T Temperature in centigrade
#' @param P Sample pressure in db
#' @param P.ref Reference pressure, typically 0 db (surface).
#' @export
calc.ptemp = function(S, T, P, P.ref) {
    ## Calculates the potential temperature

    del.P  = P.ref - P
    del.th = del.P * calc.adiabatic.temp.grad(S, T, P)
    th     = T * 1.00024 + 0.5 * del.th
    q      = del.th

    ## theta2
    del.th = del.P * calc.adiabatic.temp.grad(S, th/1.00024, P + 0.5 * del.P)
    th     = th + (1 - 1 / sqrt(2)) * (del.th - q)
    q      = (2-sqrt(2)) * del.th + (-2 + 3/sqrt(2)) * q

    ## theta3
    del.th = del.P * calc.adiabatic.temp.grad(S, th / 1.00024, P + 0.5 * del.P);
    th     = th + (1 + 1/sqrt(2)) * (del.th - q)
    q      = (2 + sqrt(2)) * del.th + (-2 - 3/sqrt(2)) * q

    ## theta4
    del.th = del.P * calc.adiabatic.temp.grad(S, th/1.00024, P + del.P)

    ## Potential Temperature
    (th + (del.th - 2 * q)/6) / 1.00024
}


#' @title Calculate Potential Density Anomaly
#' @author Thomas Bryce Kelly
#' @description Function to calculate potential density anomaly from T, S and pressure. Comonly called sigma-theta, potential density anomaly is a useful water mass tracer due to conservation of potential density due to adiabatic processes.
#' @param S Salinity in PSU
#' @param T Temperature in centigrade
#' @param P Sample pressure in db
#' @param P.ref Reference pressure, typically 0 db (surface).
#' @examples
#' calc.sigma.theta(S = 35, T = 8, P = 350)
#' @export
calc.sigma.theta = function(S, T, P, P.ref = 0) {
    ## Potential density anomaly
    calc.rho(S = S, T = calc.ptemp(S = S, T = T, P = P, P.ref= P.ref), P = P.ref) - 1000
}


#' @title Convert Oxygen Units
#' @export
#' @author Thomas Bryce Kelly
#' @references https://ocean.ices.dk/tools/unitconversion.aspx
conv.uM.to.ml = function(x, rev = FALSE, verbose = FALSE) {
  if (rev) {
    if (verbose) {print('Converting oxygen from ml/l to uM')}
    return(x/0.022391)
  }
  if (verbose) {print('Converting oxygen from uM to ml/l')}
  return(x * 0.022391)
}


#' @title Convert Conductivity to Salinity
#' @author Thomas Bryce Kelly
#' @description Convert from conductivity observations to the seawater salinity.
#' @param C Conductivity in uS cm-1
#' @param T Temperature in centigrade
#' @export
conv.cond.to.sal = function(C = 10000, T = 10) {
    ## Convert conductivity in uS/cm to ppt salinity.

    a0 = 0.008
    a1 = -0.1692
    a2=25.3851
    a3=14.0941
    a4=-7.0261
    a5=2.7081

    b0=0.0005
    b1=-0.0056
    b2=-0.0066
    b3=-0.0375
    b4=0.0636
    b5=-0.0144

    c0=0.6766097
    c1=0.0200564
    c2=0.0001104259
    c3=-0.00000069698
    c4=0.0000000010031
    r = C / 42914
    r = r / (c0 + T*(c1 + T*(c2 +T *(c3 +T*c4))))

    r2 = r ^ 0.5
    ds = b0 + r2*(b1 + r2*(b2 + r2*(b3 + r2*(b4 + r2*b5))))
    ds = ds * (T - 15) / (1 + 0.0162*(T - 15))

    sal = a0 + r2*(a1 + r2*(a2 + r2*(a3 + r2*(a4 + r2*a5)))) + ds

    ## Check value = 8.122 ppt
    sal
}

#' @title Convert Pressure to Depth
#' @author Thomas Bryce Kelly
#' @description Convert between pressure and depth using UNESCO 1983 coefficients.
#' @references UNESCO 1983
#' @param p Pressure in db
#' @param latitude Latitude in degrees N
#' @param geo.anom Geopotential anomaly in meters
#' @param d Depth in meters, NULL unless using the function in reverse
#' @param rev Run the function in reverse? (i.e. Depth -> Pressure)
#' @import polynom
#' @export
conv.p.to.d = function(p=1e4, latitude = 30, geo.anom = 0, d = NULL, rev = FALSE) {
    ## UNESCO 1983
    ## geo.anom in J/Kg
    ##
    ## Check Value:
    ## DEPTH = 9712.653 M for P=10000 DECIBARS, LAT=30 DEG, DEL=0

    ## Specific volume
    a1 = 9.72659
    a2 = -2.2512e-5
    a3 = 2.279e-10
    a4 = -1.82e-15

    V = a1 * p + a2 * p^2 + a3 * p^3 + a4 * p^4

    ## Gravitational Variation
    b2 = 5.2788e-3
    b4 = 2.36e-5
    r = latitude * 3.1415926535 / 180

    g = 9.880318 * (1 + b2 * sin(r)^2 + b4 * sin(r)^4)

    if (rev) {
        ans = rep(NA, length(d))

        for (i in 1:length(d)) {
            V = (d - geo.anom/9.8) * g
            p = solve(polynom::polynomial(c(-V[i],a1, a2, a3, a4)))
            p = Re(p)
            ans[i] = p[which(p > 0 & p < 2 * d[i])]
        }
        return(ans)
    }

    ## Return
    V/g + geo.anom/9.8
}



#' @title Calculate AOU
#' @export
#' @param T temperature
#' @param S salinity
#' @param P pressure
#' @param oxy oxygen in umol kg-1
calc.AOU = function(T, S, P, oxy) { #T = Celcius, P = dbar, P, oxy = umol kg-1
  Osat = calc.O2sol(S = S, T = T) ##uM
  Osat = (Osat*1000)/calc.rho(S=S, T=T, p=P) #convert to umol kg-1
  Osat - oxy
}

#' @title Calculate Phosphate Star (PO4Star)
#' @export
#' @references broeker and peng 1998
calc.POstar = function(Oxy, PO4, R = 1/175){ ## using defined redfield. #umol kg
  PO4 + oxy*R - 1.95
}

#' @title Calculate Preformed Phosphate
#' @export
#' @references broeker and peng 1998
calc.pref.PO4 = function(AOU, PO4, R = 1/175){ #umol kg (calculate preformed PO4)
  PO4 - R * AOU
}



#####################################
#####################################
##### Air Sea Gas Exchange  #########
#####################################
#####################################

#' @title Calculate Schmidt Number
#' @author Thomas Bryce Kelly
#' @description Calcualtes the schmidt number for seawater at a given temperature.
#' @references Wanninkhof, R. (2014). Relationship between wind speed and gas exchange over the ocean revisited: Gas exchange and wind speed over the ocean. Limnology and Oceanography: Methods 12, 351–362. doi:10.4319/lom.2014.12.351.
#' @param SST Water temperature in degree centigrade
#' @export
calc.schmidt.number = function(SST) {

    ## Parameters of fit from Wannikhof 2014.
    # Coefficients for Least Squares Forth Order Polynomial fits of Schmidt Number vs Temperature for
    # Seawater (35ppt) ranging from -2 to 40C. Below is for O2 in seawater:
    a = 1920.4
    b = -135.6
    c = 5.2122
    d = -0.10939
    e = 0.00093777

    ## Calculate the Schmidt Number
    a + b*SST + c*(SST^2) + d*(SST^3) + e*(SST^4)
}


#' @title Calculate k
#' @description Calculate k (aka piston velocity) for a given windspeed and water temperature.
#' @references Wanninkhof, R. (2014). Relationship between wind speed and gas exchange over the ocean revisited: Gas exchange and wind speed over the ocean. Limnology and Oceanography: Methods 12, 351–362. doi:10.4319/lom.2014.12.351.
#' @references Wanninkhof, R. (1992). Relationship between wind speed and gas exchange over the ocean. Journal of Geophysical Research 97, 7373. doi:10.1029/92JC00188.
#' @references Sweeney, C., Gloor, E., Jacobson, A. R., Key, R. M., McKinley, G., Sarmiento, J. L., et al. (2007). Constraining global air-sea gas exchange for CO 2 with recent bomb 14 C measurements: BOMB 14 C AND AIR-SEA GAS EXCHANGE. Global Biogeochemical Cycles 21, n/a-n/a. doi:10.1029/2006GB002784.
#' @export
#### Calcualtes k
calc.k = function(u, SST) { # m/d
    ## Based on scaling estimate using 14C bomb data, 0.27
    #0.27 * u^2 / sqrt(schmidt.number(SST)/660) * (24/100)  ## Sweeny et al

    ## "Relationship should be applicable to deduce gas transfer velocities at steady winds
    # using shipboard anemometers"
    #0.31 * u^2 / sqrt(schmidt.number(SST)/660) * (24/100) # Wanninkov 1992

    ## Wanninkov 2014
    0.251 * u^2 / sqrt(calc.schmidt.number(SST) / 660) * (24/100) # Improvement over Sweeney 2007.
}


#' @title Calculate k (Reuer)
#' @description This function calculates the temporal length scale of the NCP measurement, aka a piston velocity (m/d). Essentially it calculates the ventilation based on water column and wind speed over a sequence of wind speed measurements (at 10 m height).
#' @param time the time at which you are calculating k for
#' @param wtime a vector of times at which a wind speed is recorded
#' @param wspeed a vector of prior windspeeds ordered from most recent to oldest
#' @param mld mixed layer depth in meters
#' @param SST mixed layer temperature in C
#' @param teeter.mod a boolean flag to turn on or off the modified version of Reuer's weighting scheme as outlined in Teeter et al. 2018.
#' @author Thomas Bryce Kelly
#' @references Reuer, M. K., Barnett, B. A., Bender, M. L., Falkowski, P. G., and Hendricks, M. B. (2007). New estimates of Southern Ocean biological production rates from O2/Ar ratios and the triple isotope composition of O2. Deep Sea Research Part I: Oceanographic Research Papers 54, 951–974. doi:10.1016/j.dsr.2007.02.007.
#' @references Teeter, L., Hamme, R. C., Ianson, D., and Bianucci, L. (2018). Accurate Estimation of Net Community Production From O 2 /Ar Measurements. Global Biogeochem. Cycles. doi:10.1029/2017GB005874.
#' @export
#' @keywords Air-sea Oceanography NCP
#' @return It will return a list with the individual weights and weighted sum (for dianostics), the k value (that you want, m/d), the instantaneous k (also useful, m/d), k.hist showing the history of k values (diagnostic), and time is time before measurement that each diagnostic is referenced against (for diagnostics, d).
calc.k.reuer = function(time, wtime, wspeed, mld, SST, teeter.mod = TRUE) { ## m d-1
  if (is.na(mld)) {
    return(list(dt = NA,
                weights = NA,
                w.sum = NA,
                k = NA,
                k.inst = NA,
                k.hist = NA,
                time = 0,
                integration.time = 0)
           )
  }
  ## Start
  weights = rep(1, length(wspeed)) # initialize to weights of 1
  k.result = calc.k(wspeed, SST) # k's for each time point (m/d)

  dt = as.numeric(difftime(time, wtime, units = 'days')) # positive values are before
  dt = c(dt[1], diff(dt)) # delta time between each point

  if (length(wspeed) > 1) {
    ## Loop through and sequentially calculate the transfer velocities and weighting of ventilation
    for (i in 2:length(wspeed)) {
      ## Calculate the weighting per Reuer et al. (2007)
      f = max(0, min(k.result[i-1] * dt[i] / mld, 1))
      weights[i] = weights[i-1] * (1 - f)
    }
  }
  w = sum(weights)
  weights = weights / w ## Normalize

  if (any(is.na(weights))) { warning('Weights contain an NA, check data for problems!')}
  #### Final Calculations
  if (teeter.mod) { k.final = sum(weights * k.result, na.rm = TRUE) }
  else { k.final = sum(weights * k.result, na.rm = TRUE) / (1 - weigths[length(weights)]) }

  tau = NA
  tau = approx(cumsum(weights), cumsum(dt), 0.9, rule = 2)$y

  list(dt = dt,
       weights = weights,
       w.sum = w,
       k = k.final,
       k.inst = k.result[1],
       k.hist = k.result,
       time = cumsum(dt),
       integration.time = tau)
}


#' @title Make TS Plot
#' @export
#' @author Thomas Bryce Kelly
plot.TS = function(x, y, xlim = c(25,36), ylim = c(-5, 15), levels = seq(0, 40, by = 2),
                   cex.lab = 1, drawlabels = TRUE, labels = NULL, freezing.line = T, freezing.col = '#00000030',
                   col = 'grey', lwd = 1, lty = 1, xlab = 'Practical Salinity', ylab = 'Potential Temperature',
                   pch = 1, col.pch = 'black', cex.pch = 1, main = NULL) {

  T = seq(ylim[1], ylim[2], length.out = 100)
  S = seq(xlim[1], xlim[2], length.out = 100)
  grid = expand.grid(S = S, T = T)
  sigma = calc.sigma.theta(S = grid$S, T = grid$T, P = 0, P.ref = 0)

  contour(x = S, y = T, z = matrix(sigma, nrow = length(S)), levels = levels, lwd = lwd, lty = lty, col = col,
          xaxs = 'i', yaxs = 'i', xlab = xlab, ylab = ylab, main = main, xlim = xlim, ylim = ylim,
          labcex = cex.lab, drawlabels = drawlabels, labels = labels)
  points(x, y, pch = pch, col = col.pch, cex = cex.pch)

  if (freezing.line) {
    x = seq(xlim[1], xlim[2], length.out = 1000)
    y = calc.freezing.point(x, 0)
    polygon(x = c(x, rev(x)), c(y, rep(-100, length(x))), col = freezing.col, border = NA, density = 20)
  }
}


#####################################
#####################################
######  Light Formulas   ##########
#####################################
#####################################

#' @description (6.02e23 quanta/Ein) / (86,400 seconds/day * 2.77e18 quanta/s/W)
#' @title Convert Light Units
#' @export
#' @author Thomas Bryce Kelly
#' @references Unit conversions as listed by MBARI: https://www3.mbari.org/bog/nopp/par.html
conv.einstein.to.watt = function(x) {
    x * 6.02e23 / 2.77e18
}

#' @description (6.02e23 quanta/Ein) / (86,400 seconds/day * 2.77e18 quanta/s/W)
#' @title Convert Light Units
#' @author Thomas Bryce Kelly
#' @export
#' @references Unit conversions as listed by MBARI: https://www3.mbari.org/bog/nopp/par.html
conv.watt.to.einstein = function(x) {
    x * 2.77e18 / 6.02e23
}



#####################################
#####################################
#####   Helpful Functions   #########
#####################################
#####################################

#' @title Calculate MLD
#' @author Thomas Bryce Kelly
#' @param depth depth
#' @param density the density of the seawater. Alternatively, if using a temperature based MLD definition, the temperature of the seawater.
#' @param set.depth the depth from which the mixed layer character is defined
#' @param delta the difference from the mixed layer signal that defines the end of the mixed layer.
#' @param resolution the resolution of the vertical interpolation for when calculating the MLD.
#' @return returns the depth where the signal has changed from set.depth by >= delta.
#' @export
calc.mld = function(depth, density, set.depth = 10, delta = 0.1, resolution = 0.1) {

    depth.new = seq(set.depth, 120, by = 0.1) ## 0.1 meter resolution
    density.new = approx(depth, density, xout = depth.new, rule = 2)$y

    l = which(density.new >= density.new[1] + delta)

    if (length(l) < 1) {
        warning('No MLD found! Returning NA.')
        return(NA)
    }

    min(depth.new[l])
}


calc.freezing.point = function(S, p, f = 0) {
  ## Constants
  c0 = 0.002519
  c1 = -5.946302841607319
  c2 =  4.136051661346983
  c3 = -1.115150523403847e1
  c4 =  1.476878746184548e1
  c5 = -1.088873263630961e1
  c6 =  2.96101883964073
  c7 = -7.433320943962606
  c8 = -1.561578562479883
  c9 =  4.073774363480365e-2
  c10 =  1.158414435887717e-2
  c11 = -4.122639292422863e-1
  c12 = -1.123186915628260e-1
  c13 =  5.715012685553502e-1
  c14 =  2.021682115652684e-1
  c15 =  4.140574258089767e-2
  c16 = -6.034228641903586e-1
  c17 = -1.205825928146808e-2
  c18 = -2.812172968619369e-1
  c19 =  1.877244474023750e-2
  c20 = -1.204395563789007e-1
  c21 =  2.349147739749606e-1
  c22 =  2.748444541144219e-3

  ## Rescale
  SA_r = S * 1e-2
  x = sqrt(SA_r)
  p_r = p * 1e-4

  tf.1 = c0 + SA_r * (c1 + x * (c2 + x * (c3 + x * (c4 + x * (c5 + c6 * x)))))
  tf.2 =  p_r * (c7 + p_r * (c8 + c9 * p_r))
  tf.3 =  SA_r * p_r * (c10 + p_r * (c12 + p_r * (c15 + c21 * SA_r)) + SA_r * (c13 + c17 * p_r + c19 * SA_r) +
                          x * (c11 + p_r * (c14 + c18 * p_r)  + SA_r * (c16 + c20 * p_r + c22 * SA_r)))
  tf = tf.1 + tf.2 + tf.3

  t_freezing = tf - f * (1e-3) * (2.4 - S / 70.33008) ## Adjust for fraction of air

  t_freezing
}

