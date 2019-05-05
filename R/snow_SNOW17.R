#' SNOW-17 Snow Module
#'
#' Calculates snow via SNOW17 Model.  Current version is configured for Anderson 1973 version with some 2006 elements
#'
#' @param Param snow module parameter list
#' @param Prcp Precipitation continuous time-series (mm) (xts object) (same duration as Tavg)
#' @param Tavg Average Temperature continuous time-series (deg C) (xts object) (same duration as Prcp)
#' @param Elevation Elevation (meters)
#' @param InitialState Initial state vector (W_i, ATI, W_q, Deficit)
#' @param dtt Constant time interval of temperature data (hours) (Default: 24 Hours)
#' @param dtp Constant time interval of precipitation data (hours) (Default: 24 Hours)
#' @param calcA_v If Latitude is above 60 degrees then calculate A_v for each time-step
#' @param meltFlag Output a vector time-series melt flag
#' @param preserveInput Flag to preserve the input as part of the output
#' @param verbose Verbose Flag
#'
#' @return A list with various time-series outputs (xts objects)
#' @details For details see Anderson (2006) and Anderson (1973)
#' @export

snow_SNOW17 <- function(Param, Prcp, Tavg, Elevation, InitialState = c(0, 0, 0, 0), dtt = 24, dtp = 24, calcA_v = FALSE, meltFlag = FALSE, preserveInput = FALSE, verbose = FALSE) {

  JDate <- format(zoo::index(Prcp), '%j')

  Dates.ts <- zoo::index(Prcp)
  Prcp.vector <- as.vector(Prcp)
  Tavg.vector <- as.vector(Tavg)

  # Preprocessing ----
  verbose.startTime <- Sys.time()
  verbose.timeStepTotal <- length(Prcp.vector)

  if (verbose) {
    print('Initializing SNOW17...')
  }

  # Environment Preparation ----
  Output <- list()
  Output$Output <- list()
  Outflow <- rep(NA, times = length(Prcp.vector))
  melt <- rep(NA, times = length(Prcp.vector))
  SWEO <- rep(NA, times = length(Prcp.vector))
  Depth <- rep(NA, times = length(Prcp.vector))
  Output$FinalState <- list()

  # If preserving the input
  if (preserveInput) {

    if (verbose) {
      print('Saving Inputs with Outputs...')
    }

    Output$Input <- list()
    Output$Input$TimeSeries <- list()
    Output$Input$TimeSeries$Prcp <- Prcp
    Output$Input$TimeSeries$Tavg <- Tavg
    Output$Input$Param <- Param
    Output$Input$Elevation <- Elevation
    Output$Input$InitialState <- InitialState
    Output$Input$dtt <- dtt
    Output$Input$dtp <- dtp
    Output$Input$calcA_v <- calcA_v
  }

  # Parameter Defining ----
  if (is.list(Param)){
    SCF <- Param$SCF # Snow cover factor (multiplying factor)
    f_s <- Param$f_s # Fraction of precipitation in the form of snow
    PXTEMP <- Param$PXTEMP # Threshold Temp (degC)
    MFMAX <- Param$MFMAX # Maximum melt factor during non-rain periods (mm * degC^-1 * 6hr^-1)
    MFMIN <- Param$MFMIN # Minimum melt factor during non-rain periods (mm * degC^-1 * 6hr^-1)
    UADJ <- Param$UADJ # Average wind function during rain-on-snow periods (mm * mb^-1)
    MBASE <- Param$MBASE # Base temperature for snowmelt during non-rain periods (degC)
    TIPM <- Param$TIPM # Antecedent temperature index parameter
    PLWHC <- Param$PLWHC # Percent liquid water holding capacity (decimal fraction)
    NMF <- Param$NMF # Maximum negative melt factor (mm * degC^-1 * 6hr^-1)
    DAYGM <- Param$DAYGM # Daily amount of melt (mm * day^-1)
    A_v <- Param$A_v # Seasonal variation Adjustment (decimal fraction) (Refer to Anderson 2006 for value equation)
  } else {
    stop('Param variable is not in an accepted format')
  }

  # Parameter Range Testing ====
  # This is a TODO


  # Initial State----
  W_i <- InitialState[1] # Water equivalent of the ice portion of the snow cover (mm)
  ATI <- InitialState[2] # Antecedent temperature index (degC)
  W_q <- InitialState[3] # Liquid water held by the snow (mm)
  Deficit <- InitialState[4] # Heat deficit (mm)
  E <- 0

  # Defining Constants ----
  SBConst <- 6.12e-10 # Stefan-Boltzman Constant (mm/K/hr)
  L_f <- 80 # Latent Heat of Fusion (cal * gm^-1)
  c_i <- 0.5 # Specific Heat of Ice (cal * gm^-1 * degC^-1)
  f_r <- 1 - f_s # Fraction of precipitation in the form of rain

  if (f_s == 1) f_r <- 1 # A bit of a loophole to coincide with MATLAB version

  if (verbose) {
    print('Running the model...')
    pb <- utils::txtProgressBar(min = 0, max = verbose.timeStepTotal, style = 3)
  }

  # Primary Loop ----
  for(i in 1:length(Prcp.vector)) {

    if (verbose) setTxtProgressBar(pb, i)

    T_a <- Tavg.vector[i]
    P_r <- Prcp.vector[i]
    Date <- Dates.ts[i]


    T_rain <- max(T_a, 0) # Temperature of the Rain
    T_n <- min(T_a, 0) # Temperature of new snow

    # Precipitation Forming ====
    if (T_a <= PXTEMP) {
      # Temp is cold for snow
      NEW_SNOW <- P_r
      RAIN <- 0
    } else {
      # No snow, all rain
      NEW_SNOW <- 0
      RAIN <- P_r
    }

    P_n <- NEW_SNOW * f_s * SCF # Water equivalent of new snowfall (mm)
    W_i <- W_i + P_n # Water equivalent of the ice portion of the snow cover (mm)
    E <- 0 # Excess liquid water in the snow cover

    # Temperature of new snow
    deltaD_p <- -((T_n * P_n) / (L_f / c_i)) # Change in heat deficit due to snowfall (mm)

    # Snow Cover Accumulation ====

    # Density of New Snow
    if(T_a <= -15) {
      rho_n <- 0.05
    } else {
      rho_n <- 0.05 + 0.0017 * (T_a ^ 1.5)
    }

    H_n <- (0.1 * P_n)/rho_n # Depth of new snowfall (cm)


    # Energy Exchange at Snow/Air Interface ====

    # Seasonal Variation and Melt Factor ####
    if (utility_isLeapYear(as.numeric(format(Date, '%Y')))) { # If Leap year
      N <- as.numeric(JDate[i]) - 81
      Days <- 366
    } else { # Not a leap year
      N <- as.numeric(JDate[i]) - 80
      Days <- 365
    }

    S_v <- 0.5 * sin((N * 2 * pi)/Days) + 0.5

    # If latitude is above 60 A_v should be recalculated per time-step
    if (calcA_v) {
      A_v <- utility_A_v(Date)
    }


    M_f <- (dtt/6) * (S_v * A_v * (MFMAX - MFMIN) + MFMIN) # Seasional varying non-rain melt feactor

    # Energy Exchange when No Surface Melt ####
    # Handling Antecedent Temperature Index
    if (P_n > (1.5 * dtt)) {
      ATI <- T_n
    } else {
      TIPM_deltat_t <- 1.0 - ((1.0 - TIPM)^(dtt/6))
      ATI <- ATI + TIPM_deltat_t * (T_a - ATI)  #Antecedent Temperature Index
    }
    ATI <- min(ATI, 0)


    if (RAIN > (0.25 * dtp)) {
      # Rain-on-Snow Melt ####
      P_a <- 33.86 * (29.9 - (0.335 * (Elevation/100)) + (0.00022 * (Elevation/100)^(2.4))) # Atmospheric Pressure (mb) (This equation is incorrectly written in Anderson 2006)
      e_sat <- 2.7489e8 * exp((-4278.63) / (T_a + 242.792)) # Saturation vapor pressure calculation
      # Melt due to rain (M_r)
      Melt <- SBConst * dtp * (((T_a + 273)^4) - (273^4)) + (0.0125 * RAIN * f_r * T_rain) + (8.5 * UADJ * (dtp/6) * ((0.9 * e_sat - 6.11) + 0.00057 * P_a * T_a))
      Melt <- max(Melt, 0)

    } else if (RAIN <= (0.25 * dtp) & (T_a > MBASE)) {
      # Non-Rain Melt ####
      # Melt in no rain (M_nr)
      Melt <- (M_f * (T_a - MBASE) * (dtp/dtt)) + 0.0125 * RAIN * T_rain
      Melt <- max(Melt, 0)
    } else {
      Melt <- 0
    }

    # Internal State of Snow Cover/Ripeness of the Snow Cover ====
    delta_D_t <- NMF * (dtp/6) * (M_f/MFMAX) * (ATI - T_n)

    Deficit <- max((Deficit + deltaD_p + delta_D_t), 0) # Heat Deficit (mm)

    # Limits of heat deficit
    if (Deficit > (0.33 * W_i)) Deficit <- (0.33 * W_i)

    if (Melt < W_i) {
      W_i <- W_i - Melt
      Q_w <- Melt + RAIN * f_r
      W_qx <- PLWHC * W_i # Water equivalent of the ice portion of the snow cover (mm)

      if ((Q_w + W_q) > (Deficit + W_qx + (PLWHC * Deficit))) {
        # Snow is Ripe
        E <- Q_w + W_q - W_qx - Deficit - (PLWHC * Deficit) # Excess liquid water (mm)
        W_i <- W_i + Deficit # W_i increases because water refreezes as heat deficit decreases
        W_q <- W_qx + PLWHC * Deficit # Fills liquid water capacity
        Deficit <- 0

      } else if ((Q_w + W_q) >= Deficit) {
        # Snow is not Ripe but ice is being melted
        E <- 0
        W_i <- W_i + Deficit # W_i increases as water refreezes as heat deficit is decreased
        W_q <- W_q + Q_w - Deficit
        Deficit <- 0

      } else if ((Q_w + W_q) < Deficit) {
        # Snow is not Ripe
        E <- 0
        W_i <- W_i + Q_w + W_q # W_i increases as water refreezes as heat deficit is decreased
        Deficit <- Deficit - Q_w - W_q

      }

    } else { # Melt >= W_i

      Melt <- W_i + W_q
      W_i <- 0
      W_q <- 0
      Q_w <- Melt + RAIN * f_r # Direct Runoff Scenario (Check f_r is supposed to be there)
      E <- Q_w

    }

    if (Deficit == 0) {
      ATI <- 0
    }

    # Density and Depth Computations ====
    # This is a TODO

    # Transmission of Water through the Snow Cover ====
    # This is a TODO

    # Heat Transfer at the Snow-Soil interface/Release ====
    # This requires an update to Anderson 2006

    if (W_i > DAYGM) {

      gmwlos <- (DAYGM/W_i) * W_q
      gmslos <- DAYGM
      gmro <- gmwlos + gmslos
      W_i <- W_i - gmslos
      W_q <- W_q - gmwlos

      E <- E + gmro
      SWE <- W_i + W_q

    } else {

      gmro <- W_i + W_q
      W_i <- 0
      W_q <- 0

      E <- E + gmro
      SWE <- 0

    }

    # Output Writing ====
    Outflow[i] <- E
    melt[i] <- Melt
    SWEO[i] <- SWE

    Output$FinalState$W_i <- W_i
    Output$FinalState$ATI <- ATI
    Output$FinalState$W_q <- W_q
    Output$FinalState$Deficit <- Deficit

  } # End of Timestep

  # Finalizing Output Writing
  Output$Output$Outflow <- xts::xts(Outflow, order.by = as.Date(zoo::index(Prcp)))
  Output$Output$melt <- xts::xts(melt, order.by = as.Date(zoo::index(Prcp)))
  Output$Output$SWE <- xts::xts(SWEO, order.by = as.Date(zoo::index(Prcp)))

  # Quick Melt Flag Calculation Afterwards
  if (meltFlag){
    meltFlag <- sapply(melt, function(ts){if (ts > 0) return(TRUE) else return(FALSE)})
    Output$Output$meltFlag <- meltFlag
  }

  verbose.endTime <- Sys.time()

  if (verbose) {
    print(paste0('SNOW17 Module Run Time: ', format(verbose.endTime - verbose.startTime)))
  }

  return(Output)
}

# SNOW17 Ranges ----
snow_SNOW17.ranges <- function(){
  list(
    SCF <- c(1, 5), # Figure this out
    f_s <- c(0, 1), # Decimal fraction from Anderson 2006
    PXTEMP <- c(-10, 10), # degC from Anderson 2006
    MFMAX <- c(NA, NA),
    MFMIN <- c(NA, MFMAX), #
    UADJ <- c(0, NA),
    MBASE <- c(NA, NA),
    TIPM <- c(0.01, 1), # Decimal fraction from Anderson 2006, should not include 1 or 0
    PLWHC <- c(0, 0.4), # Decimal fraction from Anderson 2006
    NMF <- c(NA, NA), #
    DAYGM <- c(NA, NA), # Should be between
    A_v <- c(1, 1) # Do not calibrate.  Use 1 unless Latitude is above 54 degrees then use function utility_A_v


  )
}

# SNOW17 Default parameters ----
snow_SNOW17.defaults <- function() {
  list(
    SCF <- 1.2, # Find this Out
    f_s <- 1, # Decimal fraction from Wi code
    PXTEMP <- 0, # From Default Constant
    MFMAX <- NA,
    MFMIN <- NA,
    UADJ <- 1, # Find this Out
    MBASE <- 0, # degC from Anderson 2006
    TIPM <- 0.5, # Find this Out
    PLWHC <- 0.2, # Find this Out (This is middle of the range)
    NMF <- 1, # Find this Out
    DAYGM <- 1, # Find this Out
    A_v <- 1 # Decimal Fraction from Anderson 2006 (This should be used)

  )

}
