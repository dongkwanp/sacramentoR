#' Sacramento Soil Moisture Accounting Module (SAC-SMA)
#'
#' Calculates the Sacramento Soil Moisture Accounting (SAC-SMA) Model
#'
#' @param Param hydrology module parameter list
#' @param Prcpts Precipitation/Snow Melt continuous time-series (mm) (xts object) (same duration as PETts)
#' @param PETts Potential Evapotranspiration continuous time-series (deg C) (xts object) (same duration as Prcpts)
#' @param InitialState Initial state vector (uztwc, uzfwc, lztwc, lzfsc, lzfpc, adimc)
#' @param thresholdZero Threshold for Zero
#' @param ninc_min Minimum number of time increments that interval is divided into (Default: 20)
#' @param preserveInput Flag to preserve the input as part of the output
#' @param preserveState Flag to preserve state of the model at each time-step
#' @param verbose Verbose Flag
#' @param verbose.ts Verbose per time-step Flag
#'
#' @return A list of various time-series as an xts object
#' @details For details see NOAA rfs:23sacsma.wpd (2002) and Blasone, et. al. (2008)
#' @export

hydrology_sacsma <- function(Param, Prcpts, PETts, InitialState = c(0, 0, 500, 500, 500, 0), thresholdZero = 0.00001, ninc_min = 20, preserveInput = FALSE, preserveState = FALSE, verbose = FALSE, verbose.ts = FALSE) {

  # Preprocessing ----
  verbose.startTime <- Sys.time()
  verbose.timeStepTotal <- length(Prcpts)

  if (verbose) {
    print('Initializing SAC-SMA...')
  }

  # Environment Preparation ----
  Output <- list()
  Output$State <- list()
  #Output$TotalFlow <- rep(NA, times = length(Prcp))
  #Output$SurfaceFlow <- rep(NA, times = length(Prcp))
  #Output$BaseFlow <- rep(NA, times = length(Prcp))

  if (preserveInput) {

    if (verbose) {
      print('Saving Inputs with Outputs...')
    }

      Output$Input <- list()
      Output$Input$TimeSeries <- list()
      Output$Input$TimeSeries$Prcpts <- Prcpts
      Output$Input$TimeSeries$PETts <- PETts
      Output$Input$Param <- Param
      Output$Input$InitialState <- InitialState
  }

  # Converting to vectors from xts object (but preserving xts)
  TimeSteps <- zoo::index(Prcpts)
  Prcpts.xts <- Prcpts
  Prcpts <- as.vector(Prcpts.xts)
  PETts.xts <- PETts
  PETts <- as.vector(PETts.xts)

  # Parameter Defining ----
  if (is.list(Param)) {
    # Capacity Thresholds
    uztwm <- Param$uztwm # Upper zone tension water capacity (mm)
    uzfwm <- Param$uzfwm # Upper zone free water capacity (mm)
    lztwm <- Param$lztwm # Lower zone tension water capacity (mm)
    lzfpm <- Param$lzfpm # Lower zone primary free water capacity (mm)
    lzfsm <- Param$lzfsm # Lower zone supplementary free water capacity (mm)

    # Recession Parameters
    uzk <- Param$uzk # Upper zone free water lateral depletion rate (1/day)
    lzpk <- Param$lzpk # Lower zone primary free water depletion rate (1/day)
    lzsk <- Param$lzsk # Lower zone supplementary free water depletion rate (1/day)

    # Percolation
    zperc <- Param$zperc # Percolation demand scale parameter (unitless)
    rexp <- Param$rexp # Percolation demand shape parameter (unitless): exponent of the percolation equation
    pfree <- Param$pfree # Percolating water split parameter (decimal fraction): Fraction of water percolating from upper zone directly to lower zone free water storage

    # Impervious Area
    pctim <- Param$pctim # Impervious fraction of the watershed area (decmial fraction)
    adimp <- Param$adimp # Additional impervious areas (decimal fraction)

    # Others
    riva <- Param$riva # Riparian vegetation area (decimal fraction)
    side <- Param$side # The ratio of deep recharge to channel base flow (unitless)
    rserv <- Param$rserv # Fraction of lower zone free water not transferrable to lower zone tension water (decimal fraction)

  } else {
    stop('Param variable is not in an acceptable format')
  }

  # Parameter Range Testing ====
  # This is a TODO

  # Initial State ----
  uztwc <- InitialState[1] # Upper zone tension water storage
  uzfwc <- InitialState[2] # Upper zone free water storage
  lztwc <- InitialState[3] # Lower zone tensio nwater storage
  lzfsc <- InitialState[4] # Lower zone supplementary free water storage
  lzfpc <- InitialState[5] # Lower zone primary free water storage
  adimc <- InitialState[6] # Additional impervious area storage

  # Defining Constants ----
  parea <- 1 - adimp - pctim # Permeable area = 1 - additional impervious area fraction - impervious fraction

  # Allocating Memory for Outputs and States of the Model ----
  if (preserveState) {
    # Upper Zone States
    uztwc_state <- rep(NA, times = length(Prcpts)) # State of Upper zone tension water storage (mm)
    uzfwc_state <- rep(NA, times = length(Prcpts)) # State of Upper zone free water storage (mm)

    # Lower Zone States
    lztwc_state <- rep(NA, times = length(Prcpts)) # State of Lower zone tension water storage (mm)
    lzfsc_state <- rep(NA, times = length(Prcpts)) # State of Lower zone free water supplementary storage (mm)
    lzfpc_state <- rep(NA, times = length(Prcpts)) # State of Lower zone free water primary storage (mm)
    adimc_state <- rep(NA, times = length(Prcpts)) # State of additional impervious area storages (mm)
  }


  # Model Output Array Initialization
  streamflow.ts <- rep(NA, times = length(Prcpts)) # Modelled Streamflow
  aet.ts <- rep(NA, times = length(Prcpts)) # Modelled Actual Evapotranspiration
  baseflow.ts <- rep(NA, times = length(Prcpts)) # Modelled Baseflow
  surfaceflow.ts <- rep(NA, times = length(Prcpts)) # Modelled Surface and Subsurface Flow
  interflow.ts <- rep(NA, times = length(Prcpts)) # Modelled Interflow Flow

  if (verbose) print('Running the model...')
  # Calculating through the time-steps ----
  for (i in 1:length(Prcpts)) {

    if (verbose && verbose.ts) print(paste0('Running Time-Step: ', i, ' out of ', verbose.timeStepTotal))

    Prcp <- Prcpts[i] # Adjusted precipitation from snow module (assumed)

    # Computing Evapotranspiration Loss ====
    evapDemand <- PETts[i]

    # ET Module 1 - ET from Upper Zone Tension Water Storage ####
    ET1 <- evapDemand * (uztwc/uztwm)
    red <- evapDemand - ET1 # Residual ET Demand
    uztwc <- uztwc - ET1 # Updating uztwc

    # ET Module 2 - ET from Upper Zone Free Water Storage ####
    ET2 <- 0 # Initially set

    if (uztwc <= 0) { # Means that no water in uztws
      ET1 <- ET1 + uztwc # Removing all the water from uztwc
      uztwc <- 0
      red <- evapDemand - ET1

      if (uzfwc < red) { # If Upper free water is less than residual ET
        ET2 <- uzfwc
        uzfwc <- 0
        red <- red - ET2

        # Checking if below threshold and resetting
        if (uztwc < thresholdZero) uztwc <- 0
        if (uzfwc < thresholdZero) uzfwc <- 0

      } else { # If Upper free water is more than residual ET
        ET2 <- red
        uzfwc <- uzfwc - ET2
        red <- 0
      }

    } else { # All ET is handled at uztwc so none taken from ET2 so transfer free water to tension water storage

      if ((uztwc/uztwm) < (uzfwc/uzfwm)) {
        uzrat <- ((uztwc + uzfwc) / (uztwm + uzfwm))
        uztwc <- uztwm * uzrat
        uzfwc <- uzfwm * uzrat
      }

      # Checking if below threshold and resetting
      if (uztwc < thresholdZero) uztwc <- 0
      if (uzfwc < thresholdZero) uzfwc <- 0
    }

    # ET Module 3 - ET from Lower tension water storage ####
    ET3 <- red * (lztwc/(uztwm + lztwm))
    lztwc <- lztwc - ET3

    if (ET3 < 0) warning(paste0('Timestep: ', i, ' has ET3 < 0 where ET3 = ', ET5))

    if (lztwc < 0) {
      ET3 <- ET3 + lztwc
      lztwc <- 0
    }

    # Resupplying Water from Lower Free Water Storage to Lower Tension Storage ####
    saved <- rserv * (lzfpm + lzfsm)
    ratlzt <- (lztwc/lztwm)
    ratlz <- ((lztwc + lzfpc + lzfsc - saved) / (lztwm + lzfpm + lzfsm - saved))

    if (ratlzt < ratlz) { # Resupplying from supplementary storage
      del <- (ratlz - ratlzt) * lztwm
      lztwc <- lztwc + del
      lzfsc <- lzfsc - del
      if (lzfsc < 0) { # If it uses all of lzfsc, it brings more in from lzfps
        lzfpc <- lzfpc + lzfsc
        lzfsc <- 0
      }
    }

    # Checking if below threshold and resetting
    if (lztwc < thresholdZero) lztwc <- 0

    # ET Module 5 - ET from additional impervious (ADIMP) area ####
    ET5 <- ET1 + (((red + ET2) * (adimc - ET1 - uztwc)) / (uztwm + lztwm))

    if (ET5 < 0) warning(paste0('Timestep: ', i, ' has ET5 < 0 where ET5 = ', ET5))

    adimc <- adimc - ET5

    if (adimc < 0) { # ET5 cannot exceed adimc
      ET5 <- ET5 + adimc
      admic <- 0
    }

    ET5 <- ET5 * adimp

    # Time Interval Available Moisture in Excess of uztw Requirements ####
    twx <- Prcp + uztwc - uztwm

    if (twx < 0) { # All moisture in uztw (no excess)
      uztwc <- uztwc + Prcp
      twx <- 0
    } else { # Moisture available in excess of uztw storage
      uztwc <- uztwm
    } # twx is excess rainfall after fillin uztwc
    adimc <- adimc + Prcp - twx


    # Computing Impervious Area Runoff ====
    roimp <- Prcp * pctim


    # Initializing Time Interval Sums ====
    sbf <- 0 # Sum of Total Baseflow (from Primary and Supplemental Storages)
    ssur <- 0 # Sum of Surface Runoff
    sif <- 0 # Sum of Interflow
    sperc <- 0 # Time Interval Summation of Percolation
    sdro <- 0 # Sum of Direct Runoff from the Additional Impervious Area

    # Determine computational time increments for the basic time interval
    ninc <- max(floor(1 + (0.2 * (uzfwc + twx))), ninc_min) # Number of time increments that interval is divided into for SMA

    dinc <- 1.0 / ninc # Length of each increment in days
    pinc <- twx / ninc # Amount of available moisture for each increment

    # Compute free water depletion fractions for the time increment
    duz <- 1 - (1 - uzk) ^ dinc
    dlzp <- 1 - (1 - lzpk) ^ dinc
    dlzs <- 1 - (1 - lzsk) ^ dinc

    # Incremental loop for time interval ====
    for (n in 1:ninc) {

      adsur <- 0 # Amount of surface runoff

      # Computing Direct Runoff from adimp area ####
      ratio <- ((adimc - uztwc) / lztwm)

      if (ratio < 0) ratio <- 0

      addro <- pinc * (ratio ^ 2) # Amount of direct runoff from the additional impervious area

      # Computing Baseflow & keeping track of time interval sum ####
      bf_p <- lzfpc * dlzp # Baseflow from free water primary storage
      lzfpc <- lzfpc - bf_p

      # Checking if below threshold and resetting
      if (lzfpc <= 0.0001) {
        bf_p <- bf_p + lzfpc
        lzfpc <- 0
      }

      sbf <- sbf + bf_p

      bf_s <- lzfsc * dlzs # Baseflow from free water supplemental storage
      lzfsc <- lzfsc - bf_p

      # Checking if below threshold and resetting
      if (lzfsc <= 0.0001) {
        bf_s <- lzfsc + bf_s
        lzfsc <- 0
      }

      sbf <- sbf + bf_s # Total Baseflow from primary and supplemental storages

      # Computing Percolation ####
      if ((pinc + uzfwc) <= 0.01) {
        uzfwc <- uzfwc + pinc
      } else {
        percm <- (lzfpm * dlzp) + (lzfsm * dlzs) # Limiting drainage rate from the combined saturated lower zone storages
        perc <- percm * (uzfwc / uzfwm)

        defr <- 1 - ((lztwc + lzfpc + lzfsc) / (lztwm + lzfpm + lzfsm)) # defr is the lower zone moisture deficiency ratio

        if (defr < 0) defr <- 0

        perc <- perc * (1 + (zperc * (defr ^ rexp)))

        # Percolation occurs from uzfws before pav is added

        if (perc >= uzfwc) perc <- uzfwc # Perclation rate exceeds uzfws

        uzfwc <- uzfwc - perc # Percolation rate is less than uzfws

        check <- lztwc + lzfpc + lzfsc + perc - lztwm - lzfpm - lzfsm

        if (check > 0) { # Check to see if percolation exceeds lower zone deficiency
          perc <- perc - check
          uzfwc <- uzfwc + check
        }

        sperc <- sperc + perc # sperc is time interval sum of perc

        # Computing Interflow and keeping track of time interval sum
        del <- uzfwc * duz # Amount of Interflow
        sif <- sif + del # Summing Interflow
        uzfwc <- uzfwc - del

        # Distribute percolated water to lower zones
        # Tension water filled first except for PFREE area
        # PERCT is percolation to tension water and PERCF is percolation going to free water
        perct <- perc * (1 - pfree) # Percolation going to the tension water storage

        if ((perct + lztwc) <= lztwm) {
          lztwc <- lztwc + perct
          percf <- 0 # Percolation going to the lower zone free water storages
        } else {
          percf <- lztwc + perct - lztwm
          lztwc <- lztwm
        }

        # Distribute percolation in excess of tension requirements among free water storage
        percf <- percf + (perc * pfree)

        if (percf != 0) {
          hpl <- (lzfpm / (lzfpm + lzfsm)) # Relative size of primary storage as compared with total lower zone free water storages

          # Relative fullness of each storage
          ratlp <- (lzfpc / lzfpm)
          ratls <- (lzfsc / lzfsm)

          fracp <- (hpl * 2 * (1 - ratlp)) / (1 - ratlp + 1 - ratls) # Fraction going to primary

          if (fracp > 1) fracp <- 1

          percp <- percf * fracp # Amount of excess percolation to primary
          percs <- percf - percp # Amount of excess percolation to supplemental

          lzfsc <- lzfsc + percs

          if (lzfsc > lzfsm) {
            percs <- percs - lzfsc + lzfsm
            lzfsc <- lzfsm
          }

          lzfpc <- lzfpc + percf - percs

          if (lzfpc >= lzfpm) { # Check to make sure lzfps does not exceed lzfpm
            excess <- lzfpc - lzfpm
            lztwc <- lztwc + excess
            lzfpc <- lzfpm
          }

        }

        # Distribute pinc between uzfws and surface runoff ####
        if (pinc != 0) {
          if ((pinc + uzfwc) <= uzfwm) { # If pinc exceeds uzfwm
            uzfwc <- uzfwc + pinc # No surface runoff
          } else {
            sur <- pinc + uzfwc - uzfwm # Surface Runoff
            uzfwc <- uzfwm

            ssur <- ssur + (sur * parea)

            # adsur is the amount of surface runoff from that portion of adimp that is not direct runoff
            # addro/pinc is the fraction of adimp currently generating direct runoff
            adsur <- sur * (1 - (addro / pinc))
            ssur <- ssur + adsur * adimp
          }
        }
      }

      adimc <- adimc + pinc - addro - adsur

      if (adimc > (uztwm + lztwm)) {
        addro <- addro + adimc - (uztwm + lztwm)
        adimc <- uztwm + lztwm
      }

      sdro <- sdro + (addro * adimp) # Direct runoff from the additional impervious area

      if (adimc < thresholdZero) adimc <- 0

    } # Loop ending

    # Compute sum and adjust runoff amounts by area over which they are generated ====

    # eused is the ET from PAREA which is 1 - adimp - pctim
    eused <- ET1 + ET2 + ET3
    sif <- sif * parea

    # Separate channel components of baseflow from non-channel component
    tbf <- sbf * parea # tbf is the total baseflow
    bfcc <- tbf * (1 / (1 + side)) # bfcc is baseflow, channel component

    # Groundflow and Surfaceflow
    base <- bfcc # ZBaseflow and interflow are considered ground inflow to the channel
    surf <- roimp + sdro + ssur + sif # Surface flow = direct runoff + surface inflow to the channel

    # ET Module 4 - ET from riparian vegetation ####
    ET4 <- (evapDemand - eused) * riva

    # Check that adimc >= uztws
    if (adimc < uztwc) adimc <- uztwc

    # Total Inflow to Channel for a Timestep ====
    ch_inflow <- surf + base - ET4

    if (ch_inflow <= 0) {
      ET4 <- surf + base
      ch_inflow <- 0
      surf <- 0
      base <- 0
    } else {
      # surface runoff is reduced first
      surf_remainder <- surf - ET4
      surf <- max(0, surf_remainder)

      if (surf_remainder < 0) base <- base + surf_remainder # If surf_remainder is left over then baseflow is reduced
    }

    eused <- eused * parea

    tet <- eused + ET4 + ET5

    # Outputting for timesteps ====
    streamflow.ts[i] <- surf + base
    surfaceflow.ts[i] <- surf
    baseflow.ts[i] <- base
    aet.ts[i] <- tet
    interflow.ts[i] <- sif # Calculated into surface flow, so subtract from surface flow to find surface flow without interflow

    if (preserveState) {
      uztwc_state[i] <- uztwc
      uzfwc_state[i] <- uzfwc
      lztwc_state[i] <- lztwc
      lzfpc_state[i] <- lzfpc
      lzfsc_state[i] <- lzfsc
      adimc_state[i] <- adimc
    }
  }

  if (preserveState) {
    Output$State$uztwc_state <- xts::xts(uztwc_state, order.by = TimeSteps)
    Output$State$uzfwc_state <- xts::xts(uzfwc_state, order.by = TimeSteps)
    Output$State$lztwc_state <- xts::xts(lztwc_state, order.by = TimeSteps)
    Output$State$lzfpc_state <- xts::xts(lzfpc_state, order.by = TimeSteps)
    Output$State$lzfsc_state <- xts::xts(lzfsc_state, order.by = TimeSteps)
    Output$State$adimc_state <- xts::xts(adimc_state, order.by = TimeSteps)
  }

  Output$Output$TotalStreamflow <- xts::xts(streamflow.ts, order.by = TimeSteps)
  Output$Output$Surfaceflow <- xts::xts(surfaceflow.ts, order.by = TimeSteps) # Surfaceflow includes Interflow, so direct runoff is surfaceflow - interflow
  Output$Output$Baseflow <- xts::xts(baseflow.ts, order.by = TimeSteps)
  Output$Output$AET <- xts::xts(aet.ts, order.by = TimeSteps)
  Output$Output$Interflow <- xts::xts(interflow.ts, order.by = TimeSteps)

  if (verbose) {
    verbose.endTime <- Sys.time()
    print(paste0('SAC-SMA Module Run Time: ', format(verbose.endTime - verbose.startTime)))
  }

  return(Output)
}


# SACSMA Ranges ----
hydrology_sacsma.ranges <- function() {
  list(
    # Capacity Threshold
    uztwm <- c(1, 150), # From Blasone, et. al. (2008)
    uzfwm <- c(1, 150), # From Blasone, et. al. (2008)
    lztwm <- c(1, 500), # From Blasone, et. al. (2008)
    lzfpm <- c(1, 1000), # From Blasone, et. al. (2008)
    lzfsm <- c(1, 1000), # From Blasone, et. al. (2008)

    # Recession Parameters
    uzk <- c(0.1, 0.5), # From Blasone, et. al. (2008)
    lzpk <- c(0.0001, 0.25), # From Blasone, et. al. (2008)
    lzsk <- c(0.01, 0.25), # From Blasone, et. al. (2008)

    # Percolation
    zperc <- c(1, 250), # From Blasone, et. al. (2008)
    rexp <- c(0, 5), # From Blasone, et. al. (2008)
    pfree <- c(0, 0.6), # From Blasone, et. al. (2008)

    # Impervious Area
    pctim <- c(0.000001, 0.1), # From Blasone, et. al. (2008)
    adimp <- c(0, 0.4), # From Blasone, et. al. (2008) (max is basically 1 - pctim)

    # Others
    riva <- c(0, 1), # Depends on land use and how deep the roots are (deeper roots = higher riva)
    side <- c(0, 0.5), # Decimal Fraction (needs confirming)
    rserv <- c(0, 0.5) # Decimal Fraction (needs confirming)
  )
}

# SAC-SMA Defaults ----
hydrology_sacsma.defaults <- function() {
  list(
    # Capacity Threshold
    uztwm <- NA,
    uzfwm <- NA,
    lztwm <- NA,
    lzfpm <- NA,
    lzfsm <- NA,

    # Recession Parameters
    uzk <- NA,
    lzpk <- NA,
    lzsk <- NA,

    # Percolation
    zperc <- NA,
    rexp <- NA,
    pfree <- NA,

    # Impervious Area
    pctim <- NA,
    adimp <- NA,

    # Others
    riva <- 0.3,
    side <- 0,
    rserv <- 0.3 # From UTexas GradHydro2001 Presentation file
  )
}




