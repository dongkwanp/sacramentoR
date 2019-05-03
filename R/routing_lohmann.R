#' Lohmann Routing Module
#'
#' Calculates the Lohmann Routing model coupled with land surface parameterization
#'
#' @param Param routing module parameter list
#' @param directInflow Timeseries of Direct Inflow (mm) (xts object)
#' @param baseInflow Timeseries of Base Inflow (mm) (xts object)
#' @param flowLength Travel distance from catchment outlet to watershed outlet (meters)
#' @param Outlet Flag to determine if this is the outlet
#' @param WatershedCharacteristics Watershed Characteristics Vector c(KE, UH_DAY, DT) (not recommended to change this unless you know what you're doing)
#' @param preserveInput Flag to preserve the input as part of the output
#' @param verbose Verbose Flag
#'
#' @return A list of various time-series as an xts object
#' @details For details see Lohmann, Nolte-Holube, and Raschke (1996)
#' @export

routing_lohmann <- function(Param, directInflow, baseInflow, flowLength, Outlet = FALSE, WatershedCharacteristics = c(12, 96, 3600), preserveInput = FALSE, verbose = FALSE) {

  # Preprocessing ----
  Dates <- zoo::index(directInflow)
  verbose.startTime <- Sys.time()
  verbose.timeStepTotal <- length(directInflow)

  if (verbose) {
    print('Initializing Lohmann Routing Module...')
  }

# Environment Preparation ----
  Output <- list()
  Output$Output <- list()

  if (preserveInput) {

    if (verbose) {
      verbose.subStartTime <- Sys.time()
      print('Saving Inputs with Outputs...')
    }

    Output$Input <- list()
    Output$Input$TimeSeries <- list()
    Output$Input$TimeSeries$directInflow <- directInflow
    Output$Input$TimeSeries$baseInflow <- baseInflow
    Output$Input$flowLength <- flowLength
    Output$Input$Param <- Param
  }


  # Parameter Defining ----
  if (verbose) {
    print('Assigning Parameter Values...')
  }

  if (is.list(Param)) {
    N <- Param$N # Catchment's UH Shape Parameter (N)
    K <- Param$K # Catchment's UH Scale Parameter (K)
    VELO <- Param$VELO # Wave Velocity in the Linearized Saint-Venant Equation (m/s)
    DIFF <- Param$DIFF # Diffusivity in the Linearized Saint-Venant Equation (m2/s)
  } else stop('Param variable is not in an acceptable format')

  if (is.vector(WatershedCharacteristics)) {
    KE <- WatershedCharacteristics[1] # Base Time for HRU UH (day)
    UH_DAY <- WatershedCharacteristics[2] # Base time for river routing UH (day)
    DT <- WatershedCharacteristics[3] # Time step in second for solving Saint-Venant equation (This will affect TMAX)
  } else if (is.list(WatershedCharacteristics)) {
    KE <- WatershedCharacteristics$KE
    UH_DAY <- WatershedCharacteristics$UH_DAY
    DT <- WatershedCharacteristics$DT
  } else stop ('WatershedCharacteristics is not in an acceptable format')

  # Parameter Range Testing ====
  # This is a TODO

  # Defining Constants ----
  TMAX <- UH_DAY * 24 # Base time of river routing UH in hour because DT is for an 1 hour
  LE <- 48 * 50 # Base time (hr) for Green function values

  directInflow.vector <- as.vector(directInflow)
  baseInflow.vector <- as.vector(baseInflow)

  # Defining Functions ----
  if (verbose) {
    print('Loading internal functions...')
  }

  hruh_fun <- function(x) return(stats::dgamma(x, shape = N, scale = 1/K))

  # Allocating Memory for Outputs and States of the Model ----
  # Initializing memory for model logics
  UHriver <- rep(0, UH_DAY)

  # Model Logic ----
  # If watershed outlet ====
  if (verbose) {
    print('Beginning model logic...')
  }
  if (Outlet) UHriver[1] <- 1 else {
    # if not...
    # Calculate Green's function to solve Saint-Venant Equation ====
    if (verbose) {
      verbose.subStartTime <- Sys.time()
      print('Calculating Green\'s function...')
    }

    t <- 0
    uhm_grid <- rep(NA, LE)

    if (verbose) {
      verbose.subforStartTime <- Sys.time()
      print('Calculating Green\'s function: Loop 1...')
    }

    for (k in 1:LE) {

      t <- t + DT

      pot <- (((VELO * t) - flowLength) ^ 2) / (4 * DIFF * t)

      if (pot <= 69) H <- ((flowLength / (2 * t * sqrt(pi * t * DIFF))) * exp(-pot)) else H <- 0

      uhm_grid[k] <- H
    }

    if (verbose) {
      verbose.subforEndTime <- Sys.time()
      print(paste0('Green\'s Function: Loop 1 took: ', format(verbose.subforEndTime - verbose.subforStartTime)))
    }

    if (sum(uhm_grid) == 0) uhm_grid[1] <- 1 else uhm_grid <- (uhm_grid/sum(uhm_grid))

    UHM <- uhm_grid

    FR <- matrix(data = 0, nrow = TMAX, ncol = 2)
    FR[1:24,1] <- 1/24

    if (verbose) {
      verbose.subforStartTime <- Sys.time()
      print('Calculating Green\'s function: Loop 2...')
    }

    for (t in 1:TMAX) {
      L <- 1:(TMAX + 24)
      L <- L[(t - L) > 0]

      FR[t, 2]<- FR[t,2] + sum(FR[t - L, 1] * UHM[L])

      # for (L in 1:(TMAX + 24)) if ((t - L) > 0) FR$var2[t] <- FR$var2[t] + FR$var1[t-L] * UHM[L]
    }

    if (verbose) {
      verbose.subforEndTime <- Sys.time()
      print(paste0('Green\'s Function: Loop 2 took: ', format(verbose.subforEndTime - verbose.subforStartTime)))
    }

    if (verbose) {
      verbose.subforStartTime <- Sys.time()
      print('Calculating Green\'s function: Loop 3...')
    }

    UHriver <- sapply(1:UH_DAY, function(t) sum(FR[(24 * t - 23):(24 * t), 2]))

    #for (t in 1:UH_DAY) UHriver[t] <- sum(FR$var2[((24 * t) - 23):(24 * t)])

    if (verbose) {
      verbose.subforEndTime <- Sys.time()
      print(paste0('Green\'s Function: Loop 3 took: ', format(verbose.subforEndTime - verbose.subforStartTime)))
    }

  }
  if (verbose) {
    verbose.subEndTime <- Sys.time()
    print(paste0('Entire Green\'s function took: ', format(verbose.subEndTime - verbose.subStartTime)))
  }

  # HRU's UH Represented by Gamma Distribution ====
  if (verbose) {
    verbose.subStartTime <- Sys.time()
    print('Calculating HRU\'s UH through Gamma Distribution...')
  }

  UH_HRU_direct <- rep(0, KE)

  for (i in 1:KE) {
    UH_HRU_direct[i] <- integrate(hruh_fun, (24 * (i - 1)), (24 * i))$value
  }

  UH_HRU_base <- rep(0, KE)
  UH_HRU_base[1] <- 1

  if (verbose) {
    verbose.subEndTime <- Sys.time()
    print(paste0('That took: ', format(verbose.subEndTime - verbose.subStartTime)))
  }

  # Combined UH for HRU's response at the watershed outlet ====
  if (verbose) {
    verbose.subStartTime <- Sys.time()
    print('Combining UH for HRU\'s response...')
  }

  UH_direct <- rep(0, (KE + UH_DAY - 1))
  UH_base <- rep(0, (KE + UH_DAY - 1))

  for (k in 1:KE) {
    for (u in 1:UH_DAY) {
      UH_direct[(k + u - 1)] <- UH_direct[(k + u - 1)] + UH_HRU_direct[k] * UHriver[u]
      UH_base[(k + u - 1)] <- UH_base[(k + u - 1)] + UH_HRU_base[k] * UHriver[u]
    }
  }

  UH_direct <- (UH_direct/sum(UH_direct))
  UH_base <- (UH_base/sum(UH_base))

  if (verbose) {
    verbose.subEndTime <- Sys.time()
    print(paste0('That took: ', format(verbose.subEndTime - verbose.subStartTime)))
  }


  # Make convolution for watershed outlet total flow ====
  if (verbose) {
    verbose.subStartTime <- Sys.time()
    print('Calculating watershed total flow...')
  }

  directFlow <- rep(0, length(directInflow.vector))
  baseFlow <- rep(0, length(directInflow.vector))

  for (i in 1:length(directInflow)) {

    j <- 1:(KE + UH_DAY - 1)
    j <- j[((i - j) + 1) >= 1]

    directFlow[i] <- directFlow[i] + sum(UH_direct[j] * directInflow.vector[(i - j + 1)])
    baseFlow[i] <- baseFlow[i] + sum(UH_base[j] * baseInflow.vector[(i - j + 1)])

  }



  TotalStreamflow <- directFlow + baseFlow

  if (verbose) {
    verbose.subEndTime <- Sys.time()
    print(paste0('That took: ', format(verbose.subEndTime - verbose.subStartTime)))
  }



  # Writing to Output ----
  if (verbose) {
    print('Writing to output...')
  }

  Output$Output$surfaceFlow <- xts::xts(directFlow, order.by = Dates)
  Output$Output$baseFlow <- xts::xts(baseFlow, order.by = Dates)
  Output$Output$TotalStreamflow <- xts::xts(TotalStreamflow, order.by = Dates)


  if (verbose){
    verbose.endTime <- Sys.time()
    print(paste0('Lohmann Module Run Time: ', format(verbose.endTime - verbose.startTime)))
  }

  return(Output)
}

# Lohmann Ranges ----
routing_lohmann.ranges <- function() {
  list(
    N <- c(NA, NA),
    K <- c(0, 1),
    VELO <- c(1, 3), # VIC Documentation for a basin in Germany
    DIFF <- c(200, 4000) # VIC Documentation for a basin in Germany
  )
}

# Lohmann Defaults ----
routing_lohmann.ranges <- function() {
  list(
    N <- NA,
    K <- NA,
    VELO <- 2, # Just a number
    DIFF <- 2000 # Just a number
  )
}
