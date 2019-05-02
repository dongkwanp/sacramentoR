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
#' @param debug Debug Flag
#'
#' @return A list of various time-series as an xts object
#' @details For details see Lohmann, Nolte-Holube, and Raschke (1996)
#' @export

routing_lohmann <- function(Param, directInflow, baseInflow, flowLength, Outlet = FALSE, WatershedCharacteristics = c(12, 96, 3600), preserveInput = FALSE, verbose = FALSE) {

  # Preprocessing ----
  verbose.startTime <- Sys.time()
  verbose.timeStepTotal <- length(Prcpts)

  if (verbose) {
    print('Initializing Lohmann Routing Module...')
  }

# Environment Preparation ----
  Output <- list()
  Output$Output <- list()

  if (preserveInput) {

    if (verbose) {
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

  # Defining Functions ----
  hruh_fun <- function(x) {
    shape <- get('N', envir = parent.frame())
    scale <- (1/get('K', envir = parent.frame()))
    return(stats::dgamma(x, shape = shape, scale = scale))
  }

  # Allocating Memory for Outputs and States of the Model ----
  # Initializing memory for model logics
  UHriver <- rep(0, UH_DAY)

  # Model Logic ----
  # If watershed outlet ====
  if (Outlet) UHriver[1] <- 1 else {
    # if not...
    # Calculate Green's function to solve Saint-Venant Equation ====
    t <- 0
    uhm_grid <- rep(NA, LE)

    for (k in 1:LE) {
      t <- t + DT

      pot <- (((VELO * t) - flowLength) ^ 2) / (4 * DIFF * t)

      if (pot <= 69) H <- (flowLength / (2 * t * sqrt(pi * t * DIFF) * exp(-pot))) else H <- 0

      uhm_grid[k] <- H
    }

    if (sum(uhm_grid) == 0) uhm_grid[1] <- 1 else uhm_grid <- (uhm_grid/sum(uhm_grid))

    UHM <- uhm_grid

    FR <- data.frame(var1 = rep(0, TMAX), var2 = rep(0, TMAX))
    FR[1:24,1] <- 1/24

    for (t in 1:TMAX) for (L in 1:(TMAX + 24)) if ((t - L) > 0) FR$var2[t] <- FR$var2[t] + FR$var1[t-L] * UHM[L]

    for (t in 1:UH_DAY) UHriver[t] <- sum(FR$var2[((24 * t) - 23):(24 * t)])

  }

  # HRU's UH Represented by Gamma Distribution ====
  UH_HRU_direct <- rep(0, KE)

  for (i in 1:KE) {
    UH_HRU_direct[i] <- integrate(hruh_fun, (24 * (i - 1)), (24 * i))$value
  }
  UH_HRU_base <- rep(0, KE)
  UH_HRU_base[1] <- 1

  # Combined UH for HRU's response at the watershed outlet ====
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

  # Make convolution for watershed outlet total flow ====
  directFlow <- rep(0, length(directInflow))
  baseFlow <- rep(0, length(directInflow))

  for (i in 1:length(directInflow)) {
    for (j in 1:(KE + UH_DAY - 1)) {
      if ((i - j + 1) >= 1) {
        directFlow[i] <- directFlow[i] + UH_direct[j] * directInflow[(i - j + 1)]
        baseFlow[i] <- baseFlow[i] + UH_base[j] * baseInflow[(i - j + 1)]
      }

    }
  }

  TotalStreamflow <- directFlow + baseFlow

  # Writing to Output ----
  Output$Output$directFlow <- directFlow
  Output$Output$baseFlow <- baseFlow
  Output$Output$TotalStreamflow <- TotalStreamflow

}

# Lohmann Ranges ----
routing_lohmann.ranges <- function() {
  list(
    N <- c(NA, NA),
    K <- c(NA, NA),
    VELO <- c(NA, NA),
    DIFF <- c(NA, NA)
  )
}

# Lohmann Defaults ----
routing_lohmann.ranges <- function() {
  list(
    N <- NA,
    K <- NA,
    VELO <- NA,
    DIFF <- NA
  )
}
