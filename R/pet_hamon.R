#' Hamon Method for Potential Evapotranspiration
#'
#' Calculation of potential evapotranspiration based on the Hamon method
#'
#' @param Coeff Hamon PET Model Coefficient (Default: 1 from Lu (2005))
#' @param Date Date (Time-series vector or as a Date Object)
#' @param Tavg Average Temperature (deg C) (Time-series vector or numeric)
#' @param Latitude Latitude (Degrees)
#' @param dayFUN Daylight hours function (Defaults to CBM Model)
#'
#' @return Hamon Potential Evapotranspiration Value
#' @details For details see Lu, Sun, McNulty, and Amatya (2005) and Haith and Shoemaker (1987)
#' @export

pet_hamon <- function(par, Date, Tavg, Latitude, dayFUN = sacramentoR::daylight_CBM) {

  JDate <- format(as.Date(Date), '%j')

  daylight <- dayFUN(JDate, Latitude)

  # ESAT is the saturated vapor pressure (mb) at the given T
  ESAT <- 6.108 * base::exp(17.26939 * (Tavg / (Tavg + 237.3)))

  RhoSAT <- 216.7 * (ESAT / (Tavg + 273.3))

  # Cumulation of everything plus dividing by 12
  Hamon <- par * 0.1651 * (daylight / 12) * RhoSAT

  return(Hamon)
}


# Hamon Ranges ----
pet_hamon.ranges <- function() {
  list(
    par <- c(0.8, 3) # Lu 2005 suggests default of 1.2 but other docs suggest 1
  )
}

pet_hamon.defaults <- function() {
  list(
    par <- 1 # Default: 1.2 from Lu (2005)
  )
}



