#' Hamon Method for Potential Evapotranspiration
#'
#' Calculation of potential evapotranspiration based on the Hamon method
#'
#' @param Coeff Hamon PET Model Coefficient
#' @param JDate Julian Date
#' @param Tavg Average Temperature (deg C)
#' @param Latitude Latitude (Degrees)
#' @param dayFUN Daylight hours function (Defaults to CBM Model)
#'
#' @return Hamon Potential Evapotranspiration Value
#' @details For details see Lu, Sun, McNulty, and Amatya (2005) and Haith and Shoemaker (1987)
#' @export

pet_hamon <- function(par, JDate, Tavg, Latitude, dayFUN = sacramentoR::daylight_CBM, ...) {

  daylight <- dayFUN(JDate, Latitude)

  # ESAT is the saturated vapor pressure (mb) at the given T
  ESAT <- 6.108 * base::exp(17.26939 * (Tavg / (Tavg + 237.3)))

  RhoSAT <- 216.7 * (ESAT / (Tavg + 273.3))

  # Cumulation of everything plus dividing by 12
  Hamon <- par * 0.1651 * (daylight / 12) * RhoSAT

  return(Hamon)
}
