#' N Daylight Hours Model
#'
#' Calculating the N daynight hours.
#'
#' @param Latitude Latitude in Degrees
#' @param JDate Julian Day of the Year
#' @param LeapYear Leap Year Flag
#'
#' @return Daylight hours
#' @details For details see Allen (1998)
#' @export

daylight_N <- function(JDate, Latitude, LeapYear = FALSE) {

  # Declination
  if (LeapYear) Year <- 366 else Year <- 365
  delta <- 0.409 * base::sin(((2 * pi)/Year) * JDate - 1.39)

  # Sunet hour angle
  omega <- base::acos(-base::tan(Latitude * pi/180) * base::tan(delta))

  N = (24/pi) * omega

  return(N)
}
