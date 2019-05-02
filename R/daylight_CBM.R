#' CBM Daylight Model
#'
#' Calculating the N daynight hours.
#'
#' @param Latitude Latitude in Degrees
#' @param Date Date vector
#' @param p Daylength coefficient (Default: Constant based on Astronomical calendar)
#'
#' @return Daylight hours
#' @details For details see Forsythe (1995) and Schoolfield (1982)
#' @export

daylight_CBM <- function(Date, Latitude, p = 0.833) {

  JDate <- as.numeric(format(as.Date(Date), '%j'))

  # Predicting angle of revolution angle (theta) from the day of the year (JDate)
  theta <- 0.2163108 + 2 * base::atan(0.9671396 * base::tan(0.0086 * (JDate - 186)))

  # Predicting sun's declination angle (phi)
  phi <- base::asin(0.39795 * cos(theta))



  # Daylight Hour
  D = 24 - (24/pi) * base::acos((base::sin(((pi/180) * p)) + (base::sin((pi/180) * Latitude) * base::sin(phi)))
                                /(base::cos((pi/180) * Latitude) * base::cos(phi)))

  return(D)
}
