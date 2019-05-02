#' N Daylight Hours Model
#'
#' Calculating the N daynight hours.
#'
#' @param Latitude Latitude in Degrees
#' @param Date Date vector
#'
#' @return Daylight hours
#' @details For details see Allen (1998)
#' @export

daylight_N <- function(Date, Latitude) {

  JDate <- as.numeric(format(as.Date(Date), '%j'))
  LeapYear <- utility_isLeapYear(as.numeric(format(Date, '%Y')))
  Year <- sapply(LeapYear, function(x) if(x) return(366) else return(365))

  delta <- 0.409 * base::sin(((2 * pi)/Year) * JDate - 1.39)

  # Sunet hour angle
  omega <- base::acos(-base::tan(Latitude * pi/180) * base::tan(delta))

  N = (24/pi) * omega

  return(N)
}
