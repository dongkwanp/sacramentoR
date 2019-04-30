#' A_v Calculator
#'
#' Calculates the seasonal variation adjustment (A_v) based on Anderson 2006 equation 7
#'
#' @param Date Date (as Date)
#' @param Latitude Latitude (Degrees) (Default value triggers calculation)
#'
#' @return A_v
#' @details For details see Anderson (2006)
#' @export

utility_A_v <- function(Date, Latitude = 60) {
  if (Latitude < 54) A_v <- 1.0
  else {
    JDate <- as.numeric(format(Date, '%j'))
    Sept24 <- as.numeric(format(as.Date(paste0(format(Date, '%Y'), '-09-24')), '%j'))
    Mar18 <- as.numeric(format(as.Date(paste0(format(Date, '%Y'), '-03-18')), '%j'))
    Apr27 <- as.numeric(format(as.Date(paste0(format(Date, '%Y'), '-04-27')), '%j'))
    Aug15 <- as.numeric(format(as.Date(paste0(format(Date, '%Y'), '-08-15')), '%j'))

    if (JDate <= Mar18) {
      A_v <- 0
    } else if (JDate < Apr27 && JDate > Mar18) {
      A_v <- (JDate - Mar18)/(Apr27 - Mar18)
    } else if (JDate >= Apr27 && JDate <= Aug15) {
      A_v <- 1
    } else if (JDate < Sept24 && JDate > Aug15) {
      A_v <- 1 - (JDate - Aug15)/(Sept24 - Aug15)
    } else if (JDate >= Sept24) {
      A_v <- 0
    } else {
      stop('Error in Date')
    }
  }
  return(A_v)
}
