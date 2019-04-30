#' Leap Year Checker
#'
#' Checks to see if it's a leap year
#'
#' @param Year Year
#'
#' @return True if leap year
#' @export

utility_isLeapYear <- function(Year) {
  return(((Year %% 4 == 0) && (Year %% 100 != 0)) || (Year %% 400 == 0))
}
