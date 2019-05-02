#' Area absolute to ratio of total area
#'
#' Converts area from absolute values to a ratio of the total area
#'
#' @param Area A vector list of absolute area values
#'
#' @return A vector of ratio that all adds up to 1
#' @export

utility_areaRatio <- function(Area) {
  return(Area/sum(Area))
}
