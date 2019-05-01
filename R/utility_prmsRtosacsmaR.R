#' Translating PRMS Configuration Files to SAC-SMA-DS Format
#'
#' These translates prmsR's data format and PRMS's configuration files to SAC-SMA-DS use.
#'
#' Note: This function is used specifically for California and could
#'
#' @param prms.param prmsR parameter variable (assumes units are in English, except for Elevation)
#' @param projection.from PROJ4 string for CRS (Default: CONUS Albers Projection)
#' @param projection.to PROJ4 string for CRS (Default: WGS84 Google Maps)
#' @param elev_convert Converting elevation from feet to meters
#' @param area_convert Converting area from acres to square meters
#'
#' @return output.configuration Outputs translated configurations
#' @export

# Note to self: This requires the following additional libraries: sp

utility_prsmRtosacsmaR <- function(prms.param,
                           projection.from = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=1,1,-1,0,0,0,0 +units=m +no_defs",
                           projection.to = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs",
                           elev_convert = TRUE,
                           area_convert = TRUE) {

  library('sp')
  library('rgdal')

  RawImport <- data.frame(1:length(prms.param$param$hru_x[[5]]), prms.param$param$hru_x[[5]], prms.param$param$hru_y[[5]],
                          prms.param$param$hru_area[[5]], prms.param$param$hru_elev[[5]])

  colnames(RawImport) <- c('HRUID', 'xalbers', 'yalbers', 'area', 'elevation')

  Coordinates_Albers <- cbind(RawImport$xalbers, RawImport$yalbers)
  colnames(Coordinates_Albers) <- c('xalbers', 'yalbers')

  Coordinates_Albers <- sp::coordinates(Coordinates_Albers)
  PointLabels <- data.frame(HRUID=RawImport$HRUID, area=RawImport$area, elevation=RawImport$elevation)
  # Assuming input is EPSG:5070 (What I built PRMS Models in)
  SPoint_Albers <- sp::SpatialPointsDataFrame(Coordinates_Albers, PointLabels, proj4string=sp::CRS(projection.from))
  # Converting to EPSG:3857 - Google Maps
  SPoint_latlon <- sp::spTransform(SPoint_Albers, CRSobj = sp::CRS(projection.to))

  Coordinates_latlon <- sp::coordinates(SPoint_latlon)
  colnames(Coordinates_latlon) <- c('X_lon', 'Y_lat')

  Output <- data.frame(Coordinates_latlon, SPoint_latlon)
  Output <- Output[,1:(length(colnames(Output)) - 3)]

  if(elev_convert) {
    Output$elevation <- Output$elevation * 0.3048
  }

  if(area_convert) {
    Output$area <- Output$area * 4046.86
  }

  return(Output)
}
