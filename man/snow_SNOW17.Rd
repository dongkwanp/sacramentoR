% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/snow_SNOW17.R
\name{snow_SNOW17}
\alias{snow_SNOW17}
\title{SNOW-17 Snow Module}
\usage{
snow_SNOW17(Param, Prcp, Tavg, Elevation, InitialState = c(0, 0, 0, 0),
  dtt = 24, dtp = 24, calcA_v = FALSE, meltFlag = FALSE,
  preserveInput = FALSE, verbose = FALSE)
}
\arguments{
\item{Param}{snow module parameter list}

\item{Prcp}{Precipitation continuous time-series (mm) (xts object) (same duration as Tavg)}

\item{Tavg}{Average Temperature continuous time-series (deg C) (xts object) (same duration as Prcp)}

\item{Elevation}{Elevation (meters)}

\item{InitialState}{Initial state vector (W_i, ATI, W_q, Deficit)}

\item{dtt}{Constant time interval of temperature data (hours) (Default: 24 Hours)}

\item{dtp}{Constant time interval of precipitation data (hours) (Default: 24 Hours)}

\item{calcA_v}{If Latitude is above 60 degrees then calculate A_v for each time-step}

\item{meltFlag}{Output a vector time-series melt flag}

\item{preserveInput}{Flag to preserve the input as part of the output}

\item{verbose}{Verbose Flag}
}
\value{
A list with various time-series outputs (xts objects)
}
\description{
Calculates snow via SNOW17 Model.  Current version is configured for Anderson 1973 version with some 2006 elements
}
\details{
For details see Anderson (2006) and Anderson (1973)
}
