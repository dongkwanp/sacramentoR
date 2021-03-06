% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hydrology_sacsma.R
\name{hydrology_sacsma}
\alias{hydrology_sacsma}
\title{Sacramento Soil Moisture Accounting Module (SAC-SMA)}
\usage{
hydrology_sacsma(Param, Prcpts, PETts, InitialState = c(0, 0, 500, 500,
  500, 0), thresholdZero = 1e-05, ninc_min = 20,
  preserveInput = FALSE, preserveState = FALSE, verbose = FALSE)
}
\arguments{
\item{Param}{hydrology module parameter list}

\item{Prcpts}{Precipitation/Snow Melt continuous time-series (mm) (xts object) (same duration as PETts)}

\item{PETts}{Potential Evapotranspiration continuous time-series (deg C) (xts object) (same duration as Prcpts)}

\item{InitialState}{Initial state vector (uztwc, uzfwc, lztwc, lzfsc, lzfpc, adimc)}

\item{thresholdZero}{Threshold for Zero}

\item{ninc_min}{Minimum number of time increments that interval is divided into (Default: 20)}

\item{preserveInput}{Flag to preserve the input as part of the output}

\item{preserveState}{Flag to preserve state of the model at each time-step}

\item{verbose}{Verbose Flag}
}
\value{
A list of various time-series as an xts object
}
\description{
Calculates the Sacramento Soil Moisture Accounting (SAC-SMA) Model
}
\details{
For details see NOAA rfs:23sacsma.wpd (2002) and Blasone, et. al. (2008)
}
