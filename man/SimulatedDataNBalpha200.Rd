% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{SimulatedDataNBalpha200}
\alias{SimulatedDataNBalpha200}
\title{Simulated data with Negative Binomial noise and alpha = 200}
\format{
A list of 100 elements where each element is a matrix with 100
rows and 96 columns.
}
\source{
Details on the simulation setup are available at \url{ourpaper}
}
\usage{
SimulatedDataNBalpha200
}
\description{
100 simulated mutational counts data sets.
Each data set contains 100 patients, 96 different mutation types.
}
\details{
The true number of underlying signatures for each data set equals to 5.

We added Negative Binomial noise to the simulated mutational counts
data with alpha  = 200 (i.e. intermediate overdispersion).
}
\keyword{datasets}
