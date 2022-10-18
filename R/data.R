#' Mutational counts for 21 breast cancer patients
#'
#' A dataset containing the mutational counts of 96 different mutation types
#' for 21 breast cancer patients.
#'
#'
#' @format A data matrix with 21 rows and 96 columns.
#'
#' @source \url{https://doi.org/10.1016/j.cell.2012.04.023}
"BRCA21"

#' Simulated data with Poisson noise
#'
#' 100 simulated mutational counts data sets.
#' Each data set contains 100 patients, 96 different mutation types.
#'
#' The true number of underlying signatures for each data set equals to 5.
#'
#' We added Poisson noise to the simulated mutational counts data.
#'
#' @format A list of 100 elements where each element is a matrix with 100
#' rows and 96 columns.
#'
#' @source Details on the simulation setup are available at \url{ourpaper}
"SimulatedDataPoisson"

#' Simulated data with Negative Binomial noise and alpha = 10
#'
#' 100 simulated mutational counts data sets.
#' Each data set contains 100 patients, 96 different mutation types.
#'
#' The true number of underlying signatures for each data set equals to 5.
#'
#' We added Negative Binomial noise to the simulated mutational counts
#' data with alpha  = 10 (i.e. high overdispersion).
#'
#' @format A list of 100 elements where each element is a matrix with 100
#' rows and 96 columns.
#'
#' @source Details on the simulation setup are available at \url{ourpaper}
"SimulatedDataNBalpha10"

#' Simulated data with Negative Binomial noise and alpha = 200
#'
#' 100 simulated mutational counts data sets.
#' Each data set contains 100 patients, 96 different mutation types.
#'
#' The true number of underlying signatures for each data set equals to 5.
#'
#' We added Negative Binomial noise to the simulated mutational counts
#' data with alpha  = 200 (i.e. intermediate overdispersion).
#'
#' @format A list of 100 elements where each element is a matrix with 100
#' rows and 96 columns.
#'
#' @source Details on the simulation setup are available at \url{ourpaper}
"SimulatedDataNBalpha200"

#' Simulated data with patient specific Negative Binomial noise
#'
#' 100 simulated mutational counts data sets. Each data set contains 100 patients, 96 different mutation types.
#'
#' The true number of underlying signatures for each data set equals to 5.
#'
#' We added Negative Binomial noise to the simulated mutational counts
#' data with alpha values inspired from the ones of the patients in
#' \url{https://doi.org/10.1016/j.cell.2012.04.023} (i.e. patient
#' specific overdispersion).
#'
#' @format A list of 100 elements where each element is a matrix with 100
#' rows and 96 columns.
#'
#' @source Details on the simulation setup are available at \url{ourpaper}
"SimulatedDataNBn"
