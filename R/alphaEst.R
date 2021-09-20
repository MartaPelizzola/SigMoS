#' Estimation of the overdispersion parameter for Negative Binomial non-negative matrix factorization
#'
#' Empirical Bayes estimation of the overdispersion parameter for Negative Binomial non-negative matrix factorization.
#' The overdispersion parameter can be either vector of length one or vector of length no. of patients if patient-specific overdispersion is used.
#'
#'
#' @param data Numeric matrix of mutational counts data. Matrix size: no. of mutation types x no. of patients.
#' @param k Number of signatures to be used for the non-negative matrix factorization
#' @param patient_specific Logical. If TRUE patient-specific overdispersion is used in the Negative Binomial model.
#'
#'
#' @return Overdispersion parameter. Either vector of length one or vector of length no. of patients if patient-specific overdispersion is used.
#'
#' @export
#'
#'
alphaEst <- function(data, k, patient_specific = FALSE){
  if (k!=round(k)){
    stop("The number of signatures must be an integer.")
  }
  if(is.null(data)){
    stop("The data set of the mutational counts is missing.")
  }
  if(is.null(k)){
    stop("A value for the number of signatures to be estimated is missing.")
  }
  res_p <- NMFPoisEM(data,k,tol = 1e-5,
                     seed = sample(100000,1))
  h_p <- res_p$P
  w_p <- res_p$E

  a_empbayes <- data/(h_p%*%w_p)
  if (patient_specific){
    alpha <- 1/apply(a_empbayes, 2, var)
  } else {
    alpha <- 1/var(as.vector(a_empbayes))
  }
  return(alpha)
}
