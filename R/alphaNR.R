#' @title Estimation of the overdispersion parameter for Negative Binomial non-negative matrix factorization
#'
#' @description Likehood estimation of the dispersion parameter in negative binomial using Newton-Raphson.
#' The overdispersion parameter can be either vector of length one or vector of length no. of patients if patient-specific overdispersion is used.
#'
#'
#' @param data Numeric matrix of mutational counts data. Matrix size: no. of patients x no. of mutation types.
#' @param k Number of signatures to be used for the non-negative matrix factorization
#' @param patient_specific Logical. If TRUE patient-specific overdispersion is used in the Negative Binomial model.
#'
#'
#' @return Overdispersion parameter. Either vector of length one or vector of length no. of patients if patient-specific overdispersion is used.
#'
#' @export
#'
#'
alphaNR <- function(data, k=NULL, patient_specific = FALSE){
  if (k!=round(k)){
    stop("The number of signatures must be an integer.")
  }
  if(is.null(data)){
    stop("The data set of the mutational counts is missing.")
  }
  data <- t(data)
  if(is.null(k)){
    stop("A value for the number of signatures to be estimated is missing.")
  }
  if(length(k)!=1){
    stop("'k' has length larger than 1.")
  }
  res_p <- NMFPois(data,k,tol = 1e-2,
                     seed = sample(100000,1))
  h_p <- res_p$P
  w_p <- res_p$E

  # differentiated once
  neglikdiff1 = function(alpha, data, estimate){
    sum(digamma(data + alpha) - digamma(alpha) - data/(alpha+estimate) - alpha/(alpha+estimate) + log(alpha/(alpha+estimate)) + 1)
  }

  # differentiated twice
  neglikdiff2 = function(alpha, data, estimate){
    sum(trigamma(data + alpha) - trigamma(alpha) + data/(alpha+estimate)^2 + 1/alpha - 2/(alpha+estimate) + alpha/(alpha+estimate)^2)
  }

  NR_alpha = function(data,estimate){
    alpha <- 1/var(data/estimate)
    alphaold = alpha + 5
    for(i in 1:10){
      alpha = alpha - neglikdiff1(alpha, data = data, estimate = estimate)/neglikdiff2(alpha, data = data, estimate = estimate)
      if(!(alpha > 0)){ alpha = runif(1,1,10)}

      if(abs(alpha - alphaold) < 0.01) break
      alphaold = alpha
    }
    return(alpha)
  }

  if(patient_specific){
    alpha = numeric(ncol(data))
    estimate = h_p%*%w_p
    for(i in 1:ncol(data)){
      alpha[i] = NR_alpha(data[,i],estimate[,i])
    }
    }else{
    data = as.vector(data)
    estimate = as.vector(h_p%*%w_p)

    alpha = NR_alpha(data,estimate)
  }

  return(alpha)
}

alphaNR2 <- function(data, k=NULL, patient_specific = FALSE){
  if (k!=round(k)){
    stop("The number of signatures must be an integer.")
  }
  if(is.null(data)){
    stop("The data set of the mutational counts is missing.")
  }
  if(is.null(k)){
    stop("A value for the number of signatures to be estimated is missing.")
  }
  if(length(k)!=1){
    stop("'k' has length larger than 1.")
  }
  res_p <- NMFPois(data,k,tol = 1e-2,
                     seed = sample(100000,1))
  h_p <- res_p$P
  w_p <- res_p$E

  # differentiated once
  neglikdiff1 = function(alpha, data, estimate){
    sum(digamma(data + alpha) - digamma(alpha) - data/(alpha+estimate) + estimate/(alpha+estimate) + log(alpha/(alpha+estimate)))
  }

  # differentiated twice
  neglikdiff2 = function(alpha, data, estimate){
    sum(trigamma(data + alpha) - trigamma(alpha) + data/(alpha+estimate)^2 + (estimate/alpha)/(alpha+estimate) - estimate/(alpha+estimate)^2)
  }

  NR_alpha = function(data,estimate){
    alpha <- 1/var(data/estimate)
    alphaold = alpha + 5
    for(i in 1:20){
      alpha = alpha - neglikdiff1(alpha, data = data, estimate = estimate)/neglikdiff2(alpha, data = data, estimate = estimate)
      if(!(alpha > 0)){ alpha = runif(1,1,10)}

      if(abs(alpha - alphaold) < 0.01) break
      alphaold = alpha
    }
    return(alpha)
  }

  if(patient_specific){
    alpha = numeric(ncol(data))
    estimate = h_p%*%w_p
    for(i in 1:ncol(data)){
      alpha[i] = NR_alpha(data[,i],estimate[,i])
    }
  }else{
    data = as.vector(data)
    estimate = as.vector(h_p%*%w_p)

    alpha = NR_alpha(data,estimate)
  }

  return(alpha)
}
