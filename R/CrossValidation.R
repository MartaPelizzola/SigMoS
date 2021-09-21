#' @import foreach doParallel


#' @title Cross-validation for the number of mutational signatures
#'
#' @description Cross-validation algorithm for a given number of signatures k using either Poisson or Negative Binomial model.
#'
#'
#' @param data Numeric matrix of mutational counts data. Matrix size: no. of mutation types x no. of patients.
#' @param k Number of signatures to be used in the cross-validation.
#' @param n_iterations Number of iterations for cross-validation procedure. Default is 100.
#' @param method Non-negative matrix factorization model: 'Poisson' or 'NB' can be used. Default is 'NB' corresponding to the Negative Binomial model.
#' @param cost_f Cost function to be used for the cross-validation error. Default is GKL (Generalized Kullback-Leibler). Alternatives are: Frobenius norm ("Frobenius") and Itakura-Saito divergence ("IS").
#' @param size_train Size of the train set. Default is 0.9, i.e. 90% of the patients are used for training and 10% for testing in the cross-validation.
#' @param patient_specific Logical. If TRUE patient-specific overdispersion is used in the Negative Binomial model.
#' @param tol Threshold for convergence of the EM and MM algorithms for estimating mutational signatures. Default is 1e-5.
#'
#'
#' @return List of length two. 'cost' is a vector of length n_iterations where the cost for each iteration is stored.
#' cost_k is the median of the values in 'cost' and it is used for comparison to choose the best number of signatures for a given data set.
#'
#' @export

CrossValidation <- function(data,k=NULL,n_iterations=100,method = "NB", cost_f="GKL",size_train=0.9,patient_specific=FALSE,tol = 1e-5){
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
    stop("'k' has length larger than 1. Cross-validation is performed for one value of the number of signatures 'k' at a time.")
  }

  if (!(method %in% c("NB", "Poisson"))){
    method <- "NB"
    warning("Method for non-negative matrix factorization is set to Negative Binomial 'NB'.")
  }

  if(method=="NB"){
    alpha <- alphaEst(data,k,patient_specific)
    res_nb <- NMFNBMM(M=data, N=k, alpha = alpha)

    cores=detectCores()
    cl <- makeCluster(cores[1]-1)
    registerDoParallel(cl)

    cost <- foreach(i=1:n_iterations, .combine=c, .packages=c('SQUAREM', 'SigModeling'), .export = ls(globalenv())) %dopar% {
      n <- ncol(data)
      train_set = sample(1:n, size_train*n)
      res_train <- NMFNBMM(M=data[,train_set], N=k, alpha = alpha)
      h_train <- res_train$P

      ord <- corrSignatures(h_train,res_nb$P)

      if (cost_f=="GKL"){
        tmp_cost <- gklDiv(as.vector(data[,-train_set]), as.vector(h_train %*% res_nb$E[ord,-train_set]))
      } else if (cost_f=="Frobenius"){
        tmp_cost = sqrt(sum((data[,-train_set]-(h_train %*% res_nb$E[ord,-train_set]))^2))
      } else if (cost_f=="IS"){
        a_div = as.vector(data[,-train_set])
        b_div = as.vector(h_train %*% res_nb$E[ord,-train_set])

        zeros <- which(a_div>0)
        tmp_costIS <- b_div
        tmp_costIS[zeros] <- (a_div/b_div - log(a_div/b_div))[zeros]
        tmp_cost <- sum(tmp_costIS)
      } else {
        warning(paste0("The tmp_cost function ", cost_f, "is not implemented. Generalized Kullback-Leibler will be used."))
        tmp_cost <- gklDiv(as.vector(data[,-train_set]), as.vector(h_train %*% res_nb$E[ord,-train_set]))
      }
      tmp_cost
    }
    stopCluster(cl)
    cv_results <- list()
    cv_results$cost <- cost
    cv_results$cost_k <- median(cost)
    return(cv_results)
  } else{
    res_nb <- NMFPoisEM(M=data, N=k)

    cores=detectCores()
    cl <- makeCluster(cores[1]-1)
    registerDoParallel(cl)

    cost <- foreach(i=1:n_iterations, .combine=c, .packages=c('SQUAREM', 'SigModeling'), .export = ls(globalenv())) %dopar% {
      n <- ncol(data)
      train_set = sample(1:n, size_train*n)
      res_train <- NMFPoisEM(M=data[,train_set], N=k)
      h_train <- res_train$P

      ord <- corrSignatures(h_train,res_nb$P)

      if (cost_f=="GKL"){
        tmp_cost <- gklDiv(as.vector(data[,-train_set]), as.vector(h_train %*% res_nb$E[ord,-train_set]))
      } else if (cost_f=="Frobenius"){
        tmp_cost = sqrt(sum((data[,-train_set]-(h_train %*% res_nb$E[ord,-train_set]))^2))
      } else if (cost_f=="IS"){
        a_div = as.vector(data[,-train_set])
        b_div = as.vector(h_train %*% res_nb$E[ord,-train_set])

        zeros <- which(a_div>0)
        tmp_costIS <- b_div
        tmp_costIS[zeros] <- (a_div/b_div - log(a_div/b_div))[zeros]
        tmp_cost <- sum(tmp_costIS)
      } else {
        warning(paste0("The tmp_cost function ", cost_f, "is not implemented. Generalized Kullback-Leibler will be used."))
        tmp_cost <- gklDiv(as.vector(data[,-train_set]), as.vector(h_train %*% res_nb$E[ord,-train_set]))
      }
      tmp_cost
    }
    stopCluster(cl)
    cv_results <- list()
    cv_results$cost <- cost
    cv_results$cost_k <- median(cost)
    return(cv_results)
  }

}
