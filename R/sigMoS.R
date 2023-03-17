#' @import foreach
#' @import doParallel
#' @import parallel
#'
#' @title Model selection algorithm for the number of mutational signatures
#'
#' @description Model selection algorithm for a given number of signatures k using either Poisson or Negative Binomial model.
#'
#'
#' @param data Numeric matrix of mutational counts data. Matrix size: no. of patients x no. of mutation types.
#' @param k Number of signatures to be used in the cross-validation.
#' @param n_iterations Number of iterations for cross-validation procedure. Default is 100.
#' @param method Non-negative matrix factorization model: 'Poisson' or 'NB' can be used. Default is 'NB' corresponding to the Negative Binomial model.
#' @param cost_f Cost function to be used for the cross-validation error. Default is GKL (Generalized Kullback-Leibler). Alternatives are: Frobenius norm ("Frobenius") and Itakura-Saito divergence ("IS").
#' @param size_train Size of the train set. Default is 0.9, i.e. 90% of the patients are used for training and 10% for testing in the cross-validation.
#' @param patient_specific Logical. If TRUE patient-specific overdispersion is used in the Negative Binomial model.
#' @param tol Threshold for convergence of the EM and MM algorithms for estimating mutational signatures. Default is 1e-3.
#'
#'
#' @return List of length two including the costs.
#'  - **cost** is a vector of length n_iterations where the cost for each iteration is stored.
#'  - **cost_k** is the median of the values in 'cost' and it is used for comparison to choose the best number of signatures for a given data set.
#'  - **FullDataCost** is the likelihood value for the fit of the full data
#'  - **Signatures** The signature matrix for the full data
#'  - **Exposures** The exposure matrix for the full data
#'
#'
#' @examples
#' # Use SigMoS with the Negative Binomial distribution:
#' res <- sigmos(BRCA21,k=3,patient_specific = TRUE)
#' str(res$Signatures)
#' str(res$Exposures)
#' # Use SigMoS with the Poisson distribution:
#' res <- sigmos(BRCA21,k=3,method="Poisson")
#' str(res$Signatures)
#' str(res$Exposures)
#' #### Not Run
#' ## Find estimated number of signatures for the given example using the Negative
#' ## Binomial distribution. Evaluate no. of signatures between 2 and 7:
#' # res_cv <- list()
#' # CVcost <- rep(0,6)
#' # for (i in 2:7){
#' #  res_cv[[i]] <- sigmos(BRCA21,k=i,method="NB")
#' #  CVcost[i-1] = res_cv[[i]]$cost_k
#' # }
#' # which.min(CVcost)+1 #estimated number of signatures
#'
#' @export
sigmos <- function(data,k=NULL,n_iterations=100,method = "NB", cost_f="GKL",size_train=0.9,patient_specific=FALSE,tol = 1e-5, beta = 2){
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
    alpha <- alphaNR(data,k,patient_specific)
    res_nb <- NMFNB(M=data, N=k, alpha = alpha, tol = tol, seed = sample(1:1000,5))

    cores=detectCores()
    cl <- makeCluster(cores[1]-1)
    registerDoParallel(cl)

    cost <- foreach(i=1:n_iterations, .combine=c, .packages=c('SQUAREM'), .export = ls(globalenv())) %dopar% {
      set.seed(i)
      n <- nrow(data)
      train_set = sample(1:n, size_train*n)
      if(patient_specific){
        res_train <- NMFNB(M=data[train_set,], N=k, alpha = alpha[train_set], tol = tol)
      }else{
        res_train <- NMFNB(M=data[train_set,], N=k, alpha = alpha, tol = tol)
      }

      h_train <- res_train$P

      ord = corSig(h_train,res_nb$P)$match

      if (cost_f=="GKL"){
        tmp_cost <- gklDiv(as.vector(data[-train_set,]), as.vector(res_nb$E[-train_set,ord] %*% h_train))
      } else if (cost_f=="Frobenius"){
        tmp_cost = sqrt(sum((data[-train_set,]-(res_nb$E[-train_set,ord] %*% h_train))^2))
      } else if (cost_f=="IS"){
        a_div = as.vector(data[-train_set,])
        b_div = as.vector(res_nb$E[-train_set,ord] %*% h_train)

        zeros <- which(a_div>0)
        tmp_costIS <- b_div
        tmp_costIS[zeros] <- (a_div/b_div - log(a_div/b_div))[zeros]
        tmp_cost <- sum(tmp_costIS)
      } else if (cost_f=="Beta"){
        a_div = as.vector(data[-train_set,])
        b_div = as.vector(res_nb$E[-train_set,ord] %*% h_train)

        zeros <- which(a_div>0)
        tmp_costBETA <- b_div
        tmp_costBETA[zeros] <- (1/(beta*(beta+1)))*(a_div^(beta+1)-b_div^(beta+1) - (beta+1)*(b_div^beta)*(a_div-b_div))[zeros]
        tmp_cost <- sum(tmp_costBETA)
      } else {
        warning(paste0("The tmp_cost function ", cost_f, "is not implemented. Generalized Kullback-Leibler will be used."))
        tmp_cost <- gklDiv(as.vector(data[-train_set,]), as.vector(res_nb$E[-train_set,ord] %*% h_train))
      }
      tmp_cost
    }
    stopCluster(cl)
    cv_results <- list()
    cv_results$cost <- cost
    cv_results$cost_k <- median(cost)
    cv_results$FullDataCost = res_nb$div
    cv_results$Signatures = res_nb$P
    cv_results$Exposures = res_nb$E
    return(cv_results)
  } else{
    res_nb <- NMFPois(M=data, N=k, tol = tol, seed = sample(1:1000,5))

    cores=detectCores()
    cl <- makeCluster(cores[1]-1)
    registerDoParallel(cl)

    cost <- foreach(i=1:n_iterations, .combine=c, .packages=c('SQUAREM'), .export = ls(globalenv())) %dopar% {
      set.seed(i)
      n <- nrow(data)
      train_set = sample(1:n, size_train*n)
      res_train <- NMFPois(M=data[train_set,], N=k, tol = tol)
      h_train <- res_train$P

      ord <- corSig(h_train,res_nb$P)$match

      if (cost_f=="GKL"){
        tmp_cost <- gklDiv(as.vector(data[-train_set,]), as.vector(res_nb$E[-train_set,ord] %*% h_train))
      } else if (cost_f=="Frobenius"){
        tmp_cost = sqrt(sum((data[-train_set,]-(res_nb$E[-train_set,ord] %*% h_train))^2))
      } else if (cost_f=="IS"){
        a_div = as.vector(data[-train_set,])
        b_div = as.vector(res_nb$E[-train_set,ord] %*% h_train)

        zeros <- which(a_div>0)
        tmp_costIS <- b_div
        tmp_costIS[zeros] <- (a_div/b_div - log(a_div/b_div))[zeros]
        tmp_cost <- sum(tmp_costIS)
      } else if (cost_f=="Beta"){
        a_div = as.vector(data[-train_set,])
        b_div = as.vector(res_nb$E[-train_set,ord] %*% h_train)

        zeros <- which(a_div>0)
        tmp_costBETA <- b_div
        tmp_costBETA[zeros] <- (1/(beta*(beta+1)))*(a_div^(beta+1)-b_div^(beta+1) - (beta+1)*(b_div^beta)*(a_div-b_div))[zeros]
        tmp_cost <- sum(tmp_costBETA)
      } else {
        warning(paste0("The tmp_cost function ", cost_f, "is not implemented. Generalized Kullback-Leibler will be used."))
        tmp_cost <- gklDiv(as.vector(data[-train_set,]), as.vector(res_nb$E[-train_set,ord] %*% h_train))
      }
      tmp_cost
    }
    stopCluster(cl)
    cv_results <- list()
    cv_results$cost <- cost
    cv_results$cost_k <- median(cost)
    cv_results$FullDataCost = -res_nb$gkl
    cv_results$Signatures = res_nb$P
    cv_results$Exposures = res_nb$E
    return(cv_results)
  }
}
