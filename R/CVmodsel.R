#' @import foreach
#' @import doParallel
#' @import parallel
#'
#'
#' @title Cross-validation and model selection pipeline for learning of mutational signatures
#'
#' @description Analysis pipeline with model selection and cross-validation for extracting mutational signatures from
#' mutational counts data using non-negative matrix factorization. A Poisson model is tested first and, if
#' overdispersion is observed, then a Negative Binomial model is tested.
#'
#'
#' @param data Numeric matrix of mutational counts data. Matrix size: no. of mutation types x no. of patients.
#' @param k Vector of integers for the number of signatures values to be tested in cross-validation.
#' @param n_iterations Number of iterations for cross-validation procedure.
#' @param cost_f Cost function to be used for the cross-validation error. Default is GKL (Generalized Kullback-Leibler). Alternatives are: Frobenius norm ("Frobenius") and Itakura-Saito divergence ("IS").
#' @param size_train Size of the train set. Default is 0.9, i.e. 90% of the patients are used for training and 10% for testing in the cross-validation.
#' @param patient_specific Logical. If TRUE patient-specific overdispersion is used in the Negative Binomial model.
#' @param tol Threshold for convergence of the EM and MM algorithms for estimating mutational signatures. Default is 1e-5.
#'
#'
#' @return List of non-negative matrix factorization (with number of signatures chosen with the sigmos function) and residuals results for Poisson model if only the Poisson model is used.
#' If the Negative Binomial model is used too, the function returns a list of length two with one element per model (Poisson and NegBin). Each of them is a list of two elements for the results of the matrix factorization (with number of signatures chosen with the sigmos) and the residuals.
#'
#' @export

CVmodsel <- function(data,k=2:8,n_iterations=100,cost_f="GKL",size_train=0.9,patient_specific=FALSE,tol = 1e-5, seed = sample(1:1000,3)){
  if(is.null(data)){
    stop("The data set of the mutational counts is missing.")
  }

  if(!all.equal(k,round(k))){
    stop("The values for the number of signatures to be tested must be integers.")
  }

  #### check Poisson
  cost <- c()
  for (i in k){
    tmp_cost <- sigmos(data,i,n_iterations,method = "Poisson", cost_f,size_train,patient_specific,tol)
    cost <- c(cost, tmp_cost$cost_k)
  }
  k_poisson <- k[which.min(cost)]
  res <- NMFPois(data, k_poisson, tol = tol, seed = seed)

  W <- res$E
  H <- res$P
  #### check residuals
  rsd = (data - H%*%W)
  pp <- length(which(rsd<2*sqrt(data) & rsd>-2*sqrt(data)))/length(rsd)

  if(pp < 0.95){
    #### Poisson fit is not good. Ask if NB check should be done
    x <- menu(c("Yes", "No"), title=paste0("The percentage of residual points within the expected variance lines is ", round(pp,2), " showing overdispersion. Do you want to apply the negative binomial model?"))
    if (x==1){
      #### Build the result variable for the negative binomial case
      result <- list()
      result$Poisson <- list()
      result$Poisson$res_nmf <- res
      result$Poisson$residuals <- rsd
      result$NegBin <- list()

      #### check Negative Binomial
      cost_nb <- c()
      for (i in k){
        tmp_cost <- sigmos(data,i,n_iterations,method = "NB", cost_f,size_train,patient_specific,tol)
        cost_nb <- c(cost_nb,tmp_cost$cost_k)
      }

      k_nb <- k[which.min(cost_nb)]
      alpha <- alphaNR(data,k_nb,patient_specific)
      res_nb <- NMFNB(data, k_nb,alpha = alpha, tol = tol, seed = seed)

      W_nb <- res_nb$E
      H_nb <- res_nb$P
      #### check residuals
      rsd_nb = (data - H_nb%*%W_nb)
      pp_nb <- length(which(rsd_nb<2*sqrt(data*(1+data/alpha)) & rsd_nb>-2*sqrt(data*(1+data/alpha))))/length(rsd_nb)


      if (pp_nb>=0.95){
        print(paste0("The percentage of residual points within the expected variance lines is ", round(pp_nb,2), " showing that the Negative Binomial model is appropriate for this data set."))
      } else {
        if(patient_specific==FALSE){
          print(paste0("The percentage of residual points within the expected variance lines is ", round(pp_nb,2), " showing overdispersion under the Negative Binomial model. Using patient specific overdispersion may improve the results"))
        } else {
          print(paste0("The percentage of residual points within the expected variance lines is ", round(pp_nb,2), " showing overdispersion under the Negative Binomial model."))
        }
      }

      #### Return result Poisson and Negative Binomial
      result$NegBin$res_nmf <- res_nb
      result$NegBin$residuals <- rsd_nb
    }
    if (x == 2){
      #### Do not check for negative binomial & return Poisson results
      print(paste0("The percentage of residual points within the expected variance lines is ", round(pp,2), " showing overdispersion under the Poisson model."))

      result <- list()
      result$res_nmf <- res
      result$residuals <- rsd
    }
  } else {
    #### Poisson model is good, return Poisson results
    print(paste0("The percentage of residual points within the expected variance lines is ", round(pp,2), " showing that the Poisson model is appropriate for this data set."))

    result <- list()
    result$res_nmf <- res
    result$residuals <- rsd
  }
  return(result)
}
