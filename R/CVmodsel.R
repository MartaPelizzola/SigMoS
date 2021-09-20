#' Cross-validation and model selection pipeline for learning of mutational signatures
#'
#' Analysis pipeline with model selection and cross-validation for extracting mutational signatures from 
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
#' @return List of non-negative matrix factorization and residuals results for Poisson model if only the Poisson model is used. 
#' If the Negative Binomial model is used too returns a list of length two with one element per model (Poisson and NegBin). Each of them is a list of two elements for the results of the matrix factorization and the residuals. 
#'
#' @export

CVmodsel <- function(data,k=2:3,n_iterations=100,cost_f="GKL",size_train=0.9,patient_specific=FALSE,tol = 1e-5){
  #### check Poisson 
  cost <- foreach(i=k, .combine=c, .packages='SQUAREM', .export = ls(globalenv())) %dopar% {
    tmp_cost <- CrossValidation(data,i,n_iterations,method = "Poisson", cost_f,size_train,patient_specific,tol)
    tmp_cost$cost_k
  }
  k_poisson <- k[which.min(cost)]
  res <- NMFPoisEM(data, k_poisson)
  
  W <- res$E
  H <- res$P
  #### check residuals 
  rsd = (data - H%*%W)
  pp <- length(which(rsd<2*sqrt(data) & rsd>-2*sqrt(data)))/length(rsd)
  
  result <- list()
  result$Poisson <- list()
  result$Poisson$res_nmf <- res
  result$Poisson$residuals <- rsd
  
  if(pp < 0.95){
    #### Poisson fit is not good. Ask if NB check should be done
    x <- menu(c("Yes", "No"), title=paste0("The percentage of residual points within the expected variance lines is ", pp, " showing overdispersion. Do you want to apply the negative binomial model?"))
    if (x==1){
      #### check Negative Binomial
      cost_nb <- foreach(i=k, .combine=c, .packages='SQUAREM', .export = ls(globalenv())) %dopar% {
        tmp_cost <- CrossValidation(data,i,n_iterations,method = "NB", cost_f,size_train,patient_specific,tol)
        tmp_cost
      }
      k_nb <- k[which.min(cost_nb)]
      res_nb <- NMFNBMM(data, k_nb)
      
      W_nb <- res_nb$E
      H_nb <- res_nb$P
      #### check residuals
      rsd_nb = (data - H_nb%*%W_nb)
      pp_nb <- length(which(rsd_nb<2*sqrt(data) & rsd_nb>-2*sqrt(data)))/length(rsd_nb)
      
      result$NegBin <- list()
      result$NegBin$res_nmf <- res_nb
      result$NegBin$residuals <- rsd_nb
      
      if (pp>=0.95){
        print(paste0("The percentage of residual points within the expected variance lines is ", pp, " showing that the Negative Binomial model is appropriate for this data set."))
      } else {
        print(paste0("The percentage of residual points within the expected variance lines is ", pp, " showing overdispersion."))
      }
      #### Return result Poisson and Negative Binomial
      return(result)
    } else {
      #### Return result Poisson only (when no Negative Binomial analysis)
      return(result)
    }
  } else {
    print(paste0("The percentage of residual points within the expected variance lines is ", pp, " showing that the Poisson model is appropriate for this data set."))
    #### Return result Poisson only (when the Poisson fit is good)
    return(result)
  }
}