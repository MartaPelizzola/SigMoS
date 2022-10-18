#' @import SQUAREM


#' @title Negative Binomial non-negative matrix factorization
#'
#' @description MM algorithm for negative binomial non-negative matrix factorization
#'
#'
#' @param M Numeric matrix of mutational counts data. Matrix size: no. of mutation types x no. of patients.
#' @param N Integer. Number of signatures to be used for matrix factorization.
#' @param alpha Overdispersion parameter. Either vector of length one or vector of length no. of patients if patient-specific overdispersion is used.
#' @param tol Threshold for convergence of the MM algorithm for estimating mutational signatures. Default is 1e-5.
#' @param seed Vector of seeds for initializing the values of the signatures and exposures matrices. If default only one initialization is performed.
#'
#'
#' @return List of non-negative matrix factorization results. P is the signature matrix and E is the exposures matrix
#'
#' @export
#'
#'
NMFNB = function(M, N=NULL, alpha, tol = 1e-3, seed = sample(1:1000,1)){
  if (N!=round(N)){
    stop("The number of signatures must be an integer.")
  }
  if(is.null(M)){
    stop("The data set of the mutational counts is missing.")
  }
  M <- t(M)
  if(is.null(N)){
    stop("A value for the number of signatures to be estimated is missing.")
  }
  if(length(N)!=1){
    stop("More than one value for the number of signatures is used as input. NB-NMF can only be performed for one value of 'N' at a time.")
  }

  K <- dim(M)[1]
  G <- dim(M)[2]

  if(!(length(alpha)==1 || length(alpha)==G)){
    stop(paste0("The length of alpha is ", length(alpha),". The overdispersion parameter 'alpha' should have length 1 or length equal to the number of patients: ",G))
  }

  div <- rep(0,length(seed)) # vector of different GKLD values
  Plist <- list()            # list of P matrices
  Elist <- list()            # list of E matrices
  reslist <- list()

  alphamat = matrix(alpha, nrow = K, ncol = G, byrow = T)

  NB.em = function(x){
    x = exp(x)
    P = matrix(x[1:(K*N)], nrow = K, ncol = N)
    E = matrix(x[-c(1:(K*N))], nrow = N, ncol = G)

    PE = P%*%E
    P <- P * ( ( ( M / (PE) ) %*% t(E) ) / ( ( (alphamat + M) / (alphamat + PE) ) %*% t(E) ) )     # update of signatures
    PE = P%*%E
    E <- E * ( (t(P) %*% ( M / (PE) ) ) / (t(P) %*% ( (alphamat + M) / (alphamat + PE) ) ) )     # update of exposures

    par = c(as.vector(P),as.vector(E))
    par[par <= 0] = 1e-50
    return(log(par))
  }

  # divergence
  # NBobj = function(x){
  #   x = exp(x)
  #   P = matrix(x[1:(K*N)], nrow = K, ncol = N)
  #   E = matrix(x[-c(1:(K*N))], nrow = N, ncol = G)
  #
  #   PE = as.vector(t(P%*%E))
  #   M = as.vector(t(M))
  #   r <- -(M+alpha)*(log(alpha + M) - log(alpha + PE))
  #   p <- which(M > 0)
  #   r[p] <- (M * (log(M)- log(PE)) - (M+alpha)*(log(alpha + M) - log(alpha + PE)))[p]
  #
  #   obj = sum(r) # euclidean distance
  #
  #   return(obj)
  # }

  # likelihood function
  NBlik = function(x){
    x = exp(x)
    P = matrix(x[1:(K*N)], nrow = K, ncol = N)
    E = matrix(x[-c(1:(K*N))], nrow = N, ncol = G)

    y = as.vector(t(M))
    mu = as.vector(t(P%*%E))

    prob = 1-(mu/(alpha + mu))

    r = dnbinom(y, alpha, ifelse(prob != 0,prob, 1e-16), log = T)

    return(sum(r))
  }

  for(i in 1:length(seed)){
    set.seed(seed[i])

    P <- matrix(runif(K*N), nrow = K, ncol = N)  # Initialize P
    E <- matrix(runif(N*G), nrow = N, ncol = G)  # Initialize E

    init = log(c(as.vector(P),as.vector(E)))
    sres = squarem(init, fixptfn = NB.em, objfn = NBlik, control = list(tol = tol, minimize = F))

    P = matrix(exp(sres$par[1:(K*N)]), nrow = K, ncol = N)
    E = matrix(exp(sres$par[-c(1:(K*N))]), nrow = N, ncol = G)
    E = diag(colSums(P)) %*% E # normalizing
    P = P %*% diag(1/colSums(P))

    Plist[[i]] <- P # signatures
    Elist[[i]] <- E # exposures
    div[i] <- NBlik(sres$par) # final likelihood value
    reslist[[i]] = sres
  }

  best <- which.max(div) # Largest Likelihood value
  P = Plist[[best]]
  E = Elist[[best]]


  Output <- list()
  Output$P <-  t(P)
  Output$E <-  t(E)
  Output$div <- div[best]
  Output$results <- reslist

  return(Output)
}
