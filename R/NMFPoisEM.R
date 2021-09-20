#' Poisson non-negative matrix factorization
#'
#' EM algorithm for poisson non-negative matrix factorization
#'
#'
#' @param M Numeric matrix of mutational counts data. Matrix size: no. of mutation types x no. of patients.
#' @param N Integer. Number of signatures to be used for matrix factorization.
#' @param tol Threshold for convergence of the EM algorithm for estimating mutational signatures. Default is 1e-5.
#' @param seed Vector of seeds for initializing the values of the signatures and exposures matrices. If default only one initialization is performed.
#'
#'
#' @return List of non-negative matrix factorization results. P is the signature matrix and E is the exposures matrix
#'
#' @export
#'
#'

NMFPoisEM = function(M,N=NULL, tol = 1e-5, seed = sample(1:1000,1)){
  if (N!=round(N)){
    stop("The number of signatures must be an integer.")
  }
  if(is.null(M)){
    stop("The data set of the mutational counts is missing.")
  }
  if(is.null(N)){
    stop("A value for the number of signatures to be estimated is missing.")
  }
  if(length(N)!=1){
    stop("More than one value for the number of signatures is used as input. NB-NMF can only be performed for one value of 'N' at a time.")
  }

  K <- dim(M)[1]  # mutations
  G <- dim(M)[2]  # patients

  div <- rep(0,length(seed)) # vector of different GKLD values
  Plist <- list()            # list of P matrices
  Elist <- list()            # list of E matrices
  reslist <- list()

  poisson.em = function(x){
    x = exp(x)
    P = matrix(x[1:(K*N)], nrow = K, ncol = N)
    E = matrix(x[-c(1:(K*N))], nrow = N, ncol = G)

    PE <- P%*%E
    P <- P * ((M/PE) %*% t(E))      # update of signatures
    P <- P %*% diag(1/colSums(P))   # make sure the columns sum to one

    PE <- P%*%E
    E <- E * (t(P) %*% (M/PE))      # update of exposures

    par = c(as.vector(P),as.vector(E))
    par[par <= 0] = 1e-10
    return(log(par))
  }

  gklobj = function(x){
    x = exp(x)
    P = matrix(x[1:(K*N)], nrow = K, ncol = N)
    E = matrix(x[-c(1:(K*N))], nrow = N, ncol = G)

    GKL <- gklDiv(as.vector(M),as.vector(P%*%E)) # GKLD value

    return(GKL)
  }

  for(i in 1:length(seed)){
    set.seed(seed[i])

    P <- matrix(runif(K*N), nrow = K, ncol = N)  # Initialize P
    E <- matrix(runif(N*G), nrow = N, ncol = G)  # Initialize E

    init = log(c(as.vector(P),as.vector(E)))
    sres = squarem(init, fixptfn = poisson.em, objfn = gklobj, control = list(tol = tol))

    P = matrix(exp(sres$par[1:(K*N)]), nrow = K, ncol = N)
    E = matrix(exp(sres$par[-c(1:(K*N))]), nrow = N, ncol = G)

    Plist[[i]] <- P # signatures
    Elist[[i]] <- E # exposures
    div[i] <- gklobj(sres$par)   # final generalized Kullback-Leibler divergence
    reslist[[i]] = sres
  }

  best <- which.min(div) # Smallest GKLD value
  P = Plist[[best]]
  E = Elist[[best]]

  Output <- list()
  Output$P <-  P
  Output$E <-  E
  Output$gkl <- div[best]
  Output$results <- reslist

  return(Output)
}
