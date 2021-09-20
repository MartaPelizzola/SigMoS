#' Negative Binomial non-negative matrix factorization
#'
#' MM algorithm for negative binomial non-negative matrix factorization
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
NMFNBMM = function(M, N, alpha, tol = 1e-5, seed = sample(1:1000,1)){
  if (N!=round(N)){
    stop("The number of signatures must be an integer.")
  }
  if(is.null(M)){
    stop("The data set of the mutational counts is missing.")
  }
  if(is.null(N)){
    stop("A value for the number of signatures to be estimated is missing.")
  }

  K <- dim(M)[1]
  G <- dim(M)[2]

  if(!(length(alpha)==1 | length(alpha)==G)){
    stop("The overdispersion parameter 'alpha' should have length 1 or length equal to the number of patients.")
  }

  div <- rep(0,length(seed)) # vector of different GKLD values
  Plist <- list()            # list of P matrices
  Elist <- list()            # list of E matrices
  reslist <- list()

  NB.em = function(x){
    x = exp(x)
    P = matrix(x[1:(K*N)], nrow = K, ncol = N)
    E = matrix(x[-c(1:(K*N))], nrow = N, ncol = G)

    P <- P * ( ( ( M / (P%*%E) ) %*% t(E) ) / ( ( (alpha + M) / (alpha + P%*%E) ) %*% t(E) ) )     # update of signatures

    E <- E * ( (t(P) %*% ( M / (P%*%E) ) ) / (t(P) %*% ( (alpha + M) / (alpha + P%*%E) ) ) )     # update of exposures

    par = c(as.vector(P),as.vector(E))
    par[par <= 0] = 1e-10
    return(log(par))
  }

  NBobj = function(x){
    x = exp(x)
    P = matrix(x[1:(K*N)], nrow = K, ncol = N)
    E = matrix(x[-c(1:(K*N))], nrow = N, ncol = G)

    y = as.vector(M)
    mu = as.vector(P%*%E)
    r = mu
    p = which(y > 0)
    r[p] = (y * (log(y)- log(mu)) - (alpha + y) * log((alpha + y)/(alpha + mu)))[p]

    obj = sum(r) # euclidean distance

    return(obj)
  }

  for(i in 1:length(seed)){
    set.seed(seed[i])

    P <- matrix(runif(K*N), nrow = K, ncol = N)  # Initialize P
    E <- matrix(runif(N*G), nrow = N, ncol = G)  # Initialize E

    init = log(c(as.vector(P),as.vector(E)))
    sres = squarem(init, fixptfn = NB.em, objfn = NBobj, control = list(tol = tol))

    P = matrix(exp(sres$par[1:(K*N)]), nrow = K, ncol = N)
    E = matrix(exp(sres$par[-c(1:(K*N))]), nrow = N, ncol = G)
    E = diag(colSums(P)) %*% E # normalizing
    P = P %*% diag(1/colSums(P))

    Plist[[i]] <- P # signatures
    Elist[[i]] <- E # exposures
    div[i] <- NBobj(sres$par) # final generalized Kullback-Leibler divergence
    reslist[[i]] = sres
  }

  best <- which.min(div) # Smallest GKLD value
  P = Plist[[best]]
  E = Elist[[best]]


  Output <- list()
  Output$P <-  P
  Output$E <-  E
  Output$obj <- div[best]
  Output$results <- reslist

  return(Output)
}
