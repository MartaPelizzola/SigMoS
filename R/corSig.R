#' @importFrom stats cor
#'
#' @title Correlation of mutational signatures.
#'
#'
#' @param H1 Numeric matrix of mutational signatures. Each column should represent a signature.
#' @param H2 Numeric matrix of mutational signatures. Each column should represent a signature.
#'
#'
#' @return Vector of indexes to reorder the second matrix to match the first one.
#' Also returns the full correlation matrix between the signatures.
#'
#' @examples
#' # Find cosine similarity between signatures estimated using SigMoS
#' # with the Negative Binomial distribution and the Poisson distribution:
#' res_nb <- sigmos(BRCA21,k=3,patient_specific = TRUE)
#' res_pois <- sigmos(BRCA21,k=3,method="Poisson")
#' corSig(res_pois$Signatures, res_nb$Signatures)
#'
#' @export
corSig <- function(H1, H2){
  if (!all.equal(dim(H1),dim(H2))){
    stop("The two signature matrices need to have the same dimensions")
  }
  dist <- cor(t(H1),t(H2))
  K <- nrow(H1)
  m <- numeric(K)
  distmat = dist
  for(s in 1:K){
    max.dist <- max(dist)
    remove = which(dist == max.dist, arr.ind = TRUE)
    dist[remove[1,1],] <- -2
    dist[,remove[1,2]] <- -2
    m[remove[1,1]] <- remove[1,2]

  }
  Output <- list()
  Output$distmat <- distmat      # the distance matrix of correlations
  Output$match <- as.integer(m)  # the best matched signatures
  return(Output)
}

