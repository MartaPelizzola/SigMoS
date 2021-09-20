#' Correlation of mutational signatures
#'
#'
#' @param H1 Numeric matrix of mutational signatures.
#' @param H2 Numeric matrix of mutational signatures.
#'
#'
#' @return Vector of indexes to reorder the second matrix to match the first one.
#' @export


corrSignatures <- function(H1, H2){
  if (all.equal(dim(H1),dim(H2))){
    stop("The two signature matrices need to have the same dimensions")
  }
  cormat <- cor(H1,H2)
  nsig    <- ncol(H1)
  ord   = t(sapply(1:nsig, function (x) order(cormat[x,], decreasing = TRUE)[1:nsig]))
  while(length(unique(ord[,1]))<length(ord[,1])){
    dups <- duplicated(ord[,1])
    for (i in which(dups)){
      tmp <- ord[i,1]
      ord[i,1] <- ord[i,2]
      ord[i,2:ncol(ord)] <- c(ord[i,3:ncol(ord)], tmp)
    }
  }
  return(ord[1:nsig])
}
