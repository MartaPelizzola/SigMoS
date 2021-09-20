gklDiv <- function(y, mu){
  r <- mu
  p <- which(y > 0)
  r[p] <- (y * (log(y)- log(mu)) - y + mu)[p]
  return(sum(r))
}