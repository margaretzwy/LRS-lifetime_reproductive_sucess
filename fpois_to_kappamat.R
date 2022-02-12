# Convert fertility to number of offspring distribution by Poisson distribution, truncated up if used.
# It can be replaced if any other number of offspring distribution is used.

fpois_to_kappamat <- function(Fmatrix, max.offspring = 12){ # for Tsuga case only
  # Since there is only one state for offspring in Tsuga case, kappamat is a matrix here
  # otherwise it should be an array or a list
  n.stages <- dim(Fmatrix)[2]
  kappamat <- matrix(NA, nrow = max.offspring + 1, ncol = n.stages)
  for(i in seq_len(n.stages)){
    kappamat[, i] <- Fmatrix[1,i]^c(0:max.offspring)*exp(-Fmatrix[1,i])/factorial(c(0:max.offspring))
  }
  #put the cut off pobablity back to the last number of offspring to make sure the colSums(kappamat) = 1s
  kappamat[dim(kappamat)[1], ] <- kappamat[dim(kappamat)[1], ] + 1 - colSums(kappamat)
  return(kappamat)
}
