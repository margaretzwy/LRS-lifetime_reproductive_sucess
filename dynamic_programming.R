# Calculate offspring distribution for individuals end up at different stages.
# It is used for the hybrid method for age + stage cases which calculates gamma part1 in hybrid method.



dynamic_programming <- function(kappamat, G, pvec, max.offspring, initial.stage){
  omega <- length(kappamat)
  n.stages <- dim(pvec)[1]
  V <- R <- list() # V is new transition list, R is a new reproduction list. All matrices in the list have an extra stage called dead.
  for(i in 1:omega){
    U <- G[[i]] %*% diag(pvec[, i])
    U <- rbind(U, 1-pvec[, i]) # add a stage of death
    V[[i]] <- cbind(U, c(rep(0, n.stages), 1)) # add the transition from the stage of death to all stages
    R[[i]] <- cbind(kappamat[[i]], c(1, rep(0, (dim(kappamat[[i]])[1] - 1))))# add the reproduction of the stage of death
  }
  ks_join <- matrix(0, nrow = max.offspring, ncol = n.stages+1) #number of offspring and stage joined distribution
  ks_join[initial.stage, 1] <- 1
  for (t in 1:omega) {
    ks_new <- matrix(0, nrow = max.offspring, ncol = n.stages+1)
    for (s in 1:(n.stages + 1)) {
      z.stage <- convolve(ks_join[, s], rev(R[[t]][, s]), type = "o")
      if(length(z.stage) < max.offspring){
        z.stage <- c(z.stage, rep(0, max.offspring - length(z.stage)))
      }else{
        z.stage <- z.stage[1:max.offspring]
      }
      ks_new = ks_new + z.stage %*% t(V[[t]][, s])
    }
    ks_join <- ks_new
  }
  if(round(sum(ks_join), 10) < 1){
    print(round(sum(ks_join), 10))
    stop("The max.offspring is too small so that cut off causes that the offspring proboblities of all stages are not sum to 1. Please increase your max.offspring and run again.")
  }else if(round(sum(ks_join), 10) == 1){
    return(ks_join)
  }else{
    round(sum(ks_join), 10)
    stop("The sum of the offspring proboblities of all stages is larger than 1. Check your data!")
  }
}
