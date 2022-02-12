# Calculate the mean, variance, standard deviation, and first 3 moments of lifetime reproductive success
# based on van Daalen and Caswell daalen (2013, Henceforth it will be referred as vDC).
# All the parameters in this function use the same symbols in the vDC paper.


daalen_caswell_equations <- function(Fmatrix,
                                     Umatrix,
                                     fertility.distr = c("bernoulli", "poisson"),
                                     initial.stage = 1,
                                     last.stage = initial.stage){
  # Use van Daalen and Caswell 2013 6(3) Theoretical Ecology formular to calculat mean of LRO
  # for alpha = 1 situation only
  U <- Umatrix
  tau <- dim(U)[1]
  qvec <- 1 - colSums(U) # the probability of dying
  P <- rbind(cbind(U, rep(0, tau)),
             c(qvec, 1))
  R1 <- matrix(rep(c(Fmatrix[1,],0), tau + 1), nrow = tau + 1, byrow = T)
  Z <- cbind(diag(tau), rep(0, tau))
  N <- solve(diag(tau) - U)
  Rhat1 <- Z %*% R1 %*% t(Z)

  if(fertility.distr == "bernoulli"){
    R2 <- R3 <- R4 <- R1
  }else if(fertility.distr == "poisson"){
    R2 = R1 + R1 * R1
    R3 = R1 + 3 * R1 * R1 + R1 * R1 * R1
    #R4 <- R1 + 3 * R1 * R1 + R1 * R1 * R1 # need to figure out from tulja'a note
  }
  Rhat2 <- Z %*% R2 %*% t(Z)
  Rhat3 <- Z %*% R3 %*% t(Z)
  rho1 <- t(N) %*% Z %*% t(P * R1) %*% rep(1, tau + 1)
  rho2 <- t(N) %*% (Z %*% t(P * R2) %*% rep(1, tau + 1) + 2 * t(U * Rhat1) %*% rho1)
  rho3 <- t(N) %*% (Z %*% t(P * R3) %*% rep(1, tau + 1) + 3 * t(U * Rhat2) %*% rho1 + 3 * t(U * Rhat1) %*% rho2)
  #rho4 <- t(N) %*% (Z %*% t(P * R4) %*% rep(1, tau + 1) + 4 * t(U * Rhat3) %*% rho1 + 4 * t(U * Rhat2) %*% rho2 + 4 * t(U * Rhat1) %*% rho3)
  vdc.var <- rho2 - rho1*rho1
  vdc.mean <- rho1
  return(data.frame(mean = vdc.mean[initial.stage:last.stage], var = vdc.var[initial.stage:last.stage], sd = sqrt(vdc.var[initial.stage:last.stage]),
                    mom1 = rho1[initial.stage:last.stage], mom2 = rho2[initial.stage:last.stage], mom3 = rho3[initial.stage:last.stage]))
}
