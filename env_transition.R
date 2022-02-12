# obtain the environment transition matrix, D, according to equilibrium environmental ratio (evn1/evn2)
# and, rho,the autocorrelation of evn1 and evn2
env_transition <- function(ratio_pi, rho){
  b <- 1 - ratio_pi*(1-rho)
  a <- rho + 1 - b

  D <- matrix(c(c(a, 1 - b),  # r value will impact the frequent changes between evn1 and evn2
                c(1 - a, b)),
              byrow = T, nrow = 2)
  # make sure D is a transition matrix, but of course pi and rho need to be calculated
  D[which(D > 1)] <- 1
  D[which(D < 0)] <- 0
  return(D)
}
