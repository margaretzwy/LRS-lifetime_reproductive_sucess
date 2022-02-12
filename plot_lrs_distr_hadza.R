#  Plot the lifetime reproductive success of Hadza

plot_lrs_distr_hadza <- function(distr, theoretic_distri, n= 15){
  lrs.mean <- round(distr_mean(distr),2)
  lrs.sd <- round(sqrt(distr_var(distr)), 2)
  plot(distr$mids[1:n], distr$density[1:n], xlab = "Number of offspring (LRS)", ylab = "Probablity",
       pch = 20, lwd = 3)
  polygon(x = c(max(0, lrs.mean - lrs.sd), lrs.mean + lrs.sd, lrs.mean + lrs.sd, max(0, lrs.mean - lrs.sd), max(0, lrs.mean - lrs.sd)),
          y = c(0, 0, max(distr$density)*0.05,max(distr$density)*0.05,0),
          density = -1, col = "lightgray", border = "lightgray")
  arrows(lrs.mean, max(distr$density)*0.075, x1 = lrs.mean, y1 = 0, length = 0.05, col = "darkgray", lwd = 1.5)
  lines(theoretic_distri$mids, theoretic_distri$density, lty = 3)
  text(lrs.mean, max(distr$density)*0.09, lrs.mean, col = "darkgray", cex = 0.85)
  text(lrs.mean + lrs.sd, max(distr$density)*0.07, lrs.mean + lrs.sd, col = "lightgray", cex = 0.75)
  text(max(0, lrs.mean - lrs.sd), max(distr$density)*0.07, max(0, lrs.mean - lrs.sd), col = "lightgray", cex = 0.75)
  points(distr$mids[1:n], distr$density[1:n], pch = 20)
}
