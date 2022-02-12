#  Plot the lifetime reproductive success of Tsuga

plot_lrs_distr_tsuga <- function(distr, n= 800){
  lrs.mean <- round(distr_mean(distr),2)
  lrs.sd <- round(sqrt(distr_var(distr)), 2)
  plot(distr$mids[1:n], log10(distr$density[1:n]), xlab = "Number of offspring (LRS)", ylab = expression(log[10](Probablity)),#"log10(Probablity)",
       pch = 20)
  y.range <- c(max(distr$density[1:n]), min(distr$density[1:n][!is.na(distr$density[1:n])]))
  polygon(x = c(max(0, lrs.mean - lrs.sd), lrs.mean + lrs.sd, lrs.mean + lrs.sd, max(0, lrs.mean - lrs.sd), max(0, lrs.mean - lrs.sd)),
          y = c(log10(min(y.range)),log10(min(y.range)),log10(min(y.range))*0.95,log10(min(y.range))*0.95,log10(min(y.range))),
          density = -1, col = "lightgray", border = "lightgray")
  arrows(lrs.mean, log10(min(y.range))*0.925, x1 = lrs.mean, y1 = log10(min(y.range)), length = 0.05, col = "darkgray")
  #lines(theoretic_distri$mids, theoretic_distri$density, lty = 3)
  text(lrs.mean + 15, log10(min(y.range))*0.9, lrs.mean, col = "darkgray", cex = 0.85)
  text(lrs.mean + lrs.sd + 40, log10(min(y.range))*0.935, lrs.mean + lrs.sd, col = "lightgray", cex = 0.75)
  text(-9, log10(min(y.range))*0.935, max(0, lrs.mean - lrs.sd), col = "lightgray", cex = 0.75)
  arrows(distr$mids[1] + 50, log10(distr$density[1]), x1 = distr$mids[1]+ 10, y1 = log10(distr$density[1]), length = 0.05, col = "darkgray")
  text(distr$mids[1] + 180, log10(distr$density[1]), paste("Pr(LRS = 0) = ",round(distr$density[1],2), sep = ""), col = "darkgray", cex = 0.75)
}
