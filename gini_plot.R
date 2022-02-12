# To plot the Lorenz Curve and find the Gini index.


gini_plot <- function(population.name, gamma.distr){
  gini_val <- round(reldist::gini(gamma.distr$mids, gamma.distr$density), 3)
  y <- c(0, cumsum(gamma.distr$mids * gamma.distr$density)/distr_mean(gamma.distr))
  x <- c(0, cumsum(gamma.distr$density))

  plot(c(0, x),c(0, y), xaxs = "i", yaxs = "i",
       type = "l",
       xlab = "Cumulative share of individual",
       ylab = "Cumulative share of offspring reproduced")
  polygon(x = c(x, 0), y = c(y, 0), density = -1, col = "lightgray", border = F)
  abline(0,1, lty = 2)
  lines(x,y)
  legend(0, 0.9, paste(population.name, ": GI = ", gini_val, sep = ""), bty = "n")
}
