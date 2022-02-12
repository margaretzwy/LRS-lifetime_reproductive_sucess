# Calculate the lifetime reproductive success distribution for stage only cases with F and U matrix provided,
# as in examples of van Daalen and Caswell daalen (2013)s.

stage_only_fu <- function(population.name,
                          max.offspring = 2^14,
                          initial.stage = 1,
                          distri_type = "poisson"){
  ### For cases has F and U matrice file provided
  ### max.offspring <- 2^14 is for "Tusga" only, other population need to reconsider the maximum offspring for poisson distribution
  Fmatrix <- as.matrix(read.csv(paste(population.name, "_F.csv", sep = ""), header = F)) #fertility
  Umatrix <-  as.matrix(read.csv(paste(population.name, "_U.csv", sep = ""), header = F))

  gamma <- lrs_fft_stage_only(Fmatrix, Umatrix,
                              max.offspring,
                              initial.stage,
                              distri_type = distri_type)
  gamma.distr <- data.frame(density = gamma,
                            mids = seq_along(gamma) - 1)
  ##====plot the LRS distribution and the theoretic distribution
  plot_lrs_distr_allpop(population.name, gamma.distr)
  mtext(population.name)
  ##====Gini index
  gini_plot(population.name, gamma.distr)
  ##====Return the LRS distribution
  return(gamma.distr)
}
