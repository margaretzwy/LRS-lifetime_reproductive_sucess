# Calculate the lifetime reproductive success distribution for age only cases with R1 and P matrix provided,
# as in examples of van Daalen and Caswell daalen (2013)s.

age_only_r1p <- function(population.name){
  ### for cases provide R1 and P matrice, such as Hutt, Hadza and Ache
  ##==== Get the LRS distribution
  R1 <- as.matrix(read.csv(paste(population.name, "_R1.csv", sep = ""), header = F)) #fertility
  P <-  as.matrix(read.csv(paste(population.name, "_P.csv", sep = ""), header = F))
  Umatrix <- P[-dim(P)[1], -dim(P)[2]]
  Fmatrix <- R1[-dim(R1)[1], -dim(R1)[2]]
  gamma <- age_only_repro_distr(Fmatrix, Umatrix,
                                max.kids = 1,
                                distri_type = "bernoulli")# because they are human populations we assume each birth only one kid
  gamma.distr <- data.frame(density = gamma, mids = (seq_along(gamma) - 1)) # Write distribution and its corresponse kids number into a data frame
  ##==== Get the theoretic distribution from first 3 moments
  caswell_result <- daalen_caswell_equations(Fmatrix, Umatrix, "bernoulli") # Calculate dVC results based on their paper
  raw_moments <- as.numeric(caswell_result[4:6]) # Get the moments calculated by dVC formular
  theoretic_distri <- distri_from_moments(raw_moments) # Get a best fitting distribution based on the first 3 moments
  ##====plot the LRS distribution and the theoretic distribution
  plot_lrs_distr_allpop(population.name, gamma.distr, theoretic_distri)
  mtext(population.name)
  ##====Gini index
  gini_plot(population.name, gamma.distr) # plot Gini Index figure
  ##====Return the LRS distribution
  return(gamma.distr)
}
