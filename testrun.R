testrun <- function(max.offspring){
  #########=======     normal and early spring cases for reo deer   ==========#########
  pi <- 0.5 # longterm frequency
  rho <- 0.5 # environment autocorrelation
  ratio_pi <- pi/(1-pi)
  D <- env_transition(ratio_pi, rho)
  ### Read the data file ####
  matrices <- roedeer_matrices_list()
  Klist <- r_to_kappamat(matrices) # list of n.stage matrix, each of them has 12 columns for 12 ages and 2 rows for 0 or 1 fawns
  ## size transition matrix, it is age dependent
  # Four age classes: yearlings, 2-7 years old, 8-11 years old and >11 years
  Glist <- g_to_gmat(matrices)
  ## size conditional survival, n.stage stages, row is corresponded to stage
  Pmatrix <- s_to_pvec(matrices)
  omega <- length(Klist) # total number of ages
  ## arbitary environment affect on survival pvec, two environment
  #env_effect_p <- diag(c(0.5, 0.75, rep(1, omega - 2))) # only on age unified on stage
  env_effect_p <- diag(c(rep(0.819/0.911, 7), rep(1, omega - 7))) # 0.819/0.911 = 0.90
  # only on age unified on stage, 0.819/0.911 from 2013 Jean-M (I goofed, the early/normal rate in paper are 0.955/0.911 at Trois; 0.895/0.958 at Chize)
  pvec_env1 <- Pmatrix # normal spring
  pvec_env2 <- Pmatrix %*% env_effect_p # early spring
  Pmatrix_env <- rbind(pvec_env1, pvec_env2)
  Klist_env <- Klist # set list structure
  for (i in seq_len(omega)) {
    Klist_env[[i]] <- cbind(Klist[[i]], Klist[[i]]) ### fertility is the same for both environments
  }
  gamma <- lrs_with_diff_evntrans(D,
                                  Glist,
                                  Klist_env,
                                  Pmatrix_env,
                                  max.offspring)
  return(gamma)
}
