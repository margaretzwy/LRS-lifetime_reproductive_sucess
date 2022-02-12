#########=======     Age-stage cases   ==========#########
### an example for block method ####

agestage_roedeer <- function(newborn_stage = 50){
  ###==== get Roe deer matrices
  matrices <- roedeer_matrices_list()
  n.stages <- dim(matrices[[1]]$R)[2] # obtain number of stages
  kappamat <- r_to_kappamat(matrices) # list of n.stages matrix, each of them has 12 columns for 12 ages and 2 rows for 0 or 1 fawns
  G <- g_to_gmat(matrices) # list of n.stages matrix, each of them is n.stages X n.stages
  pvec <- s_to_pvec(matrices) #size conditional survival, n.stages stages, row is corresponded to stage
  max.offspring <- 16
  #===== use block matrices method
  ptm <- proc.time()
  gamma <- block_matrices_solve(Klist = kappamat, Glist = G, Pmatrix = pvec,
                                max.offspring = max.offspring,
                                initial.age = 1,
                                end.age = 1,
                                stage_per_age = n.stages)
  proc.time() - ptm
  #plot the LRS distribution and the theoretic distribution
  initial.stage <- newborn_stage
  gamma.distr <- data.frame(density = gamma[, initial.stage],
                            mids = seq_along(gamma[,1]) - 1)
  plot_lrs_distr_allpop("Roe deer", gamma.distr)
  mtext(paste0("Roe deer born at stage ", initial.stage))
  return(gamma)
}
