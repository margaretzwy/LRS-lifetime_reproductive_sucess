# D is the environment transition matrix
lrs_with_diff_evntrans <- function(D,
                                   Glist_env,
                                   Klist_env,
                                   Pmatrix_env,
                                   max.offspring){
  n.age <- length(Glist_env) # total number of ages
  n.stage <- dim(Glist_env[[1]])[2] # number of stages
  n.env <- dim(D)[2]
  Ins <- diag(1, n.stage/n.env) # identity with dimention N_env which is the number of environment

  Phat <- kronecker(D, Ins) # put environment transition matrices into a big matrix
  if(n.age > 1){
    Glist_env_tran <- Klist_env # set G_env the same list structure as kappamat_env
    for (i in seq_len(n.age)) {
      Glist_env_tran[[i]] <- Glist_env[[i]] %*% Phat
    }
  }else{
    Glist_env_tran <- list(Glist_env[[1]] %*% Phat)
  }
  ### Calculate LRS distribution ####
  gamma_evn <- block_matrices_solve(Klist = Klist_env,
                                    Glist = Glist_env_tran,
                                    Pmatrix = Pmatrix_env,
                                    stage_per_age = n.stage,
                                    initial.age = 1,
                                    end.age = 1,
                                    max.offspring = max.offspring)

  return(gamma_evn)
}
