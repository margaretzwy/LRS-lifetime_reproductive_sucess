# Using fast Fourier transform (FFT) methods to obtain the lifetime reproductive success distribution
# for stage only cases.


lrs_fft_stage_only <- function(Fmatrix, Umatrix,
                               max.offspring,
                               initial.stage = 1,
                               end.stage = 1,
                               distri_type){# = "poisson"){
  ## For reproduction, there are two type of distribution considered, poisson and bernoulli
  ####===== FFT method to calculate the LRS distribution
  U <- Umatrix
  if(distri_type == "poisson"){
    kappamat <- fpois_to_kappamat(Fmatrix, max.offspring - 1)
  }else if(distri_type == "bernoulli"){
    kappamat <- fbern_to_kappamat(Fmatrix, max.offspring - 1)
  }
  kappamat[is.na(kappamat)] <- 0
  n.kids <- dim(kappamat)[1]
  n.stages <- dim(kappamat)[2]
  p <- fftw::planFFT(n.kids)
  alphahat.mat <- betahat.mat <- matrix(0, nrow = n.kids, ncol = n.stages)
  for (i in 1:n.stages) {
    alphahat.mat[,i] <- fftw::FFT(kappamat[,i], plan=p)
  }
  I <- diag(1, n.stages)
  e <- rep(1, n.stages)
  for (i in 1:n.kids) {
    W <- diag(alphahat.mat[i,])
    betahat.mat[i,] <- W %*% solve(I - t(U) %*% W) %*% (I - t(U)) %*% e
  }
  gamma <- Re(fftw::IFFT(betahat.mat[, initial.stage], plan = p))
  if(end.stage < initial.stage){
    end.stage <- initial.stage
  }
  if(end.stage > initial.stage){
    for (i in (initial.stage + 1):end.stage) {
      gamma <- cbind(gamma, Re(fftw::IFFT(betahat.mat[,i], plan = p)))
    }
    colnames(gamma) <- paste("stage", initial.stage:end.stage, sep = "")
  }
  return(gamma)
}
