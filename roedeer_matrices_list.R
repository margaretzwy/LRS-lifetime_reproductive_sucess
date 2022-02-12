# This function uses code from the IPM Roedeer.R.
# It produces all relevant matrices for Roe deer.
# Output is a list with components r, g, and s for Roe deer reproduction, transition and survival;
# each is also a list for the 4 different age classes.


roedeer_matrices_list <- function(){
  # from the IPM_Roedeer.R
  # Four age classes: yearlings, 2-7 years old, 8-11 years old >11 years
  # Return 4 list of parameter matrices
  minsize <- 1 # minimum possible size of a roe deer
  maxsize <- 40 # maximum possible size of a roe deer
  nn <- 200 # nn of classes within the model
  n.age <- 12 # number of ages -- senescents pooled into final class
  max.age <- 50 # maximum possible age a roe deer can live until

  # Part (II) ##############################################################
  # Compute the kernel component functions from the fitted models
  # Survival function -- read in intercept (a) and slope (b) as well as mid-point trait values
  sx<-function(x,a,b) {
    u<-exp(a+b*x)
    return((u/(1+u)))
  }

  # growth kernel -- read in intercepts (a,c) and slopes (b,d) for the mean (a,b) and for the variance (c,d) components of the kernel
  gxy<-function(x,y,a,b,c,d) {
    sigmax2 <- c + d*x
    sigmax2[sigmax2<=0] <- 0.001 # added it arbitrily to make sure G has non zero colunm sums
    sigmax <- sqrt(sigmax2)

    mux <- a+b*x
    fac1 <- sqrt(2*pi)*sigmax
    fac2 <- ((y-mux)^2)/(2*sigmax2)
    return(exp(-fac2)/fac1)
  }

  # a and b -- intercept and slope of reproduce or not; c and d -- intercept and slope of # of offspring.
  rx<-function(x,a,b,c,d){
    u1 <- exp(a+b*x)
    u1 <- u1/(1+u1)
    u2 <- exp(c+d*x)
    u2 <- u2/(1+u2)
    nkids <- u1 * u2
    return(nkids)
  }

  # # inheritance kernel -- read in intercepts (a,c) and slopes (b,d) for the mean (a,b) and for the variance (c,d) components of the kernel
  # dxy<-function(x,y,a,b,c,d) {
  #   mux <- a+b*x
  #   sigmax2 <- c+d*x
  #   sigmax <- sqrt(sigmax2)
  #   fac1<-sqrt(2*pi)*sigmax
  #   fac2<-((y-mux)^2)/(2*sigmax2)
  #   f<-exp(-fac2)/fac1
  #   return(f)
  # }

  ############## The 'big matrix' M of size n x n
  bigmatrix<-function(n,p) {
    # upper and lower integration limits
    L<-0.9*minsize; U<-1.1*maxsize
    # boundary points b and mesh points y
    b<-L+c(0:n)*(U-L)/n
    y<-0.5*(b[1:n]+b[2:(n+1)])
    # create S, R, G and D matrices
    S <- diag(sx(y,p[1],p[2]))
    R <- matrix(c(1 - rx(y,p[7],p[8],p[9],p[10]), rx(y,p[7],p[8],p[9],p[10])), nrow = 2, byrow = T) # the reproduction distribution, kappamat
    G <- t(outer(y,y,gxy,p[3],p[4],p[5],p[6]))# stage transit from colunm to row, sum(colunm)=1, sum(row)!=1
    #D <- t(outer(y,y,dxy,p[11],p[12],p[13],p[14]))
    # scale D and G so columns sum to 1
    G <- G/matrix(as.vector(apply(G,2,sum)),nrow=n,ncol=n,byrow=TRUE)
    #D <- D/matrix(as.vector(apply(D,2,sum)),nrow=n,ncol=n,byrow=TRUE)
    G <- ifelse(is.nan(G),0,G)
    return(list(S=S,R=R,G=G,meshpts=y))
  }

  # Part (IV) ##############################################################

  # parameter list
  # Four age classes: yearlings, 2-7 years old, 8-11 years old >11 years
  # Parameter lists are always the same
  # intercepts (yearlings, adults, 8-11, > 11), body mass slope, random effect for id, random effect for year
  # dummy values put in for yearling fecundity, reproduction and inheritance.  These do not make it into the final model

  surv.params <- c(1.42959,1.42959-0.46893,1.42959-2.07102,1.42959-2.83337,0.06384,0,0.000189)
  fec.params <- c(-1000,-4.1559,-4.1559-0.7671,-4.1559-1.2029,0.1828,0.44628,0.51036) # yearlings don't produce
  rep.params <- c(-1000,-2.09521,-2.09521-0.16431,-2.09521-0.18806,0.13401,0.25776,0.23816)
  growth.params.mean <- c(10.52787,10.52787-2.58298,10.52787-3.20964,10.52787-2.46930,0.68205,0,0.08039)
  growth.params.var <- c(1.54171,1.54171-1.97124,1.54171-1.43275,1.54171-1.85049,0.11441,0.16755,1.41896)
  # inheritance.params.mean <- c(0,3.9785-0.816,3.9785-0.755-0.816,3.9785-0.755-0.816,0.5331,1.9822,0.7319)
  # inheritance.params.var <- c(0,1.29627+0.07541,1.29627-0.11433+0.07541,1.29627-0.11433+0.07541,0.01777,2.4520,0)

  sub.mats.list <- NULL
  for (i in 1:4){
    # params = (SI,SS,GMI,GMS,GVI,GVS,FI,FS,RI,RS,DMI,DMS,DVI,DVS)
    params <- c(surv.params[i],surv.params[5],
                growth.params.mean[i],growth.params.mean[5],
                growth.params.var[i],growth.params.var[5],
                fec.params[i],fec.params[5],
                rep.params[i],rep.params[5])
    sub.mats.list[[i]] <- bigmatrix(nn,params) # return(list(S=S,R=R,G=G,meshpts=y))
  }
  return(sub.mats.list)
}
