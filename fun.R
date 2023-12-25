#####################################################
# fun.R
# Jiahui Feng
# update Dec. 24, 2023
# 
# support functions for data simulation and mixture cure rate model
# 
# Functions:
#
# simul_cure    # simulate cure rate model survival outcomes
# auc_cure      # auc based on estimation
# prederr_cure  # prediction error based on estimation
# bernstein     # bernstein polynomial basis functions over triangulation
# fpca.img      # functional principal componenet analysis for imaging data
# sfpca.img     # supervised functional principal componenet analysis for imaging data
# 
#####################################################


simul_cure <- function(N, kappa, rho, covs){
  
  # N: sample size
  # kappa: shape parameter of Weibull distribution
  # rho: scale parameter of Weibull distribution
  # covs: observed covariates 
  
  X0 <- 1 # intercept for the logistic regression, to control the probability of cure
  p <- 1 - exp(covs + X0)/(1+exp(covs + X0)) ## probability of cured
  ind_cure <- sapply(p, function(x) rbinom(1,1,x)) 
  v <- runif(n = N)
  temp <- -log(v) / (exp(covs))
  Tlat <- temp^(1/kappa) / rho
  
  # censoring times
  C <- runif(n = N, min = 0, max = 30)
  # follow-up times and event indicators
  time <- pmin(Tlat, C)
  status <- as.numeric(Tlat <= C)
  status[which(ind_cure==1)] <- 0
  time[ind_cure==1] <- C[ind_cure==1]
  # data set
  data.frame(id=1:N,
             time=time,
             event=status,
             cure=ind_cure)
}

## auc of the mixture cure rate model
auc_cure <- function(dat, uncureprob, tau){
  
  # dat: input data set
  # uncureprob: estimated uncured probability
  # tau: pre-specified maximum survival time that an uncured subject can have
  
  ind_uncure <- which(dat$event == 1)
  ind_cure <- which(dat$event == 0 & dat$time > tau)
  n1 <- length(ind_uncure)
  n0 <- length(ind_cure)
  w <- outer(uncureprob[ind_cure], uncureprob[ind_uncure], "<") + 0.5 * outer(uncureprob[ind_cure], uncureprob[ind_uncure], "==")
  sum(w) / (n0*n1)
}

## prediction error ##
prederr_cure <- function(dat, uncureprob, tau){
  
  # dat: input data set
  # uncureprob: estimated uncured probability
  # tau: pre-specified maximum survival time that an uncured subject can have
  
  dat$uncure <- c(uncureprob)
  dat$minTime <- pmin(dat$time, tau)
  fitC <- survfit(Surv(minTime, (1-event))~1, dat,  type=c("kaplan-meier"))
  KM_c <- summary(fitC,times=dat$time,extend=TRUE) 
  dat$indi <- 1- ifelse(dat$event==0&dat$time<tau, 1, 0)
  dat <- dat[order(dat$time),]
  dat$ipcw <- dat$indi/KM_c$surv
  dat$ipcw[is.nan(dat$ipcw)] <- 0
  dat$ipcw[is.infinite(dat$ipcw)] <- 0
  dat <- dat[dat$ipcw<10,]
  error <- dat$ipcw * (1-dat$cure - dat$uncure)^2
  mean(error)
}



bernstein <- function(Y, V.est, Tr.est, d.est, r, Z, lambda){
  
  # Y: imaging data
  # V.est: vertical points of triangles
  # Tr.est: triangles
  # d.est: degree of piecewise polynomials
  # r: smoothness parameter
  # Z: expanded grid points
  # lambda: tuning parameter
  
  n <- nrow(Y)
  Bfull.est <- basis(V.est,Tr.est,d.est,r,Z)
  B <- Bfull.est$B
  ind.inside <- Bfull.est$Ind.inside
  Q2 <- Bfull.est$Q2
  K <- Bfull.est$K
  Y <- matrix(Y[,ind.inside],nrow=n)
  lambda <- as.matrix(lambda)
  t.area <- Bfull.est$tria.all 
  
  this.call <- match.call()
  n <- nrow(Y)
  npix <- ncol(Y)
  J <- ncol(Q2)
  
  W <- as.matrix(B%*%Q2)
  WW <- crossprod(W,W)
  rhs <- crossprod(W,t(Y))
  D <- crossprod(t(crossprod(Q2,as.matrix(K))),Q2)
  D <- as.matrix(D)
  
  flag <- (rankMatrix(WW)<J)
  if(!flag){
    Ainv <- chol(WW,pivot=TRUE)
    A <- solve(t(Ainv))
    ADA <- A%*%D%*%t(A)
    eigs <- eigen(ADA)
    Cval <- eigs$values
  }
  
  nl <- length(lambda)
  
  gcv_all <- sapply(lambda,FUN=function(Lam){  
    Dlam <- Lam*D
    lhs <- WW+Dlam
    lhs.inv <- chol2inv(chol(lhs));
    theta <- crossprod(t(lhs.inv),rhs)
    gamma <- crossprod(t(Q2),theta)
    Yhat <- crossprod(t(W),theta)
    res <- t(Y)-Yhat
    sse <- apply(res^2,2,sum)
    if(!flag){
      df <- sum(1/(1+Cval*Lam))
    }
    if(flag){
      Hmtx <- crossprod(t(crossprod(t(W),lhs.inv)),t(W))
      df <- sum(diag(Hmtx))
    }
    gcv <- npix*sse/(npix-df)^2
  })
  gcv_all <- matrix(gcv_all,nrow=n)
  lam.ind <- apply(gcv_all,1,which.min)
  lambdac <- lambda[lam.ind]
  
  theta <- c()
  gamma <- c()
  Yhat <- c()
  df <- c()
  for (i in 1:n){
    lamc.tmp <- lambdac[i]
    Dlam <- lamc.tmp*D
    lhs <- WW+Dlam
    lhs.inv <- chol2inv(chol(lhs));
    rhs.tmp <- as.matrix(rhs[,i],ncol=1)
    theta.tmp <- crossprod(t(lhs.inv),rhs.tmp)
    theta <- cbind(theta,theta.tmp)
    gamma.tmp <- crossprod(t(Q2),theta.tmp) 
    gamma <- cbind(gamma,gamma.tmp)
    Yhat.tmp <- crossprod(t(W),theta.tmp)
    Yhat <- cbind(Yhat,Yhat.tmp)
    if(!flag){
      df.tmp <- sum(1/(1+Cval*lamc.tmp))
    }
    if(flag){
      Hmtx <- crossprod(t(crossprod(t(W),lhs.inv)),t(W))
      df.tmp <- sum(diag(Hmtx))
    }
    df <- c(df,df.tmp)
  }
  
  return(list(B, gamma, lambdac, t.area, ind.inside))
}



fpca_img <- function(type = "bernstein", Y, lambda, npc, V.est, Tr.est, d.est, r, Z, ncr){
  
  # type: basis function type
  # Y: imaging data
  # lambda: tuning parameter
  # npc: number of PCs
  # V.est: vertical points of triangles
  # Tr.est: triangles
  # d.est: degree of piecewise polynomials
  # r: smoothness parameter
  # Z: expanded grid points
  # ncr: dimension of the input images
  
  if(type == 'bernstein'){
    est <- bernstein(Y= Y, V.est = V.est, Tr.est = Tr.est, d.est = d.est, 
                    r = r, Z = Z, lambda = lambda)
    basis <- est[[1]];scores <- t(est[[2]]);t.area <- est[[4]];ind.inside <- est[[5]]
  }
  
  org.idx <- c(1:ncr*ncr)
  desMat <- matrix(0, nrow = ncr*ncr, ncol = ncol(scores))
  for(zz in 1:ncol(desMat)){
    desMat[ind.inside,zz] = basis[,zz]
  }
  
  
  B <- aperm(array(desMat, c(ncr, ncr,ncol(scores))), c(3, 1, 2)) 
  mya <- array(dim = c(1, ncr, ncr))
  for(i in 1){
    mya[i,,] <- matrix(Y[i,], ncr, ncr)
  }
  g <- funData(list(c(1:ncr), c(1:ncr)), mya) 
  
  W <- MFPCA:::calcBasisIntegrals(B, 2,g@argvals) 
  
  S <- t(scores)
  sqrM <- function (X) 
  {
    EX <- eigen(X)
    VX <- EX$values
    QX <- EX$vectors
    YX <- QX %*% diag(1/sqrt(VX)) %*% t(QX)
    return(YX)
  }
  
  
  U <-1/nrow(Y)*W%*%S%*%t(S)%*%W; G <- W 
  halfG_inv  <- sqrM(G)
  tM <- t(halfG_inv)%*%U%*%halfG_inv
  eigen_res <- eigen(tM) 
  
  fd_list <- lapply(1:npc, function(ipc){
    coef_pc<- halfG_inv%*%as.matrix(Real(eigen_res$vectors[,ipc])) 
  })
  
  sup_basis <- NULL
  for(k in 1:npc){
    kth <- matrix(fd_list[[k]], nrow = 1) %*% t(basis)
    sup_basis <- cbind(sup_basis, t(kth))
  }
  
  
  return(list(sup_basis, eigen_res$values[1:npc], basis, scores, t.area))
}


sfpca_img <- function(type = "bernstein", Y , train_dat.id, theta, lambda, npc, V.est, Tr.est, d.est, r, Z, ncr, tau){
  
  # type: basis function type
  # Y: imaging data
  # train_dat.id: 
  # theta: tuning parameter for sFPCA
  # lambda: tuning parameter for bernstein splines
  # npc: number of PCs
  # V.est: vertical points of triangles
  # Tr.est: triangles
  # d.est: degree of piecewise polynomials
  # r: smoothness parameter
  # Z: expanded grid points
  # ncr: dimension of the input images
  # tau: pre-specified maximum survival time that an uncured subject can have
  
  
  N <- nrow(train_dat.id)
  surv_dat <- train_dat.id;surv_log = surv_dat
  surv_log$logTime <- log(surv_dat$time) 
  surv_log$minTime <- pmin(surv_log$logTime, log(tau))
  fitC <- survfit(Surv(minTime, (1-event))~1, surv_log,  type=c("kaplan-meier"))
  KM_c <- summary(fitC,times=surv_log$logTime,extend=TRUE) 
  surv_log$indi <- 1- ifelse(surv_log$event==0&surv_log$time<tau, 1, 0)
  surv_dat <- surv_log[order(surv_log$logTime),]
  surv_dat$ipcw <- surv_dat$indi/KM_c$surv
  surv_dat$ipcw[is.nan(surv_dat$ipcw)] <- 0
  surv_dat$ipcw[is.infinite(surv_dat$ipcw)] <- 0
  Y_bar <- mean(surv_dat$ipcw * surv_dat$logTime)
  
  surv_dat$Y_de_mean <- surv_dat$logTime - Y_bar
  surv_dat <- surv_dat[match(surv_log$id, surv_dat$id),]
  train_y <- surv_dat$Y_de_mean * surv_dat$ipcw # here it's already ipcw-ed
  
  if(type == 'bernstein'){
    est <- bernstein(Y= Y, V.est = V.est, Tr.est = Tr.est, d.est = d.est, 
                    r = r, Z = Z, lambda = lambda)
    basis <- est[[1]];scores <- t(est[[2]]);t.area <- est[[4]];ind.inside <- est[[5]]
  }
  
  org.idx <- c(1:ncr*ncr)
  desMat <- matrix(0, nrow = ncr*ncr, ncol = ncol(scores))
  for(zz in 1:ncol(desMat)){
    desMat[ind.inside,zz] = basis[,zz]
  }
  
  B <- aperm(array(desMat, c(ncr, ncr,ncol(scores))), c(3, 1, 2)) 
  mya <- array(dim = c(1, ncr, ncr))
  for(i in 1){
    mya[i,,] <- matrix(Y[i,], ncr, ncr)
  }
  g <- funData(list(c(1:ncr), c(1:ncr)), mya) 
  
  W <- MFPCA:::calcBasisIntegrals(B, 2,g@argvals) 
  
  S <- t(scores)
  maty <- matrix(rep(train_y,each=nrow(S)),nrow=nrow(S))
  M <- rowSums(maty*W%*%S)
  MM <- as.matrix(M)%*%t(as.matrix(M))
  sqrM <- function (X) 
  {
    EX <- eigen(X)
    VX <- EX$values
    QX <- EX$vectors
    YX <- QX %*% diag(1/sqrt(VX)) %*% t(QX)
    return(YX)
  }
  
  
  U <-theta/length(train_y)*W%*%S%*%t(S)%*%W+(1-theta)*MM/(length(train_y)^2); G <- W 
  halfG_inv  <- sqrM(G)
  tM <- t(halfG_inv)%*%U%*%halfG_inv
  eigen_res <- eigen(tM) 
  
  fd_list <- lapply(1:npc, function(ipc){
    coef_pc<- halfG_inv%*%as.matrix(Real(eigen_res$vectors[,ipc])) 
  })
  
  sup_basis <- NULL
  for(k in 1:npc){
    kth <- matrix(fd_list[[k]], nrow = 1) %*% t(basis)
    sup_basis <- cbind(sup_basis, t(kth))
  }
  
  
  return(list(sup_basis, eigen_res$values[1:npc], basis, scores, t.area, surv_dat$ipcw))
}




quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
}