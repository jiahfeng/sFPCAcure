#####################################################
# sim_curerate.R
# Jiahui Feng
# update Dec. 24, 2023
# 
# simulation study to illustrate the implementation 
# of mixture cure rate model
#
#####################################################


library(funData)
library(MASS)
library(caret)
library(readr)
library(survival)
library(BPST)
library(Triangulation)
library(smcure)
library(tidyverse)
library(MFPCA)

source("fun.R")

V.est <- read_rds('V.est') # vertices for triangles for simulation setting
Tr.est <- read_rds('Tr.est') # triangles
Z <- read_rds('Z') # expanded grid points
d.est <- 2 # degree of piecewise polynomials
r <- 1 # smoothness parameter
ind.inside <- inVT(V.est, Tr.est, Z[,1], Z[,2])$ind.inside
#TriPlot(V.est,Tr.est)

### calculate eigenfunctions
# univariate Legendre polynomials basis functions
basis <- eFun(argvals <- seq(0, 1, length.out=40), M = 2,  type = "Poly")@X
# tensor product (ind.inside adopted from ImageInB)
basis1 <- basis[1,] %*% t(basis[2,]); dim(basis1) <- c(1,1600)
basis1[1, -ind.inside] <- NA
basis2 <- basis[2,] %*% t(basis[1,]); dim(basis2) <- c(1,1600)
basis2[1, -ind.inside] <- NA
basis3 <- basis[2,] %*% t(basis[2,]); dim(basis3) <- c(1,1600)
basis3[1, -ind.inside] <- NA
eigenfun <- rbind(basis1, basis2, basis3)


theta <- seq(0.001, 1, length.out=10) # pre-specified theta 
npc <- 1 # number of FPC/sFPC
ncr <- 40 # dimension of the simulation images
lambda <- c(0, 1, 10, 10^2, 10^3, 10^6) # tuning parameter for bernstein splines
trueVals <- c(10, 8, 4)  # true coefficient values for imaging data

N <- 400 #sample size
argvals <- list(seq(0, 1, length.out=40), seq(0,1,length.out=40))
dis <- seq(0, 1, length.out=40)[2]


auc.sFPCA.all <- c()
auc.FPCA.all <- c()
prederr.sFPCA.all <- c()
prederr.FPCA.all <- c()
theta.all <- c()

t <- 100 # number of replicates

for (i in 1:t) {
  set.seed(i)
  ### generate scores
  scores <- t(replicate(N, rnorm(3, sd = sqrt(trueVals))))
  resX <- scores %*% eigenfun; Y <- resX
  
  # integral of C(s)*Z(s) over Omega
  ## Setting 1 ##
  integ_omega <- apply(Y, 1, function(x)
    sum(dis^2*x*(10*basis3), na.rm = T))
  
  # # setting 2 ##
  # integ_omega <- apply(Y, 1, function(x)
  #   sum(dis^2*x*(basis1 + 2*basis2 + 8*basis3), na.rm = T))
  
  
  # generate data
  surv_data <- simul_cure(N=N, kappa = 2, rho = 0.16, covs = integ_omega)
  
  tau <- max(surv_data[surv_data$event==1,]$time) # set tau as largest survival time among the subjects with events
  Y_train <- Y[1:300,]
  Y_test <- Y[301:400,]
  surv_data_train <- surv_data[1:300,]
  surv_data_test <- surv_data[301:400,]
  
  #----------- for training ---------------#
  
  ## 5-fold corss-validation for selecting theta
  folds <- createFolds(1:nrow(Y_train), k = 5)
  auc_all <- c()
  for (j in 1:length(theta)) {
    the <- theta[j]
    auc <- sapply(folds, function(x){
      cv_fun <- sfpca_img(type = "bernstein", Y = Y_train[-x,], train_dat.id = surv_data_train[-x,],
                          theta = the, lambda = lambda, npc = npc, V.est = V.est, Tr.est = Tr.est, d.est = d.est,
                          r = r, Z = Z, ncr = ncr, tau = tau)
      o_tfs_cv <- as.matrix(cv_fun[[1]]) ## supervised FPC
      b.basis_cv <- as.matrix(cv_fun[[3]]); b.scores_cv <- as.matrix(cv_fun[[4]])
      score_sup_cv <- t(o_tfs_cv) %*% (b.basis_cv%*%t(b.scores_cv))
      if(npc == 2){
        score_sup_cv <- t(score_sup_cv)
      }else{
        score_sup_cv <- cbind(score_sup_cv[1,])
      }
      
      score_names <- c()
      for(q in 1:(ncol(score_sup_cv))){
        tname <- paste("score", as.character(q), sep = "")
        score_names <- c(score_names, tname)
      }
      
      tdat_cv <- surv_data_train[-x,]
      tdat_cv <- cbind(tdat_cv, as.matrix(score_sup_cv))
      
      colnames(tdat_cv)[(ncol(tdat_cv) - ncol(score_sup_cv) + 1) : ncol(tdat_cv)] <- score_names
      
      fmla_fpc <- as.formula(paste("Surv(time,event) ~ ", paste(score_names, collapse= "+")))
      fcure_fpc <- as.formula(paste("~ ", paste(score_names, collapse= "+")))
      
      fitted_obj_cv <- quiet(smcure(formula = fmla_fpc,
                                    cureform = fcure_fpc, model = "ph",
                                    data = tdat_cv, Var = F, nboot = 100))
      
      est_cv <- bernstein(Y = Y_train[x,], V.est = V.est, Tr.est = Tr.est, d.est = d.est, r = r, Z = Z,
                         lambda = lambda)
      b.basis.test_cv <- as.matrix(est_cv[[1]]); b.scores.test_cv = as.matrix(t(est_cv[[2]]))
      
      score_sup_test_cv <- t(o_tfs_cv) %*% ((b.basis.test_cv%*%t(b.scores.test_cv)))
      
      if(npc == 2){
        score_sup_test_cv <- t(score_sup_test_cv)
      }else{
        score_sup_test_cv <- cbind(score_sup_test_cv[1,])
      }
      
      test_dat_cv <- surv_data_train[x,]
      pred_cv <- predictsmcure(object = fitted_obj_cv, newX = score_sup_test_cv, newZ = score_sup_test_cv, model = "ph")
      auc_cure(dat <- test_dat_cv, uncureprob = pred_cv$newuncureprob, tau = tau)
    }
    )
    auc_all <- rbind(auc_all, auc)
  }
  
  index_max <- which.max(apply(auc_all, 1, function(x) sum(x, na.rm = TRUE)))
  thetac <- theta[index_max]
  theta.all <- c(theta.all, thetac)
  
  
  ### sFPCA ###
  sfpc_fun <- sfpca_img(type = "bernstein", Y = Y_train, train_dat.id = surv_data_train,
                       theta = thetac, lambda = lambda, npc = npc, V.est = V.est, Tr.est = Tr.est, d.est = d.est,
                       r = r, Z = Z, ncr = ncr, tau = tau)
  o_tfs <- as.matrix(sfpc_fun[[1]])
  
  b.basis <- as.matrix(sfpc_fun[[3]]); b.scores <- as.matrix(sfpc_fun[[4]])
  
  score_sup <- t(o_tfs) %*% (b.basis%*%t(b.scores))
  
  
  if(npc == 2){
    score_sup <- t(score_sup)
  }else{
    score_sup <- cbind(score_sup[1,])
  }
  
  score_names <- c()
  for(q in 1:(ncol(score_sup))){
    tname <- paste("score", as.character(q), sep = "")
    score_names <- c(score_names, tname)
  }
  
  
  tdat.id <- surv_data_train
  tdat.id <- cbind(tdat.id, score_sup)
  
  colnames(tdat.id)[(ncol(tdat.id) - ncol(score_sup) + 1) : ncol(tdat.id)] <- score_names
  
  fmla_fpc <- as.formula(paste("Surv(time,event) ~ ", paste(score_names, collapse= "+")))
  fcure_fpc <- as.formula(paste("~ ", paste(score_names, collapse= "+")))
  
  fitted_obj_sup <- quiet(smcure(formula = fmla_fpc,
                                cureform = fcure_fpc, model = "ph",
                                data = tdat.id, Var = F, nboot = 100))
  
  
  
  ### FPCA ###
  fpc_fun <- fpca_img(type = "bernstein", Y = Y_train, lambda = lambda, npc = npc, V.est = V.est, Tr.est = Tr.est, d.est = d.est, 
                     r = r, Z = Z, ncr = ncr)
  
  o_tfs_fpc <- as.matrix(fpc_fun[[1]]) ## supervised FPC
  
  b.basis_fpc <- as.matrix(fpc_fun[[3]]); b.scores_fpc <- as.matrix(fpc_fun[[4]])
  
  score_fpc <- t(o_tfs_fpc) %*% (b.basis_fpc%*%t(b.scores_fpc)) ## supervised score
  
  if(npc == 2){
    score_fpc <- t(score_fpc)
  }else{
    score_fpc <- cbind(score_fpc[1,])
  }
  
  
  tdat.id_fpc <- surv_data_train
  tdat.id_fpc<- cbind(tdat.id_fpc, score_fpc)
  
  colnames(tdat.id_fpc)[(ncol(tdat.id_fpc) - ncol(score_fpc)+ 1) : ncol(tdat.id_fpc)] <- score_names
  
  fitted_obj_fpc <- quiet(smcure(formula = fmla_fpc,
                                cureform = fcure_fpc, model = "ph",
                                data = tdat.id_fpc, Var = F, nboot = 100))
  
  
  #### TEST ####
  
  est <- bernstein(Y = Y_test, V.est = V.est, Tr.est = Tr.est, d.est = d.est, r = r, Z = Z,
                  lambda = lambda)
  b.basis.test <- as.matrix(est[[1]]); b.scores.test <- as.matrix(t(est[[2]]))
  
  
  #sFPCA
  score_sup_test <- t(o_tfs) %*% ((b.basis.test%*%t(b.scores.test)))
  if(npc == 2){
    score_sup_test <- t(score_sup_test)
  }else{
    score_sup_test <- cbind(score_sup_test[1,])
  }
  
  test_dat.id <- surv_data_test
  
  pred_sup <- predictsmcure(object = fitted_obj_sup, newX = score_sup_test, newZ = score_sup_test, model = "ph")
  auc.sup <- auc_cure(dat = test_dat.id, uncureprob = pred_sup$newuncureprob, tau = tau)
  prederr.sup <- prederr_cure(dat = test_dat.id, uncureprob = pred_sup$newuncureprob, tau = tau)
  
  auc.sFPCA.all <- c(auc.sFPCA.all, auc.sup)
  prederr.sFPCA.all <- c(prederr.sFPCA.all, prederr.sup)
  
  
  #FPCA
  score_fpc_test <- t(o_tfs_fpc) %*% ((b.basis.test%*%t(b.scores.test)))
  if(npc == 2){
    score_fpc_test <- t(score_fpc_test)
  }else{
    score_fpc_test <- cbind(score_fpc_test[1,])
  }
  
  
  pred_fpc <- predictsmcure(object = fitted_obj_fpc, newX = score_fpc_test, newZ = score_fpc_test, model = "ph")
  auc.fpc <- auc_cure(dat = test_dat.id, uncureprob = pred_fpc$newuncureprob, tau = tau)
  prederr.fpc <- prederr_cure(dat = test_dat.id, uncureprob = pred_fpc$newuncureprob, tau = tau)
  
  auc.FPCA.all <- c(auc.FPCA.all, auc.fpc)
  prederr.FPCA.all <- c(prederr.FPCA.all, prederr.fpc)
}



