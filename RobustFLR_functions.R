
###############################################
########                              #########
########          Robust FLR          #########
########       Required Functions     #########
###############################################


#### Functions ####################

# Applies Least Squares
# coefX : the FPCA scores for predictor X
# coefY : the FPCA scores for predictor y
MFLR = function(coefX,coefY)  {
  # coefX = basis coefficients of predictors X
  # coefY = basis coefficients of predictors Y
  # Bestimate = regression matrix
  Bestimate = ginv(t(coefX)%*%coefX)%*%t(coefX) %*% coefY
  return(Bestimate)
}

# gives estimated functions using FPCA coefficents
Pred_MFLR = function(new_coefX, Bestimate, harmonics)  {
  # new_coefX = basis coefficents of X
  # Bestimate = regression matrix
  # harmonics = FPCs of Y
  # num_comp = number of FPCs
  num_comp = dim(harmonics$coefs)[2]
  est = fit_function( t(new_coefX%*% Bestimate), harmonics) 
  
  # returns a prediction of corresponding Y function (centred)
  return(est)
}

# used within Pred_MFLR
fit_function = function(scores,harmonics)  {
  # scores = principal coefficents
  # harmonics = FPCs
  # num_comp = number of FPCs
  num_comp = dim(harmonics$coefs)[2]
  h = scores[1]*harmonics[1,]
  for (j in 2:num_comp){
    h = h + scores[j]*harmonics[j,]  
  }
  
  # return FPC approximation
  return(h)
}

# Applies Robust FPCA (Bali et al 2011) with time series directly
# data : time series matrix
# harmonics : FPCs as FDA object
# gaitbasis : Bspline basis used to represent FPCs
# num_comp : number of FPCs to estimate
# N : length of the time series 
robust_FPCA_Bali = function(data, harmonics_tgt, gaitbasis, num_comp, N){
  
  time <- seq(0,1,len=N)
  
  pc <- PcaGrid(data, num_comp)
  #PcaGrid(t(t(Y)-mu.hat.y), num_comp)
  
  plot(getCenter(pc))
  #plot(pc)
  #summary(pc)
  rob_mu = getCenter(pc)
  
  #getEigenvalues(pc)
  
  rob_coef = getScores(pc)
  loadings = getLoadings(pc)
  
  # Converts FPCs into FDA objects
  rob_fd_y = smooth.basisPar(time, loadings, gaitbasis, Lfdobj=NULL, lambda=0)$fd
  
  med_mu_y = smooth.basisPar(time, rob_mu, gaitbasis, Lfdobj=NULL, lambda=0)$fd
  robust_FPCA_y = harmonics_tgt
  for (ii in 1:num_comp){
    robust_FPCA_y$coefs[,ii] = rob_fd_y$coefs[,ii]
  }
  
  # return robust FPC objects
  rob= list()
  rob$harmonics= robust_FPCA_y
  rob$med_mu_y = med_mu_y
  rob$rob_coef = rob_coef
  
  return(rob)
  
}


# Applies Robust FPCA (Bali et al 2011) with B-spline coefficents of time series 
# temp_harmonics : an FPCA FDA object which we will use as a template 
robust_FPCA_Bali_Bspline = function(data, harmonics_tgt, gaitbasis, num_comp, N){
  
  #data = Y
  #num_comp = 3
  gaittime <- seq(0,1,len=N)
  gaitrange <- c(0,1)
  nord=50
  gaitbasis <- create.bspline.basis(c(0,1), nbasis = nord, norder = 3)
  
  #data = t(X.all)
  #num_comp=3
  fdParobj = fdPar(gaitbasis,int2Lfd(2),lambda=0) #Applies smoothness to the functions
  gaitfd <- smooth.basisPar(gaittime, t(data), gaitbasis, Lfdobj=NULL, lambda=0)$fd
  
  resp_pca = pca.fd(gaitfd,nharm=num_comp,fdParobj)
  
  # FPCA components of Y
  harmonics_tgt = resp_pca$harmonics
  
  pc <- PcaGrid(t(gaitfd$coefs), num_comp)
  #PcaGrid(t(t(Y)-mu.hat.y), num_comp)
  
  #plot(getCenter(pc))
  #plot(pc)
  #summary(pc)
  
  
  #getEigenvalues(pc)
  
  rob_coef = getScores(pc)
  loadings = getLoadings(pc)
  
  rob_mu = as.vector(getCenter(pc))
  med_mu_y = fd(rob_mu, gaitbasis)
  
  #mu.hat.y = l1median(X=data,trace=-1) # pcaPP::
  #med_mu_y = smooth.basisPar(gaittime, mu.hat.y, gaitbasis, Lfdobj=NULL, lambda=0)$fd
  
  
  #plot(fd(gaitfd$coefs[,10], gaitbasis))
  robust_FPCA_y = harmonics_tgt
  for (ii in 1:num_comp){
    robust_FPCA_y$coefs[,ii] = loadings[,ii]
  }
  
  # return robust FPC objects
  rob= list()
  rob$harmonics= robust_FPCA_y
  rob$med_mu_y = med_mu_y
  rob$rob_coef = rob_coef
  
  return(rob)
  
}


# Applies MLTS
MLTS_iterated = function(coefX, coefY,nsamples){
  
  library(rrcov)
  library(MASS)
  
  # coefX = basis coefficients of predictors X
  # coefY = basis coefficients of predictors Y
  # Bestimate = regression matrix
  
  # breakdown_p = number of samples used in MLTS
  breakdown_p = round(nsamples*0.8) #round( (nsamples)/2+1)
  
  # make m=1000 draws
  m = 1000
  
  sets = matrix(nrow=m, ncol = breakdown_p)
  min_dist = rep(0,m)
  
  for (ii in 1:m){
    index = sample.int(nsamples, size = breakdown_p)
    
    B = solve(t(coefX[index,])%*%coefX[index,])%*%t(coefX[index,]) %*% coefY[index,]  
    sigma = t(coefY[index,] - coefX[index,]%*%B)%*%(coefY[index,] - coefX[index,]%*%B)#/(nsamples-num_comp)
    det_sigma = det(sigma)
    #if (det_sigma=0){
    #  index = c( index, sample.int(size=1, setdiff(index,1:nsamples)) ) 
    #}
    residuals_LTS = coefY - coefX%*%B
    
    LTS_dist = rep(0,nsamples)
    
    for (jj in 1:nsamples){
      LTS_dist[jj] = t(residuals_LTS[jj,])%*%solve(sigma)%*%residuals_LTS[jj,]
    }
    
    # uses scores from MCD identify samples with smallest error
    order_samples = sort.int(LTS_dist, decreasing = FALSE, index.return = TRUE)
    
    # using samples with smallest error calculate regression matrix 
    new_ind = order_samples$ix[1:breakdown_p]
    Xtrain = coefX[new_ind,]
    Ytrain = coefY[new_ind,]
    
    sets[ii,] = new_ind
    min_dist[ii] = sum(LTS_dist[new_ind])
  }
  
  # best samples
  best_samples = sort.int(min_dist, decreasing = FALSE, index.return = TRUE)
  new_ind = best_samples$ix[1:10]
  
  new_sets = sets[new_ind,]
  best_sets = sets[new_ind,]
  top10 = rep(0,10)
  # apply C-steps on these samples
  for (kk in 1:10){
    det_sigma = 0
    det_sigma_new = 1
    index = new_sets[kk,]
    while(abs(det_sigma_new-det_sigma)>0.1){
      
      B = solve(t(coefX[index,])%*%coefX[index,])%*%t(coefX[index,]) %*% coefY[index,]  
      sigma = t(coefY[index,] - coefX[index,]%*%B)%*%(coefY[index,] - coefX[index,]%*%B)#/(nsamples-num_comp)
      det_sigma = det(sigma)
      
      residuals_LTS = coefY - coefX%*%B
      
      LTS_dist = rep(0,nsamples)
      
      for (jj in 1:nsamples){
        LTS_dist[jj] = t(residuals_LTS[jj,])%*%solve(sigma)%*%residuals_LTS[jj,]
      }
      
      # uses scores from MCD identify samples with smallest error
      order_samples = sort.int(LTS_dist, decreasing = FALSE, index.return = TRUE)
      
      # using samples with smallest error calculate regression matrix 
      index = order_samples$ix[1:breakdown_p]
      B = solve(t(coefX[index,])%*%coefX[index,])%*%t(coefX[index,]) %*% coefY[index,]  
      sigma = t(coefY[index,] - coefX[index,]%*%B)%*%(coefY[index,] - coefX[index,]%*%B)#/(nsamples-num_comp)
      det_sigma_new = det(sigma)
      
      residuals_LTS = coefY - coefX%*%B
      
      LTS_dist = rep(0,nsamples)
      
      for (jj in 1:nsamples){
        LTS_dist[jj] = t(residuals_LTS[jj,])%*%solve(sigma)%*%residuals_LTS[jj,]
      }
      best_dist = sum(LTS_dist[index])
      best_sets[kk,] = index
    }
    
    top10[kk]=best_dist
    
  }
  
  index = best_sets[which.min(top10),]
  Bestimate = solve(t(coefX[index,])%*%coefX[index,])%*%t(coefX[index,]) %*% coefY[index,] 
  return(list(trimmed_samples = index, beta = Bestimate))
}

# Applies MLTS with one-step reweighting to improve efficiency
MLTS_iterated_reweighted = function(coefX, coefY,nsamples){
  
  library(rrcov)
  library(MASS)
  
  # coefX = basis coefficients of predictors X
  # coefY = basis coefficients of predictors Y
  # Bestimate = regression matrix
  M=dim(coefX)[2]
  L=dim(coefY)[2]
  
  # breakdown_p = number of samples used in MLTS
  breakdown_p = round(nsamples*0.8) #round( (nsamples)/2+1)
  MLTS_solution = MLTS_iterated(coefX, coefY,nsamples)
  Bestimate = MLTS_solution$beta
  trimmed_samples = MLTS_solution$trimmed_samples
  residuals_LTS = coefY - coefX%*%Bestimate
  cov_ep =  t(residuals_LTS[trimmed_samples,]) %*% residuals_LTS[trimmed_samples,] #MLTS covariance 
  
  LTS_dist = rep(0,nsamples)
  
  num_comp = L
  # reweighting 
  new_ind = c()
  threshold = qchisq(0.01,num_comp)
  for (jj in 1:nsamples){
    LTS_dist[jj] = t(residuals_LTS[jj,])%*%solve(cov_ep)%*%residuals_LTS[jj,]
    if (LTS_dist[jj]<threshold){
      new_ind = c(new_ind,jj)
    }
  }
  Xtrain = coefX[new_ind,]
  Ytrain = coefY[new_ind,]
  Bestimate = MFLR(Xtrain, Ytrain) #MLTS beta
  
  return(list(trimmed_samples = new_ind, beta = Bestimate))
}


# Gives Fitting error using RFLR with components for the RFPCA: rob_comp_X,rob_comp_Y and trimming given by breakdown_p 
# Sim_output : Object containing information from simulation 
BIC_fit_bali_trim = function(Sim_output, rob_comp_X,rob_comp_Y, breakdown_p){
  
  N = Sim_output$N 
  nsamples  = Sim_output$nsamples
  num_out = Sim_output$num_out
  
  normal_samples = nsamples-num_out
  outlier_samples = nsamples - normal_samples
  
  X.all = Sim_output$Predictor
  Y = Sim_output$Response
  
  set.seed(123)
  gaittime <- seq(0,1,len=N)
  gaitrange <- c(0,1)
  nord=200
  gaitbasis <- create.bspline.basis(c(0,1), nbasis = nord, norder = 3)
  
  oneWipe = array( dim = c(N,nsamples,2) )
  
  oneWipe[,,1] = t(X.all)
  
  
  oneWipe[,,2] = t(Y)
  
  fdParobj = fdPar(gaitbasis,int2Lfd(2),lambda=0) #Applies smoothness to the functions
  gaitfd <- smooth.basisPar(gaittime, oneWipe, gaitbasis, Lfdobj=NULL, lambda=0)$fd
  
  # FPCA on Y
  resp_fd<- gaitfd[,2]
  resp_pca = pca.fd(resp_fd,nharm=rob_comp_Y,fdParobj)
  resp_scores=resp_pca$scores
  
  # FPCA components of Y
  harmonics_tgt = resp_pca$harmonics
  mu_y = resp_pca$meanfd
  print(sum(resp_pca$varprop[1:rob_comp_Y]))
  
  # FPCA on X
  pcs = pca.fd(gaitfd[,1],nharm=rob_comp_X,fdParobj)
  harmonics_x = pcs$harmonics
  pred_scores=pcs$scores
  #print(sum(pcs$varprop[1:comp_X]))
  
  ### Robust FPCA 
  robY = robust_FPCA_Bali(Y, harmonics_tgt, gaitbasis, rob_comp_Y, N)
  
  robust_FPCA_y = robY$harmonics
  med_mu_y = robY$med_mu_y
  rob_coefY = robY$rob_coef
  
  # X
  robX = robust_FPCA_Bali(X.all, harmonics_x, gaitbasis, rob_comp_X, N)
  rob_coefX = robX$rob_coef
  
  beta_est  = MLTS_iterated_trim(rob_coefX, rob_coefY,nsamples,breakdown_p)$beta
  
  Y_hat_rob = matrix(ncol= N, nrow = nsamples)
  PE = rep(0,nsamples)
  rob_residual = matrix(ncol=nsamples, nrow=N)
  for (dd in 1:nsamples){
    Xtest <- rob_coefX[dd,]
    truth=Y[dd,] 
    rec_one_out= Pred_MFLR(Xtest, beta_est, robust_FPCA_y) + med_mu_y
    
    Y_hat_rob[dd,] = eval.fd(rec_one_out, gaittime)
    PE[dd] = inprod(resp_fd[dd,]-rec_one_out,resp_fd[dd,]-rec_one_out)
    rob_residual[,dd] = truth-eval.fd(rec_one_out, gaittime)
  }
  
  dist_FLR = rep(0,normal_samples)
  for (pp in 1:normal_samples){
    dist_FLR[pp] = t(Y[pp,]-Y_hat_rob[pp,])%*%(Y[pp,]-Y_hat_rob[pp,])
  }
  
  rob_fit = sum(dist_FLR)/(normal_samples)
  
  return(rob_fit)
}

# Gives Fitting error using classical FLR with components for the FPCA: comp_X, comp_Y
class_BIC_fit = function(Sim_output, rob_comp_X,rob_comp_Y){
  
  N = Sim_output$N 
  nsamples  = Sim_output$nsamples
  num_out = Sim_output$num_out
  
  normal_samples = nsamples-num_out
  outlier_samples = nsamples - normal_samples
  
  gaittime <- seq(0,1,len=N)
  gaitrange <- c(0,1)
  nord=200
  gaitbasis <- create.bspline.basis(c(0,1), nbasis = nord, norder = 3)
  
  oneWipe = array( dim = c(N,nsamples,2) )
  
  oneWipe[,,1] = t(X.all)
  
  
  oneWipe[,,2] = t(Y)
  
  fdParobj = fdPar(gaitbasis,int2Lfd(2),lambda=0) #Applies smoothness to the functions
  gaitfd <- smooth.basisPar(gaittime, oneWipe, gaitbasis, Lfdobj=NULL, lambda=0)$fd
  
  # FPCA on Y
  resp_fd<- gaitfd[,2]
  resp_pca = pca.fd(resp_fd,nharm=rob_comp_Y,fdParobj)
  resp_scores=resp_pca$scores
  
  # FPCA components of Y
  harmonics_tgt = resp_pca$harmonics
  mu_y = resp_pca$meanfd
  #print(sum(resp_pca$varprop[1:rob_comp_Y]))
  
  # FPCA on X
  pcs = pca.fd(gaitfd[,1],nharm=rob_comp_X,fdParobj)
  harmonics_x = pcs$harmonics
  pred_scores=pcs$scores
  
  coefX=pred_scores
  coefY=resp_scores
  
  # Estimate regression matrix
  classical_beta = MFLR(coefX, coefY) #MLTS_iterated(coefX, coefY, nsamples)$beta #
  #classical_beta  = MLTS_iterated(coefX, coefY,nsamples)$beta
  #### residuals using classical FLR ####  
  
  Y_hat = matrix(ncol= N, nrow = nsamples) # estimates 
  classical_residuals = matrix(ncol=nsamples, nrow=N) #residuals
  
  for (dd in 1:nsamples){
    Xtest <- coefX[dd,]
    truth= Y[dd,] 
    rec_one_out= Pred_MFLR(Xtest, classical_beta, harmonics_tgt) + mu_y
    Y_hat[dd,] = eval.fd(rec_one_out, gaittime)
    
    classical_residuals[,dd] = truth-eval.fd(rec_one_out, gaittime)
  }
  
  #### Classical BIC and RBIC ####
  
  dist_FLR = rep(0,normal_samples)
  for (pp in 1:normal_samples){
    dist_FLR[pp] = t(Y[pp,]-Y_hat[pp,])%*%(Y[pp,]-Y_hat[pp,])
  }
  
  class_fit = sum(dist_FLR)/(normal_samples)
  
  return(class_fit)
}


# Depth values for residuals using RFLR with components: rob_comp_X,rob_comp_Y
rob_depth = function(Sim_output, rob_comp_X,rob_comp_Y){
  
  N = Sim_output$N 
  nsamples  = Sim_output$nsamples
  num_out = Sim_output$num_out
  
  normal_samples = nsamples-num_out
  outlier_samples = nsamples - normal_samples
  
  X.all = Sim_output$Predictor
  Y = Sim_output$Response
  set.seed(123)
  gaittime <- seq(0,1,len=N)
  gaitrange <- c(0,1)
  nord=50
  gaitbasis <- create.bspline.basis(c(0,1), nbasis = nord, norder = 3)
  
  oneWipe = array( dim = c(N,nsamples,2) )
  
  oneWipe[,,1] = t(X.all)
  
  
  oneWipe[,,2] = t(Y)
  
  fdParobj = fdPar(gaitbasis,int2Lfd(2),lambda=0) #Applies smoothness to the functions
  gaitfd <- smooth.basisPar(gaittime, oneWipe, gaitbasis, Lfdobj=NULL, lambda=0)$fd
  
  # FPCA on Y
  resp_fd<- gaitfd[,2]
  resp_pca = pca.fd(resp_fd,nharm=rob_comp_Y,fdParobj)
  resp_scores=resp_pca$scores
  
  # FPCA components of Y
  harmonics_tgt = resp_pca$harmonics
  mu_y = resp_pca$meanfd
  print(sum(resp_pca$varprop[1:rob_comp_Y]))
  
  # FPCA on X
  pcs = pca.fd(gaitfd[,1],nharm=rob_comp_X,fdParobj)
  harmonics_x = pcs$harmonics
  pred_scores=pcs$scores
  print(sum(pcs$varprop[1:rob_comp_X]))
  
  robY = robust_FPCA_Bali(Y, harmonics_tgt, gaitbasis, rob_comp_Y, N)
  
  robust_FPCA_y = robY$harmonics
  med_mu_y = robY$med_mu_y
  rob_coefY = robY$rob_coef
  
  # X
  robX = robust_FPCA_Bali(X.all, harmonics_x, gaitbasis, rob_comp_X, N)
  rob_coefX = robX$rob_coef
  
  beta_est  = MLTS_iterated_reweighted(rob_coefX, rob_coefY,nsamples)$beta
  
  Y_hat_rob = matrix(ncol= N, nrow = nsamples)
  PE = rep(0,nsamples)
  rob_residual = matrix(ncol=nsamples, nrow=N)
  for (dd in 1:nsamples){
    Xtest <- rob_coefX[dd,]
    truth=Y[dd,] 
    rec_one_out= Pred_MFLR(Xtest, beta_est, robust_FPCA_y) + med_mu_y
    
    Y_hat_rob[dd,] = eval.fd(rec_one_out, gaittime)
    PE[dd] = inprod(resp_fd[dd,]-rec_one_out,resp_fd[dd,]-rec_one_out)
    rob_residual[,dd] = truth-eval.fd(rec_one_out, gaittime)
  }
  
  dist_FLR = rep(0,normal_samples)
  for (pp in 1:normal_samples){
    dist_FLR[pp] = t(Y[pp,]-Y_hat_rob[pp,])%*%(Y[pp,]-Y_hat_rob[pp,])
  }
  
  rob_fit = sum(dist_FLR)/(normal_samples)
  
  func_data = fdata(t(rob_residual),argvals=NULL,rangeval=NULL,names=NULL,fdata2d=FALSE)
  plot(func_data, main='Residuals using RFLR')
  
  rob_depth = depth.mode(func_data)
  
  return(rob_depth$dep)
}


# Depth values for residuals using classical FLR with components: comp_X,comp_Y
class_depth = function(Sim_output, comp_X, comp_Y){
  
  N = Sim_output$N 
  nsamples  = Sim_output$nsamples
  num_out = Sim_output$num_out
  
  normal_samples = nsamples-num_out
  outlier_samples = nsamples - normal_samples
  
  X.all = Sim_output$Predictor
  Y = Sim_output$Response
  
  gaittime <- seq(0,1,len=N)
  gaitrange <- c(0,1)
  nord=200
  gaitbasis <- create.bspline.basis(c(0,1), nbasis = nord, norder = 3)
  
  oneWipe = array( dim = c(N,nsamples,2) )
  
  oneWipe[,,1] = t(X.all)
  
  
  oneWipe[,,2] = t(Y)
  
  fdParobj = fdPar(gaitbasis,int2Lfd(2),lambda=0) #Applies smoothness to the functions
  gaitfd <- smooth.basisPar(gaittime, oneWipe, gaitbasis, Lfdobj=NULL, lambda=0)$fd
  
  # FPCA on Y
  resp_fd<- gaitfd[,2]
  resp_pca = pca.fd(resp_fd,nharm=comp_Y,fdParobj)
  resp_scores=resp_pca$scores
  
  # FPCA components of Y
  harmonics_tgt = resp_pca$harmonics
  mu_y = resp_pca$meanfd
  #print(sum(resp_pca$varprop[1:comp_Y]))
  
  # FPCA on X
  pcs = pca.fd(gaitfd[,1],nharm=comp_X,fdParobj)
  harmonics_x = pcs$harmonics
  pred_scores=pcs$scores
  #print(sum(pcs$varprop[1:comp_X]))
  
  coefX=pred_scores
  coefY=resp_scores
  
  # Estimate regression matrix
  classical_beta = MFLR(coefX, coefY)
  
  Y_hat = matrix(ncol= N, nrow = nsamples) # estimates 
  classical_residuals = matrix(ncol=nsamples, nrow=N) #residuals
  
  for (dd in 1:nsamples){
    Xtest <- coefX[dd,]
    truth= Y[dd,] 
    rec_one_out= Pred_MFLR(Xtest, classical_beta, harmonics_tgt) + mu_y
    Y_hat[dd,] = eval.fd(rec_one_out, gaittime)
    
    classical_residuals[,dd] = truth-eval.fd(rec_one_out, gaittime)
  }
  
  func_data_est2 = fdata(t(classical_residuals),argvals=NULL,rangeval=NULL,names=NULL,fdata2d=FALSE)
  #plot(func_data_est2, main='Residuals using LS regression')
  
  
  dist_FLR = rep(0,normal_samples)
  for (pp in 1:normal_samples){
    dist_FLR[pp] = t(Y[pp,]-Y_hat[pp,])%*%(Y[pp,]-Y_hat[pp,])
  }
  
  class_fit = sum(dist_FLR)/(normal_samples)
  
  func_data = fdata(t(classical_residuals),argvals=NULL,rangeval=NULL,names=NULL,fdata2d=FALSE)
  class_depth = depth.mode(func_data)
  
  return(class_depth$dep)
}
