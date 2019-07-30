###############################################
########                              #########
########          Robust FLR          #########
########       Technometrics paper    #########
###############################################

#breakdown_p = round(nsamples*(1-alpha))

# This code chooses parameters for Robust FLR model 
#### Models selection using BIC and RBIC ####

Model_select_FLR = function(Sim_output, max_number){
  
  N = Sim_output$N 
  nsamples  = Sim_output$nsamples
  num_out = Sim_output$num_out
  
  normal_samples = nsamples-num_out
  outlier_samples = nsamples - normal_samples
  
  X.all = Sim_output$Predictor
  Y = Sim_output$Response
  
  breakdown_p = round(nsamples*(1-alpha))
  
  # BIC and RBIC are estimated using classical FPCA  
  BIC = matrix(nrow = max_number,ncol = max_number)
  RBIC = matrix(nrow = max_number,ncol = max_number)
  
  # BIC_rob and RBIC_rob are estimated using robust FPCA
  BIC_rob = matrix(nrow = max_number,ncol = max_number)
  RBIC_rob = matrix(nrow = max_number,ncol = max_number)
  
  
  for(hhh in 2:(max_number+1)){
    for (ggg in 2:(max_number+1)){
      
      #### classical FLR ####
      
      # number of eigenfunctions
      num_comp_X = hhh
      num_comp_Y = ggg
      
      #### setup B-spline basis ####
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
      resp_pca = pca.fd(resp_fd,nharm=num_comp_Y,fdParobj)
      resp_scores=resp_pca$scores
      
      # FPCA components of Y
      harmonics_tgt = resp_pca$harmonics
      mu_y = resp_pca$meanfd
      #print(sum(resp_pca$varprop[1:num_comp_Y]))
      
      # FPCA on X
      pcs = pca.fd(gaitfd[,1],nharm=num_comp_X,fdParobj)
      harmonics_x = pcs$harmonics
      pred_scores=pcs$scores
      #print(sum(pcs$varprop[1:num_comp_X]))
      
      coefX=pred_scores
      coefY=resp_scores
      
      # Estimate regression matrix
      classical_beta = MFLR(coefX, coefY) 
      
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
      
      func_data_class = fdata(t(classical_residuals),argvals=NULL,rangeval=NULL,names=NULL,fdata2d=FALSE)
      #plot(func_data_class, main='Residuals using classical FLR')
      #lines(func_data_class[(normal_samples+1):nsamples,], col=3)
      
      #### Classical BIC and RBIC ####
      
      dist_FLR = rep(0,nsamples)
      for (pp in 1:nsamples){
        dist_FLR[pp] = t(Y[pp,]-Y_hat[pp,])%*%(Y[pp,]-Y_hat[pp,])
      }
      
      scale_FLR = sd(dist_FLR)
      
      index = sort.int(dist_FLR, decreasing = FALSE, index.return = TRUE)
      num_trim = round(0.8*nsamples)
      new_ind = index$ix[1:num_trim]
      
      BIC[hhh-1,ggg-1] = (nsamples)*log( sum(dist_FLR)/(nsamples) ) + (num_comp_X*num_comp_Y+1)*log(nsamples)  
      RBIC[hhh-1,ggg-1] = (num_trim)*log( sum(dist_FLR[new_ind])/(num_trim) ) + (num_comp_X*num_comp_Y+1)*log(num_trim)
      
      
      #### Residuals using robust FPCA ####
      
      robY = robust_FPCA_Bali_Bspline(Y, harmonics_tgt, gaitbasis, num_comp_Y, N)
      
      robust_FPCA_y = robY$harmonics
      med_mu_y = robY$med_mu_y
      rob_coefY = robY$rob_coef
      
      # X
      robX = robust_FPCA_Bali_Bspline(X.all, harmonics_x, gaitbasis, num_comp_X, N)
      
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
        rob_residual[,dd] = truth-eval.fd(rec_one_out, gaittime)
      }
      
      
      func_data_rob = fdata(t(rob_residual),argvals=NULL,rangeval=NULL,names=NULL,fdata2d=FALSE)
      #plot(func_data_rob, main='Residuals using robust FLR', col=1)
      #lines(func_data_rob[(normal_samples+1):nsamples,], col=3)
      
      #### Robust BIC and RBIC ####
      
      dist_FLR = rep(0,nsamples)
      for (pp in 1:nsamples){
        dist_FLR[pp] = t(Y[pp,]-Y_hat_rob[pp,])%*%(Y[pp,]-Y_hat_rob[pp,])
      }
      
      scale_FLR_rob = sd(dist_FLR)
      
      index = sort.int(dist_FLR, decreasing = FALSE, index.return = TRUE)
      num_trim = round(0.8*nsamples)
      new_ind = index$ix[1:num_trim]
      
      BIC_rob[hhh-1,ggg-1] = (nsamples)*log( sum(dist_FLR)/(nsamples) ) + (num_comp_X*num_comp_Y+1)*log(nsamples)  
      RBIC_rob[hhh-1,ggg-1] = (0.8*nsamples)*log( sum(dist_FLR[new_ind])/(nsamples*0.8) ) + (num_comp_X*num_comp_Y+1)*log(nsamples*0.8)
      
    }
    #print ("Model iteration")
  }
  
  #### Model choice ####
  RBIC_model = which(RBIC == min(RBIC), arr.ind = TRUE)
  
  RBIC_rob_model = which(RBIC_rob == min(RBIC_rob), arr.ind = TRUE)
  BIC_model = which(BIC_rob == min(BIC_rob), arr.ind = TRUE)
  
  BIC_rob_model = which(BIC_rob == min(BIC_rob), arr.ind = TRUE)
  
  RBIC_comp_X = RBIC_rob_model[1]+1
  RBIC_comp_Y = RBIC_rob_model[2]+1
  
  BIC_comp_X = BIC_rob_model[1]+1
  BIC_comp_Y = BIC_rob_model[2]+1
  
  comp_X = BIC_model[1]+1
  comp_Y = BIC_model[2]+1
  
  RBIC_rob_val = BIC_fit_bali_trim(Sim_output, RBIC_comp_X,RBIC_comp_Y,breakdown_p)
  BIC_rob_val = BIC_fit_bali_trim(Sim_output, BIC_comp_X,BIC_comp_Y,breakdown_p)
  
  BIC_val = class_BIC_fit(Sim_output, comp_X,comp_Y)
  
return(Models = list(RBIC_comp_X = RBIC_comp_X, RBIC_comp_Y = RBIC_comp_Y, BIC_comp_X = BIC_comp_X, BIC_comp_Y = BIC_comp_Y, comp_X = comp_X, comp_Y = comp_Y,
                     RBIC_val = RBIC_rob_val, BIC_val = BIC_rob_val, BIC_class_val = BIC_val) )

}
