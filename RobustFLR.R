
###############################################
########                              #########
########          Robust FLR          #########
########       Technometrics paper    #########
###############################################


.libPaths("C:/Users/hullait/Documents/R/win-library/3.2")

# Required packages 
library('MASS')
library('fda')
library('fda.usc')
library('rrcov')
library('robustbase')
library('ROCR')

# read functions
source('robustFLR_functions_all.R')

source('Simulations_FLR.R')


scenario = 1 # two scenarios
nsamples = 400 # number of samples

a = 0.2 # the proportion of samples that are outliers
num_out = 400*a # number of outliers

alpha = 0.2 # proportion trimmed in MLTS

# Number of simulations 
nsim = 20 # using 100 simulations is slow

RBIC_vec = rep(0,nsim) # fitting error using RBIC with robust FLR
BIC_vec = rep(0,nsim) # fitting error using BIC with robust FLR
BIC_class_vec = rep(0,nsim) # fitting error using BIC with classical FLR

auc_class = rep(0,nsim) # AUC values using classical FLR
auc_rob = rep(0,nsim) # AUC values using robust FLR
auc_y = rep(0,nsim) # AUC values using DIRECT method

for (jj in 1:nsim){
  
  set.seed(jj)
# Simulation
Output = Simulations(scenario, num_out)

X.all = Output$Predictor
Y = Output$Response
rob_fit = Output$non_out_fit #minimum average fitting error using true eigenfunctions and regression matrix 
out_fit = Output$out_fit

# Model Selection 
source('Model_select_FLR.R')

# max_number : the maximum number of FPCs to be estimated 
FLR_models = Model_select_FLR(Output, max_number=5)
print ("Model selected")

  RBIC_comp_X = FLR_models$RBIC_comp_X
  RBIC_comp_Y = FLR_models$RBIC_comp_Y
  
  comp_X = FLR_models$comp_X
  comp_Y = FLR_models$comp_Y
  
  # Model fit 
  RBIC_val = FLR_models$RBIC_val #model fit using robust model with RBIC
  BIC_val = FLR_models$BIC_val #model fit using robust model with BIC
  BIC_class_val = FLR_models$BIC_class_val #model fit using classical model with BIC
  
  RBIC_vec[jj] = RBIC_val
  BIC_vec[jj] = BIC_val
  BIC_class_vec[jj] = BIC_class_val
  
  #### Depth values ####
  class_depth_val = class_depth(Output, comp_X,comp_Y)
  rob_depth_val = rob_depth(Output, RBIC_comp_X,RBIC_comp_Y)
  
  func_data = fdata(Y,argvals=NULL,rangeval=NULL,names=NULL,fdata2d=FALSE)
  Y_depth = depth.mode(func_data)$dep
  
  
  #### ROC curves ####

  
  labels = c(rep(1, nsamples-num_out), rep(0,num_out))
  pred_class <- prediction( class_depth_val, labels )
  pred_rob <- prediction( rob_depth_val, labels )
  pred_y <- prediction(Y_depth, labels )
  
  perf_class = performance(pred_class, 'auc')
  perf_rob = performance(pred_rob, 'auc')
  perf_y = performance(pred_y, 'auc')
  
  auc_class[jj] = as.numeric(perf_class@y.values)
  auc_rob[jj] = as.numeric(perf_rob@y.values)
  auc_y[jj] = as.numeric(perf_y@y.values)
  
  }


# Plot of ROC curves given in Figure 3 - using set.seed(1)
perf_class <- performance( pred_class, "tpr", "fpr" )
perf_rob <- performance(pred_rob, "tpr", "fpr")
perf_y <- performance(pred_y, "tpr", "fpr")
plot( perf_class, colorize = FALSE, col=2, cex.lab=1.25)
plot(perf_rob, add = TRUE, colorize = FALSE, col=1)
plot(perf_y, add = TRUE, colorize = FALSE, col=3)
legend(0.6, 0.5, c( "Robust", "Classical", "Direct"), col = c(1, 2, 3),
       text.col = "black", lty = c(1,1,1), bg = "gray90") 


# Table 4: AUC values a=0.2, alpha=0.2
direct_auc = sum(auc_y)/nsim
classic_auc = sum(auc_class)/nsim 
robust_auc = sum(auc_rob)/nsim


# Table 2: a=0.2, alpha=0.2
av_RBIC = sum(RBIC_vec)/nsim # RFLR with RBIC
av_BIC = sum(BIC_vec)/nsim # RFLR with BIC
av_BIC_class = sum(BIC_class_vec)/nsim # FLR with BIC
