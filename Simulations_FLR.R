
###############################################
########                              #########
########          Robust FLR          #########
########       Technometrics paper    #########
###############################################


###############################################
########                              #########
########        Simulation Study      #########
########                              #########
########                              #########
###############################################

# 3 eigenfunctions in X and 3 in Y which come from Fourier bases
# B generated from unif(-3,3) 
# Outliers generated using R~N(0,0.2)

# num_out = number of outliers to introduce

# There are two scenarios: Scenario 1 gives global outliers whereas Scenario 2 gives localised outliers.

Simulations = function(scenario, num_out){

  
choice = 3 # B generated from unif(-choice,choice) 
out_var = 0.5 # Outliers from B'=B+R, where R~N(0,out_var)
alpha=0.3 #proportion of trimming 


#### predictors ####
N = 500
tt=seq(0,1, by=1/(N-1)) 
tt_other = seq(1:N)
t.y = seq(-1,1, by=2/(N-1))
nsamples = 400
normal_samples = nsamples-num_out
outlier_samples = nsamples - normal_samples


mu=-10*(t.y)^2 + 2

tt_other = seq(1:N)
phi1=sqrt(2/N)*sin(pi*tt_other/N)
phi2=sqrt(2/N)*sin(7*pi*tt_other/N)
phi3=sqrt(2/N)*cos(3*pi*tt_other/N)


#plot(phi1, main="First eigenfunction", xlab="time", ylab="", type="l")
#plot(phi2, main="Second eigenfunction", xlab="time", ylab="", type="l")
#plot(phi3, main="Third eigenfunction", xlab="time", ylab="", type="l")


#plot(mu, main="Mean function", ylab="", xlab="time", type="l")

lambda1=40 #first eigenvalue
lambda2=10 #second eigenvalue
lambda3=1


# normal case
score1=rnorm(nsamples, 0,lambda1)
score2=rnorm(nsamples, 0,lambda2)
score3=rnorm(nsamples, 0,lambda3)

score_x = cbind(score1,score2,score3)

x_values = matrix(nrow= nsamples, ncol = N)
for (i in 1:nsamples){
  x_values[i,]=score1[i]*phi1+score2[i]*phi2+score3[i]*phi3 +mu
}

#### response ####


yphi1=sqrt(2/N)*cos(9*pi*tt_other/N)
yphi2=sqrt(2/N)*sin(5*pi*tt_other/N)
yphi3=sqrt(2/N)*cos(2*pi*tt_other/N)

truemean_y = 60*exp(-(tt-1)^2)

#plot(yphi1, main="First eigenfunction", xlab="time", ylab="", type="l")
#plot(yphi2, main="Second eigenfunction", xlab="time", ylab="", type="l")
#plot(yphi3, main="Third eigenfunction", xlab="time", ylab="", type="l")

# normal case
scorey1=rnorm(nsamples, 0,lambda1)
scorey2=rnorm(nsamples, 0,lambda2)
scorey3=rnorm(nsamples, 0,lambda3)

score_y = cbind(scorey1,scorey2,scorey3)

#variance of the noise
sigma=0.1

num_eig_Y = 3
num_eig_X = 3

# regression matrix  
B = matrix(sample(seq(-choice, choice, 0.01),num_eig_X*num_eig_Y),num_eig_Y,num_eig_X)


epsilon = matrix( rnorm(num_eig_Y*nsamples,mean=0,sd=0.01), num_eig_Y, nsamples) 
y_values = matrix(nrow= nsamples, ncol = N)
for (i in 1:normal_samples){
  y_values[i,] = t(B%*%c(score1[i],score2[i],score3[i]))%*%t(cbind(yphi1,yphi2,yphi3)) +truemean_y + t( (cbind(yphi1,yphi2,yphi3))%*%epsilon[,i])
}

if (scenario==1){
  B_out = B + matrix( rnorm(num_eig_X*num_eig_Y,mean=0,sd=out_var),num_eig_Y,num_eig_X)
  for (i in (normal_samples+1):nsamples){
    y_values[i,] = t(B_out%*%c(score1[i],score2[i],score3[i]))%*%t(cbind(yphi1,yphi2,yphi3)) + truemean_y + t((cbind(yphi1,yphi2,yphi3))%*%epsilon[,i]) 
  }
  
}
if (scenario==2){
  gaittime <- seq(0,1,len=N)
  gaitrange <- c(0,1)
  nord=10
  gaitbasis <- create.bspline.basis(c(0,1), nbasis = nord, norder = 3)
  evaluated_basis = eval.basis(gaitbasis, gaittime)
  normalised_basis = eval.basis(gaitbasis, gaittime)
  
  values = diag(inprod(gaitbasis,gaitbasis)*N)
  intu = sample(nord,1)
  out_basis = normalised_basis[,intu]
  B_out = rbind(B, rnorm(3,3,1))
  for (i in (normal_samples+1):nsamples){
    y_values[i,] = t(cbind(yphi1,yphi2,yphi3,out_basis/values[intu])%*%(B_out%*%c(score1[i],score2[i],score3[i])) +truemean_y) + t( (cbind(yphi1,yphi2,yphi3))%*%epsilon[,i])
    
  }
}

E <- matrix( rnorm(N*nsamples,mean=0,sd=sigma), nsamples, N) 

Y <- y_values + E
X.all = x_values + E

#### plot Predictior and Response curves ####

#matplot(t(X.all), type='l', ylab='X(t)', xlab='time', main='Plot of predictor curves', col=rgb(0,0,0,alpha=0.4))
#matlines(t(X.all[(normal_samples+1):nsamples,]), type='l', lwd=3, lty=1)

#matplot(t(Y), type='l', ylab='Y(t)', xlab='time', main='Plot of response curves', col=rgb(0,0,0,alpha=0.6))
#matlines(t(Y[(normal_samples+1):nsamples,]), type='l', lwd=3, lty=1)

#### Residuals using Truth ####

true_residuals = matrix(nrow=nsamples, ncol=N)
dist_FLR = rep(0,nsamples)
for (pp in 1:nsamples){
  Y_hat = cbind(yphi1,yphi2,yphi3)%*%(B%*%c(score1[pp],score2[pp],score3[pp])) +truemean_y
  true_residuals[pp,] = (Y[pp,]-Y_hat)
  dist_FLR[pp] = t(Y[pp,]-Y_hat)%*%(Y[pp,]-Y_hat)
}

rob_fit = sum(dist_FLR[1:normal_samples])/(normal_samples)
out_fit = sum(dist_FLR[(normal_samples+1):nsamples])/(outlier_samples)

#matplot(t(true_residuals), type='l', ylab='r(t)', xlab='time', main='Plot of residual curves', col=rgb(0,0,0,alpha=0.7))
#matlines(t(true_residuals[(normal_samples+1):nsamples,]), type='l', col=2:6, lwd=1, lty=1)


return(list(Predictor = X.all, Response = Y, non_out_fit = rob_fit, out_fit = out_fit, N=N, nsamples = nsamples, num_out = num_out))

}