#!/usr/bin/env Rscript
options(repos = c(CRAN = "https://cloud.r-project.org"))

args = commandArgs(trailingOnly=TRUE)
# problem setting
# reps <- 100

# R <- 20
# n = 100
# alpha=.05
# ep=2
shape1=1
shape2=1
beta0=1/2
beta1=beta1_true=2
tol <- 10^-4

nSIM_range <- as.integer(args[1])
nSIM_start <- as.integer(args[2])
epsilon <- as.double(args[3])
n <- as.integer(args[4])
R_synthetic <- as.integer(args[5])


nSIM_list <- c((1 + (nSIM_start-1)*nSIM_range):(nSIM_start*nSIM_range))
nsim_sub <- 200 # number of paramboot samples
# R_synthetic <- 40
confidence_level_list <- c(0.9)



theta = c(beta1,beta0,log(shape1),log(shape2))
nuisance = c(beta0,log(shape1),log(shape2))


set.seed(42)
dir.create(file.path('.', 'results/logistic'), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path('.', 'results/summary'), showWarnings = FALSE, recursive = TRUE)




###################################################
# automatic install of packages if they are not installed already
list.of.packages <- c(
  "foreach",
  "doSNOW",
  "ddalpha",
  "tictoc"
)
new.packages <- list.of.packages[!(list.of.packages %in% 
                                     installed.packages()[,"Package"])]
if(length(new.packages) > 0){
  install.packages(new.packages, dep=TRUE)
}

# loading packages
for(package.i in list.of.packages){
  suppressPackageStartupMessages(
    library(
      package.i,
      character.only = TRUE
    )
  )
}

# use doSNOW for parallel computing
n.cores <- min(124, parallel::detectCores() - 1)
cl <- makeSOCKcluster(n.cores)
registerDoSNOW(cl)

pb <- txtProgressBar(max = nSIM_range, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

###################################################

# lower and upper bound for perturbation
lower = function(x){
  bottom = (x+1)^2-1
  middle = x-1/4
  top = x^2
  return(bottom*(x<=-1/2)+middle*(x>-1/2 & x<1/2)+top*(x>=1/2))
}

upper = function(x){
  bottom = -x^2
  middle = x+1/4
  top = -(x-1)^2+1
  return(bottom*(x<=-1/2)+middle*(x>-1/2 & x<1/2)+top*(x>=1/2))
}

expit = function(x){
  return(exp(x)/(1+exp(x)))
}

ObjPertLog = function(ep,norm,q,X,Y,N){
  m = dim(X)[2]
  lambda = (1/4)*(m)
  xi=1
  b=rep(0,m)
  if(ep==0){
    Delta=0
  } else {
    if(norm==1){
      xi = (m)*2
    }
    if(norm==2){
      xi = sqrt(m)*2
    }
    if(norm==3){# 3 represents infinity.
      xi = 2
    }
    b = (xi/(q*ep))*N
    
    Delta = lambda/(exp((1-q)*ep)-1)
  }
  
  
  obj = function(beta){
    return( -(1/n)*sum(Y*(X%*%beta) - log(1+exp(X%*%beta)))+Delta/(2*n)*sum(beta^2) + sum(b*beta)/n)
  }
  grad = function(beta){
    logLike = (apply((X)*matrix(Y-exp(X%*%beta)/(1+exp(X%*%beta)),nrow=n,ncol=m),2,function(x) sum(x)))
    # 2 means cols
    return( -(1/n)*logLike    + (Delta/n)*beta + b/n)
  }
  #min = optim(par = rep(0,m),fn=obj,gr = grad,method = "L-BFGS-B")
  min = optim(par = rep(0,m),fn=obj,gr = grad,method = "BFGS")
  
  if(min$convergence!=0){
    print("Obj-Pert did not converge") 
  }
  return(min$par)
}###   END OBJ  PERT

# DP mechanism
sdp_fun = function(seed,U,N1,N2,ep,theta,n){
  beta1=theta[1]
  beta0=theta[2]
  shape1=exp(theta[3])
  shape2=exp(theta[4])
  
  ep1 = ep*.9
  ep2 = (ep-ep1)
  n = length(seed)
  x = qbeta(seed,shape1,shape2)
  x2 = 2*x-1
  Y = U<=expit(beta0+beta1*x2)
  X = matrix(c(rep(1,n),x2),nrow=n,ncol=2)
  
  s12 = ObjPertLog(ep=ep1,norm=3,q=.85,X=X,Y=Y,N=N1)
  s3 = mean(x)+(1/(ep2*n))*N2[1]
  s4 = mean(x^2)+(1/(ep2*n))*N2[2]
  
  return(c(s12,s3,s4))
}


theta_hat_naive <- function(s_dp){

    # get intial values for theta by plugging in obj pert values for beta_0, beta_1
    # and backsolving for a,b in terms of privatized T(x)
    z_bar_private <- s_dp[3]
    s2_private <- s_dp[4] - z_bar_private^2
    a_initial <- max((z_bar_private^2 * (1-z_bar_private)-z_bar_private*s2_private)/s2_private, 0.001)
    b_initial <- max(((1-z_bar_private)^2 * z_bar_private - (1-z_bar_private) * s2_private) / s2_private, 0.001)

    return(c(s_dp[1], s_dp[2], log(a_initial), log(b_initial)))
}

# #### SIM #####
solved_theta_list <- foreach(
  sim_idx = nSIM_list,
  .combine = 'rbind',
  .packages = c("ddalpha", "ald", "cbinom", "extraDistr"),
  .options.snow=opts
) %dopar% {
    set.seed(sim_idx) 
  #   browser() 
    # set.seed(4)
    # generate simulation data
    rejection_num = 100
    seed = runif(n,min=0,max=1)
    U = runif(n,min=0,max=1)

    U1 = runif(2,min=-1,max=1)
    G1 = rgamma(1,shape=3,rate = 1)
    N1 = G1*U1

    W1 = matrix(runif(rejection_num*2*1,min=-1,max=1),nrow=rejection_num*1,ncol=2)
    index = which(W1[,2]>=lower(W1[,1]) & W1[,2]<=upper(W1[,1]))

    U2 = as.numeric(W1[index[1],])
    G2 = rgamma(1,shape=3,rate = 1)
    N2 = G2*U2
  
    s_dp <- sdp_fun(seed,U,N1,N2,epsilon,theta,n)
  
    solved_theta <- theta_hat_naive(s_dp)
    solved_theta
}

txt_name <- paste("results/logistic_naive_compare_solved_theta-gdp_ep=", epsilon, "-n=", n, 
                  "-R_synthetic=", R_synthetic, "-nsim_sub=", nsim_sub, 
                  "-nSIM=1000.rdata", sep='')
file.create(txt_name)

save(solved_theta_list, file = txt_name)

# load(txt_name)

pb2 <- txtProgressBar(max = nSIM_range*nsim_sub, style = 3)
progress2 <- function(n) setTxtProgressBar(pb2, n)
opts2 <- list(progress=progress2)
paramboot_solved_theta_list <- foreach(
  i = 1:(nSIM_range*nsim_sub),
  .combine = 'rbind',
  .packages = c("ddalpha", "ald", "cbinom", "extraDistr"),
  .options.snow=opts2
) %dopar% {
    set.seed(i + (nSIM_start-1)*nSIM_range*nsim_sub)  
    nSIM_idx <- as.integer((i-1) / nsim_sub) + 1
    solved_theta <- solved_theta_list[nSIM_idx,1:length(solved_theta_list[1,])]
    ## generate data noise 
    rejection_num = 100
    seed = runif(n,min=0,max=1)
    U = runif(n,min=0,max=1)

    U1 = runif(2,min=-1,max=1)
    G1 = rgamma(1,shape=3,rate = 1)
    N1 = G1*U1

    W1 = matrix(runif(rejection_num*2*1,min=-1,max=1),nrow=rejection_num*1,ncol=2)
    index = which(W1[,2]>=lower(W1[,1]) & W1[,2]<=upper(W1[,1]))

    U2 = as.numeric(W1[index[1],])
    G2 = rgamma(1,shape=3,rate = 1)
    N2 = G2*U2
    s_dp <- sdp_fun(seed,U,N1,N2,epsilon,solved_theta,n)
  
    paramboot_solved_theta <- theta_hat_naive(s_dp)
    
    paramboot_solved_theta
}
txt_name <- paste("results/logistic_naive_compare_pb_solved_theta-gdp_ep=", epsilon, "-n=", n, 
                  "-R_synthetic=", R_synthetic, "-nsim_sub=", nsim_sub, "-nSIM=", nSIM_range, "_", nSIM_start, 
                  ".rdata", sep='')
# file.create(txt_name)

# save(paramboot_solved_theta_list, file = txt_name)


# CI_results <- c()
# for (i in c(1:nSIM_range)){
#   current_result <- c()
#   paramboot_solved_theta <- paramboot_solved_theta_list[((1+(i-1)*nsim_sub):(i*nsim_sub)),1]
#   for (confidence_level in confidence_level_list) {
#     CI_beta1_ends <- quantile(paramboot_solved_theta,
#                               probs=c((1-confidence_level)/2, (1+confidence_level)/2),
#                               names=FALSE)
#     current_result <- cbind(current_result, 
#                             CI_beta1_ends[1], CI_beta1_ends[2])
#   }
#   CI_results <- rbind(CI_results, current_result)
# }
# CI_results

CI_results <- c()
for (i in c(1:nSIM_range)){
  current_result <- c()
  # Get all parameters for current simulation
  paramboot_solved_theta <- paramboot_solved_theta_list[((1+(i-1)*nsim_sub):(i*nsim_sub)),]
  
  # Calculate CIs for each parameter
  for (param_idx in 1:ncol(paramboot_solved_theta)) {
    for (confidence_level in confidence_level_list) {
      CI_ends <- quantile(paramboot_solved_theta[,param_idx],
                         probs=c((1-confidence_level)/2, (1+confidence_level)/2),
                         names=FALSE)
      current_result <- cbind(current_result, 
                            CI_ends[1], CI_ends[2])
    }
  }
  CI_results <- rbind(CI_results, current_result)
}
CI_results


txt_name <- paste("results/logistic_naive_compare_CI_results-gdp_ep=", epsilon, "-n=", n, 
                  "-R_synthetic=", R_synthetic, "-nsim_sub=", nsim_sub, "-nSIM=", nSIM_range, "_", nSIM_start, 
                  ".rdata", sep='')
file.create(txt_name)

save(CI_results, file = txt_name)

