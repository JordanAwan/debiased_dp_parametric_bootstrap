# part1_run.R
#!/usr/bin/env Rscript
options(repos = c(CRAN = "https://cloud.r-project.org"))

args = commandArgs(trailingOnly=TRUE)

# Parameters from command line
eps_cur = as.double(args[1])
n_cur = as.integer(args[2])
rep_idx = as.integer(args[3])

# Fixed configuration
reps <- 1000
R <- 200
alpha = .1
ep = eps_cur
n = n_cur
shape1 = 1
shape2 = 1
beta0 = 1/2
beta1 = beta1_true = 2
tol <- 10^-4

theta = c(beta1,beta0,log(shape1),log(shape2))
nuisance = c(beta0,log(shape1),log(shape2))

set.seed(42)
dir.create(file.path('.', 'results/logistic'), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path('.', 'results/summary'), showWarnings = FALSE, recursive = TRUE)

filename_prefix <- paste0("Repro_Logistic-shape1_", shape1,
                          "-shape2_", shape2,
                          "-beta0=", beta0,
                          "-beta1=", beta1,
                          "-R_", R,
                          "-n_", n,
                          "-alpha_", alpha,
                          "-ep_", ep,
                          "-reps_", reps,
                          "-", rep_idx)
state_file <- paste0("./results/logistic/", filename_prefix, "_ci_state.RData")

###################################################
# Load required packages
list.of.packages <- c("foreach", "doSNOW", "ddalpha", "tictoc")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages) > 0){ install.packages(new.packages, dep=TRUE) }
for(package.i in list.of.packages){
  suppressPackageStartupMessages(library(package.i, character.only = TRUE))
}

# Parallel cluster
n.cores <- min(124, parallel::detectCores() - 1)
cl <- makeSOCKcluster(n.cores)
registerDoSNOW(cl)
pb <- txtProgressBar(max = 1, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

# Utility functions (lower, upper, expit, ObjPertLog, sdp_fun, etc.)
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

# objective perturbation
ObjPertLog = function(ep,norm,q,X,Y,N){
  m = dim(X)[2]-1
  lambda = (1/4)*(m+1)
  xi=1
  b=rep(0,m+1)
  if(ep==0){
    Delta=0
  } else {
    if(norm==1){
      xi = (m+1)*2
    }
    if(norm==2){
      xi = sqrt(m+1)*2
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
    logLike = (apply((X)*matrix(Y-exp(X%*%beta)/(1+exp(X%*%beta)),nrow=n,ncol=m+1),2,function(x) sum(x)))
    # 2 means cols
    return( -(1/n)*logLike    + (Delta/n)*beta + b/n)
  }
  #min = optim(par = rep(0,m+1),fn=obj,gr = grad,method = "L-BFGS-B")
  min = optim(par = rep(0,m+1),fn=obj,gr = grad,method = "BFGS")
  
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

# DP mechanism(vectorized version)
sdp_vec = function(seed,U,N1,N2,ep,theta,n){
  s = matrix(rep(0,R*4),nrow=R,ncol=4)
  ## below is tricky to vectorize
  for(i in 1:R){
    s[i,] = sdp_fun(seed[i,],U[i,],N1[i,],N2[i,],ep,theta,n)
  }
  return(s)
}

# repro optimization objection function
scoreMCE_theta = function(theta,seed,U,N1,N2,s_dp,ep,n){
  synth = sdp_vec(seed,U,N1,N2,ep,theta,n)
  synth = rbind(synth,s_dp)  
  D_synth = depth.Mahalanobis(synth,synth)
  r = rank(D_synth, ties.method="max")[R+1]
  s = r + D_synth[R+1]
  return(-s)
}


# check whether the proposal (optim_par) is acceptable
accept = function(optim_par,seed,U,N1,N2,s_dp,ep,n, search_left_bound, search_right_bound){
  proposed_result <- scoreMCE_theta(optim_par, seed,U,N1,N2,s_dp,ep,n)
  if ((-proposed_result) >= (floor((alpha)*(R+1))+1)){
    return(optim_par)
  }

  opt <-
    optim(
      par = optim_par,
      fn = scoreMCE_theta,
      lower = c(search_left_bound, -10, -5, -5),
      upper = c(search_right_bound, 10, 5, 5),
      method = "L-BFGS-B",
      seed = seed,
      U = U,
      N1 = N1,
      N2 = N2,
      s_dp = s_dp,
      ep = ep,
      n = n
    )
  if ((-opt$value) >= (floor((alpha)*(R+1))+1)){
    return(opt$par)
  } else {
    # cannot find any param that is acceptable
    return(c(NA, NA, NA, NA, NA))
  }
}

# --- Modified getCI: Only complete left CI bound ---
getCI = function(seed,U,N1,N2,s_dp,ep,n){
  theta_lower_bound <- -10
  theta_upper_bound <- 10

  optim_par <- c(beta1_true, nuisance)
  search_left_bound <- theta_lower_bound
  search_right_bound <- theta_upper_bound
  optim_par <- accept(optim_par, seed,U,N1,N2,s_dp,ep,n, search_left_bound, search_right_bound)
  if (is.na(optim_par[1])){
    print("failed to find a point in CI")
    return(c(0, 0))
  }
  theta_accepted_middle_val <- optim_par[1]

  ## Binary search for left end
  left_lower <- theta_lower_bound
  left_upper <- theta_accepted_middle_val - tol
  optim_par_new <- optim_par
  while (left_upper - left_lower > tol){
    left_middle <- (left_lower + left_upper) / 2
    search_left_bound <- left_lower
    search_right_bound <- left_middle
    optim_par_new[1] <- (search_left_bound + search_right_bound) / 2
    optim_par_new <- accept(optim_par_new, seed,U,N1,N2,s_dp,ep,n, search_left_bound, search_right_bound)
    if (is.na(optim_par_new[1])){
      left_lower <- search_right_bound
      optim_par_new <- optim_par
    } else {
      left_upper <- optim_par_new[1] - tol
    }
  }

  # Save state for part 2 to compute right bound
  save(seed, U, N1, N2, s_dp, ep, n, theta_accepted_middle_val, optim_par,
       left_lower, left_upper, file = state_file)
  cat("[Checkpoint] Saved left CI state to", state_file, "\n")

  return(c(left_lower, NA))
}

# --- Single repetition ---
CIs <- foreach(
  rep = rep_idx:rep_idx,
  .combine = 'rbind',
  .packages=c('ddalpha', 'tictoc'),
  .options.snow=opts
) %dopar% {
  result_file <- paste0("./results/logistic/", filename_prefix, ".csv")
  if (!file.exists(result_file)) {
    set.seed(rep)
    tic()

    # Generate data for DP mechanism
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
    s_dp = sdp_fun(seed,U,N1,N2,ep,theta,n)

    seed = matrix(runif(n*R,min=0,max=1),nrow=R)
    U = matrix(runif(n*R,min=0,max=1),nrow=R)
    U1 = matrix(runif(2*R,min=-1,max=1),nrow=R)
    G1 = rgamma(R,shape=3,rate = 1)
    N1 = G1*U1
    W1 = matrix(runif(rejection_num*R*2,min=-1,max=1),nrow=rejection_num*R,ncol=2)
    index = which(W1[,2]>=lower(W1[,1]) & W1[,2]<=upper(W1[,1]))
    U2 = W1[index[1:R],]
    G2 = rgamma(R,shape=3,rate = 1)
    N2 = G2*U2

    CI = getCI(seed,U,N1,N2,s_dp,ep,n)
    time_elapsed = toc()
    
    # load all variables again in order to save the elapsed time from the first script 
    load(state_file)
    time_1 = time_elapsed$toc - time_elapsed$tic
    save(seed, U, N1, N2, s_dp, ep, n, theta_accepted_middle_val, optim_par,
       left_lower, left_upper, time_1, file = state_file)
    # result = data.frame(c(CI, time_elapsed$toc - time_elapsed$tic))
    # colnames(result) = "1"
    # rownames(result) = c("l", "NA", "t")
    # write.csv(result, result_file)
    # t(result)
  }
}
stopCluster(cl)
