# part2_resume.R
#!/usr/bin/env Rscript
options(repos = c(CRAN = "https://cloud.r-project.org"))

args = commandArgs(trailingOnly=TRUE)
eps_cur = as.double(args[1])
n_cur = as.integer(args[2])
rep_idx = as.integer(args[3])

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
result_file <- paste0("./results/logistic/", filename_prefix, ".csv")

# Load required packages
list.of.packages <- c("ddalpha", "tictoc")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages) > 0){ install.packages(new.packages, dep=TRUE) }
for(package.i in list.of.packages){
  suppressPackageStartupMessages(library(package.i, character.only = TRUE))
}

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

if (!file.exists(result_file)) {
    tic()
    # Load state from part 1
    if (!file.exists(state_file)) {
      stop("State file not found: ", state_file)
    }
    load(state_file)
    cat("[Resume] Loaded state from", state_file, "\n")
    
    # Binary search for right bound only
    right_lower <- theta_accepted_middle_val + tol
    right_upper <- 10
    optim_par_new <- optim_par
    
    while (right_upper - right_lower > tol){
      right_middle <- (right_lower + right_upper) / 2
      search_left_bound <- right_middle
      search_right_bound <- right_upper
      optim_par_new[1] <- (search_left_bound + search_right_bound) / 2
      optim_par_new <- accept(optim_par_new, seed,U,N1,N2,s_dp,ep,n, search_left_bound, search_right_bound)
      if (is.na(optim_par_new[1])){
        right_upper <- search_left_bound
        optim_par_new <- optim_par
      } else {
        right_lower <- optim_par_new[1] + tol
      }
    }
    
    CI_final <- c(left_lower, right_upper)
    cat("[Result] Final CI:", paste(round(CI_final, 4), collapse = ", "), "\n")
    
    time_elapsed = toc()
    
    # Write back to CSV
    result = data.frame(c(CI_final, time_1 + time_elapsed$toc - time_elapsed$tic)) # add the time the first script took to the time the second script took 
    colnames(result) = "1"
    rownames(result) = c("l", "h", "t")
    write.csv(result, result_file)
    t(result)
    
}

