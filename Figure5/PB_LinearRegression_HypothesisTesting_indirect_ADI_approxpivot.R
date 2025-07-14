Delta = 2
ep = 1
mu = .5
tau = 1
# X ~ N(mu, tau^2)
beta = c(-.5,1) # first entry is intercept, second is slope.
sa = .5
# Y ~ N(X * beta[2] + beta[1], sa^2)
seed_num <- 0
theta = c(beta[2],beta[1],mu,tau,sa)

n <- 100  # sample size

alpha=.05
B_paramboot=200 # number of paramboot samples
R_indirect_est=50 # number of indirect estimator generated samples
reps=1000 # num of simulation for computing coverage
dir.create(file.path('.', 'results/PB_LR_HT'), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path('.', 'results/compare_clamp'), showWarnings = FALSE, recursive = TRUE)

###################################################
# automatic install of packages if they are not installed already
list.of.packages <- c(
  "foreach",
  "doSNOW",
  "ddalpha",
  "MASS"
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

# show progress bar
pb <- txtProgressBar(max = reps, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

###################################################

clamp = function(x,a,b){
  return(pmin(pmax(x,a),b))
}

truncMean = function(mu,sa,a,b){
  alpha = (a-mu)/sa
  beta = (b-mu)/sa
  Z = pnorm(beta)-pnorm(alpha)
  return(mu+(dnorm(alpha)-dnorm(beta))*sa/Z)
}

truncVar = function(mu,sa,a,b){
  alpha = (a-mu)/sa
  beta = (b-mu)/sa
  Z = pnorm(beta)-pnorm(alpha)
  return((sa^2)*(1+(alpha*dnorm(alpha)-beta*dnorm(beta))/Z-((dnorm(alpha)-dnorm(beta))/Z)^2))
}

clampMean = function(mu,sa,a,b){
  alpha = (a-mu)/sa
  beta = (b-mu)/sa
  Z = pnorm(beta)-pnorm(alpha)
  if(Z==0)
    return(a*pnorm(alpha)+b*(1-pnorm(beta)))
  else
    return(Z*truncMean(mu,sa,a,b)+a*pnorm(alpha)+b*(1-pnorm(beta)))
}

clampE2 = function(mu,sa,a,b){
  alpha = (a-mu)/sa
  beta = (b-mu)/sa
  Z = pnorm(beta)-pnorm(alpha)
  if(Z==0)
    EX2 = (a^2*pnorm(alpha)+b^2*(1-pnorm(beta)))
  else
    EX2 = (a^2*pnorm(alpha)+b^2*(1-pnorm(beta))+
             Z*(truncVar(mu,sa,a,b)+(truncMean(mu,sa,a,b))^2))
  return(EX2)
}

clampVar = function(mu,sa,a,b){
  EX2 = clampE2(mu,sa,a,b)
  return(EX2-(clampMean(mu,sa,a,b))^2)
}

clampXY = function(beta,mu,sa,tau,Delta){
  integrand = function(x){
    base = dnorm(x,mu,tau)
    # factor1 = pmax(pmin(x,Delta),-Delta)
    # factor2 = clampMean(beta[1]+beta[2]*x,sa,-Delta,Delta)
    # result = factor1*factor2*base
    factor12 = clampMean(x*(beta[1]+beta[2]*x),abs(x)*sa,-Delta^2,Delta^2)
    result = factor12*base
    if(is.nan(result)){
      print("Nan in ClampXY")
      return(0)
    }
    else
      return(result)
  }
  integrand_par = function(x){
    sapply(x,integrand)
  }
  result <- tryCatch({
    integrate(integrand_par,lower=-Inf,upper=Inf, rel.tol = .Machine$double.eps^.5)$value
    }, error=function(err) {
      return(.Machine$integer.max)
    }
  )
  return(result)
}
  
sdp_vec = function(ux,uy,N,Delta,theta,ep,n){
  ### ux should be r*n dimensional N(0,1).
  ### uy should be r*n dimensional N(0,1)
  ### N should be r*5 dimensional N(0,1)
  ### theta should be dataframe with mu, tau, beta0, beta1, sa
  beta1=as.numeric(theta[1])
  beta0=as.numeric(theta[2])
  mu=as.numeric(theta[3])
  tau=as.numeric(theta[4])
  sa=as.numeric(theta[5])
  
  x = ux*tau+mu
  y = beta0+beta1*x+uy*sa
  xc = clamp(x,-Delta,Delta)
  yc = clamp(y,-Delta,Delta)
  x2c = clamp(x^2,0,Delta^2) #TODO: the clamping is in this way
  y2c = clamp(y^2,0,Delta^2)
  xyc = clamp(x*y,-Delta^2,Delta^2)
  
  # ep = ep/5
  ep = ep / sqrt(5)
  xbar = apply((xc),1,mean)+2*Delta/(n*ep)*N[,1]
  ybar = apply((yc),1,mean)+2*Delta/(n*ep)*N[,2]
  x2bar = apply((x2c),1,mean)+Delta^2/(n*ep)*N[,3] 
  y2bar = apply((y2c),1,mean)+Delta^2/(n*ep)*N[,4]
  ### last changed to reflect Alabi and Vadhan
  xybar = apply((xyc),1,mean)+2*Delta^2/(n*ep)*N[,5]
  
  return(data.frame(xbar,ybar,x2bar,y2bar,xybar))
}


obj_ADI = function(theta,ux,uy,N,s_dp,Delta,ep,n){
  synth <- sdp_vec(ux,uy,N,Delta,theta,ep,n)
  D_synth <- depth.Mahalanobis(s_dp, synth)
  return(-D_synth) 
}

theta_hat_ADI_func = function(s_dp,Delta,ep,n,R_indirect_est,init=NULL){
  if (is.null(init)) init=c(0,0,0,1,1)

  ux = matrix(rnorm(R_indirect_est*n),ncol=n,nrow=R_indirect_est)
  uy = matrix(rnorm(R_indirect_est*n),ncol=n,nrow=R_indirect_est)
  N = matrix(rnorm(R_indirect_est*5),ncol=5,nrow=R_indirect_est)
  opt = optim(par=init, obj_ADI,
              lower=c(-5,-5,-5,0.001,0.001),
              upper=c(5,5,5,5,5),
              method="L-BFGS-B",
              ux=ux,uy=uy,N=N,
              s_dp=s_dp,Delta=Delta,ep=ep,n=n)
                
  return(opt$par)
}

partial_s_partial_theta_func = function(theta,Delta){
  diff_val = 1e-6
  all_theta_diff = list(theta)
  for (i in c(1:length(theta))){
    new_theta = theta
    new_theta[i] = new_theta[i] + diff_val
    all_theta_diff = append(all_theta_diff, list(new_theta))
  }
  all_beta_diff = list()
  # print("all_theta_diff")
  # print(all_theta_diff)
  for (params in all_theta_diff){
    beta1 = params[1]
    beta0 = params[2]
    mu = params[3]
    tau = params[4]
    sa = params[5]
    
    Ex = clampMean(mu, tau, -Delta, Delta)
    Ey = clampMean(beta0 + beta1*mu, sqrt(beta1^2*tau^2+sa^2), -Delta, Delta)
    Ex2 = clampE2(mu, tau, -Delta, Delta)
    Ey2 = clampE2(beta0+beta1*mu, sqrt(beta1^2*tau^2+sa^2), -Delta, Delta)
    Exy = clampXY(c(beta0,beta1),mu,sa,tau,Delta)

    all_beta_diff = append(all_beta_diff, list(c(Ex, Ey, Ex2, Ey2, Exy)))
  }
  # print("all_beta_diff")
  # print(all_beta_diff)
  partial_s_partial_theta = c()
  original_beta = unlist(all_beta_diff[1])
  for (i in c(1:length(theta))){
    curr_beta = unlist(all_beta_diff[i+1])
    partial_s_partial_theta = c(partial_s_partial_theta, (curr_beta - original_beta) / diff_val)
  }

  return(matrix(partial_s_partial_theta, nrow=length(theta)))
}

test = function(s_dp,beta1,Delta,ep,n, R_indirect_est, init=NULL){
  if (is.null(init)) init=c(1,1,1,1,1)
  theta_hat = theta_hat_ADI_func(s_dp, Delta, ep, n, R_indirect_est, init=init)
  beta1_hat = theta_hat[1]

  R_indirect_est = 200
  ux = matrix(rnorm(R_indirect_est*n),ncol=n,nrow=R_indirect_est)
  uy = matrix(rnorm(R_indirect_est*n),ncol=n,nrow=R_indirect_est)
  N = matrix(rnorm(R_indirect_est*5),ncol=5,nrow=R_indirect_est)
  s_dp_hat = sdp_vec(ux,uy,N,Delta,theta_hat,ep,n)
  s_dp_cov = cov(s_dp_hat)
  partial_s_partial_theta = partial_s_partial_theta_func(theta_hat, Delta)
  # print("partial_s_partial_theta")
  # print(partial_s_partial_theta)
  # print("s_dp_cov")
  # print(s_dp_cov)
  param_hat_cov = ginv(t(partial_s_partial_theta) %*% ginv(s_dp_cov) %*% (partial_s_partial_theta))
  beta1_hat_var = param_hat_cov[1,1]

  return(as.numeric((beta1_hat-beta1)^2 / beta1_hat_var))
}

tests = function(s_dp,beta1,Delta,ep,n,R_indirect_est,init=NULL){
  return(apply(X=s_dp,1,FUN=test,beta1=beta1,Delta=Delta,ep=ep,n=n, R_indirect_est=R_indirect_est,init=init))
}

theta_init_func = function(s_dp){
  xbar=s_dp[1]
  ybar=s_dp[2]
  x2bar=s_dp[3]
  y2bar=s_dp[4]
  xybar=s_dp[5]
  tau2_hat = pmax(1e-12, (x2bar-(xbar)^2))
  tau2_unb = tau2_hat*(n/(n-1))
  beta1_hat = (xybar-xbar*ybar)/tau2_hat
  beta0_hat = ((ybar*x2bar)-(xbar*xybar))/tau2_hat
  S2_hat = (y2bar-2*beta0_hat*ybar-2*beta1_hat*xybar+beta0_hat^2+
              2*beta0_hat*beta1_hat*xbar+beta1_hat^2*x2bar)*(n/(n-2))
  S2_hat = pmax(1e-12, S2_hat)            
  return(c(beta1_hat, beta0_hat, xbar, sqrt(tau2_unb), sqrt(S2_hat)))
}

PB_accept = function(s_dp,beta1,Delta,ep,n,B_paramboot,R_indirect_est){
  theta_init = theta_init_func(s_dp)
  # theta_init = c(0,0,0,1,1)
  T = test(s_dp,beta1,Delta,ep,n,R_indirect_est, init=theta_init)
  print(paste("T", T))
  if(is.nan(T))
    return(TRUE)
  
  theta_hat = theta_hat_ADI_func(s_dp,Delta,ep,n,R_indirect_est, init=theta_init) # here beta1_hat is not 0

  
  thresh = ceiling((B_paramboot+1)*(1-alpha))
  
  ux = matrix(rnorm(B_paramboot*n),ncol=n,nrow=B_paramboot)
  uy = matrix(rnorm(B_paramboot*n),ncol=n,nrow=B_paramboot)
  N = matrix(rnorm(B_paramboot*5),ncol=5,nrow=B_paramboot)
  synth = sdp_vec(ux,uy,N,Delta,theta_hat,ep,n)
  synth_t = tests(synth,theta_hat[1],Delta,ep,n,R_indirect_est, init=theta_hat) # here beta1 becomes beta1_hat
  synth_sort = sort(synth_t)
  ### make sure that it is correct (according to Alabi & Vadhan) 
  ### to "fail to reject" if there are  not enough synthetic samples 
  ### to get the needed quantile
  print(paste("length(synth_sort)", length(synth_sort)))
  print(paste("thresh", thresh))
  print(paste("synth_sort", synth_sort))
  print(paste("synth_sort[thresh]", synth_sort[thresh]))
  if(length(synth_sort)<thresh)
    return(TRUE)
  if(T>synth_sort[thresh])
    return(FALSE)
  else
    return(TRUE)
}

######################

n_list <- c(100,200,300,400,500,1000,2000,5000)
beta1_list <- c(1, 0.8, 0.6, 0.5, 0.4, 0.2, 0.1, 0)
ep_list <- c(1)
# args = commandArgs(trailingOnly=TRUE)

for (ep in ep_list) {
  
  nrow <- length(n_list)
  ncol <- length(beta1_list)
  power_table <- matrix(nrow = nrow, ncol = ncol)
  
  # for (j in c(as.integer(args[1]))){
  for (j in c(1:ncol)) {
    for (i in c(1:(nrow))) {
      n <- n_list[i]
      beta1_true <- beta1_list[j]
      
      beta[2] = beta1_true
      theta = c(beta1_true,beta[1],mu,tau,sa)
      nuisance = c(beta[1],mu,tau,sa)
      rep_list <- c(1:reps) 
      results <- foreach(
        rep = rep_list,
        .combine = 'cbind',
        .packages=c('ddalpha','MASS'),
        .options.snow=opts
      ) %dopar% {
      # for(rep in 1:1){
        set.seed(rep)
        ux = matrix(rnorm(1*n),ncol=n,nrow=1)
        uy = matrix(rnorm(1*n),ncol=n,nrow=1)
        N = matrix(rnorm(1*5),ncol=5,nrow=1)
        s_dp = sdp_vec(ux,uy,N,Delta,theta,ep,n)
        
        if(PB_accept(s_dp,beta1=0,Delta,ep,n,B_paramboot,R_indirect_est))
          rejected = 0
        else
          rejected = 1
        rejected
      }
      power_table[i,j] <- mean(results)
      power_table2 <- data.frame(power_table)
      rownames(power_table2) <- n_list
      colnames(power_table2) <- beta1_list
      
      # txt_name <- paste("./results/PB_LR_HT/PB_LR_HT_ada_indbeta1_adap", "-clamp_", Delta, 
      #                   "-gdp_ep=", ep, "-alpha=", alpha,
      #                   "-X_mu=", mu, "-tau=", tau,
      #                   "-beta=(", beta[1], ",", beta1_true, ")",
      #                   "-sa=", sa,
      #                   "-B=", B_paramboot, "-R=", R_indirect_est, "-reps=", reps, "-", as.integer(args[2]), ".csv", sep='')
      txt_name <- "./results/INT.csv"
      print("")
      print(txt_name)
      print(t(power_table2))
      print(Sys.time())
      write.table(t(power_table2), txt_name, sep=",", row.names=TRUE, col.names=TRUE)
    }
  }
}

