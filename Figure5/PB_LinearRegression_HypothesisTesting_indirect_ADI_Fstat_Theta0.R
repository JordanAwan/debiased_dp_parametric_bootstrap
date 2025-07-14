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
  "ddalpha"
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


obj_ADI = function(nuisance,beta1,ux,uy,N,s_dp,Delta,ep,n){
  theta <- c(beta1, nuisance)
  synth <- sdp_vec(ux,uy,N,Delta,theta,ep,n)
  D_synth <- depth.Mahalanobis(s_dp, synth)
  return(-D_synth) 
}

theta_hat_ADI_func = function(s_dp,Delta,ep,n,R_indirect_est,beta1,init=NULL){
  if (is.null(init)) init=c(0,0,1,1)

  ux = matrix(rnorm(R_indirect_est*n),ncol=n,nrow=R_indirect_est)
  uy = matrix(rnorm(R_indirect_est*n),ncol=n,nrow=R_indirect_est)
  N = matrix(rnorm(R_indirect_est*5),ncol=5,nrow=R_indirect_est)
  opt = optim(par=init, obj_ADI,
              lower=c(-5,-5,0.001,0.001),
              upper=c(5,5,5,5),
              method="L-BFGS-B",
              beta1=beta1,
              ux=ux,uy=uy,N=N,
              s_dp=s_dp,Delta=Delta,ep=ep,n=n)
                
  return(opt$par)
}

test = function(s_dp,beta1,Delta,ep,n){
  r=2
  xbar=s_dp[1]
  ybar=s_dp[2]
  x2bar=s_dp[3]
  y2bar=s_dp[4]
  xybar=s_dp[5]
  tau2_hat = (x2bar-(xbar)^2)
  tau2_unb = tau2_hat*(n/(n-1))
  
  if(tau2_hat<=0)
    return(NaN)
  beta1_hat = (xybar-xbar*ybar)/tau2_hat
  beta0_hat = ((ybar*x2bar)-(xbar*xybar))/tau2_hat
  # (wrong formula on Alabi's paper)
  # S2_hat = (y2bar-2*beta0_hat*ybar-2*beta1_hat*xybar+beta0_hat^2+
  #             2*beta0_hat*beta1_hat*xbar+beta1_hat^2*xybar)*(n/(n-r))
  S2_hat = (y2bar-2*beta0_hat*ybar-2*beta1_hat*xybar+beta0_hat^2+
              2*beta0_hat*beta1_hat*xbar+beta1_hat^2*x2bar)*(n/(n-r))
  # (wrong formula on Alabi's paper)
  # S2_hat0 = (y2bar-2*beta0_hat*ybar+beta0_hat^2)*(n/(n-r))
  S2_hat0 = (y2bar-ybar^2)*(n/(n-1))
  if(S2_hat<=0)
    return(NaN)
  if(S2_hat0<=0)
    return(NaN)
  return(as.numeric((beta1_hat-beta1)^2*(n*tau2_hat/S2_hat)))
}

tests = function(s_dp,beta1,Delta,ep,n){
  return(apply(X=s_dp,1,FUN=test,beta1=beta1,Delta=Delta,ep=ep,n=n))
}

theta_init_func = function(s_dp){
  xbar=s_dp[1]
  ybar=s_dp[2]
  x2bar=s_dp[3]
  y2bar=s_dp[4]
  xybar=s_dp[5]
  tau2_hat = pmax(1e-12, (x2bar-(xbar)^2))
  tau2_unb = tau2_hat*(n/(n-1))
  S2_hat0 = (y2bar-ybar^2)*(n/(n-1))
  S2_hat0 = pmax(1e-12, S2_hat0)            
  return(c(0, ybar, xbar, sqrt(tau2_unb), sqrt(S2_hat0)))
}

PB_accept = function(s_dp,beta1,Delta,ep,n,B_paramboot,R_indirect_est){
  T = test(s_dp,beta1,Delta,ep,n)
  print(paste("T", T))
  if(is.nan(T))
    return(TRUE)
  
  nuisance_hat = theta_hat_ADI_func(s_dp,Delta,ep,n,R_indirect_est,beta1=beta1, init=theta_init_func(s_dp)[2:5]) # or theta_init_func(s_dp)[2:5]
  theta_hat = c(beta1, nuisance_hat)
  
  thresh = ceiling((B_paramboot+1)*(1-alpha))
  
  ux = matrix(rnorm(B_paramboot*n),ncol=n,nrow=B_paramboot)
  uy = matrix(rnorm(B_paramboot*n),ncol=n,nrow=B_paramboot)
  N = matrix(rnorm(B_paramboot*5),ncol=5,nrow=B_paramboot)
  synth = sdp_vec(ux,uy,N,Delta,theta_hat,ep,n)
  synth_t = tests(synth,theta_hat[1],Delta,ep,n) # here beta1 becomes beta1_hat
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
args = commandArgs(trailingOnly=TRUE)

for (ep in ep_list) {
  
  nrow <- length(n_list)
  ncol <- length(beta1_list)
  power_table <- matrix(nrow = nrow, ncol = ncol)
  
  for (j in c(as.integer(args[1]))){
    for (i in c(1:nrow)) {
      n <- n_list[i]
      beta1_true <- beta1_list[j]
      
      beta[2] = beta1_true
      theta = c(beta1_true,beta[1],mu,tau,sa)
      nuisance = c(beta[1],mu,tau,sa)
      results <- foreach(
        rep = 1:reps,
        .combine = 'cbind',
        .packages=c('ddalpha'),
        .options.snow=opts
      ) %dopar% {
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
      
      # txt_name <- paste("./results/PB_LR_HT/PB_LR_HT_ADI_Fstat_Theta0", "-clamp_", Delta, 
      #                   "-gdp_ep=", ep, "-alpha=", alpha,
      #                   "-X_mu=", mu, "-tau=", tau,
      #                   "-beta=(", beta[1], ",", beta1_true, ")",
      #                   "-sa=", sa,
      #                   "-B=", B_paramboot, "-R=", R_indirect_est, "-reps=", reps, ".csv", sep='')
      txt_name <- "./results/IN.csv"
      print("")
      print(txt_name)
      print(t(power_table2))
      write.table(t(power_table2), txt_name, sep=",", row.names=TRUE, col.names=TRUE)
    }
  }
}

