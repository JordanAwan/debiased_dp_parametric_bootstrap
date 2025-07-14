
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
R=200 # number of paramboot samples
R_synthetic=200
reps=10000 # num of simulation for computing coverage
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

obj_nuisance = function(nuisance,beta1,s_dp,Delta){
  beta0 = nuisance[1]
  mu = nuisance[2]
  tau = nuisance[3]
  sa = nuisance[4]
  
  Ex = clampMean(mu, tau, -Delta, Delta)
  Ey = clampMean(beta0 + beta1*mu, sqrt(beta1^2*tau^2+sa^2), -Delta, Delta)
  Ex2 = clampE2(mu, tau, -Delta, Delta)
  Ey2 = clampE2(beta0+beta1*mu, sqrt(beta1^2*tau^2+sa^2), -Delta, Delta)
  return((s_dp[1]-Ex)^2+
           (s_dp[2]-Ey)^2+
           (s_dp[3]-Ex2)^2+
           (s_dp[4]-Ey2)^2)
}

obj_nuisance_ada = function(nuisance,ux,uy,N,beta1,s_dp,Delta){
  theta <- c(beta1, nuisance)
  synth <- sdp_vec(ux,uy,N,Delta,theta,ep,n)
  D_synth <- depth.Mahalanobis(s_dp, synth)
  return(-D_synth) 
}

nuisance_hat = function(beta1,s_dp,Delta,init=NULL, method = "adaptive"){
  if (is.null(init)) init=c(1,1,1,1)
  if (method == "analytical") {
    opt = optim(par=init,obj_nuisance,
                lower=c(-5,-5,0.001,0.001),
                upper=c(5,5,5,5),
                method="L-BFGS-B",
                beta1=beta1,s_dp=s_dp,Delta=Delta)
  }
  else {
    if (method == "adaptive") {
      ux = matrix(rnorm(R_synthetic*n),ncol=n,nrow=R_synthetic)
      uy = matrix(rnorm(R_synthetic*n),ncol=n,nrow=R_synthetic)
      N = matrix(rnorm(R_synthetic*5),ncol=5,nrow=R_synthetic)
      opt = optim(par=init, obj_nuisance_ada,
                  lower=c(-5,-5,0.001,0.001),
                  upper=c(5,5,5,5),
                  method="L-BFGS-B",
                  ux=ux,uy=uy,N=N,
                  beta1=beta1,s_dp=s_dp,Delta=Delta)
    } else{
      print("not implemented")
      return(init)
    }
  }
  return(opt$par)
}


theta_hat = function(s_dp,Delta,ep,n,beta1){
  return(c(0, nuisance_hat(beta1, s_dp, Delta, init=nuisance))) 
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
  ### technically above should should be n/(n-r)
  #return((beta1_hat-beta1)^2*tau2_hat/S2_hat)
  return(as.numeric((beta1_hat-beta1)^2*(n*tau2_hat/S2_hat)))
}

tests = function(s_dp,beta1,Delta,ep,n){
  return(apply(X=s_dp,1,FUN=test,beta1=beta1,Delta=Delta,ep=ep,n=n))
}



PB_accept = function(s_dp,beta1,Delta,ep,n,R){
  T = test(s_dp,beta1,Delta,ep,n)
  print(T)
  
  theta0 = theta_hat(s_dp,Delta,ep,n,beta1) # TODO: this only applies to beta1=0

  if(is.nan(T))
    return(c(1, theta0))

  
  thresh = ceiling((R+1)*(1-alpha))
  
  ux = matrix(rnorm(R*n),ncol=n,nrow=R)
  uy = matrix(rnorm(R*n),ncol=n,nrow=R)
  N = matrix(rnorm(R*5),ncol=5,nrow=R)
  synth = sdp_vec(ux,uy,N,Delta,theta0,ep,n)
  synth_t = tests(synth,beta1,Delta,ep,n)
  synth_sort = sort(synth_t)
  ### make sure that it is correct (according to Alabi & Vadhan) 
  ### to "fail to reject" if there are  not enough synthetic samples 
  ### to get the needed quantile
  print(length(synth_sort))
  print(thresh)
  print(synth_sort[thresh])
  if(length(synth_sort)<thresh)
    return(c(2, theta0))
  if(T>synth_sort[thresh])
    return(c(0, theta0))
  else
    return(c(3, theta0))
}

######################

clamp_list <- c(0.5, 0.8, 1.0, 1.5, 2.0, 5.0, 10.0)
n_list <- c(100,200,300,400,500,1000,2000)
ep_list <- c(1)
args = commandArgs(trailingOnly=TRUE)

for (ep in ep_list) {
  
  nrow <- length(n_list)
  ncol <- length(clamp_list)
  power_table <- matrix(nrow = nrow, ncol = ncol)
  
  for (j in c(as.integer(args[1]))){
    for (i in c(1:nrow)) {
      n <- n_list[i]
      Delta <- clamp_list[j]
      beta1_true <- 0
      
      beta[2] = beta1_true
      theta = c(beta1_true,beta[1],mu,tau,sa)
      nuisance = c(beta[1],mu,tau,sa)
      results <- foreach(
        rep = 1:reps,
        .combine = 'rbind',
        .packages=c('ddalpha'),
        .options.snow=opts
      ) %dopar% {
        set.seed(rep)
        ux = matrix(rnorm(1*n),ncol=n,nrow=1)
        uy = matrix(rnorm(1*n),ncol=n,nrow=1)
        N = matrix(rnorm(1*5),ncol=5,nrow=1)
        s_dp = sdp_vec(ux,uy,N,Delta,theta,ep,n)
        
        # if(PB_accept(s_dp,beta1=0,Delta,ep,n,R))
        #   rejected = 0
        # else
        #   rejected = 1
        # rejected
        PB_accept(s_dp,beta1=0,Delta,ep,n,R)
      }      
      txt_name <- paste("./results/PB_LR_HT/PB_LR_HT_ada_o4", "-clamp_", Delta, "-n_", n,
                        "-gdp_ep=", ep, "-alpha=", alpha,
                        "-X_mu=", mu, "-tau=", tau,
                        "-beta=(", beta[1], ", b1)",
                        "-sa=", sa,
                        "-R=", R, "-reps=", reps, "_checkclamp_estden.csv", sep='')
      write.table((results), txt_name, sep=",", row.names=TRUE, col.names=TRUE)
    }
  }
}

