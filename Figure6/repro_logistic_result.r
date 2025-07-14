library(ggplot2)
library(tidyverse)
####################################################################
# problem setting
reps <- 1000

R <- 200
n = 500
alpha=.1
ep=2
shape1=1
shape2=1
beta0=1/2
beta1=beta1_true=2

set.seed(42)
dir.create(file.path('.', 'results/logistic'), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path('.', 'results/summary'), showWarnings = FALSE, recursive = TRUE)

###
# load the computed data from each random seed (rep) and summarize them into a file
n_list = c(100, 200, 500, 1000, 2000)
ep_list = c(0.1, 0.3, 1, 3, 10)

pivot = TRUE

for(ep in ep_list){
  for(n in n_list){
    filename2 = paste("./results/summary/Repro_Logistic",
                      "-shape1_", shape1, "-shape2_", shape2,
                      "-beta0=", beta0, "-beta1=", beta1, "-R_", R, "-n_", n,
                      "-alpha_", alpha, "-ep_", ep, "-reps_", reps, ".csv", sep="")
    print(filename2)
    
    CIs = c()
    for (rep in c(1:reps)){
      filename = paste("./results/logistic/Repro_Logistic",
                      "-shape1_", shape1, "-shape2_", shape2,
                      "-beta0=", beta0, "-beta1=", beta1, "-R_", R, "-n_", n,
                      "-alpha_", alpha, "-ep_", ep, "-reps_", reps, "-", rep, ".csv", sep="")
      tryCatch({
        csv_file = read.csv(filename, header=TRUE)
        CI = csv_file$X1
        CIs = rbind(CIs, CI)
      },
      error=function(cond){print(filename)})
    }
    write.csv(CIs, filename2)
  }
}



reps_ori <- 1000

R <- 200
n = 500
alpha=.1
ep=1
shape1=1
shape2=1
beta0=1/2
beta1=beta1_true=2

# load the results to print the figure
n_list = c(100, 200, 500, 1000, 2000)
ep_list = c(0.1, 0.3, 1, 3, 10)

df_tab_coverage = data.frame(matrix(rep(0, length(n_list) * length(ep_list)), nrow=length(n_list), ncol=length(ep_list)))
colnames(df_tab_coverage) = ep_list
rownames(df_tab_coverage) = n_list
df_tab_width = df_tab_coverage
df_tab_count = df_tab_coverage

coverage_list = c()
width_list = c()
for(n in n_list){
  for(ep in ep_list){
    filename2 = paste("./results/summary/Repro_Logistic",
                      "-shape1_", shape1, "-shape2_", shape2,
                      "-beta0=", beta0, "-beta1=", beta1, "-R_", R, "-n_", n,
                      "-alpha_", alpha, "-ep_", ep, "-reps_", reps_ori, ".csv", sep="")
    # print(filename2)
    if (!file.exists(filename2)) {
      next
    }
    CIs = read.csv(filename2)
    # reps = length(CIs$V1)
    reps = min(1000, length(CIs$V1))
    covered = rep(0, reps)
    width = rep(0, reps)
    time_list = rep(0, reps)
    for (rep in c(1:reps)){
      if(CIs$V1[rep]<=beta1_true & CIs$V2[rep]>=beta1_true)
        covered[rep] = TRUE
      width[rep] = CIs$V2[rep] - CIs$V1[rep]
    }
    
    df_tab_coverage[as.character(n), as.character(ep)] = mean(covered[1:rep])
    df_tab_count[as.character(n), as.character(ep)] = reps
    df_tab_width[as.character(n), as.character(ep)] = mean(width[1:rep])
    coverage_list = rbind(coverage_list, c(ep, as.character(n), mean(covered[1:rep]), sd(covered[1:rep])/sqrt(reps)))
    width_list = rbind(width_list, c(ep, n, mean(width[1:rep]), sd(width[1:rep])))
  }
}

# write the coverage of logictic CI by repro
df_tab_coverage
write.csv(t(df_tab_coverage), "./results/Repro_logistic.csv")
df_tab_width
df_tab_count
df_coverage_list = data.frame(coverage_list)
ep_mapping = c('solid','dashed','dotted','longdash')
colnames(df_coverage_list) = c('ep', 'n', 'coverage', 'coverage_std')
df_width_list2 = data.frame(width_list)
colnames(df_width_list2) = c('ep', 'n', 'width', 'width_std')
df_width_list2$ep = factor(as.character(df_width_list2$ep), levels = c("0.1", "0.3", "1", "3", "10"))
df_width_list2$lt = df_width_list2$ep

df_coverage_list
df_width_list2
