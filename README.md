This code accompanies the paper "Optimal Debiased Inference on Privatized Data \\via Indirect Estimation and Parametric Bootstrap" by Zhanyu Wang, Arin Chang, and Jordan Awan

# Introduction
The code is for the paper 'Optimal Debiased Inference on Privatized Data
via Indirect Estimation and Parametric Bootstrap' by Zhanyu Wang, Arin Chang, and Jordan Awan. 

We denote our parametric bootstrap method as PB.

## Table 1 and Figure 2 are not using PB.
Table 1 is from https://arxiv.org/pdf/2210.06140 Figure 7. See the `README.md` file about 'Figure 7' in https://github.com/Zhanyu-Wang/Differentially_Private_Bootstrap (the data are in `results/simulation_mean_N=10000_mu=*_nSIM=2000.txt`).

Figure 2 is from https://arxiv.org/pdf/2303.05328 Figure 6. See the `README.md` file about 'Figure 6' in https://github.com/Zhanyu-Wang/Simulation-based_Finite-sample_Inference_for_Privatized_Data.

## Figure 4: Comparison of the sampling distribution of different estimates
Go to the folder `./Figure4_Table2`.

Run `PB_CI_Normal_comparison.r` for the naive estimates and simplified-t estimates. 

Run `PB_CI_Normal_ADI.r` for our estimates. 

Then run `compare_bias.r` and see the pdf in `results/bias_comparison.pdf`.

## Table 2: Comparison of the CIs of different methods
Go to the folder `./Figure4_Table2`.

Run `PB_CI_Normal_comparison.r` to generate the compared CI results in `estimates/paramboot_normal_comparison_results.csv`. 

Run `PB_CI_Normal_ADI.r` for our CI results in `estimates/paramboot_normal_meansd-clamp_0_3-gdp_mu=1-conf=0.95-N=100-nSIM=1000-R=50-B=200-seed=0.adaptive_indirect.csv`.

## Figure 5: Comparison of the rejection probability on linear regression coefficient
Go to the folder `./Figure5`.

For 'method = Monte Carlo (Naive estimator + F-statistic)', run `PB_LinearRegression_HypothesisTesting_naive.r` for the result using the method by (Alabi and Vadhan, 2022). 

For 'method = Repro', the results for Repro are from https://arxiv.org/pdf/2303.05328 Figure 5. See the `README.md` file about 'Figure 5' in https://github.com/Zhanyu-Wang/Simulation-based_Finite-sample_Inference_for_Privatized_Data (the data are in `linear_logistic/results/Repro_LR_HT/Repro_LR_HT-clamp_2-gdp_ep=1-alpha=0.05-X_mu=0.5-tau=1-beta=(-0.5, b1)-sa=0.5-R=200-reps=1000.csv`). 

For 'method = Monte Carlo (Indirect estimator + F-statistic)', run  `PB_LinearRegression_HypothesisTesting_indirect_ADI_Fstat_Theta0.r` for the result using the adaptive indirect estimator for PB and using privatized F-statistic for testing. 

For 'method = Monte Carlo (Indirect estimator + approximate-pivot)', run  `PB_LinearRegression_HypothesisTesting_indirect_ADI_approxpivot.r` for the result using the adaptive indirect estimator for PB and using approximate pivot statistic for testing. 

With the outputs ready in the folder `./results`, run `plot_paramboot_Heatmaps.ipynb` to generate the figure `LR_HT_gdp=1.pdf`.

## Figure 6: Comparison of the coverage and width of confidence intervals for logistic regression parameter
Go to the folder ./Figure6.

For the result of the debiased estimator, run `logistic_load_1000.r` with epsilon and n arguments varying over the appropirate range (n=100, 200, 500, 1000, 2000; epsilon=0.1, 0.3, 1.0, 3.0, 10.0) followed by `check_coverage_objpert_1000.r` to generate matrices containing the coverage and width measurements across 1000 replicates. The bash script `create_job_logistic_1000.sh` will submit jobs that run the `logistic_load_1000.r` file with appropriate parameters and `create_job_check_CI_objpert_1000.R` runs `check_coverage_objpert_1000.r` approriately. 

For the result of the naive estimator, run `logistic_load_naive.r` with appropriate epsilon and n arguments followed by `check_coverage_objpert_naive.r` to generate matrices containing the coverage and width measurements across 1000 replicates. The bash script `create_job_logistic_naive.sh` will similarly submit jobs that run `logistic_load_naive.r` with appropriate parameters and `create_job_check_CI_objpert_naive.sh` runs `check_coverage_objpert_naive.r`. 

For the result of the repro method, run `repro_logistic_part1.r` followed by `repro_logistic_part2.r`. The two files split up the computation so that each file does not take too long to run. Next run `repro_logistic_result.r` to get the coverage and width results. The file `create_job_repro_logistic.sh` can be used to submit jobs that run `repro_logistic_part1.r` and `repro_logistic_part2.r`. The same bash script can be used to run both files with a single modification to the name of the file being run and appropriately named output files. 

Finally, update the data matrices in `plot_objpert.py` to hold the values in the results of the above three experiments and run `plot_objpert.py` to generate the figure of Figure 6. 


## Appendix Figures 7, 8, 9, 10
Go to the folder `./appendix_Figure7_8_9_10`. The data are saved in `./estimates/` which are generated by running `Normal_estimates_ADI_tune.r` and `Normal_estimates_IND_tune.r`. 

Run `compare_bias.r` to generate the pdfs: 

`bias_comparison_indirectR.pdf` for Figure 7, 

`bias_comparison_gdp.pdf` for Figure 8, 

`bias_comparison_clamp.pdf` for Figure 9, 

`bias_comparison_clamp_IND.pdf` for Figure 10.

## Appendix Figures 11
Go to the folder `./appendix_Figure11`. The data are saved in `./results/` which are generated by running `PB_LinearRegression_HypothesisTesting_indirect_ada_only4_checkclamp_estden_10k.R`. Then run `PB_LR_HT_ind_bias.R` to generate the pdf `bias_comparison_HT_LR_null.pdf` for Figure 11.

