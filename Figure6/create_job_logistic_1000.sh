#!/bin/sh -l
for ((i = 0; i < 200; i++))
{
for epsilon in 0.1 0.3 1.0 3.0 10
{
for n in 100 200 500 1000 2000
{
echo "#!/bin/sh -l
#SBATCH --output=/home/chan1074/dp_naivebayes/job_outputs/logistic_n=$n-nsim=5_nsimsub=200_ep=$epsilon-R=200-job_idx=$i.out
#SBATCH --error=/home/chan1074/dp_naivebayes/job_outputs/logistic_n=$n-nsim=5_nsimsub=200_ep=$epsilon-R=200-job_idx=$i.err
module load r/4.4
cd \$SLURM_SUBMIT_DIR
run_file='logistic_load_1000.r 5 1 $epsilon $n 200 200 $i'
echo \$run_file
Rscript \$run_file" > _generated_job_compare_CI_$n-$epsilon.slurmjob
# sbatch -A statdept -t 4:00:00 -N1 -n123 _generated_job_adaptive.slurmjob
# sbatch -A statdept -t 4:00:00 -N1 -n123 _generated_job_analytical.slurmjob
# sbatch -A statdept -t 4:00:00 -N1 -n123 _generated_job_dpboot.slurmjob
sbatch -A standby -t 4:00:00 -N1 -n124 _generated_job_compare_CI_$n-$epsilon.slurmjob
}
}
}

