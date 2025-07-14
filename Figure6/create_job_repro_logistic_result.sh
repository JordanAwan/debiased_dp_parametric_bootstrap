#!/bin/sh -l
echo "#!/bin/sh -l
#SBATCH --output=/home/chan1074/dp_naivebayes/job_outputs/repro_logistic_result.out
#SBATCH --error=/home/chan1074/dp_naivebayes/job_outputs/repro_logistic_result.err
module load r/4.4
cd \$SLURM_SUBMIT_DIR
run_file='repro_logistic_result.r'
echo \$run_file
Rscript \$run_file" > _generated_job_compare_CI_$n-$epsilon.slurmjob
# sbatch -A statdept -t 4:00:00 -N1 -n123 _generated_job_adaptive.slurmjob
# sbatch -A statdept -t 4:00:00 -N1 -n123 _generated_job_analytical.slurmjob
# sbatch -A statdept -t 4:00:00 -N1 -n123 _generated_job_dpboot.slurmjob
sbatch -A standby -t 4:00:00 -N1 -n124 _generated_job_compare_CI_$n-$epsilon.slurmjob
