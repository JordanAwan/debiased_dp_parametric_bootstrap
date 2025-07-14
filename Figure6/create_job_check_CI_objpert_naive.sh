#!/bin/sh -l

echo "#!/bin/sh -l
#SBATCH --output=/home/chan1074/dp_naivebayes/job_outputs/check_CI_simulated_objpert_naive_1000_pb200_R200.out
#SBATCH --error=/home/chan1074/dp_naivebayes/job_outputs/check_CI_simulated_objpert_naive_1000_pb200_R200.err
module load r/4.4
cd \$SLURM_SUBMIT_DIR
run_file='check_coverage_objpert_naive.r'
echo \$run_file
Rscript \$run_file" > check_coverage_simulated.slurmjob
# sbatch -A statdept -t 4:00:00 -N1 -n123 _generated_job_adaptive.slurmjob
# sbatch -A statdept -t 4:00:00 -N1 -n123 _generated_job_analytical.slurmjob
# sbatch -A statdept -t 4:00:00 -N1 -n123 _generated_job_dpboot.slurmjob
sbatch -A standby -N1 -n1 check_coverage_simulated.slurmjob
