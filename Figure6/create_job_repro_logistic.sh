#!/bin/sh -l
for ((i = 1; i < 1001; i++))
{
for epsilon in 0.1 0.3 1.0 3.0 10
{
for n in 100 200 500 1000 2000
{
echo "#!/bin/sh -l
#SBATCH --output=/home/chan1074/dp_naivebayes/job_outputs/repro_logistic_$n-$epsilon-$i-part1.out
#SBATCH --error=/home/chan1074/dp_naivebayes/job_outputs/repro_logistic_$n-$epsilon-$i-part1.err
module load r/4.4
cd \$SLURM_SUBMIT_DIR
run_file='repro_logistic_part1.r $epsilon $n $i'
echo \$run_file
Rscript \$run_file" > _generated_job_compare_CI_$n-$epsilon-$i.slurmjob
# sbatch -A statdept -t 4:00:00 -N1 -n123 _generated_job_adaptive.slurmjob
# sbatch -A statdept -t 4:00:00 -N1 -n123 _generated_job_analytical.slurmjob
# sbatch -A statdept -t 4:00:00 -N1 -n123 _generated_job_dpboot.slurmjob
sbatch -A standby -t 4:00:00 -N1 -n124 _generated_job_compare_CI_$n-$epsilon-$i.slurmjob
}
}
}