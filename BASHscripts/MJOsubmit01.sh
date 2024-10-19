#!/bin/bash
# Submission script for Nic5 (uLiege cluster) with Slurm usage
# Job parameters
#SBATCH --job-name=joiretJOB
#SBATCH --output=outRibo.txt

# Needed resources
#SBATCH --time=00:05:00 # hh:mm:ss
#
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1000 # megabytes
#SBATCH --partition=batch
#
#SBATCH --mail-user=marc.joiret@uliege.be
#SBATCH --mail-type=ALL
#
#SBATCH --account=Ribosomer
#
# Operations
#module purge
#module load LIST_THE_MODULES_YOU_NEED_HERE

echo "Job started at $(date)"
# Job steps:
#srun ~/bin/myprog < mydata1
srun hostname
echo "Job ended at $(date)"
