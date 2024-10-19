#!/bin/bash
# Submission script for Nic5 (uLiege cluster) with Slurm usage
# Job parameters
#SBATCH --job-name=joiretJOB
#SBATCH --output=joiretJOB.log

# Needed resources
#SBATCH --time=00:07:00 # hh:mm:ss
#
#SBATCH --ntasks=6
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100 # megabytes
#SBATCH --array=0-1
#SBATCH --partition=batch
#
#SBATCH --mail-user=marc.joiret@uliege.be
#SBATCH --mail-type=ALL
#
# Operations
module --force purge
#module load LIST_THE_MODULES_YOU_NEED_HERE
#module load releases/2021b
#module load Python/3.9.6-GCCcore-11.2.0
#module load releases/2021b
#module load Scipy-bundle/2021.10-foss-2021b
module load releases/2020b
module load Python/3.8.6-GCCcore-10.2.0
module load Biopython/1.78-foss-2020b
# Paths for relevant directories
WORKDIR=/home/users/m/j/mjoiret
SOFTDIR=/home/users/m/j/mjoiret/myPyPrograms
DATADIR=/home/users/m/j/mjoiret/myData
INPUTDIR=/home/users/m/j/mjoiret/myInput
OUTPUTDIR=/home/users/m/j/mjoiret/myOutput
MYENVDIR=/home/users/m/j/mjoiret/myENV
FILES=($INPUTDIR/*)
# RATIO=(0.5 1.0 1.5 2.0 2.5 3.0 5.0 10.0)
RATIO=(0.5 1.0 2.0)
# initRATE=(4.3e-6 6e-6 11e-6 16.7e-6 25e-6 50e-6 100e-6 200e-6)
initRATE=(4.3e-6 16.7e-6) # initiation rate (average reference for all transcripts)
# Activate virtual environment
source $MYENVDIR/bin/activate
# Job steps:
# recall that a backslash allows to continue on the next line
echo "Job started at $(date)"
for j in {0..1}
	do
	for i in {0..2}
		do
		srun -N1 -n1 -c1 --exact python ${SOFTDIR}/myProg03.py ${DATADIR}/dataJSONyeast.json ${DATADIR}/tunnelElectrostaticsAF.json \
		${DATADIR}/CDSfasta01.txt ${FILES[$SLURM_ARRAY_TASK_ID]} ${OUTPUTDIR}/ratio${i}init${j}res$SLURM_ARRAY_TASK_ID.txt ${RATIO[$i]} ${initRATE[$j]} \
		${OUTPUTDIR}/debug_ratio${i}init${j}res$SLURM_ARRAY_TASK_ID.txt &
		done
		wait
	done
echo "Job ended at $(date)"
