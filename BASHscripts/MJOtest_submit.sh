#!/bin/bash
# Submission script for Nic5 (uLiege cluster) with Slurm usage
# Job parameters
#SBATCH --job-name=joiretJOB
#SBATCH --output=joiretJOB.log

# Needed resources
#SBATCH --time=00:10:00 # hh:mm:ss
#
#SBATCH --ntasks=1
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
#module load releases/2020b
#module load Biopython/1.78-foss-2020b

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
# Activate virtual environment
source $MYENVDIR/bin/activate
# Job steps:
# recall that a backslash allows to continue on the next line
echo "Job started at $(date)"
srun python ${SOFTDIR}/myProgram.py ${DATADIR}/dataJSONyeast.json ${DATADIR}/tunnelElectrostaticsAF.json \
${DATADIR}/CDSfasta01.txt ${FILES[$SLURM_ARRAY_TASK_ID]} ${OUTPUTDIR}/result$SLURM_ARRAY_TASK_ID.txt
echo "Job ended at $(date)"
