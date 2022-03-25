#!/bin/bash -l

#SBATCH --nodes=8
#SBATCH --ntasks=1
#SBATCH --output=cnv_%A_%a.txt
#SBATCH --job-name=infercnvpy
#SBATCH --mem=32G  # memory requested, units available: K,M,G,T


# Keep track of information related to the current job
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "Job ID : $JOB_ID"  ${SLURM_ARRAY_TASK_ID}
echo "=========================================================="

file=/athena/elementolab/scratch/yah2014/Projects/scRNAseq-MCL/bash/infercnvpy.py
echo $(ls -l $file)
srun python $file $SLURM_ARRAY_TASK_ID
