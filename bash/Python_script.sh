#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=velocyto_%A_%a.txt
#SBATCH --job-name=velocyto
#SBATCH --mem=20G  # memory requested, units available: K,M,G,T


# Keep track of information related to the current job
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "Job ID : $JOB_ID"  $SGE_TASK_ID"
echo "=========================================================="

file=/athena/elementolab/scratch/yah2014/Projects/scRNAseq-Glioma/bash/velocyto_single.py
echo $(ls -l $file)
python $file $SGE_TASK_ID
