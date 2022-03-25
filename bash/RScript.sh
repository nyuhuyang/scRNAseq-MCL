#!/bin/bash -l

#SBATCH --partition=panda_physbio   # cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=MCLDE_%A_%a.txt
#SBATCH --job-name=MCLDE
#SBATCH --mem=16G   # memory requested, units available: K,M,G,T

echo We are now running an R script.
echo "Job ID : $JOB_ID"  ${SLURM_ARRAY_TASK_ID}

conda activate r4.0.3
path=/athena/elementolab/scratch/yah2014/Projects/scRNAseq-MCL
cd $path
file=R/Seurat4/PALIBR_51/Rscripts/Differential_analysis_B.R
echo $(ls -l $path/$file)
Rscript $path/$file ${SLURM_ARRAY_TASK_ID}
