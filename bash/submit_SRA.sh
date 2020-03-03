#!/bin/bash -l

#SBATCH --partition=panda   # cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=submit_SRA_%A_%a.txt
#SBATCH --job-name=submit_SRA
#SBATCH --mem=8G   # memory requested, units available: K,M,G,T

#---------------------Variables to be set-------------------------#
fastq_path=/athena/elementolab/scratch/yah2014/Projects/scRNAseq-MCL/data/geo_submission_fastq/
key_file=/athena/elementolab/scratch/yah2014/Projects/scRNAseq-MCL/data/Counts/aspera.openssh

ascp -i $key_file -QT -l100m -k1 -d $fastq_path subasp@upload.ncbi.nlm.nih.gov:uploads/yah2014_med.cornell.edu_wmmoomCu
