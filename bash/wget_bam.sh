#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=wget_bam_%A_%a.txt
#SBATCH --job-name=wget_bam
#SBATCH --mem=16G  # memory requested, units available: K,M,G,T

#---------------------Variables to be set-------------------------#
path=/athena/elementolab/scratch/yah2014/Projects/scRNAseq-MCL/data/bam
echo $path
cd $path
echo -e "\n\n*** Downloading files for cart: MCL_Pt25_bam...\n";
wget --no-check-certificate --content-disposition --continue --ask-password --user=yah2014@med.cornell.edu "https://abc.med.cornell.edu/pubshare/data/0/100080/2_Pt_25_PB_C1D1-5GEX.bam" "https://abc.med.cornell.edu/pubshare/data/0/100079/3_Pt_25_PB_C1D8-5GEX.bam.bai" "https://abc.med.cornell.edu/pubshare/data/0/100077/5_Pt_25_PB_C25D1.bam.bai" "https://abc.med.cornell.edu/pubshare/data/0/100076/6_Pt_25_AMB_C25D1.bam" "https://abc.med.cornell.edu/pubshare/data/0/100076/6_Pt_25_AMB_C25D1.bam.bai" "https://abc.med.cornell.edu/pubshare/data/0/100081/1_Pt_25_SB_C1D1-5GEX.bam" "https://abc.med.cornell.edu/pubshare/data/0/100079/3_Pt_25_PB_C1D8-5GEX.bam" "https://abc.med.cornell.edu/pubshare/data/0/100078/4_Pt_25_PB_C24D1.bam" "https://abc.med.cornell.edu/pubshare/data/0/100080/2_Pt_25_PB_C1D1-5GEX.bam.bai" "https://abc.med.cornell.edu/pubshare/data/0/100078/4_Pt_25_PB_C24D1.bam.bai" "https://abc.med.cornell.edu/pubshare/data/0/100077/5_Pt_25_PB_C25D1.bam" "https://abc.med.cornell.edu/pubshare/data/0/100081/1_Pt_25_SB_C1D1-5GEX.bam.bai" "https://abc.med.cornell.edu/pubshare/download_list_items/checksums.txt?download_list_id=2331"

echo -e "\n\n*** Running checksum...\n";
sha1sum -c *.sha1;

echo -e "\n\n*** Complete\n";
