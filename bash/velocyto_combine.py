args_ID="Pt-25"
print(args_ID) # args_ID
import os
import loompy

path="/athena/elementolab/scratch/yah2014/Projects/scRNAseq-MCL/data/velocyto"
os.chdir(path) # change current path
print(os.getcwd())
# List all filer folder's names.
file_folders=os.listdir(os.getcwd())  # list files
files=[s for s in file_folders if args_ID in s]
print(files)


# on the command line do: cp file1.loom merged.loom
output_filename=args_ID+"_merged.loom"
loompy.combine(files, output_filename, key="Accession")

