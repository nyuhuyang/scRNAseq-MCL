import os
import loompy
import pandas as pd
from openpyxl import load_workbook
from openpyxl.utils.dataframe import dataframe_to_rows

wb = load_workbook("doc/20210715_scRNAseq_info.xlsx")
ws = wb["fastq"]
df = pd.DataFrame(ws.values)
df.columns = df.loc[0,:]
df = df.drop([0],axis = 0)
df =df[df.Sequence.eq("GEX")]
df =df[df.Phase.eq("PALIBR_I")]
df = df.sort_values(by=['id'])


path="/athena/elementolab/scratch/yah2014/Projects/scRNAseq-AIM/data/velocyto"
file_folders=os.listdir(path)  # list files
files=[s for s in file_folders if s in df["Sample.id"].ravel() + ".loom"]
print(files)

# on the command line do: cp file1.loom merged.loom
output_filename="MCL46_merged.loom"
os.chdir(path) # change current path
print(os.getcwd())
loompy.combine(files, os.path.join(path, output_filename), key="Accession")

