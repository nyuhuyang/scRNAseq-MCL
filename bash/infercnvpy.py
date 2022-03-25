import sys
import os
import gc
import scanpy as sc
import infercnvpy as cnv
import matplotlib.pyplot as plt

print('Hello! I am a task number: ', sys.argv[1])

path="/athena/elementolab/scratch/yah2014/Projects/scRNAseq-MCL/shinyApp/PALIBR_I_51/"
os.chdir(path) # change current path
print(os.getcwd())

sc.logging.print_header()

adata = sc.read_h5ad("sc1csr_gexpr.h5ad")
adata.obs["orig.ident"] = adata.obs.index.to_series().str.replace("-.*","",regex = True)

sample_list = list(set(adata.obs["orig.ident"]))
sample_list.sort() # "N01" is the first element

sample = sample_list[sys.argv[1]] # sys.argv[1] >= 0
adata_subset = adata[adata.obs["orig.ident"].isin(["N01", sample])]

cnv.io.genomic_position_from_gtf(gtf_file="../../../../Indexed_genome/refdata-gex-GRCh38-2020-A/genes/genes.gtf",adata = adata_subset, gtf_gene_id = 'gene_name',inplace= True)
adata_subset.var.head()
del adata
gc.collect()
cnv.tl.infercnv(adata_subset,reference_key="orig.ident",reference_cat="N01",window_size=250, n_jobs = 1) #multiprocessing error if n_jobs is not 1
