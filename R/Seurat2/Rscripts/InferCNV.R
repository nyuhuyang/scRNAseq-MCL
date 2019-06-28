invisible(lapply(c("infercnv","Seurat","magrittr","tidyr"), function(x) {
        suppressPackageStartupMessages(library(x,character.only = T))
}))
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

# SLURM_ARRAY_TASK_ID
set.seed=101
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.numeric(slurm_arrayid)
print(paste0("slurm_arrayid=",args))

# load data and prepare annotation files
(load(file = "data/B_cells_MCL_43_20190627.Rda"))
df_samples <- readxl::read_excel("doc/190626_scRNAseq_info.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
sample_n = which(df_samples$tests %in% c("control",paste0("test",2:12)))
df_samples = df_samples[sample_n,]

meta.data = B_cells_MCL@meta.data
meta.data$orig.ident = gsub("BH|DJ|MD|NZ","Normal",meta.data$orig.ident)
samples = unique(meta.data$orig.ident)
(samples = samples[!(samples %in% "Normal")])

cell.use <- rownames(meta.data)[meta.data$orig.ident %in% c("Normal",samples[args])]

# subset
counts <- B_cells_MCL@assays$RNA@counts[,cell.use]
meta.data <- meta.data[cell.use,c("Barcode","orig.ident")]

meta.data$Barcode = rownames(meta.data)
colnames(meta.data) =NULL
write.table(meta.data, file = paste0(path,samples[args],"_annotations_file.txt"), 
    row.names=FALSE,sep="\t", quote = FALSE)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix = counts,
                            annotations_file=paste0(path,samples[args],"_annotations_file.txt"),
                            delim="\t",
                            gene_order_file="data/gencode_v19_gene_pos.txt",
                            ref_group_names="Normal")
path_infercnv = paste0(path,samples[args],"_infercnv")
if(!dir.exists(path_infercnv)) dir.create(path_infercnv, recursive = T)

infercnv_obj = infercnv::run(infercnv_obj,
                     cutoff=0.1,
                     out_dir=path_infercnv, 
                     cluster_by_groups=TRUE)
