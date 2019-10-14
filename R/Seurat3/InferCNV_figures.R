invisible(lapply(c("Seurat","infercnv","magrittr","tidyr"), function(x) {
        suppressPackageStartupMessages(library(x,character.only = T))
}))
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

# load data and prepare annotation files
(load(file = "data/B_cells_MCL_43_20190713.Rda"))
Idents(B_cells_MCL) = "orig.ident"
df_samples <- readxl::read_excel("doc/190626_scRNAseq_info.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
args=2
sample_n = which(df_samples$cnv %in% args)
df_samples = df_samples[sample_n,]
df_samples

subset_B_cells_MCL = subset(B_cells_MCL, idents = c("DJ","BH","MD","NZ",df_samples$sample))

# Extracting HMM Features
path_infercnv = paste0(path,paste(df_samples$sample,collapse = "_"),"_infercnv/")
if(!dir.exists(path_infercnv)) dir.create(path_infercnv, recursive = T)
subset_B_cells_MCL = add_to_seurat(infercnv_output_path=path_infercnv,
                                     seurat_obj=subset_B_cells_MCL, # optional
                                     top_n=10)
FeaturePlot.1(subset_B_cells_MCL,reduction="tsne", features="proportion_scaled_cnv_chr6")
