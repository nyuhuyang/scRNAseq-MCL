invisible(lapply(c("infercnv","Seurat","magrittr","tidyr","dplyr"), function(x) {
        suppressPackageStartupMessages(library(x,character.only = T))
}))
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

# load data and prepare annotation files
(load(file = "data/B_cells_MCL_43_20190713.Rda"))
df_samples <- readxl::read_excel("doc/190626_scRNAseq_info.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
meta.data_list <- list()
Idents(B_cells_MCL) = "groups"
(groups <- c("Normal","Pt-10", "Pt-25", "Pt-17"))
for (i in seq_along(groups)){
        sub_B <- subset(B_cells_MCL, idents = groups[i])
        meta.data_list[[i]] = sub_B@meta.data
}

meta.data = bind_rows(meta.data_list)
unique(meta.data$orig.ident)
rownames(meta.data) = paste0(meta.data$orig.ident,"_",meta.data$Barcode)
cell.use <- rownames(meta.data)

# subset
counts <- B_cells_MCL@assays$RNA@counts[,cell.use]
meta.data <- meta.data[cell.use,c("Barcode","orig.ident")]
meta.data$orig.ident = gsub("BH|DJ|MD|NZ","Normal",meta.data$orig.ident)

meta.data$Barcode = rownames(meta.data)
colnames(meta.data) =NULL
write.table(meta.data, file = paste0(path,paste(groups,collapse = "_"),"_annotations_file.txt"), 
    row.names=FALSE,sep="\t", quote = FALSE)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix = counts,
                            annotations_file=paste0(path,paste(groups,collapse = "_"),"_annotations_file.txt"),
                            delim="\t",
                            gene_order_file="data/gencode_v19_gene_pos.txt",
                            ref_group_names="Normal",
                            chr_exclude = "chrM")
path_infercnv = paste0(path,paste(samples,collapse = "_"),"_infercnv")
if(!dir.exists(path_infercnv)) dir.create(path_infercnv, recursive = T)

infercnv_obj = infercnv::run(infercnv_obj,
                     cutoff=0.1,
                     out_dir=path_infercnv, 
                     cluster_by_groups=TRUE, 
                     plot_steps = FALSE,
                     denoise= T,
                     HMM = T,
                     analysis_mode = "subclusters",
                     tumor_subcluster_partition_method = "random_trees",
                     num_threads = 8,
                     tumor_subcluster_pval=0.05,
                     png_res=600)
