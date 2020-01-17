invisible(lapply(c("Seurat","infercnv","magrittr","tidyr"), function(x) {
        suppressPackageStartupMessages(library(x,character.only = T))
}))
source("../R/Seurat3_functions.R")
source("R/util.R")

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

# load data and prepare annotation files
(load(file = "data/B_cells_MCL_43_20190713.Rda"))
Idents(B_cells_MCL) = "orig.ident"
df_samples <- readxl::read_excel("doc/190626_scRNAseq_info.xlsx")
colnames(df_samples) %<>% tolower
(groups = grep("Normal",df_samples$group,value = T, invert = TRUE) %>% unique)
args =8
sample_n = which(df_samples$group %in% groups[args])
df_samples = df_samples[sample_n,]
df_samples
subset_B_cells_MCL = subset(B_cells_MCL, idents = c("DJ","BH","MD","NZ",df_samples$sample))

# Extracting HMM Features
path_infercnv = paste0(path,paste(df_samples$sample,collapse = "_"),"_infercnv/")
if(!dir.exists(path_infercnv)) dir.create(path_infercnv, recursive = T)
subset_B_cells_MCL = add_to_seurat.1(infercnv_output_path=path_infercnv,
                                     seurat_obj=subset_B_cells_MCL, # optional
                                     top_n=10)
subset_B_cells_MCL$orig.ident %<>% gsub("Pt-25-","",.)
subset_B_cells_MCL$orig.ident %<>% gsub("-","-\n",.)
subset_B_cells_MCL$orig.ident %<>% factor(levels = c("BH","DJ","MD","NZ","SB-\nC1",
                                                        "C1","C1D8","C24","AMB-\nC25","C25"))

colnames(subset_B_cells_MCL@meta.data) %<>% gsub("proportion_scaled_cnv_chr","",.)
FeaturePlot.1(subset_B_cells_MCL,reduction="tsne",split.by = "orig.ident",
              strip.text.size = 7,pt.size = 0.01, 
              border = T,ncol = 5,width=6, height=11,do.print = T,do.return = F,
              features=c(1:22,"X","Y"),
              title = "Proportion scaled CNV")
write.csv(subset_B_cells_MCL@meta.data,paste0(path,"meta.data.csv"))
