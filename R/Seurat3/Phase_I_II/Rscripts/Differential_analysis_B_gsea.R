########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
#SBATCH --mem=128G  # memory requested, units available: K,M,G,T

library(Seurat)
library(dplyr)
library(MAST)
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.numeric(slurm_arrayid)
print(paste0("slurm_arrayid=",args))

# load cvs
df_samples <- readxl::read_excel("doc/191030_scRNAseq_info.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
sample_n = which(df_samples$tests %in% c("control",paste0("test",2:12)))
df_samples = df_samples[sample_n,]

# load data
(load(file = "data/B_cells_MCL_43_20190917.Rda"))
B_cells_MCL@meta.data$orig.ident %<>% plyr::mapvalues(from = unique(df_samples$sample),
                                                      to = unique(df_samples$publication.id))
table(B_cells_MCL@meta.data$orig.ident)
B_cells_MCL$orig.ident %<>% as.character()
NewNames = paste0(B_cells_MCL@meta.data$orig.ident,"_",B_cells_MCL@meta.data$Barcode)
B_cells_MCL %<>% RenameCells(new.names = NewNames)
rownames(B_cells_MCL@reductions$tsne@cell.embeddings) = colnames(B_cells_MCL)

Idents(B_cells_MCL) = "groups"
B_cells_MCL %<>% subset(idents = c("AFT-03","AFT-04"),invert = T)
table(B_cells_MCL[["orig.ident"]])
# Differential anlaysis between 5 clusters
if(args == 1){
        Idents(B_cells_MCL) <- "X5_clusters"
        B_cells_MCL <- sortIdent(B_cells_MCL,numeric = T)
        table(Idents(B_cells_MCL))
        X5_clusters_markers <- FindAllMarkers.UMI(B_cells_MCL,
                                                  logfc.threshold = 0,only.pos = F, 
                                                  min.pct = 0.1,return.thresh = 1)
        write.csv(X5_clusters_markers,paste0(path,"X5_clusters_FC0_markers.csv"))
}

if(args == 2){
        B_cells_MCL@meta.data$X5_clusters_normal = as.numeric(as.character(B_cells_MCL@meta.data$X5_clusters))
        normal <- B_cells_MCL$orig.ident %in% c("N01","N02","N03","N04")
        B_cells_MCL@meta.data[normal,"X5_clusters_normal"] = "Normal"
        B_cells_MCL <- sortIdent(B_cells_MCL)
        Idents(B_cells_MCL) = "X5_clusters_normal"
        table(Idents(B_cells_MCL))
        X5_clusters_normal_markers <- FindPairMarkers(B_cells_MCL,ident.1 = 1:5, ident.2 = rep("Normal",5),
                                                      logfc.threshold = 0,only.pos = F,
                                                      min.pct = 0.1,return.thresh = 1)
        write.csv(X5_clusters_normal_markers,paste0(path,"X5_cluster_vs_Normal_FC0_markers.csv"))
}
