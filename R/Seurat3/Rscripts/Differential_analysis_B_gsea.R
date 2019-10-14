########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
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

# load data
(load(file = "data/B_cells_MCL_43_20190917.Rda"))

# Differential anlaysis between 5 clusters
if(args == 1){
        Idents(B_cells_MCL) <- "X5_clusters"
        B_cells_MCL <- sortIdent(B_cells_MCL,numeric = T)
        table(Idents(B_cells_MCL))
        X5_clusters_markers <- FindAllMarkers.UMI(B_cells_MCL,logfc.threshold = 0.01,
                                                  test.use = "MAST",
                                                  only.pos = FALSE, 
                                                  min.pct = 0.01,return.thresh = 1)
        write.csv(X5_clusters_markers,paste0(path,"X5_clusters_FC0.01_markers.csv"))
}

if(args == 2){
        B_cells_MCL@meta.data$X5_clusters_normal = as.numeric(as.character(B_cells_MCL@meta.data$X5_clusters))
        normal <- B_cells_MCL$orig.ident %in% c("BH","DJ","MD","NZ")
        B_cells_MCL@meta.data[normal,"X5_clusters_normal"] = "Normal"
        B_cells_MCL <- sortIdent(B_cells_MCL,numeric = T)
        Idents(B_cells_MCL) = "X5_clusters_normal"
        table(Idents(B_cells_MCL))
        X5_clusters_normal_markers <- FindPairMarkers(B_cells_MCL,ident.1 = 1:5, ident.2 = rep("Normal",5),
                                                      logfc.threshold = 0.01,only.pos = FALSE, 
                                                      min.pct = 0.01,return.thresh = 1,save.path = path)
        
        write.csv(X5_clusters_normal_markers,paste0(path,"X5_clusters_normal_FC0.01_markers.csv"))
}
