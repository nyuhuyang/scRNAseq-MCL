library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(magrittr)
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
marker_path <- paste0(path,"markers/")
if(!dir.exists(marker_path))dir.create(marker_path, recursive = T)

# ======== 2.1 =========== test with known markers==================
(load(file="data/MCL_41_20191205.Rda"))
DefaultAssay(object) <- "SCT"
features <- FilterGenes(object,c("CD19","CCND1","SOX11",
                                 "CD3D","CD4","CD8A",
                                 "MS4A7","CD14","FCGR1A",
                                 "GNLY","KLRC1","NCAM1"))
FeaturePlot.1(object,features = features, pt.size = 0.005, cols = c("gray90", "red"),
              alpha = 1,reduction = "tsne",
              threshold = 1, text.size = 20, border = T,do.print = T, do.return = F,ncol = 3, 
              units = "in",width=9, height=12, no.legend = T)
FeaturePlot.1(object,features = features, pt.size = 0.005, cols = c("gray90", "red"),
              alpha = 1,reduction = "umap",
              threshold = 1, text.size = 20, border = T,do.print = T, do.return = F,ncol = 3, 
              units = "in",width=9, height=12, no.legend = T)
