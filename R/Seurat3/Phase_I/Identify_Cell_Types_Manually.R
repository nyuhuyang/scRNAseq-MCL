library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(magrittr)
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

# ======== 2.1 =========== test with known markers==================
(load(file="data/MCL_41_harmony_20191231.Rda"))
DefaultAssay(object) <- "SCT"
# global
features <- FilterGenes(object,c("CD19","CCND1","SOX11",
                                 "CD3D","CD4","CD8A",
                                 "MS4A7","CD14","FCGR1A",
                                 "GNLY","KLRC1","NCAM1"))
# B and MCL
features <- c(FilterGenes(object,c("CD19","CCND1","SOX11",
                                   "IRF4","NFKB1","CD274",
                                   "PCNA","CDK1","MKI67")))
FeaturePlot.1(object,features = features, pt.size = 0.005, 
              cols = c("lightgrey", "red"),
              alpha = 1,reduction = "tsne",
              unique.name = "cell.types",label = F,
              threshold = 1, text.size = 20, border = T,do.print = T, do.return = F,ncol = 3, 
              units = "in",width=9, height=12, no.legend = T)
QC <- c("percent.mt","nCount_SCT","nFeature_SCT")
FeaturePlot.1(object,features = QC, pt.size = 0.005, 
              cols = c("lightgrey", "red"),
              alpha = 1,reduction = "tsne",
              unique.name = "cell.types",label = F,
              threshold = 1, text.size = 20, border = T,do.print = T, do.return = F,ncol = 3, 
              units = "in",width=9, height=4, no.legend = T)
