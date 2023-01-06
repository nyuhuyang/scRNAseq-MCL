library(Seurat)
library(dplyr)
library(tidyr)
library(magrittr)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")

object <- readRDS(file = "data/MCL_SCT_87_20220901.rds")
meta.data <- readRDS(file = "output/MCL_SCT_87_20220901_meta.data_v2.rds")
if(all(colnames(object) == rownames(meta.data))){
    print("all cellID match!")
    object@meta.data = meta.data
}
object@reductions <- readRDS(file = "output/20220915/3000/reductions_npcs100_dist.0.4_spread.1.rds")
object@meta.data %<>% cbind(object[["umap"]]@cell.embeddings)
object %<>% subset(UMAP_1 < 1 & UMAP_2 < 3 & label1.blue_encode %in% c("B cells","MCL"))

object@meta.data[,"SCT_snn_res.0.8"] = NULL

object %<>% FindNeighbors(reduction = "harmony",dims = 1:90)
resolutions = c(seq(0.01,0.09, by = 0.01),seq(0.1,0.9, by = 0.1),seq(1,5, by = 1))
for(i in 1:length(resolutions)){
    object %<>% FindClusters(resolution = resolutions[i], algorithm = 1,verbose = F)
    Progress(i, length(resolutions))
}
for(col in grep("SCT_snn_res",colnames(object@meta.data),value =T)){
    object@meta.data[,col] %<>% factor(levels = 0:max(as.integer(as.character(object@meta.data[,col]))))
}

saveRDS(object, file = "data/MCL_B_only_SCT_87_20220901.rds")
saveRDS(object@meta.data, file = "output/MCL_B_only_SCT_87_20220901_meta.data_v2.rds")
