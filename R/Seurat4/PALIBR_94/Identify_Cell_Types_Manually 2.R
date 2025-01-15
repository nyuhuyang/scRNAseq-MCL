library(Seurat)
library(dplyr)
library(tidyr)
library(magrittr)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

# ================== cell cycle ==================
object <- readRDS(file = "data/MCL_AIM_55_20240823.rds")
meta.data <- readRDS(file = "output/MCL_AIM_55_20240823_metadata_v2.rds")
if(all(colnames(object) == rownames(meta.data))){
    print("all cellID match!")
    object@meta.data = meta.data
}



# ================== copy old X4 cluster ==================
meta_data = readRDS("../scRNAseq-AIM/output/MCL_AIM_93_20220519_metadata_v3.rds")
meta.data$barcode = rownames(meta.data)
meta_data$barcode = rownames(meta_data)
meta.data %<>% dplyr::left_join(meta_data[,c("barcode",
                                           "X6cluster_MCL61",
                                           "X6cluster_AIM74",
                                           "X6cluster_MCL61.colors",
                                           "X6cluster_AIM74.colors"
                                           )], by = "barcode")
meta.data %<>% tibble::column_to_rownames("barcode")
for(cluster in c("X6cluster_MCL61","X6cluster_AIM74")){
    if(anyNA(meta.data[,cluster])){
        cellID = rownames(meta.data)[is.na(meta.data[,cluster])]
        meta.data[cellID,cluster] = meta.data[cellID,"cell.types"]
    }
}

for(cluster in c("X6cluster_MCL61.colors","X6cluster_AIM74.colors")){
    if(anyNA(meta.data[,cluster])){
        cellID = rownames(meta.data)[is.na(meta.data[,cluster])]
        meta.data[cellID,cluster] = meta.data[cellID,"cell.types.colors"]
    }
}
if(all(colnames(object) == rownames(meta.data))){
    print("all cellID match!")
    object@meta.data = meta.data
}

saveRDS(meta.data, file = "output/MCL_AIM_55_20240823_metadata_v3.rds")

object$X6cluster <- object$X6cluster_AIM74
meta.data <- object@meta.data

meta.data[meta.data$X6cluster %in% "MCL" & meta.data$cell.types != "MCL","X6cluster"] <-
    meta.data[meta.data$X6cluster %in% "MCL" & meta.data$cell.types != "MCL","cell.types"]

table(meta.data[meta.data$X6cluster %in% "MCL","SCT_snn_res.0.2"])

meta.data[meta.data$X6cluster %in% "MCL" & meta.data$SCT_snn_res.0.2 == 1 ,"X6cluster"] <- "C2"
meta.data[meta.data$X6cluster %in% "MCL" & meta.data$SCT_snn_res.0.2 == 2 ,"X6cluster"] <- "C6"
meta.data[meta.data$X6cluster %in% "MCL" & meta.data$SCT_snn_res.0.2 == 3 ,"X6cluster"] <- "NK cells"
meta.data[meta.data$X6cluster %in% "MCL" & meta.data$SCT_snn_res.0.2 == 4 ,"X6cluster"] <- "NK cells"
meta.data[meta.data$X6cluster %in% "MCL" & meta.data$SCT_snn_res.0.2 == 6 ,"X6cluster"] <- "C6"
meta.data[meta.data$X6cluster %in% "MCL" & meta.data$SCT_snn_res.0.2 == 7 ,"X6cluster"] <- "C4"
meta.data[meta.data$X6cluster %in% "MCL" & meta.data$SCT_snn_res.0.2 == 9 ,"X6cluster"] <- "other Myeloid cells"

meta.data[meta.data$X6cluster %in% "MCL" & meta.data$SCT_snn_res.0.2 == 10 ,"X6cluster"] <- "C1"
meta.data[meta.data$X6cluster %in% "MCL" & meta.data$SCT_snn_res.0.2 == 11 ,"X6cluster"] <- "C5"
meta.data[meta.data$X6cluster %in% "MCL" & meta.data$SCT_snn_res.0.2 == 12 ,"X6cluster"] <- "C1"
meta.data[meta.data$X6cluster %in% "MCL" & meta.data$SCT_snn_res.0.2 == 13 ,"X6cluster"] <- "T cells:CD8+"
meta.data[meta.data$X6cluster %in% "MCL" & meta.data$SCT_snn_res.0.2 == 14 ,"X6cluster"] <- "C1"
meta.data[meta.data$X6cluster %in% "MCL" & meta.data$SCT_snn_res.0.2 == 16 ,"X6cluster"] <- "C1"
meta.data[meta.data$X6cluster %in% "MCL" & meta.data$SCT_snn_res.0.2 == 17 ,"X6cluster"] <- "C1"
meta.data[meta.data$X6cluster %in% "MCL" & meta.data$SCT_snn_res.0.2 == 19 ,"X6cluster"] <- "C1"

meta.data1 <- meta.data[!duplicated(meta.data$X6cluster_AIM74),]
meta.data$X6cluster.colors <- plyr::mapvalues(meta.data$X6cluster,
                                              from = meta.data1$X6cluster_AIM74,
                                              to = meta.data1$X6cluster_AIM74.colors)
saveRDS(meta.data, file = "output/MCL_AIM_55_20240823_metadata_v4.rds")
