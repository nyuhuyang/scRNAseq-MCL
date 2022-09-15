library(Seurat)
library(dplyr)
library(tidyr)
library(magrittr)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

# ================== cell cycle ==================
object <- readRDS(file = "data/MCL_SCT_87_20220901.rds")
meta.data <- readRDS(file = "output/MCL_SCT_87_20220901_meta.data.rds")

if(all(colnames(object) == rownames(cluster_df))){
    print("all cellID match!")
    object@meta.data = meta.data
}

s.genes <- cc.genes$s.genes %>% gsub("MLF1IP","CENPU",.)
g2m.genes <- cc.genes$g2m.genes %>% plyr::mapvalues(from = c("FAM64A", "HN1"),
                                                    to = c("PIMREG","JPT1"))
object %<>% CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

colnames(object@meta.data) %<>% sub("Phase","cell cycle phase",.)
object$Mean.Reads.per.Cell %<>% gsub(",","",.) %>% as.integer()
object$Number.of.Reads %<>% gsub(",","",.) %>% as.integer()
object$Sequencing.Saturation %<>% gsub("%","",.) %>% as.numeric()
saveRDS(object@meta.data, file = "output/MCL_SCT_87_20220901_meta.data.rds")
# ================== copy old X4 cluster ==================
meta.data <- readRDS(file = "output/MCL_SCT_87_20220901_meta.data.rds")
meta_data = readRDS("../scRNAseq-AIM/output/MCL_AIM_93_20220519_metadata_v3.rds")
meta.data$barcode = rownames(meta.data)
meta_data$barcode = rownames(meta_data)
meta.data %<>% dplyr::left_join(meta_data[,c("barcode",
                                           "X6cluster_MCL61",
                                           "X6cluster_AIM74",
                                           #"X6cluster_MCL61.colors",
                                           #"X6cluster_AIM74.colors",
                                           "X9cluster"
                                           #"X9cluster.colors"
                                           )], by = "barcode")
meta.data %<>% tibble::column_to_rownames("barcode")
for(cluster in c("X6cluster_MCL61","X6cluster_AIM74","X9cluster")){
    if(anyNA(meta.data[,cluster])){
        cellID = rownames(meta.data)[is.na(meta.data[,cluster])]
        meta.data[cellID,cluster] = meta.data[cellID,"label1.blue_encode"]
    }
}


if(all(colnames(object) == rownames(cluster_df))){
    print("all cellID match!")
    object@meta.data = meta.data
}


meta_data1 = meta_data[!duplicated(meta_data$X6cluster_MCL61),]
meta_data2 = meta_data[!duplicated(meta_data$X6cluster_AIM74),]
meta_data3 = meta_data[!duplicated(meta_data$X9cluster),]

meta.data$X6cluster_MCL61.colors <- plyr::mapvalues(meta.data$X6cluster_MCL61,
                                                    from = meta_data1$X6cluster_MCL61,
                                                    to = meta_data1$X6cluster_MCL61.colors)
meta.data$X6cluster_AIM74.colors <- plyr::mapvalues(meta.data$X6cluster_AIM74,
                                                    from = meta_data2$X6cluster_AIM74,
                                                    to = meta_data2$X6cluster_AIM74.colors)
meta.data$X9cluster.colors <- plyr::mapvalues(meta.data$X9cluster,
                                                    from = meta_data3$X9cluster,
                                                    to = meta_data3$X9cluster.colors)
saveRDS(meta.data, file = "output/MCL_SCT_87_20220901_meta.data_v2.rds")
