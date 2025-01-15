library(Seurat)
library(dplyr)
library(tidyr)
library(magrittr)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
#path <- paste0("output/",gsub("-","",Sys.Date()),"/")
#if(!dir.exists(path))dir.create(path, recursive = T)

# ================== cell cycle ==================
object <- readRDS(file = "data/MCL_SCT_94_20240615.rds")
meta.data <- readRDS(file = "output/MCL_SCT_94_20240615_meta.data_v1.rds")

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

meta.data$response %<>% factor(levels = c("Normal","Untreated","Baseline","CR","PR","PD"))
meta.data$X4cluster %<>% plyr::mapvalues(from = c( "1","2","3","4","B_cells","Erythrocytes",
                                                   "HSC/progenitors","MCL","Monocytes:CD14+",
                                                   "Monocytes:CD16+","NK cells","Nonhematopoietic cells",
                                                   "other Myeloid cells","Plasma cells","T_cells:CD4+",
                                                   "T_cells:CD8+","T_cells:regs","unknown"),
                                         to = c( "C1","C2","C3","C4","B cells","Erythrocytes",
                                                   "HSC,progenitors","MCL","Monocytes:CD14+",
                                                   "Monocytes:CD16+","NK cells","Nonhematopoietic cells",
                                                   "other Myeloid cells","Plasma cells","T cells:CD4+",
                                                   "T cells:CD8+","T cells:regs","unknown"))
meta_data1 = meta_data[!duplicated(meta_data$X6cluster_MCL61),]
meta.data$X4cluster.colors <- plyr::mapvalues(meta.data$X4cluster,
                                              from = meta_data1$X6cluster_MCL61,
                                              to = meta_data1$X6cluster_MCL61.colors)
meta.data$X4cluster.colors %<>% gsub("Monocytes:CD14+","#FDDAEC",.)
meta.data$X4cluster.colors %<>% gsub("Monocytes:CD16+","#ADDFEE",.)
meta.data$X4cluster[is.na(meta.data$X4cluster)] ="unknown"
meta.data$X4cluster.colors[is.na(meta.data$X4cluster.colors)] ="#B3B3B3"
saveRDS(meta.data, file = "output/MCL_SCT_94_20240615_meta.data_v2.rds")
#===========
meta.data <- readRDS(file = "output/MCL_SCT_94_20240615_meta.data_v2.rds")
meta_data = readRDS("../scRNAseq-AIM/output/MCL_AIM_55_20240823_metadata_v4.rds")

meta.data$barcode = rownames(meta.data)
meta_data$barcode = rownames(meta_data)
meta.data %<>% dplyr::left_join(meta_data[,c("barcode",
                                             "X6cluster",
                                             "X6cluster.colors"
)], by = "barcode")
meta.data %<>% tibble::column_to_rownames("barcode")
for(cluster in c("X6cluster")){
    if(anyNA(meta.data[,cluster])){
        cellID = rownames(meta.data)[is.na(meta.data[,cluster])]
        meta.data[cellID,cluster] = meta.data[cellID,"cell.types"]
    }
}

for(cluster in c("X6cluster.colors")){
    if(anyNA(meta.data[,cluster])){
        cellID = rownames(meta.data)[is.na(meta.data[,cluster])]
        meta.data[cellID,cluster] = meta.data[cellID,"cell.types.colors"]
    }
}
object <- readRDS(file = "data/MCL_SCT_94_20240615.rds")
if(all(colnames(object) == rownames(meta.data))){
    print("all cellID match!")
    object@meta.data = meta.data
}

resolutions = c(seq(0.02,0.08, by = 0.02),seq(0.2,0.8, by = 0.2))
for(i in 1:length(resolutions)){
    object %<>% FindClusters(resolution = resolutions[i], algorithm = 1,verbose = F)
    Progress(i, length(resolutions))
}
saveRDS(object@meta.data, file = "output/MCL_SCT_94_20240615_meta.data_v3.rds")
# ================== manually assign ==================
meta.data <- readRDS(file = "output/MCL_SCT_94_20240615_meta.data_v3.rds")

meta.data$X7cluster <- meta.data$X6cluster

cellID_known = rownames(meta.data)[!meta.data[,"X6cluster_MCL61"] %in% c("MCL","unknown")]
meta.data[cellID_known,"X7cluster"] = meta.data[cellID_known,"X6cluster_MCL61"]
meta.data[cellID_known,"X6cluster"] = meta.data[cellID_known,"X6cluster_MCL61"]

MCL <- grepl("MCL",meta.data$X7cluster)
all_MCL <- grepl("MCL|C1|C2|C3|C4|C5|C6",meta.data$X7cluster)
meta.data[all_MCL & grepl("B",meta.data$celltype.l1) ,"X7cluster"] ="B cells"
meta.data[all_MCL & grepl("CD4 T",meta.data$celltype.l1) ,"X7cluster"] ="T cells:CD4+"
meta.data[all_MCL & grepl("CD8 T",meta.data$celltype.l1) ,"X7cluster"] ="T cells:CD8+"
meta.data[all_MCL & grepl("DC",meta.data$celltype.l1) ,"X7cluster"] ="other Myeloid cells"
meta.data[all_MCL & grepl("Mono",meta.data$celltype.l1) ,"X7cluster"] ="Monocytes"
meta.data[all_MCL & grepl("NK",meta.data$celltype.l1) ,"X7cluster"] ="NK cells"
meta.data[all_MCL & grepl("other",meta.data$celltype.l1) ,"X7cluster"] ="C4"
meta.data[all_MCL & grepl("unknown",meta.data$celltype.l1) ,"X7cluster"] ="unknown"

meta.data[MCL & grepl("0",meta.data$SCT_snn_res.0.6) ,"X7cluster"] ="C1"
meta.data[MCL & grepl("2",meta.data$SCT_snn_res.0.6) ,"X7cluster"] ="C1"
meta.data[MCL & grepl("10",meta.data$SCT_snn_res.0.6) ,"X7cluster"] ="C1"
meta.data[MCL & grepl("24",meta.data$SCT_snn_res.0.6) ,"X7cluster"] ="C1"
meta.data[MCL & grepl("34",meta.data$SCT_snn_res.0.6) ,"X7cluster"] ="C1"
meta.data[all_MCL & grepl("35",meta.data$SCT_snn_res.0.6) ,"X7cluster"] ="B cells"
meta.data[MCL & grepl("44",meta.data$SCT_snn_res.0.6) ,"X7cluster"] ="C1"

meta.data[all_MCL & grepl("3",meta.data$SCT_snn_res.0.6) ,"X7cluster"] ="unknown"
meta.data[all_MCL & grepl("13",meta.data$SCT_snn_res.0.6) ,"X7cluster"] ="unknown"
meta.data[all_MCL & grepl("14",meta.data$SCT_snn_res.0.6) ,"X7cluster"] ="Monocytes"
meta.data[all_MCL & grepl("17",meta.data$SCT_snn_res.0.6) ,"X7cluster"] ="unknown"
meta.data[all_MCL & grepl("18",meta.data$SCT_snn_res.0.6) ,"X7cluster"] ="unknown"
meta.data[all_MCL & grepl("20",meta.data$SCT_snn_res.0.6) ,"X7cluster"] ="unknown"
meta.data[all_MCL & grepl("22",meta.data$SCT_snn_res.0.6) ,"X7cluster"] ="unknown"

meta.data[MCL & grepl("6",meta.data$SCT_snn_res.0.6) ,"X7cluster"] ="C2"
meta.data[MCL & grepl("19",meta.data$SCT_snn_res.0.6) ,"X7cluster"] ="C2"
meta.data[MCL & grepl("27",meta.data$SCT_snn_res.0.6) ,"X7cluster"] ="C2"
meta.data[MCL & grepl("33",meta.data$SCT_snn_res.0.6) ,"X7cluster"] ="C2"
meta.data[MCL & grepl("39",meta.data$SCT_snn_res.0.6) ,"X7cluster"] ="C2"
meta.data[MCL & grepl("43",meta.data$SCT_snn_res.0.6) ,"X7cluster"] ="C2"


meta.data[MCL & grepl("7",meta.data$SCT_snn_res.0.6) ,"X7cluster"] ="C3"
meta.data[MCL & grepl("8",meta.data$SCT_snn_res.0.6) ,"X7cluster"] ="C4"
meta.data[MCL & grepl("9",meta.data$SCT_snn_res.0.6) ,"X7cluster"] ="C4"
meta.data[grepl("11",meta.data$SCT_snn_res.0.6) ,"X7cluster"] ="C5"
meta.data[MCL & grepl("16",meta.data$SCT_snn_res.0.6) ,"X7cluster"] ="C6"
meta.data[MCL & grepl("47",meta.data$SCT_snn_res.0.6) ,"X7cluster"] ="C6"


meta.data[MCL & grepl("21",meta.data$SCT_snn_res.0.6) & grepl("Pt02",meta.data$patient),"X7cluster"] ="C6-Pt02"
meta.data[MCL & grepl("37",meta.data$SCT_snn_res.0.6) & grepl("Pt02",meta.data$patient),"X7cluster"] ="C6-Pt02"


for(cluster in c("X4cluster","X6cluster_MCL61","X6cluster_AIM74","X6cluster","X7cluster")){
    print("=================")
    print(as.data.frame(table(meta.data[,cluster])))
}

#meta.data$X6cluster  %<>% gsub("C5","C5-rm",.)
#meta.data$X6cluster  %<>% gsub("C6","C5",.)
#meta.data$X6cluster.colors[meta.data$X6cluster == "C5"] = "#2055da"


#meta.data1 = meta.data[!duplicated(meta.data$X6cluster),]
#meta.data1 = meta.data1[order(meta.data1$X6cluster),]

val <- c("#E6AB02","#40A635","#FE8205","#8861AC","#E83C2D","#FF80BE","#2055da","#FF80BE","#799203",
          "#ff0000", "#6A3D9A",
         "#2055da", "#ADDFEE","#ADDFEE","#FB9A99",
         "#A65628","#B3B3B3","#FDDAEC",
         "#1B9E77","#B3DE69","#F0027F",
         "#7570B3","#808080","#F2F2F2")
names(val) <- c("B cells",paste0("C",1:6),"C5-rm","C6-Pt02",
                "Erythrocytes","HSC,progenitors",
                "MCL","Monocytes","Monocytes:CD14+","Monocytes:CD16+",
                "NK cells","Nonhematopoietic cells","other Myeloid cells",
                "Plasma cells","T cells:CD4+","T cells:CD8+",
                "T cells:regs","unassign","unknown")

for(cluster in c("X4cluster","X6cluster_MCL61","X6cluster_AIM74","X6cluster","X7cluster")){
    meta.data[,paste0(cluster,".colors")]  <- plyr::mapvalues(meta.data[,cluster],
                                                              from = sort(unique(meta.data[,cluster])),
                                                              to = val[sort(unique(meta.data[,cluster]))])
}

for(cluster in c("X4cluster","X6cluster_MCL61","X6cluster_AIM74","X6cluster","X7cluster")){
    print("=================")
    print(as.data.frame(table(meta.data[,paste0(cluster,"")])))
}


saveRDS(meta.data, file = "output/MCL_SCT_94_20240615_meta.data_v4.rds")




