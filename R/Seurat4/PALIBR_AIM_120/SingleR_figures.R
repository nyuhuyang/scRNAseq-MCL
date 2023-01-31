# conda activate r4.1.3
library(Seurat)
library(magrittr)
library(pheatmap)
library(dplyr)
library(tidyr)
library(ggpubr)
library(S4Vectors)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.1 load data ==========================================
object = readRDS(file = "data/MCL_120_SCT_20230106.rds")
object@meta.data %<>% cbind(FetchData(object,"CCND1"))
meta.data = object@meta.data

rm(object);GC()
##############################
# create singleR data frame
###############################
reference = c("MCL+blue_encode","MCL+azimuth_PBMC")[1]

if(reference == "MCL+blue_encode"){
    pred = readRDS("output/MCL_120_20230106_blueEncode_singleR_pred.rds")
    singlerDF = data.frame("label2.blue_encode" = pred$pruned.labels,
                           row.names = rownames(pred))
    table(rownames(pred) == rownames(meta.data))
    table(is.na(singlerDF$label2.blue_encode))
    singlerDF$label2.blue_encode[is.na(singlerDF$label2.blue_encode)]= "unknown"
    ##############################
    # adjust cell label
    ##############################
    singlerDF$label2.blue_encode %<>% gsub("Adipocytes|Fibroblasts|.*Endothelial cells|Keratinocytes",
                                          "Nonhematopoietic cells",.)
    singlerDF$label1.blue_encode = singlerDF$label2.blue_encode

    # combine cell types
    singlerDF[grep("CD4+",singlerDF$label1.blue_encode),"label1.blue_encode"] ="T cells:CD4+"
    singlerDF[grep("CD8+",singlerDF$label1.blue_encode),"label1.blue_encode"] ="T cells:CD8+"
    singlerDF[grep("B-cells",singlerDF$label1.blue_encode),"label1.blue_encode"] ="B cells"

    singlerDF$label1.blue_encode %<>% gsub("Tregs","T cells:regs",.)
    singlerDF$label1.blue_encode %<>% gsub("MEP|CLP|HSC|CMP|GMP|MPP","HSC,progenitors",.)
    singlerDF$label1.blue_encode %<>% gsub("DC|Macrophages|Macrophages M1","other Myeloid cells",.)
    singlerDF$label1.blue_encode %<>% gsub("Eosinophils|Megakaryocytes|Mesangial cells","other Myeloid cells",.)

    # reduce false positive results (B cells are labeled as MCL in normal samples)
    # and false negative results (MCL cells are labeled as B cells in MCL samples)
    # singler1sub false negative results  =========
    singlerDF$CCND1 = meta.data$CCND1
    singlerDF[(singlerDF$CCND1 >0 & singlerDF$label1.blue_encode %in% "B cells"),"label1.blue_encode"] = "MCL"
    singlerDF[(singlerDF$CCND1 >0 & grepl("B-cells",singlerDF$label2.blue_encode)),"label2.blue_encode"] = "MCL"
    # cell.types false positive results  ========
    #table(singlerDF$cell.types, meta.data$orig.ident) %>% kable %>% kable_styling()
    normal_cells <- meta.data$orig.ident %in% c("N01","N02","N03","N04") %>% rownames(singlerDF)[.]
    singlerDF[normal_cells,"label1.blue_encode"] %<>% gsub("MCL","B cells",.)
    singlerDF[normal_cells,"label2.blue_encode"] %<>% gsub("MCL","B cells",.)

    #singlerDF$CCND1 = NULL
    ##############################
    # process color scheme
    ##############################
    table(singlerDF$label1.blue_encode)
    singlerDF$label1.blue_encode.colors = singlerDF$label1.blue_encode
    singlerDF$label1.blue_encode.colors %<>% plyr::mapvalues(from = c("B cells","Erythrocytes","HSC,progenitors",
                                                                      "MCL","Monocytes",
                                                                      "NK cells","Nonhematopoietic cells","other Myeloid cells",
                                                                      "Plasma cells","T cells:CD4+","T cells:CD8+",
                                                                      "T cells:regs","unknown"),
                                                             to = c("#E6AB02", "#ff0000", "#6A3D9A",
                                                                    "#2055da", "#ADDFEE",
                                                                    "#A65628","#B3B3B3","#FDDAEC",
                                                                    "#1B9E77","#B3DE69","#F0027F",
                                                                    "#7570B3","#F2F2F2"))
    table(rownames(meta.data) == rownames(singlerDF))
    table(singlerDF$label1.blue_encode.colors)
    meta.data %<>% cbind(singlerDF[,c("label2.blue_encode","label1.blue_encode","label1.blue_encode.colors")])

}

if(reference == "MCL+azimuth_PBMC"){
    pred = readRDS("output/MCL_120_20230106_azimuth_MCL_singleR_pred.rds")
    singlerDF = data.frame("celltype.l3" = pred$pruned.labels,
                           row.names = rownames(pred))
    table(rownames(pred) == rownames(meta.data))
    table(is.na(singlerDF$celltype.l3))
    singlerDF$celltype.l3[is.na(singlerDF$celltype.l3)]= "unknown"
    FilePath = "../seurat_resources/azimuth/PBMC/"
    meta_data = read.csv(paste0(FilePath,"GSE164378/GSE164378_sc.meta.data_5P.csv"),row.names =1)
    meta_data = meta_data[!duplicated(meta_data$celltype.l3),]

    singlerDF$celltype.l2 = plyr::mapvalues(singlerDF$celltype.l3,
                                            from = meta_data$celltype.l3, to = meta_data$celltype.l2)
    singlerDF$celltype.l1 = plyr::mapvalues(singlerDF$celltype.l3,
                                            from = meta_data$celltype.l3, to = meta_data$celltype.l1)
    # reduce false positive results (B cells are labeled as MCL in normal samples)
    # and false negative results (MCL cells are labeled as B cells in MCL samples)
    # singler1sub false negative results  =========
    singlerDF$CCND1 = meta.data$CCND1
    MCL = singlerDF$CCND1 >0 & singlerDF$celltype.l1 %in% "B"
    singlerDF[MCL,"celltype.l1"] = "MCL"
    singlerDF[MCL,"celltype.l2"] = "MCL"
    singlerDF[MCL,"celltype.l3"] = "MCL"

    # cell.types false positive results  ========
    #table(singlerDF$cell.types, meta.data$orig.ident) %>% kable %>% kable_styling()
    normal_cells <- meta.data$orig.ident %in% c("N01","N02","N03","N04")
    singlerDF[normal_cells,"celltype.l1"] %<>% gsub("MCL","B",.)
    singlerDF[normal_cells,"celltype.l2"] %<>% gsub("MCL","B",.)
    singlerDF[normal_cells,"celltype.l3"] %<>% gsub("MCL","B",.)

    singlerDF$CCND1 = NULL
    ##############################
    # process color scheme
    ##############################
    singlerDF$celltype.l1.colors = singlerDF$celltype.l1
    singlerDF$celltype.l1.colors %<>% plyr::mapvalues(from = c("B",
                                                               "MCL",
                                                               "Mono",
                                                               "NK","other","DC",
                                                               "CD4 T","CD8 T",
                                                               "other T","unknown"),
                                                      to = c("#E6AB02",
                                                             "#2055da",#"#40A635","#FE8205","#8861AC","#E83C2D",
                                                             "#ADDFEE",
                                                             "#A65628","#B3B3B3","#FDDAEC",
                                                             "#B3DE69","#F0027F",
                                                             "#7570B3","#F2F2F2"))

    table(rownames(meta.data) == rownames(singlerDF))
    table(singlerDF$celltype.l1.colors)
    meta.data %<>% cbind(singlerDF)
}

if(reference == "MCL_51+azimuth_PBMC"){
    pred = readRDS("output/MCL_87_20220901_azimuth_MCL51_singleR_pred.rds")
    singlerDF = data.frame("celltype.l3" = pred$pruned.labels,
                           row.names = rownames(pred))
    table(rownames(pred) == rownames(meta.data))
    table(is.na(singlerDF$celltype.l3))
    singlerDF$celltype.l3[is.na(singlerDF$celltype.l3)]= "unknown"
    FilePath = "../seurat_resources/azimuth/PBMC/"
    meta_data = read.csv(paste0(FilePath,"GSE164378/GSE164378_sc.meta.data_5P.csv"),row.names =1)
    meta_data = meta_data[!duplicated(meta_data$celltype.l3),]

    singlerDF$celltype.l2 = plyr::mapvalues(singlerDF$celltype.l3,
                                            from = meta_data$celltype.l3, to = meta_data$celltype.l2)
    singlerDF$celltype.l1 = plyr::mapvalues(singlerDF$celltype.l3,
                                            from = meta_data$celltype.l3, to = meta_data$celltype.l1)
    # reduce false positive results (B cells are labeled as MCL in normal samples)
    # and false negative results (MCL cells are labeled as B cells in MCL samples)
    # singler1sub false negative results  =========
    singlerDF$CCND1 = meta.data$CCND1
    MCL = singlerDF$CCND1 >0 & singlerDF$celltype.l1 %in% "B"
    singlerDF[MCL,"celltype.l1"] = "MCL_0"
    singlerDF[MCL,"celltype.l2"] = "MCL_0"
    singlerDF[MCL,"celltype.l3"] = "MCL_0"

    # cell.types false positive results  ========
    #table(singlerDF$cell.types, meta.data$orig.ident) %>% kable %>% kable_styling()
    normal_cells <- meta.data$orig.ident %in% c("N01","N02","N03","N04")
    singlerDF[normal_cells,"celltype.l1"] %<>% gsub("MCL_.*","B",.)
    singlerDF[normal_cells,"celltype.l2"] %<>% gsub("MCL_.*","B",.)
    singlerDF[normal_cells,"celltype.l3"] %<>% gsub("MCL_.*","B",.)

    singlerDF$CCND1 = NULL
    singlerDF$celltype = singlerDF$celltype.l2
    singlerDF$celltype %<>% gsub("MCL_.*","MCL",.)
    singlerDF$celltype %<>% gsub("B|B intermediate|B memory|B naive","B",.)
    singlerDF$celltype %<>% gsub("CD4.*","CD4 T",.)
    singlerDF$celltype %<>% gsub("CD8.*","CD8 T",.)
    singlerDF$celltype %<>% gsub("gdT|MAIT","other T",.)
    singlerDF$celltype %<>% gsub("cDC1|cDC2|pDC","DC",.)
    singlerDF$celltype %<>% gsub("NK.*","NK",.)


    ##############################
    # process color scheme
    ##############################
    singlerDF$celltype.l1.colors = singlerDF$celltype.l1
    singlerDF$celltype.l1.colors %<>% plyr::mapvalues(from = c("B","Mono",
                                                               "MCL_0", "MCL_C1",  "MCL_C2",  "MCL_C3",  "MCL_C4",
                                                               "NK","other","DC",
                                                               "CD4 T","CD8 T",
                                                               "other T","unknown"),
                                                      to = c("#E6AB02","#ADDFEE",
                                                             "#2055da","#40A635","#FE8205","#8861AC","#E83C2D",
                                                             "#A65628","#DDD399","#FDDAEC",
                                                             "#B3DE69","#F0027F",
                                                             "#7570B3","#F2F2F2"))

    singlerDF$celltype.colors <- singlerDF$celltype
    color_df <- data.frame("celltype" = c("B","CD14 Mono","CD16 Mono",
                                          "CD4 T","CD8 T","DC","Eryth","HSPC",
                                          "MCL",
                                          "NK","other T","Plasmablast",
                                          "Platelet","Proliferating","Treg","unknown"),
                           "celltype.colors" = c("#E6AB02","#ADDFEE","#FDDAEC",
                                                 "#B3DE69","#F0027F","#EA4A34","#FDB762","#6A3D9A",
                                                 "#2055da",
                                                 "#A65628","#ED8F47","#1B9E77",
                                                 "#DDD399","#BBA0CD","#7570B3","#F2F2F2"))
    singlerDF$celltype.colors %<>% plyr::mapvalues(from = color_df$celltype,
                                                      to = color_df$celltype.colors)
    table(rownames(meta.data) == rownames(singlerDF))
    table(singlerDF$celltype.colors)
    table(singlerDF$celltype.l1.colors)

    meta.data <- meta.data[,-which(colnames(meta.data) %in% colnames(singlerDF))]
    meta.data %<>% cbind(singlerDF)
    meta.data$CCND1 = NULL
}
saveRDS(object@meta.data, file = "output/MCL_120_SCT_20230106_meta.data_v1.rds")
