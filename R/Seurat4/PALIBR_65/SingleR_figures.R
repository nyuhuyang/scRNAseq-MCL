# conda activate r4.1.1
library(Seurat)
library(magrittr)
library(pheatmap)
library(kableExtra)
library(dplyr)
library(tidyr)
library(ggpubr)
library(S4Vectors)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.2 SingleR specifications ==========================================

##############################
# create singleR data frame
###############################
pred = readRDS("output/MCL_SCT_65_20220411_singleR_pred.rds")
object = readRDS(file = "data/MCL_SCT_65_20220411.rds")

singlerDF = data.frame("label.fine" = pred$pruned.labels,
                       row.names = rownames(pred))
table(rownames(pred) == rownames(object@meta.data))
table(is.na(singlerDF$label.fine))
singlerDF$label.fine[is.na(singlerDF$label.fine)]= "unknown"

##############################
# adjust cell label
##############################
# combine cell types
singlerDF$cell.types = singlerDF$label.fine
grep("B-cells|B cells, PB",singlerDF$cell.types, value = T) %>% table
singlerDF[grep("B-cells|B cells",singlerDF$cell.types),"cell.types"] ="B_cells"
singlerDF[grep("CD4+",singlerDF$cell.types),"cell.types"] ="T_cells:CD4+"
singlerDF[grep("CD8+",singlerDF$cell.types),"cell.types"] ="T_cells:CD8+"
singlerDF$cell.types %<>% gsub("Tregs","T_cells:regs",.)
singlerDF$cell.types %<>% gsub("MEP|CLP|HSC|CMP|GMP|MPP","HSC/progenitors",.)
singlerDF$cell.types %<>% gsub("DC|Macrophages|Macrophages M1","other Myeloid cells",.)
singlerDF$cell.types %<>% gsub("Eosinophils|Megakaryocytes","other Myeloid cells",.)
singlerDF$cell.types %<>% gsub("Adipocytes|Fibroblasts|mv Endothelial cells|Endothelial cells|Keratinocytes|Melanocytes|Mesangial cells","Nonhematopoietic cells",.)

# reduce false positive results (B cells are labeled as MCL in normal samples)
# and false negative results (MCL cells are labeled as B cells in MCL samples)
# singler1sub false negative results  =========
CCND1 = FetchData(object,"CCND1")
singlerDF$CCND1 = CCND1$CCND1
singlerDF[(singlerDF$CCND1 >0 & singlerDF$cell.types %in% "B_cells"),"cell.types"] = "MCL"
# cell.types false positive results  ========
#table(singlerDF$cell.types, object@meta.data$orig.ident) %>% kable %>% kable_styling()
normal_cells <- object$orig.ident %in% c("N01","N02","N03","N04") %>% rownames(singlerDF)[.]
singlerDF[normal_cells,"cell.types"] %<>% gsub("MCL","B_cells",.)
singlerDF$CCND1 = NULL

##############################
# process color scheme
##############################
singlerDF$cell.types.colors = singlerDF$cell.types
singlerDF$cell.types.colors %<>% plyr::mapvalues(from = c("B_cells","Erythrocytes","HSC/progenitors",
                                                          "Macrophages","MCL","Monocytes",
                                                          "NK cells","Nonhematopoietic cells","other Myeloid cells",
                                                          "Plasma cells","T_cells:CD4+","T_cells:CD8+",
                                                          "T_cells:regs","unknown"),
                                                 to = c("#E6AB02", "#ff0000", "#6A3D9A",
                                                        "#FB9A99", "#2055da", "#ADDFEE",
                                                        "#A65628","#B3B3B3","#FDDAEC",
                                                        "#1B9E77","#B3DE69","#F0027F",
                                                        "#7570B3","#F2F2F2"))
table(colnames(object) == rownames(singlerDF))
table(singlerDF$cell.types.colors)
object@meta.data %<>% cbind(singlerDF)

saveRDS(object@meta.data, file = "data/MCL_SCT_65_20220411_metadata.rds")
#==============================================
meta.data = object@meta.data
meta.data$barcode  = rownames(meta.data)

# find X4cluster
meta_data = read.csv("shinyApp/PALIBR_I_51/cnv_meta_data.csv",row.names = 1)
meta_data$barcode  = rownames(meta_data)

meta.data %<>% left_join(meta_data[,c("X4cluster","barcode")], by = "barcode")
meta.data %<>% tibble::column_to_rownames(var = "barcode")

table(colnames(object) == rownames(meta.data))
object@meta.data = meta.data



saveRDS(object@meta.data, file = "data/MCL_SCT_65_20220411_metadata.rds")
