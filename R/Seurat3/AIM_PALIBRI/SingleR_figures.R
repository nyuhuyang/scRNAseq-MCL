# conda activate r4.0
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
pred = readRDS("output/20210119_MCL_AIM_58_singleR_pred.rds")
(load(file = "data/MCL_AIM_58_20201009.Rda"))

singlerDF = data.frame("label.fine" = pred$pruned.labels)
singlerDF$label.fine %<>% plyr::mapvalues(from = c("B cells, PB",
                                                   "CD4+ T-cells",
                                                   "CD4+ Tcm",
                                                   "CD4+ Tem",
                                                   "CD8+ T-cells",
                                                   "CD8+ Tcm",
                                                   "CD8+ Tem",
                                                   "Class-switched memory B-cells",
                                                   "Memory B-cells",
                                                   "mv Endothelial cells",
                                                   "naive B-cells",
                                                   "NK cells",
                                                   "Plasma cells",
                                                   "Tregs"),
                                          to = c("B_cells, PB",
                                                 "T_cells CD4+",
                                                 "T_cells CD4+, cm",
                                                 "T_cells CD4+, em",
                                                 "T_cells CD8+",
                                                 "T_cells CD8+, cm",
                                                 "T_cells CD8+, em",
                                                 "B_cells, Class-switched memory",
                                                 "B_cells, Memory",
                                                 "Endothelial cells, mv",
                                                 "B_cells, naive",
                                                 "NK_cells",
                                                 "Plasma_cells",
                                                 "T_cells CD4+, regs"))
singlerDF$label.main = gsub(", .*$","", singlerDF$label.fine)

# reduce false positive results (B cells are labeled as MCL in normal samples)
# and false negative results (MCL cells are labeled as B cells in MCL samples)
# label.main false negative results  =========
CCND1 = FetchData(object,"CCND1")
singlerDF$CCND1 = CCND1$CCND1
singlerDF[(singlerDF$CCND1 >0 & singlerDF$label.main %in% "B_cells"),"label.main"] = "MCL"
singlerDF[(singlerDF$CCND1 >0 & singlerDF$label.fine %in% c("B_cells, Memory",
                                                            "B_cells, naive",
                                                            "B_cells, PB")),"label.fine"] = "MCL"

# cell.types false positive results  ========
rownames(singlerDF) = colnames(object)
normal_cells <- object$orig.ident %in% paste0("N0",1:4) %>% rownames(singlerDF)[.]
singlerDF[normal_cells,"label.main"] %<>% gsub("MCL","B_cells",.)
singlerDF[normal_cells,"label.fine"] %<>% gsub("MCL","B_cells, PB",.)

##############################
# adjust cell label
##############################
singlerDF$cell.types = gsub("MEP|CLP|HSC|CMP|GMP|MPP","HSC/progenitors",singlerDF$label.main)
singlerDF$cell.types = gsub("DC|Macrophages|Erythrocytes|Eosinophils|Megakaryocytes|Neutrophils|Monocytes","Myeloid_cells",singlerDF$cell.types)
singlerDF$cell.types = gsub("Adipocytes|Fibroblasts|Endothelial cells","Nonhematopoietic_cells",singlerDF$cell.types)

##############################
# process color scheme
##############################
table(colnames(object) == rownames(singlerDF))
object@meta.data %<>% cbind(singlerDF[,c("label.main","label.fine","cell.types")])

df_samples <- readxl::read_excel("doc/singler.colors.xlsx")
object$cell.types.colors = plyr::mapvalues(object$cell.types,
                                           from = na.omit(df_samples$Cell.types),
                                           to = na.omit(df_samples$singler.color2))

lapply(c("label.fine","label.main","cell.types"), function(group.by)
    TSNEPlot.1(object = object, label = T, label.repel = T,group.by = group.by,
        no.legend = T,
        pt.size = 0.1,label.size = ifelse(group.by=="label.fine",3,5),
        do.print = T,do.return = F,
        title = paste(group.by, "labeling by Blueprint + Encode + MCL")))

save(object, file = "data/MCL_AIM_58_20201009.Rda")
