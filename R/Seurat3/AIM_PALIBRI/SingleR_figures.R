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
pred = readRDS("output/20201008_MCL_AIM_27_singleR_pred.rds")
(load(file = "data/MCL_AIM_27_20201008.Rda"))

singlerDF = data.frame("label.fine" = pred$pruned.labels,
                       "label.main" = gsub(",.*","",pred$pruned.labels))

# reduce false positive results (B cells are labeled as MCL in normal samples)
# and false negative results (MCL cells are labeled as B cells in MCL samples)
# singler1sub false negative results  =========
CCND1 = FetchData(object,"CCND1")
singlerDF$CCND1 = CCND1$CCND1
singlerDF[(singlerDF$CCND1 >0 & singlerDF$main.labels %in% "B cells"),"label.main"] = "MCL"
singlerDF[(singlerDF$CCND1 >0 & singlerDF$main.labels %in% "B cells"),"label.fine"] = "MCL"

# cell.types false positive results  ========
rownames(singlerDF) = colnames(object)
normal_cells <- object$orig.ident %in% c("BH","DJ","MD","NZ") %>% rownames(singlerDF)[.]
singlerDF[normal_cells,"label.main"] %<>% gsub("MCL","B cells",.)
singlerDF[normal_cells,"label.fine"] %<>% gsub("MCL","B cells, PB",.)

table(colnames(object) == rownames(singlerDF))
object@meta.data %<>% cbind(singlerDF[,c("label.main","label.fine")])
##############################
# process color scheme
##############################
singler_colors2 = c("#E6AB02","#6A3D9A", "#2055da", "#FB9A99", "#A65628", "#B3B3B3", "#B3DE69", "#F0027F")
object <- AddMetaData(object = object,metadata = singlerDF["label.main"])
object <- AddMetaColor(object = object, label= "label.main", colors = Singler.colors)

object <- AddMetaData(object = object,metadata = singlerDF["label.fine"])
object <- AddMetaColor(object = object, label= "label.fine", colors = Singler.colors)

lapply(c("label.fine","label.main"), function(group.by)
    TSNEPlot.1(object = object, label = T, label.repel = T,group.by = group.by,
        cols = Singler.colors,no.legend = T,
        pt.size = 0.1,label.size = ifelse(group.by=="label.fine",3,5),
        do.print = T,do.return = F,
        title = "Cell type labeling by Blueprint + Encode + MCL"))

lapply(c("label.fine","label.main"), function(group.by)
    UMAPPlot.1(object = object, label = T, label.repel = T,group.by = group.by,
               cols = Singler.colors,no.legend = T,
               pt.size = 0.1,label.size = ifelse(group.by=="label.fine",3,5),
               do.print = T,do.return = F,
               title = "Cell type labeling by Blueprint + Encode + MCL"))
