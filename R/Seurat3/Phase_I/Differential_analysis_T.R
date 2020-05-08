# ######################################################################
library(Seurat)
library(dplyr)
library(kableExtra)
library(magrittr)
library(ggplot2)
library(cowplot)
library(fgsea)
library(tibble)
library(ggsci)
library(fgsea)
library(openxlsx)
source("../R/Seurat3_functions.R")

(load(file="data/MCL_41_harmony_20200225.Rda"))
Idents(object) = "Doublets"
object %<>% subset(idents = "Singlet")

DefaultAssay(object)  = "SCT"
object$orig.ident_cell.types = paste0(object$orig.ident,"_", object$cell.types)
Idents(object)= "orig.ident_cell.types"
ident.1 = "Pt25_24_T_cells:CD8+"; ident.2 = "Pt11_28_T_cells:CD8+"
sub_object <- subset(object, idents = c(ident.1, ident.2))
table(Idents(sub_object))

system.time(T_markers <- FindAllMarkers.UMI(sub_object, 
                                            logfc.threshold = 0.25, 
                                            only.pos = T,
                                            test.use = "MAST",
                                            return.thresh = 0.05, 
                                            latent.vars = "nCount_SCT"))
write.csv(T_markers,paste0(path,"T_41-FC0.25_",ident.1,"_",ident.2,".csv"))

VlnPlot