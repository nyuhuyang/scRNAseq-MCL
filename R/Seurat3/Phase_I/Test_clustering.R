########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
#
# ######################################################################
invisible(lapply(c("Seurat","dplyr","cowplot","fgsea",
                   "magrittr","data.table","future","ggplot2","tidyr"), function(x) {
                       suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

# Need 64GB
# load files
(load(file = "data/MCL_41_harmony_20200225.Rda"))
Idents(object) = "Doublets"
object %<>% subset(idents = "Singlet")

object@meta.data$SCT_snn_res.0.8 %<>% as.character %>% as.integer
lapply(c(TSNEPlot.1, UMAPPlot.1), function(fun)
    fun(object, group.by="SCT_snn_res.0.8",pt.size = 1,
     label = T, label.repel = T,alpha = 0.9,cols = Singler.colors,
     no.legend = T,label.size = 4, repel = T,
     title = paste(formals(fun)$reduction, "plot for all clusters of PALIBR phase I"),
     do.print = T, do.return = F))

lapply(c(TSNEPlot.1, UMAPPlot.1), function(fun)
    fun(object, group.by="cell.types",pt.size = 1,
        label = T, label.repel = T,alpha = 0.9,#cols = Singler.colors,
        no.legend = F,label.size = 4, repel = T,
        title = paste(formals(fun)$reduction, "plot for cell types of PALIBR phase I"),
        do.print = T, do.return = F))

Idents(object) = "SCT_snn_res.0.8"
c_11 <- subset(object, idents = 11)

lapply(c(TSNEPlot.1, UMAPPlot.1), function(fun)
    fun(c_11, group.by="orig.ident",pt.size = 0.5,label = F,
        label.repel = T,alpha = 0.9,cols = Singler.colors,
        no.legend = F,label.size = 4, repel = T,title = "Cluster 11",
        do.print = T, do.return = F))

object@meta.data %<>% cbind(object[["umap"]]@cell.embeddings)
object$clusters = as.character(object$SCT_snn_res.0.8)
# change cluster 11 Pt2
c11_pt2 <- object$SCT_snn_res.0.8 == 11 & object$UMAP_1 < 0
object@meta.data[c11_pt2,"clusters"] %<>% paste0("_Pt2")
# subset cluster 5
c_5 <- subset(object, idents = 5)
DefaultAssay(c_5)  = "SCT"
c_5 %<>% FindNeighbors(reduction = "harmony",dims = 1:85) %>%
    FindClusters(resolution = 0.05)

TSNEPlot.1(c_5,do.print = T, title = "re-cluster cluster 5 ")
features <- FilterGenes(object,c("FCN1","ITGAL","ITGAM","FCGR1A",
                                 "MS4A7","CDKN1C", "CSF1R","FCGR3A",
                                 "VCAN","S100A8","CD14","CSF3R"))
FeaturePlot.1(c_5,features = features, pt.size = 0.005, cols = c("gray90", "red"),
              alpha = 1,reduction = "tsne",
              threshold = 1, text.size = 20, border = T,do.print = T, do.return = F,ncol = 4,
              units = "in",width=12, height=9, no.legend = T)

CD14 <- colnames(c_5)[c_5$SCT_snn_res.0.05 == 0]
CD16 <- colnames(c_5)[c_5$SCT_snn_res.0.05 == 1]
object@meta.data[CD14,"clusters"] %<>% paste0("_CD14")
object@meta.data[CD16,"clusters"] %<>% paste0("_CD16")
table(object$clusters)
lapply(c(TSNEPlot.1, UMAPPlot.1), function(fun)
    fun(object, group.by="clusters",pt.size = 1,
        label = T, label.repel = T,alpha = 0.9,cols = Singler.colors,
        no.legend = T,label.size = 4, repel = T,
        title = paste(formals(fun)$reduction, "plot for re-clustering of PALIBR phase I"),
        do.print = T, do.return = F))
object$UMAP_1 = NULL
object$UMAP_2 = NULL
save(object, file = "data/MCL_41_harmony_20210424.Rda")
