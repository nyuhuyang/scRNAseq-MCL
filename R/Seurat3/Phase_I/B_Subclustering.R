########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
devtools::install_github("immunogenomics/harmony", ref= "ee0877a",force = T)
invisible(lapply(c("Seurat","dplyr","kableExtra","ggplot2","cowplot","sctransform",
                   "harmony"), function(x) {
    suppressPackageStartupMessages(library(x,character.only = T))
}))
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#======3.1 subset B and MCL =========================
(load(file = "data/MCL_41_harmony_20191231.Rda"))
Idents(object) = "Doublets"
object %<>% subset(idents = "Singlet")
Idents(object) = "cell.types"
object %<>% subset(idents = c("B_cells","MCL"))

jpeg(paste0(path,"B_MCL_subset.jpeg"), units="in", width=10, height=7,res=600)
UMAPPlot(object, cols = c("#E6AB02","#2055da"),
           group.by = "cell.types")+
    ggtitle("remove sparse B and MCL cells")+
    TitleCenter()+
    geom_segment(aes(x = -16, y = -4, xend = 4, yend = -4))+
    geom_segment(aes(x = 4, y = -4, xend = 4, yend = 2))+
    geom_segment(aes(x = 4, y = 2, xend = 10, yend = 2))
dev.off()
object@meta.data %<>% cbind(object[["umap"]]@cell.embeddings )
rm_cell <- (object$UMAP_2 < -4) | (object$UMAP_2 < 2 & object$UMAP_1 > 4)
object %<>% subset(subset = UMAP_2 > -4)
object %<>% subset(subset = UMAP_1 > 4 & UMAP_2 < 2, invert = T)

lapply(c(UMAPPlot.1,TSNEPlot.1), function(fun)
    fun(object, group.by = "cell.types",
       cols = c("#E6AB02","#2055da"),
       unique.name = "cell.types",do.print = T, 
       do.return = F)
)

#======3.2 rerun harmony =========================
object@reductions = list();GC()

DefaultAssay(object)  = "SCT"
object <- FindVariableFeatures(object, selection.method = "vst",
                               num.bin = 20, nfeatures = 2000,
                               mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))
object %<>% ScaleData
object %<>% RunPCA(verbose = T,npcs = 100)

jpeg(paste0(path,"S1_B_ElbowPlot.jpeg"), units="in", width=10, height=7,res=600)
ElbowPlot(object, ndims = 100)
dev.off()
a <- seq(1,97, by = 6)
b <- a+5
for(i in seq_along(a)){
    jpeg(paste0(path,"DimHeatmap_B_pca_",i,"_",a[i],"_",min(b[i],100),".jpeg"), units="in", width=10, height=7,res=600)
    DimHeatmap(object, dims = a[i]:min(b[i],100),
               nfeatures = 30,reduction = "pca")
    dev.off() 
}

object %<>% JackStraw(num.replicate = 20,dims = 100)
object %<>% ScoreJackStraw(dims = 1:100)
a <- seq(1,100, by = 10)
b <- a+9
for(i in seq_along(a)){
    jpeg(paste0(path,"JackStrawPlot_B_",i,"_",a[i],"_",min(b[i],100),".jpeg"), units="in", width=10, height=7,res=600)
    print(JackStrawPlot(object, dims = a[i]:min(b[i],100)))
    Progress(i,length(a))
    dev.off()
}

npcs = 70
jpeg(paste0(path,"S1_RunHarmony.jpeg"), units="in", width=10, height=7,res=600)
system.time(object %<>% RunHarmony.1(group.by = "orig.ident", dims.use = 1:npcs,
                                   theta = 2, plot_convergence = TRUE,
                                   nclust = 50, max.iter.cluster = 100))
dev.off()

object %<>% FindNeighbors(reduction = "harmony",dims = 1:npcs)
system.time(object %<>% RunTSNE(reduction = "harmony", dims = 1:npcs))
object %<>% RunUMAP(reduction = "harmony", dims = 1:npcs)
saveRDS(object, file = "data/MCL_41_B_20200204.rds")

Idents(object)="cell.types"
lapply(c(TSNEPlot.1, UMAPPlot.1), function(fun) 
    fun(object, group.by="cell.types",pt.size = 0.5,label = F,
        label.repel = T,alpha = 0.9,
        cols = ExtractMetaColor(object),
        unique.name = "cell.types",
        no.legend = T,label.size = 4, repel = T, 
        title = "rerun Harmony on MCL and B cells",
        do.print = T, do.return = F))
for(group.by in c("groups","orig.ident","conditions","tissues")){
    lapply(c(TSNEPlot.1, UMAPPlot.1), function(fun) 
        fun(object, group.by=group.by,pt.size = 0.5,label = F,
            cols = Singler.colors,
            label.repel = T,alpha = 0.9,
            unique.name = "cell.types",
            no.legend = F,label.size = 4, repel = T, 
            title = paste("Harmony Integration by",group.by),
            do.print = T, do.return = F))
}
res = c(seq(0.001,0.009, by = 0.001),seq(0.01,0.09, by = 0.01),seq(0.1,2, by = 0.1))
for(i in seq_along(res)){
    object %<>% FindClusters(resolution = res[i])
    Idents(object) = paste0("SCT_snn_res.",res[i])
    UMAPPlot.1(object, group.by=paste0("SCT_snn_res.",res[i]),pt.size = 0.3,label = T,
               label.repel = T,alpha = 0.9,
               do.return = F,
               no.legend = T,label.size = 4, repel = T, 
               title = paste("res =",res[i],"in B and MCL based on harmony"),
               do.print = T, save.path = path)
    file.rename(paste0(path,"UMAPPlot_object_SCT_snn_res.",res[i],".jpeg"),
                paste0(path,i,"-UMAPPlot_object_SCT_snn_res.",res[i],".jpeg"))
    Progress(i,length(res))
}
object@meta.data = object@meta.data[,c("orig.ident","sample","tests","percent.mt","conditions",
                                "projects","tissues","nCount_SCT","nFeature_SCT",
                                "singler1sub","singler1main","cell.types","cell.types.colors",
                                "Doublets","SCT_snn_res.0.3")]
Idents(object) = "SCT_snn_res.0.3"
object %<>% RenameIdents("0" = "C1",
                         "1" = "C2",
                         "2" = "C4",
                         "3" = "C3",
                         "4" = "C1",
                         "5" = "C5",
                         "6" = "C5",
                         "7" = "C5",
                         "8" = "C1",
                         "9" = "C5",
                         "10" = "C4",
                         "11" = "C5",
                         "12" = "C1",
                         "13" = "C5")
object[["X5clusters"]] = as.character(Idents(object))
Idents(object) = "X5clusters"
object %<>% sortIdent()
saveRDS(object, file = "data/MCL_41_B_20200204.rds")
