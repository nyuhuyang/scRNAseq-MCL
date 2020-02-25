########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
#devtools::install_github("immunogenomics/harmony", ref= "ee0877a",force = T)
invisible(lapply(c("Seurat","dplyr","kableExtra","ggplot2","cowplot","sctransform",
                   "harmony"), function(x) {
    suppressPackageStartupMessages(library(x,character.only = T))
}))
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#======3.1 subset B and MCL =========================
(load(file = "data/MCL_41_harmony_20200225.Rda"))
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
object %<>% subset(subset = UMAP_2 > -4)
object %<>% subset(subset = UMAP_1 > 4 & UMAP_2 < 2, invert = T)

TSNEPlot.1(object, group.by = "cell.types",
   cols = c("#E6AB02","#2055da"),
   unique.name = "cell.types",do.print = T, 
   do.return = F)
lapply(c("groups","orig.ident","conditions","tissues"), function(group.by)
    TSNEPlot.1(object, group.by=group.by,pt.size = 0.5,label = F,
               cols = Singler.colors,
               label.repel = T,alpha = 0.9,
               unique.name = "group.by",
               no.legend = F,label.size = 4, repel = T, 
               title = paste("Harmony Integration by",group.by),
               do.print = T, do.return = F))
# re-run cluster
npcs = 70
object %<>% FindNeighbors(reduction = "harmony",dims = 1:npcs)
res =1
object %<>% FindClusters(resolution = res)
Idents(object) = "SCT_snn_res.1"
object %<>% RenameIdents("0" = "C1",
                         "1" = "C1",
                         "2" = "C2",
                         "3" = "C2",
                         "4" = "C4",
                         "5" = "C3",
                         "6" = "C4",
                         "7" = "C2",
                         "8" = "C2",
                         "9" = "C2",
                         "10"= "C1",
                         "11"= "C3",
                         "12"= "C4",
                         "13"= "C1",
                         "14"= "C1",
                         "15"= "C1",
                         "16"= "C1",
                         "17"= "C2",
                         "18"= "C1",
                         "19"= "C1",
                         "20"= "C1",
                         "21"= "C1",
                         "22"= "C1")
object[["X4clusters"]] = as.character(Idents(object))
Idents(object) = "X4clusters"
lapply(c("SCT_snn_res.1","X4clusters"), function(group.by)
    TSNEPlot.1(object, group.by=group.by,pt.size = 0.3,label = T,
               label.repel = T,alpha = 0.9,
               do.return = F,
               no.legend = T,label.size = 4, repel = T, 
               title = paste("res =0.1"),
               do.print = T))
saveRDS(object, file = "data/MCL_41_B_20200207.rds")

# cell.types false positive results  ========
table(object$cell.types, object$orig.ident) %>% kable %>% kable_styling()
normal_cells <- object$sample %in% c("BH","DJ","MD","NZ") %>% colnames(object)[.]
object@meta.data[normal_cells,"cell.types"] %<>% gsub("MCL","B_cells",.)

# UMI
object$orig.ident_X4cluster = gsub("N01|N02|N03","Normal",object$orig.ident)
object$orig.ident_X4cluster %<>% paste0("_",object$X4clusters)
object$orig.ident_X4cluster %<>% gsub("Normal_.*","Normal",.)

Idents(object) = "orig.ident_X4cluster"
exp = AverageExpression(object,assays = "SCT")
write.csv(exp$SCT,paste0(path,"MCL_41_UMI.csv"))

cell_number = table(object$orig.ident_X4cluster) %>% 
    as.data.frame() %>% t
rownames(cell_number) = c("Sample","Cell.number")
write.csv(cell_number,paste0(path,"MCL_41_UMI_cell_number.csv"))
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
TSNEPlot.1(object, group.by="cell.types",pt.size = 0.5,label = F,
    label.repel = T,alpha = 0.9,
    cols = ExtractMetaColor(object),
    unique.name = "cell.types",
    no.legend = T,label.size = 4, repel = T, 
    title = "rerun Harmony on MCL and B cells",
    do.print = T, do.return = F)

lapply(c("groups","orig.ident","conditions","tissues"), function(group.by)
    TSNEPlot.1(B_cells_MCL, group.by=group.by,pt.size = 0.5,label = F,
               cols = Singler.colors,
               label.repel = T,alpha = 0.9,
               unique.name = "group.by",
               no.legend = F,label.size = 4, repel = T, 
               title = paste("Harmony Integration by",group.by),
               do.print = T, do.return = F))

res = c(seq(0.1,2, by = 0.1))
for(i in seq_along(res)){
    object %<>% FindClusters(resolution = res[i])
    Idents(object) = paste0("SCT_snn_res.",res[i])
    TSNEPlot.1(object, group.by=paste0("SCT_snn_res.",res[i]),pt.size = 0.3,label = T,
               label.repel = T,alpha = 0.9,
               do.return = F,
               no.legend = T,label.size = 4, repel = T, 
               title = paste("res =",res[i],"in B and MCL based on harmony"),
               do.print = T, save.path = path)
    file.rename(paste0(path,"TSNEPlot_object_SCT_snn_res.",res[i],".jpeg"),
                paste0(path,i,"-TSNEPlot_object_SCT_snn_res.",res[i],".jpeg"))
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
object[["X4clusters"]] = as.character(Idents(object))
Idents(object) = "X4clusters"
object %<>% sortIdent()
saveRDS(object, file = "data/MCL_41_B_20200204.rds")
#======3.3 fine adjust cluster =========================
X4clusters_markers = read.csv(file= paste0("Yang/Figure 2/Figure Sources/",
                                           "MCL_41-FC0.05.csv"),
                              row.names = 1, stringsAsFactors=F)
Top_n = 40
top = X4clusters_markers %>% group_by(cluster) %>%
    top_n(Top_n, cluster) %>% top_n(Top_n, avg_logFC)
C2_top <- top[top$cluster %in% "C2",]
Idents(object) = "X4clusters"
C1 <- subset(object, idents = "C1")
C1 %<>% AddModuleScore(features = list(C2_top$gene),ctrl = 5,name = "C2_top_gene")

FeaturePlot.1(C1,features = "C2_top_gene1", pt.size = 1, cols = c("gray90", "red"),
              alpha = 1,reduction = "tsne", 
             text.size = 20, border = T,do.print = T, do.return = F,
              units = "in",width=9, height=12, no.legend = T)
C1 %<>% ScaleData(features=C2_top$gene)
Idents(C1) = "conditions"
UMAPPlot(C1,label=T)
DoHeatmap.1(C1, features = C2_top$gene, Top_n = Top_n,
            do.print=T, angle = 90, group.bar = T, title.size = 0, no.legend = F,size=5,hjust = 0.5,
            group.bar.height = 0.02, label=T, cex.row= 2, legend.size = 0,width=10, height=6.5,
            pal_gsea = FALSE,
            unique.name = "cell.types",
            title = "Top 40 DE genes in C1 clusters using new DEGs",
            save.path = path)
rerunH <- readRDS(paste0(path,"rerunH_B_SCT_snn_res.0.3.rds"))