########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
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
source("../R/Seurat3_functions.R")
path <- "Yang/Figure 2/Figure Sources/"
if(!dir.exists(path)) dir.create(path, recursive = T)

# load data
(load(file="data/MCL_V3_Harmony_43_20190627.Rda"))


df_samples <- readxl::read_excel("doc/191030_scRNAseq_info.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
sample_n = which(df_samples$tests %in% c("control",paste0("test",2:12)))
df_samples = df_samples[sample_n,]

object@meta.data$orig.ident %<>% plyr::mapvalues(from = unique(df_samples$sample),
                                                      to = unique(df_samples$publication.id))
table(object@meta.data$orig.ident)
NewNames = paste0(object@meta.data$orig.ident,"_",object@meta.data$Barcode)
object %<>% RenameCells(new.names = NewNames)
rownames(object@reductions$tsne@cell.embeddings) = colnames(object)

Idents(object) = "groups"

object %<>% subset(idents = c("Pt-203","Pt-204"),invert = T)
Idents(object) = "orig.ident"
(df_samples = df_samples[df_samples$publication.id %in% object$orig.ident,])
(samples = df_samples$publication.id[df_samples$publication.id %in% object$orig.ident])
Idents(object) %<>% factor(levels = samples)
table(Idents(object))

# preprocess
#Idents(object) = "Doublets"
#object %<>% subset(idents = "Singlet")
Idents(object) = "manual"
object <- subset(object,idents = c("HSC/progenitors","Nonhematopoietic cells"), invert = TRUE)
Idents(object) = "res.0.6"
object %<>% subset(idents = 11,invert = T)
object %<>% sortIdent()
remove <- object$orig.ident %in% "Pt2_30Pd" & object$res.0.6 %in% 9
object <- object[,!remove]

#==== Figure 2-A ===========
object@meta.data$manual %<>% plyr::mapvalues(from = c("B_cells","MCL",
                                                      "Myeloid cells",
                                                      "NK_cells","T_cells:CD4+",
                                                      "T_cells:CD8+"),
                                             to = c("B","MCL",
                                                    "Monocytes",
                                                    "NK","CD4 T",
                                                    "CD8 T"))
object@meta.data$manual %<>% as.character()
Idents(object) = "manual"
object %<>% sortIdent()
TSNEPlot.1(object = object, label = F, label.repel = F, group.by = "manual",
           cols = ExtractMetaColor(object),no.legend = F,border = T,
           pt.size = 0.5, do.print = T,do.return = F,legend.size = 25,
           title.size = 20,legend.title = "tSNE plots for cell types of 36 samples",
           units= "in",width=10, height=7,hjust =0.5)
npcs =50
object %<>% RunUMAP(reduction = "harmony", dims = 1:npcs)
UMAPPlot.1(object = object, label = F, label.repel = F, group.by = "manual",
           cols = ExtractMetaColor(object),no.legend = F,border = T,
           pt.size = 0.5, do.print = T,do.return = F,legend.size = 25,
           title.size = 20,legend.title = "UMAP plots for cell types of 36 samples",
           units= "in",width=10, height=7,hjust =0.5)
#==== Figure 2-B ===========
features <- FilterGenes(object,c("CD19","CCND1","SOX11",
                                 "CD3D","CD4","CD8A",
                                 "MS4A7","CD14","FCGR1A",
                                 "GNLY","KLRC1","NCAM1"))
FeaturePlot.1(object,features = features, pt.size = 0.005, cols = c("gray90", "red"),
              alpha = 1,reduction = "tsne",
              threshold = 1, strip.text.size = 20, border = T,do.print = T, do.return = F,ncol = 3, 
              units = "in",width=9, height=12, no.legend = T)
FeaturePlot.1(object,features = features, pt.size = 0.005, cols = c("gray90", "red"),
              alpha = 1,reduction = "umap",
              threshold = 1, strip.text.size = 20, border = T,do.print = T, do.return = F,ncol = 3, 
              units = "in",width=9, height=12, no.legend = T)

#==== Figure 2-C ===========
(load(file = "data/B_cells_MCL_43_20190917.Rda"))

B_cells_MCL@meta.data$orig.ident %<>% plyr::mapvalues(from = unique(df_samples$sample),
                                                 to = unique(df_samples$publication.id))
table(B_cells_MCL@meta.data$orig.ident)
B_cells_MCL$orig.ident %<>% as.character()
NewNames = paste0(B_cells_MCL@meta.data$orig.ident,"_",B_cells_MCL@meta.data$Barcode)
B_cells_MCL %<>% RenameCells(new.names = NewNames)
rownames(B_cells_MCL@reductions$tsne@cell.embeddings) = colnames(B_cells_MCL)

Idents(B_cells_MCL) = "orig.ident"
(samples = df_samples$publication.id[df_samples$publication.id %in% B_cells_MCL$orig.ident])
B_cells_MCL %<>% subset(idents = samples)
chose <- c("Normal","X5_clusters")[1]
if(chose == "Normal"){
        B_cells_MCL$X5_clusters_normal = as.numeric(as.character(B_cells_MCL@meta.data$X5_clusters))
        normal <- grepl("N02|N01|N03|N04",B_cells_MCL$orig.ident)
        B_cells_MCL@meta.data[normal,"X5_clusters_normal"] = "Normal"
        Idents(B_cells_MCL) = "X5_clusters_normal"
        B_cells_MCL %<>% sortIdent()
        table(Idents(B_cells_MCL))
        X5_clusters_normal_markers <- FindPairMarkers(B_cells_MCL,ident.1 = 1:5, ident.2 = rep("Normal",5),
                                                      logfc.threshold = 0,only.pos = F,
                                                      min.pct = 0.1,return.thresh = 1,save.path = path)
        
        write.csv(X5_clusters_normal_markers,paste0(path,"X5_clusters_normal_FC0_markers.csv"))
        X5_clusters_markers = read.csv(file=paste0(path,"X5_clusters_normal_FC0_markers.csv"),
                                              row.names = 1, stringsAsFactors=F)
        colnames(X5_clusters_markers)[grep("cluster",colnames(X5_clusters_markers))] = "cluster"
}

if(chose == "X5_clusters"){
        Idents(B_cells_MCL) = "X5_clusters"
        B_cells_MCL <- sortIdent(B_cells_MCL,numeric = T)
        table(Idents(B_cells_MCL))
        X5_clusters_markers <- FindAllMarkers.UMI(B_cells_MCL,
                                                  logfc.threshold = 0,only.pos = F, 
                                                  min.pct = 0.1,return.thresh = 1)
        
        write.csv(X5_clusters_markers,paste0(path,"X5_clusters_FC0_markers.csv"))
        X5_clusters_markers = read.csv(file=paste0(path,"X5_clusters_FC0_markers.csv"),
                                       row.names = 1, stringsAsFactors=F)

        markers <- FilterGenes(B_cells_MCL,c("CCND1","CD19","CD5","CDK4","RB1","BTK","SOX11"))
        (MT_gene <- grep("^MT-",X5_clusters_markers$gene))
        X5_clusters_markers = X5_clusters_markers[-MT_gene,]
        Top_n = 40
        top = X5_clusters_markers %>% group_by(cluster) %>% top_n(Top_n, avg_logFC)
        features = c(as.character(top$gene),
                     tail(VariableFeatures(object = B_cells_MCL), 2),
                     markers)
        B_cells_MCL %<>% ScaleData(features=features)
        featuresNum <- make.unique(features, sep = ".")
        B_cells_MCL = MakeUniqueGenes(object = B_cells_MCL, features = features)
        
        Idents(B_cells_MCL) = ifelse(chose == "Normal","X5_clusters_normal","X5_clusters")
        Idents(B_cells_MCL) %<>% factor(levels = c(ifelse(chose == "Normal","Normal",""),
                                                   1:5))
        DoHeatmap.1(B_cells_MCL, features = featuresNum, Top_n = Top_n,
                    do.print=T, angle = 0, group.bar = F, title.size = 0, no.legend = F,size=5,hjust = 0.5,
                    group.bar.height = 0, label=F, cex.row= 2, legend.size = 0,width=10, height=6.5,
                    pal_gsea = FALSE,
                    title = NULL)
}


#==== Figure 2-D ===========
chose <- c("Normal","X5_clusters")[2]
if(chose == "Normal"){
        res = read.csv(file=paste0(path,"X5_clusters_normal_FC0.01_markers.csv"),
                       row.names = 1, stringsAsFactors=F)
        table(res$cluster1.vs.cluster2)
        res$cluster1.vs.cluster2 = gsub(" vs.Normal","",res$cluster1.vs.cluster2)
}
if(chose == "X5_clusters"){
        res = read.csv(file=paste0(path,"X5_clusters_FC0.01_markers.csv"),
                       row.names = 1, stringsAsFactors=F)
        table(res$cluster)
}

head(res)
res = res[order(res["p_val_adj"]),]
head(res, 20)
(clusters <- unique(res$cluster))
hallmark <- fgsea::gmtPathways("../seurat_resources/msigdb/h.all.v6.2.symbols.gmt")
names(hallmark) = gsub("HALLMARK_","",names(hallmark))
names(hallmark) = gsub("\\_"," ",names(hallmark))

FgseaDotPlot(stats=res, pathways=hallmark, nperm=1000,padj = 0.25,pval = 0.05,
             order.by = c(4,"NES"),decreasing = F,
             size = "-log10(pval)", fill = "NES",sample = "B_MCL clusters", 
             pathway.name = "Hallmark",rotate.x.text = F,
             font.ytickslab = 13,font.xtickslab = 15,hjust = 1,
             width=6.5, height=7)


#==== Figure 2-G ===========
table(Idents(object))
object@meta.data$orig.ident = gsub("N02|N01|N03|N04","Normal",object@meta.data$orig.ident)
Idents(object) = "orig.ident"


df_samples <- readxl::read_excel("doc/191030_scRNAseq_info.xlsx")
colnames(df_samples) <- tolower(colnames(df_samples))
tests <- paste0("test",c(9))
for(test in tests){
        sample_n = which(df_samples$tests %in% test)
        df <- as.data.frame(df_samples[sample_n,])
        (samples <- unique(df$publication.id))
        rownames(df) = samples
        
        samples <- c(df$publication.id[order(df$tsne)])
        print(samples <- samples[!is.na(samples)])
        
        subset_object <- subset(object, idents = samples)
        subset_object$orig.ident %<>% factor(levels = samples)
        
        Idents(subset_object) = "manual"
        
        subset_object %<>% sortIdent()
        TSNEPlot.1(subset_object, pt.size =0.3, 
                   strip.text.size = 14,no.legend = T,
                   group.by = "manual",split.by = "orig.ident",legend.size = 0,
                   cols = ExtractMetaColor(subset_object), ncol = length(samples),
                   unique.name = "groups", do.print = T,do.return = F,border = T,
                   width=8.5, height=2)
}


#==== Figure 3-C ===========
#(load("data/B_cells_MCL_43_20190904.Rda"))
#B_cells_MCL@meta.data$orig.ident = gsub("N02|N01|N03|N04","Normal",B_cells_MCL@meta.data$orig.ident)
B_cells_MCL@meta.data$X5_clusters_normal = as.numeric(as.character(B_cells_MCL@meta.data$X5_clusters))
normal <- B_cells_MCL$orig.ident %in% "Normal"
B_cells_MCL@meta.data[normal,"X5_clusters_normal"] = "Normal"
B_cells_MCL <- sortIdent(B_cells_MCL,numeric = T)
Idents(B_cells_MCL) = "X5_clusters_normal"
table(Idents(B_cells_MCL))

Idents(B_cells_MCL) = "orig.ident"
Pt_10 <- subset(B_cells_MCL, idents = c("Pt-10-LN-C2","Normal"))
Pt_10@meta.data$X5_clusters_normal = as.numeric(as.character(Pt_10@meta.data$X5_clusters))
normal <- Pt_10$orig.ident %in% "Normal"
Pt_10@meta.data[normal,"X5_clusters_normal"] <- 0
Idents(Pt_10) = "X5_clusters_normal"
Pt_10 <- sortIdent(Pt_10,numeric = T)
Idents(Pt_10) = "X5_clusters_normal"
table(Idents(Pt_10))
Pt_10_DE <- FindPairMarkers(Pt_10,ident.1 = 1:4, ident.2 = rep(0,4),
                            logfc.threshold = 0.1,only.pos = T,
                            min.pct = 0.1,return.thresh = 0.05,save.path = path)

write.csv(Pt_10_DE,paste0(path,"Pt_10_DE_X5_clusters_normal_FC0.1_markers.csv"))
Pt_10_DE = read.csv("output/20190904/Pt_10_DE_X5_clusters_normal_FC0.1_markers.csv",row.names = 1)

markers <- FilterGenes(Pt_10,c("CCND1","CD19","CD5","CDK4","RB1","BTK","SOX11"))
(MT_gene <- grep("^MT-",Pt_10_DE$gene))
Pt_10_DE = Pt_10_DE[-MT_gene,]
Top_n = 40
top = Pt_10_DE %>% group_by(cluster1.vs.cluster2) %>% top_n(Top_n, avg_logFC)
Idents(Pt_10) = "orig.ident"
Pt_10 <- subset(Pt_10, idents = c("Pt-10-LN-C2"))
Pt_10 %<>% ScaleData(features=unique(c(as.character(top$gene),markers)))
featuresNum <- make.unique(c(as.character(top$gene),markers), sep = ".")
Pt_10 %<>% MakeUniqueGenes(top = top, features = markers)

Idents(Pt_10) = "X5_clusters_normal"
DoHeatmap.1(Pt_10, features = featuresNum, Top_n = Top_n,
            do.print=T, angle = 0, group.bar = F, title.size = 20, no.legend = F,size=5,hjust = 0.5,
            group.bar.height = 0, label=F, cex.row= 2, legend.size = 0,width=10, height=6.5,
            pal_gsea = FALSE,
            title = paste("Top",Top_n,"differentially expressed genes in MCL vs Normal in each cluster"))


#==== Figure 2-C ===========
(load(file = "data/B_cells_MCL_43_20190713.Rda"))
B_cells_MCL <- subset(object,idents = c("B_cells","MCL"))
Idents(B_cells_MCL) <-  "res.0.6"
B_cells_MCL <- subset(B_cells_MCL, idents = c(0,1,5,6,9,13,14,18,19,20))
B_cells_MCL@meta.data$X5_clusters <- plyr::mapvalues(x = B_cells_MCL@meta.data$res.0.6,
                                                     from = c(0,1,13,14,18,19,20,5,6,9),
                                                     to =   c(1,2,2,  1, 1, 2, 1, 3,4,5))
Idents(B_cells_MCL) <- "X5_clusters"

TSNEPlot.1(B_cells_MCL, pt.size = 0.2,label = T, label.repel = T,
           do.print = T,no.legend = T,border = T,alpha = 1,
           label.size = 5, repel = T, title = "5 clusters in B/MCL cells",
           width=7, height=7)
B_cells_MCL <- PercentageFeatureSet(B_cells_MCL, pattern = "^RPL", col.name = "percent.RPL")

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
B_cells_MCL <- CellCycleScoring(B_cells_MCL, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

FeaturePlot.1(B_cells_MCL,features = c("percent.RPL","CD274","S.Score","G2M.Score"), 
              pt.size = 0.1, cols = c("gray90", "red"), alpha = 1,
              threshold = 1, strip.text.size = 20, border = T,do.print = T, 
              do.return = F,ncol = 4,
              units = "in",width=13, height=4)
