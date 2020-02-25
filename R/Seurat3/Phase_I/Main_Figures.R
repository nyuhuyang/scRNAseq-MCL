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
library(fgsea)
source("../R/Seurat3_functions.R")

# load data
(load(file="data/MCL_41_harmony_20191231.Rda"))
df_samples <- readxl::read_excel("doc/191120_scRNAseq_info.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
table(object$orig.ident)

(df_samples = df_samples[df_samples$`sample name` %in% object$orig.ident,])
(samples = df_samples$`sample name`[df_samples$`sample name` %in% object$orig.ident])

# preprocess
object$orig.ident %<>% factor(levels = samples)
Idents(object) = "orig.ident"
object %<>% subset(idents = "Pt2_30Pd", invert = T)
Idents(object) = "Doublets"
object %<>% subset(idents = "Singlet")
Idents(object) = "cell.types"
object %<>% subset(idents = c("HSC/progenitors","Nonhematopoietic cells"), invert = TRUE)
table(Idents(object))

#==== Figure 3-A ===========
path <- "Yang/Figure 3/Figure Sources/"
if(!dir.exists(path)) dir.create(path, recursive = T)

object$cell_types <- plyr::mapvalues(object@meta.data$cell.types,from = c("B_cells","MCL",
                                                      "Myeloid cells",
                                                      "NK_cells","T_cells:CD4+",
                                                      "T_cells:CD8+"),
                                             to = c("B","MCL",
                                                    "Monocytes",
                                                    "NK","CD4 T",
                                                    "CD8 T"))
object$cell_types %<>% as.character()
object$cell_types.colors = object$cell.types.colors
Idents(object) = "cell_types"
object %<>% sortIdent()
TSNEPlot.1(object = object, label = F, label.repel = F, group.by = "cell_types",
           cols = ExtractMetaColor(object),no.legend = F,border = T,
           pt.size = 0.1, do.print = T,do.return = F,legend.size = 25,
           title.size = 20,title = "tSNE plots for cell types of 41 samples",
           units= "in",width=9, height=7,hjust =0.5, save.path = path)
UMAPPlot.1(object = object, label = F, label.repel = F, group.by = "cell_types",
           cols = ExtractMetaColor(object),no.legend = F,border = T,
           pt.size = 0.5, do.print = T,do.return = F,legend.size = 25,
           title.size = 20,title = "UMAP plots for cell types of 41 samples",
           units= "in",width=9, height=7,hjust =0.5, save.path = path)
#==== Figure 2-B ===========
features <- FilterGenes(object,c("CD19","CCND1","SOX11",
                                 "CD3D","CD4","CD8A",
                                 "MS4A7","CD14","FCGR1A",
                                 "GNLY","KLRC1","NCAM1"))
FeaturePlot.1(object,features = features, pt.size = 0.005, cols = c("gray90", "red"),
              alpha = 1,reduction = "tsne",
              threshold = 1, text.size = 20, border = T,do.print = T, do.return = F,ncol = 3, 
              units = "in",width=9, height=12, no.legend = T, save.path = path)
file.rename(paste0(path,"FeaturePlot__object_cell_types_CD19-CCND1-SOX11-CD3D-CD4-CD8A-MS4A7-CD14-FCGR1A-GNLY-KLRC1-NCAM1_tsne__.jpeg"),
            paste0(path,"FeaturePlot_label.jpeg"))
FeaturePlot.1(object,features = features, pt.size = 0.005, cols = c("gray90", "red"),
              alpha = 1,reduction = "tsne",
              threshold = 1, text.size = 0, border = T,do.print = T, do.return = F,ncol = 3, 
              units = "in",width=9, height=12, no.legend = T, save.path = path)
file.rename(paste0(path,"FeaturePlot__object_cell_types_CD19-CCND1-SOX11-CD3D-CD4-CD8A-MS4A7-CD14-FCGR1A-GNLY-KLRC1-NCAM1_tsne__.jpeg"),
            paste0(path,"FeaturePlot_nolabel.jpeg"))

#==== Figure 2-C ===========
path <- "Yang/Figure 3/Figure Sources/"
if(!dir.exists(path)) dir.create(path, recursive = T)

B_cells_MCL = readRDS(file = "data/MCL_41_B_20200207.rds")
Idents(B_cells_MCL) = "orig.ident" 
B_cells_MCL %<>% subset(idents = "Pt2_30Pd", invert = T)
markers <- FilterGenes(B_cells_MCL,c("CCND1","CD19","CD5","CDK4","RB1","BTK","SOX11"))
group_colors = c("#181ea4","#5f66ec","#f46072","#e6001c")
choose = c("X4clusters","X4cluster_vs_Normal")[1]
if(choose == "X4clusters"){
        Idents(B_cells_MCL) = "X4clusters"
        
        B_cells_MCL %<>% sortIdent()
        Idents(B_cells_MCL) %<>% factor(levels = paste0("C",1:4))
        table(Idents(B_cells_MCL))
        X4clusters_markers = read.csv(file= paste0(path,choose,"/X4clusters_41-FC0.csv"),
                                      row.names = 1, stringsAsFactors=F)
        table(X4clusters_markers$cluster)
        X4clusters_markers$cluster %<>% factor(levels = paste0("C",1:4))
        markers <- FilterGenes(B_cells_MCL,c("CCND1","CD19","CD5","CDK4","RB1","BTK","SOX11"))
        (MT_gene <- grep("^MT-",X4clusters_markers$gene))
        X4clusters_markers = X4clusters_markers[-MT_gene,]
        Top_n = 40
        top = X4clusters_markers %>% group_by(cluster) %>%
                top_n(Top_n, cluster) %>% top_n(Top_n, avg_logFC)
        unique(top$cluster)
        top = top[order(top$cluster),]
        write.csv(top,paste0(path,choose,"/top40_genes_heatmap.csv"))
        features = c(as.character(top$gene),
                     tail(VariableFeatures(object = B_cells_MCL), 2),
                     markers)
        B_cells_MCL %<>% ScaleData(features=features)
        featuresNum <- make.unique(features, sep = ".")
        B_cells_MCL %<>% MakeUniqueGenes(features = features)
        
        DoHeatmap.1(B_cells_MCL, features = featuresNum, Top_n = Top_n,
                    do.print=T, angle = 0, group.bar = T, title.size = 0, no.legend = F,size=5,hjust = 0.5,
                    group.bar.height = 0.02, label=T, cex.row= 2, legend.size = 0,width=10, height=6.5,
                    pal_gsea = FALSE,
                    unique.name = "cell.types",
                    title = "Top 40 DE genes in 4 B/MCL clusters",
                    save.path = paste0(path,choose,"/"))
        file.rename(paste0(path,choose,"/Heatmap_top40_B_cells_MCL_B_cells_MCL_X4clusters_Legend.jpeg"),
                    paste0(path,choose,"/Heatmap_top40_MCL_B_cells_X4clusters.jpeg"))
}


if(choose == "X4cluster_vs_Normal"){
        B_cells_MCL$X4clusters_normal = as.character(B_cells_MCL$X4clusters)
        normal <- grepl("N01|N02|N03",B_cells_MCL$orig.ident)
        B_cells_MCL@meta.data[normal,"X4clusters_normal"] = "Normal"
        Idents(B_cells_MCL) = "X4clusters_normal"
        B_cells_MCL %<>% sortIdent()
        table(Idents(B_cells_MCL))
        X4clusters_normal_markers <- FindPairMarkers(B_cells_MCL,
                                                     ident.1 = c("Normal",paste0("C",1:4)), 
                                                     ident.2 = list(paste0("C",1:4),
                                                                    "Normal","Normal","Normal","Normal"),
                                                      logfc.threshold = 0.25,only.pos = T,
                                                      min.pct = 0.1,return.thresh = 0.05,
                                                     save.path = path,
                                                     latent.vars = "nCount_SCT")
        
        write.csv(X4clusters_normal_markers,paste0(path, choose,"/X4clusters_normal_FC0.25_markers.csv"))
        X4clusters_markers = read.csv(file=paste0(path, choose,"/X4clusters_normal_FC0.25_markers.csv"),
                                       row.names = 1, stringsAsFactors=F)
        colnames(X4clusters_markers)[grep("cluster",colnames(X4clusters_markers))] = "cluster"
        X4clusters_markers$cluster %<>% gsub(" /.*","",.)
        table(X4clusters_markers$cluster)
        markers <- FilterGenes(B_cells_MCL,c("CCND1","CD19","CD5","CDK4","RB1","BTK","SOX11"))
        (MT_gene <- grep("^MT-",X4clusters_markers$gene))
        X4clusters_markers = X4clusters_markers[-MT_gene,]
        Top_n = 40
        top = X4clusters_markers %>% group_by(cluster) %>% top_n(Top_n, avg_logFC)
        table(top$cluster)
        top = top[top$cluster %in% c("C1","C2","C3","C4"),]
        write.csv(top,paste0(path,choose,"/top40_4clusters_over_normal_genes_heatmap.csv"))
        features = c(as.character(top$gene),
                     tail(VariableFeatures(object = B_cells_MCL), 2),
                     markers)
        B_cells_MCL %<>% ScaleData(features=features)
        featuresNum <- make.unique(features, sep = ".")
        B_cells_MCL = MakeUniqueGenes(object = B_cells_MCL, features = features)
        
        Idents(B_cells_MCL) = "X4clusters_normal"
        Idents(B_cells_MCL) %<>% factor(levels = c("Normal", paste0("C",1:4)))
        table(Idents(B_cells_MCL))
        DoHeatmap.1(B_cells_MCL, features = featuresNum, Top_n = Top_n,
                    do.print=T, angle = 0, group.bar = T, title.size = 0, no.legend = F,size=5,hjust = 0.5,
                    group.bar.height = 0.02, label=T, cex.row= 2, legend.size = 0,width=10, height=6.5,
                    pal_gsea = FALSE,
                    unique.name = "cell.types",
                    title = "Top 40 DE genes in 4 B/MCL clusters vs Normal",
                    save.path = path)
        file.rename(paste0(path,choose,"/Heatmap_top40_B_cells_MCL_B_cells_MCL_X4clusters_normal_Legend.jpeg"),
                    paste0(path,choose,"/Heatmap_top40_MCL_B_cells_X4clusters_normal.jpeg"))
}

# average expression
if(choose == "X4clusters"){
        Idents(B_cells_MCL) = "X4clusters"
        
        B_cells_MCL %<>% sortIdent()
        Idents(B_cells_MCL) %<>% factor(levels = paste0("C",1:4))
        table(Idents(B_cells_MCL))
        Top_n = 40
        top = read.csv(file= paste0(path,choose,"/top40_genes_heatmap.csv"),
                       row.names = 1, stringsAsFactors=F)
        top = top %>% group_by(cluster) %>% top_n(Top_n, avg_logFC)
        table(top$cluster)
        features = c(as.character(top$gene),
                     tail(VariableFeatures(object = B_cells_MCL), 2),
                     markers)
        featuresNum <- make.unique(features, sep = ".")
        exp = AverageExpression(B_cells_MCL[features,], 
                                assays = "SCT") %>% .$SCT
        exp[tail(VariableFeatures(object = B_cells_MCL), 2),] =0
        exp = MakeUniqueGenes(object = exp, features = features)
        
        scale_exp <- exp %>% t %>% scale %>% t
        DoHeatmap.matrix(scale_exp, features = featuresNum,
                         group.by = 1:4,size = 6,angle = 0,
                         draw.lines =F, raster = FALSE,
                         pal_gsea = FALSE,
                         width=1.5, height=10,res=600,no.legend = T,
                         cex.row=5,
                         group.colors = group_colors,
                         do.print = T,
                         unique.name = "cell.types",
                         title = "40 DEGs",
                         save.path = paste0(path,choose,"/"))
        }

if(choose == "X4cluster_vs_Normal"){
        B_cells_MCL$X4clusters_normal = as.character(B_cells_MCL$X4clusters)
        normal <- grepl("N01|N02|N03",B_cells_MCL$orig.ident)
        B_cells_MCL@meta.data[normal,"X4clusters_normal"] = "Normal"
        Idents(B_cells_MCL) = "X4clusters_normal"
        table(Idents(B_cells_MCL))
        top = read.csv(file= paste0(path,choose,"/top40_4clusters_over_normal_genes_heatmap.csv"),
                       row.names = 1, stringsAsFactors=F)
        table(top$cluster)
        top = top[top$cluster %in% c("C1","C2","C3","C4"),]
        features = c(as.character(top$gene),
                     tail(VariableFeatures(object = B_cells_MCL), 2),
                     markers)
        featuresNum <- make.unique(features, sep = ".")
        exp = AverageExpression(B_cells_MCL[features,], 
                                assays = "SCT") %>% .$SCT
        exp = MakeUniqueGenes(object = exp, features = features)
        
        scale_exp <- exp %>% t %>% scale %>% t
        #colnames(scale_exp) = c("Normal",1:4)
        #scale_exp = scale_exp[,c("Normal",1:4)]
        group.by = factor(c("Normal",1:4), levels = c("Normal",1:4))
        DoHeatmap.matrix(scale_exp, features = featuresNum,
                         group.by = group.by,size = 6,angle =60,
                         draw.lines =F, raster = FALSE,
                         pal_gsea = FALSE,
                         group.colors = c("#31aa3a",group_colors),
                         width=1.8, height=11,res=600,no.legend = T,
                         cex.row=5,
                         do.print = T,
                         unique.name = "cell.types",
                         title = "40 DEGs",
                         save.path = paste0(path,choose,"/"))
}

#==== Figure 3-D ===========
choose <- c("X4clusters","X4cluster_vs_Normal")[1]
res = read.csv(file = paste0("Yang/Figure 3/Figure Sources/",
                             choose,"/",choose,"_41-FC0.csv"),
               row.names = 1, stringsAsFactors=F)
if(choose == "X4clusters") table(res$cluster)
if(choose == "X4cluster_vs_Normal"){
        table(res$cluster1.vs.cluster2)
        res$cluster1.vs.cluster2 = gsub("Normal vs. ","",res$cluster1.vs.cluster2)
}

head(res)
res = res[order(res["p_val_adj"]),]
head(res, 20)
(clusters <- unique(res$cluster))
hallmark <- fgsea::gmtPathways("../seurat_resources/msigdb/h.all.v6.2.symbols.gmt")
names(hallmark) = gsub("HALLMARK_","",names(hallmark))
names(hallmark) = gsub("\\_"," ",names(hallmark))
hallmark$`NFKB SIGNALING` =  read.delim("data/200222 NFKB pathway gene list.txt") %>% 
        pull %>% as.character()
# Now, run the fgsea algorithm with 1000 permutations:
fgseaRes = FgseaDotPlot(stats=res, pathways=hallmark,
                        padj = 0.25,pval = 0.05,
                        order.yaxis.by = c("C4","NES"),
                        order.xaxis = paste0("C",1:4),
                        decreasing = F,
                        Rowv = F,Colv = F,
                        size = " -log10(pval)", fill = "NES",
                        pathway.name = "Hallmark",rotate.x.text = T,
                        title = "in B and MCL sub-cluster",
                        font.xtickslab=10, font.main=14, font.ytickslab = 10,
                        font.legend = list(size = 12),font.label = list(size = 12),
                        do.return = T,save.path = path, do.print = T,
                        width = 5,height = 6)
write.csv(fgseaRes, file = paste0(path,choose,"_FDR0.25_pval0.05.csv"))

##########################################
#==== Figure 3C ===========
##########################################
path <- "Yang/Figure 3/Figure Sources/"
if(!dir.exists(path)) dir.create(path, recursive = T)

Idents(object) = "groups"
table(Idents(object))
sub_object <- subset(object, idents = c("Normal", "Untreated"))
sub_object$orig.ident %<>% gsub("N02|N01|N03","Normal",.)
keep_normal = which(sub_object$orig.ident %in% "Normal")
keep_normal = sample(keep_normal,size = length(keep_normal)/3, replace = F)
keep_untreated = which(!sub_object$orig.ident %in% "Normal")
sub_object = sub_object[,c(keep_normal,keep_untreated)]
Idents(sub_object) = "orig.ident"
table(Idents(sub_object))

Idents(sub_object) = "cell_types"
TSNEPlot.1(sub_object, pt.size =0.3, 
           text.size = 14,no.legend = T,
           group.by = "cell_types",split.by = "orig.ident",legend.size = 0,
           cols = ExtractMetaColor(sub_object), ncol = length(unique(sub_object$orig.ident)),
           unique.name = "groups", do.print = T,do.return = F,border = T,
           width=8.5, height=2, save.path = path)

###############################
# All pairwise heatmaps:
###############################
path <- "Yang/B_pairwise_heatmaps/"
if(!dir.exists(path)) dir.create(path, recursive = T)
group_colors = c("#181ea4","#5f66ec","#f46072","#e6001c")
# =========Doheatmap for Normal / MCL ============
B_cells_MCL = readRDS(file = "data/MCL_41_B_20200207.rds")
B_cells_MCL$orig.ident %<>% gsub("N02|N01|N03","Normal",.)
B_cells_MCL$X4_orig.ident = paste(B_cells_MCL$orig.ident,
                                   B_cells_MCL$X4clusters, sep = "_")
B_cells_MCL@meta.data$X4_orig.ident = gsub('^Normal_.*', 'Normal', B_cells_MCL@meta.data$X4_orig.ident)
table(B_cells_MCL@meta.data$X4_orig.ident)
Idents(B_cells_MCL) = "X4_orig.ident"

df_samples <- readxl::read_excel("doc/191001_scRNAseq_info.xlsx",sheet = "heatmap")
list_samples <- df2list(df_samples)
scRNAseq_info <- readxl::read_excel("doc/191120_scRNAseq_info.xlsx")
colnames(scRNAseq_info) <- colnames(scRNAseq_info) %>% tolower
list_samples <- lapply(list_samples[c("MCL","MCL.1","MCL.2")],
                       function(x) plyr::mapvalues(x,
                                               from = scRNAseq_info$sample,
                                               to = scRNAseq_info$`sample name`))
all(list_samples %>% unlist %>% as.vector %>% unique %in% 
            B_cells_MCL$orig.ident)
Idents(B_cells_MCL) = "orig.ident"
markers <- FilterGenes(B_cells_MCL,c("CCND1","CD19","CD5","CDK4","RB1","BTK","SOX11"))
(block <- VariableFeatures(B_cells_MCL) %>% tail(2))

choose = c("X4cluster_vs_Normal","X4clusters")[1]
for(sample in list_samples$MCL[1]){
        subset.MCL <- subset(B_cells_MCL, idents = c("Normal",sample))
        
        # SplitTSNEPlot======
        Idents(subset.MCL) = "X4_orig.ident"
        subset.MCL %<>% sortIdent()
        TSNEPlot.1(subset.MCL, split.by = "orig.ident",pt.size = 1,label = F,
                   do.print = F, unique.name = T)
        
        # remove cluster with less than 3 cells======
        table_subset.MCL <- table(subset.MCL$X4_orig.ident) %>% as.data.frame
        keep.MCL <- table_subset.MCL[table_subset.MCL$Freq > 2,"Var1"] %>% as.character()
        
        (X4_cluster <- keep.MCL %>% unique %>% 
                        gsub('.*\\_C',"",.) %>% 
                        sub("Normal","",.) %>% 
                        as.numeric %>% sort )
        
        print(ident.1 <- rep("Normal",length(X4_cluster)))
        print(ident.2 <- paste(sample,X4_cluster,sep="_C"))
        subset.MCL <- subset(subset.MCL, idents = c(ident.1,ident.2))
        
        # FindAllMarkers.UMI======
        
        gde.markers <- FindPairMarkers(subset.MCL, ident.1 = c(ident.1,ident.2),
                                       ident.2 = c(ident.2,ident.1), only.pos = T,
                                       logfc.threshold = 0.1,min.cells.group =3,
                                       min.pct = 0.1,return.thresh = 0.05,
                                       latent.vars = "nCount_SCT")
        write.csv(gde.markers, paste0(path,"DE_analysis_files/","Normal_vs_",sample,".csv"))
        
        gde.markers = read.csv(paste0(path,"DE_analysis_files/","Normal_vs_",sample,".csv"),row.names = 1)
        (mito.genes <- grep(pattern = "^MT-", x = gde.markers$gene))
        if(length(mito.genes)>0) gde.markers = gde.markers[-mito.genes,]
        GC()
        #DoHeatmap.1======
        Top_n = 40
        top <-  gde.markers %>% group_by(cluster1.vs.cluster2) %>% 
                top_n(Top_n, avg_logFC) %>% as.data.frame()
        write.csv(top, paste0(path,"DE_analysis_files/","top40_Normal_vs_",sample,".csv"))
        features = c(as.character(top$gene),
                     tail(VariableFeatures(object = B_cells_MCL), 2),
                     markers)
        featuresNum <- make.unique(features, sep = ".")
        exp = AverageExpression(subset.MCL[features,], 
                                assays = "SCT") %>% .$SCT
        exp = MakeUniqueGenes(object = exp, features = features)
        exp[tail(VariableFeatures(object = B_cells_MCL), 2),] =0
        scale_exp <- exp %>% t %>% scale %>% t
        colnames(scale_exp) 
        group.by = factor(c("Normal",ident.2), levels = c("Normal",ident.2))
        DoHeatmap.matrix(scale_exp, features = featuresNum,
                         group.by = group.by,size = 6,angle =90,
                         draw.lines =F, raster = FALSE,
                         pal_gsea = FALSE,
                         group.colors = c("#31aa3a",group_colors),
                         width=2, height=22,res=600,no.legend = T,
                         cex.row=5,
                         do.print = T,
                         unique.name = "cell.types",
                         title = paste("40 DEGs in",sample,"_Normal"),
                         save.path = paste0(path,choose,"/"))
}

# =========Doheatmap for MCL.1 / MCL.2 ============
df_samples <- readxl::read_excel("doc/191001_scRNAseq_info.xlsx",sheet = "heatmap")
list_samples <- df2list(df_samples)
scRNAseq_info <- readxl::read_excel("doc/191120_scRNAseq_info.xlsx")
colnames(scRNAseq_info) <- colnames(scRNAseq_info) %>% tolower
list_samples <- lapply(list_samples[c("MCL","MCL.1","MCL.2")],
                       function(x) plyr::mapvalues(x,
                                                   from = scRNAseq_info$sample,
                                                   to = scRNAseq_info$`sample name`))
all(list_samples %>% unlist %>% as.vector %>% unique %in% 
            B_cells_MCL$orig.ident)
(block <- VariableFeatures(B_cells_MCL) %>% tail(2))
markers <- FilterGenes(B_cells_MCL,c("CCND1","CD19","CD5","CDK4","RB1","BTK","SOX11"))

Idents(B_cells_MCL) = "orig.ident"
choose = c("X4cluster_vs_Normal","X4clusters")[2]
for(i in 1:9){
        
        (samples1 = list_samples$MCL.1[i])
        (samples2 = list_samples$MCL.2[i])
        
        subset.MCL <- subset(B_cells_MCL, idents = c(samples1,samples2))
        subset.MCL@meta.data$orig.ident %<>% factor(levels = c(samples1,samples2))
        # remove cluster with less than 3 cells======
        
        table_subset.MCL <- table(subset.MCL@meta.data$X4_orig.ident) %>% as.data.frame
        (keep.MCL <- table_subset.MCL[table_subset.MCL$Freq > 2,"Var1"] %>% as.character())
        (X4_cluster <- keep.MCL %>% unique %>% 
                        gsub('.*\\_C',"",.) %>% as.numeric %>% sort %>% .[duplicated(.)])
        
        print(ident.1 <- paste(samples1,X4_cluster,sep="_C"))
        print(ident.2 <- paste(samples2,X4_cluster,sep="_C"))
        
        #---SplitTSNEPlot----
        Idents(subset.MCL) = "X4_orig.ident"
        subset.MCL <- subset(subset.MCL, idents = c(ident.1,ident.2))
        
        Idents(subset.MCL) %<>% factor(levels = c(ident.1,ident.2))
        TSNEPlot.1(subset.MCL, split.by = "orig.ident",pt.size = 1,label = F,
                   do.return = F,do.print = F, unique.name = T)
        
        #gde.markers <- FindPairMarkers(subset.MCL, ident.1 = c(ident.1,ident.2),
        #                               ident.2 = c(ident.2,ident.1), only.pos = T,
        #                               logfc.threshold = 0.1,min.cells.group =3,
        #                               min.pct = 0.1,return.thresh = 0.05,
        #                               latent.vars = "nCount_SCT")
        #write.csv(gde.markers, paste0(path,"DE_analysis_files/",samples1,"_vs_",samples2,".csv"))
        gde.markers = read.csv(paste0(path,"DE_analysis_files/",samples1,"_vs_",samples2,".csv"),row.names = 1)
        print(table(gde.markers$cluster1.vs.cluster2))
        (mito.genes <- grep(pattern = "^MT-", x = gde.markers$gene))
        if(length(mito.genes)>0) gde.markers = gde.markers[-mito.genes,]
        GC()
        #DoHeatmap.1======
        Top_n = 40
        top <-  gde.markers %>% group_by(cluster1.vs.cluster2) %>% 
                top_n(Top_n, avg_logFC) %>% as.data.frame()
        features = c(as.character(top$gene),
                     tail(VariableFeatures(object = B_cells_MCL), 2),
                     markers)
        featuresNum <- make.unique(features, sep = ".")
        exp = AverageExpression(subset.MCL[features,], 
                                assays = "SCT") %>% .$SCT
        exp = MakeUniqueGenes(object = exp, features = features)
        exp[tail(VariableFeatures(object = B_cells_MCL), 2),] =0
        scale_exp <- exp %>% t %>% scale %>% t
        colnames(scale_exp) 
        (group.by = c(ident.1, ident.2))
        DoHeatmap.matrix(scale_exp, features = featuresNum,
                         group.by = group.by,size = 6,angle =90,
                         draw.lines =F, raster = FALSE,
                         pal_gsea = FALSE,
                         group.colors = rep(group_colors[X4_cluster], 2),
                         width=3, height=25,res=600,no.legend = T,
                         cex.row=5,
                         do.print = T,
                         unique.name = "cell.types",
                         title = paste0("40 DEGs in ",samples1,"_",samples2),
                         save.path = paste0(path,choose,"/"))
}

#### Fig 4E ########################################################################
path <- "Yang/Figure 4/"
if(!dir.exists(path)) dir.create(path, recursive = T)
(samples1 = list_samples$MCL.1[1])
(samples2 = list_samples$MCL.2[1])

subset.MCL <- subset(B_cells_MCL, idents = c(samples1,samples2))
subset.MCL$orig.ident %<>% factor(levels = c(samples1,samples2))
# merge C1+C2 and C3+C4
subset.MCL$X4_orig.ident %<>% gsub("C1|C2","C1+2",.)
subset.MCL$X4_orig.ident %<>% gsub("C3|C4","C3+4",.)

print(ident.1 <- paste0(samples1,c("_C1+2","_C3+4")))
print(ident.2 <- paste0(samples2,c("_C1+2","_C3+4")))

Idents(subset.MCL) = "X4_orig.ident"
TSNEPlot.1(subset.MCL, split.by = "orig.ident",pt.size = 1,label = F,
           do.return = T,do.print = F, unique.name = T)

gde.markers <- FindPairMarkers(subset.MCL, ident.1 = ident.2,
                               ident.2 = ident.1, only.pos = F,
                               logfc.threshold = 0.1,min.cells.group =3,
                               min.pct = 0.1,return.thresh = 0.05,
                               latent.vars = "nCount_SCT")
write.csv(gde.markers, paste0(path,samples2,"_vs_",samples1,".csv"))
gde.markers = read.csv(paste0(path,samples2,"_vs_",samples1,".csv"),row.names = 1)
print(table(gde.markers$cluster1.vs.cluster2))
(mito.genes <- grep(pattern = "^MT-", x = gde.markers$gene))
if(length(mito.genes)>0) gde.markers = gde.markers[-mito.genes,]
GC()
gde.markers$cluster1.vs.cluster2 %<>% as.character()
Clusters <- c("C1","C3")
for(i in seq_along(Clusters)){
        gde <- gde.markers[grepl(Clusters[i], gde.markers$cluster1.vs.cluster2),]
        p <- VolcanoPlots(data = gde, cut_off_pvalue = 0.0000001, cut_off_logFC = 0.25,
                          top = 20, cols = c("#0000ff","#d2dae2","#ff0000"),alpha=0.8, size=2,
                          legend.size = 12)+
                ggtitle(paste0(Clusters[i],"/",i*2," AMB/SB"))+
                theme(plot.title = element_text(hjust = 0.5,size=15,face = "plain"))
        jpeg(paste0(path,"VolcanoPlots_",Clusters[i],"_",i*2,"_AMB_SB.jpeg"), units="in", width=10, height=7,res=600)
        print(p)
        dev.off()
}

        
        
### Fig. 7C ==========
path <- "Yang/Figure 7/Figure Sources/"
if(!dir.exists(path)) dir.create(path, recursive = T)
group_colors = c("#181ea4","#5f66ec","#f46072","#e6001c")

B_cells_MCL = readRDS(file = "data/MCL_41_B_20200207.rds")
choose = c("X4clusters","X4cluster_vs_Normal")[1]
if(choose == "X4clusters"){
        Idents(B_cells_MCL) = "orig.ident"
        B_cells_MCL %<>% subset(idents = "Pt10_LN2Pd")
        Idents(B_cells_MCL) = "X4clusters"
        B_cells_MCL %<>% sortIdent()
        table(Idents(B_cells_MCL))
        #system.time(MCL_markers <- FindAllMarkers.UMI(B_cells_MCL, 
        #                                              only.pos = T,
        #                                              test.use = "MAST",
        #                                              logfc.threshold = 0.05,
        #                                              min.pct = 0.1,return.thresh = 0.05,
        #                                              latent.vars = "nCount_SCT"))
        #write.csv(MCL_markers,paste0(path,"Pt10_LN2Pd_X4clusters_FC0.05_markers.csv"))
        X4clusters_markers = read.csv(file= paste0(path,"Pt10_LN2Pd_X4clusters_FC0.05_markers.csv"),
                                      row.names = 1, stringsAsFactors=F)
        table(X4clusters_markers$cluster)
        X4clusters_markers$cluster %<>% factor(levels = paste0("C",1:4))
        markers <- FilterGenes(B_cells_MCL,c("CCND1","CD19","CD5","CDK4","RB1","BTK","SOX11"))
        (MT_gene <- grep("^MT-",X4clusters_markers$gene))
        X4clusters_markers = X4clusters_markers[-MT_gene,]
        Top_n = 40
        top = X4clusters_markers %>% group_by(cluster) %>%
                top_n(Top_n, cluster) %>% top_n(Top_n, avg_logFC)
        unique(top$cluster)
        top = top[order(top$cluster),]
        write.csv(top,paste0(path,"top40_genes_heatmap.csv"))
        features = c(as.character(top$gene),
                     tail(VariableFeatures(object = B_cells_MCL), 2),
                     markers)
        #DoHeatmap.1======
        featuresNum <- make.unique(features, sep = ".")
        exp = AverageExpression(B_cells_MCL[features,], 
                                assays = "SCT") %>% .$SCT
        exp = MakeUniqueGenes(object = exp, features = features)
        exp[tail(VariableFeatures(object = B_cells_MCL), 2),] =0
        scale_exp <- exp %>% t %>% scale %>% t
        colnames(scale_exp) 
        (group.by = unique(top$cluster))
        DoHeatmap.matrix(scale_exp, features = featuresNum,
                         group.by = 1:4,size = 6,angle = 0,
                         draw.lines =F, raster = FALSE,
                         pal_gsea = FALSE,
                         width=1.5, height=10,res=600,no.legend = T,
                         cex.row=5,
                         group.colors = group_colors,
                         do.print = T,
                         unique.name = "cell.types",
                         title = "40 DEGs",
                         save.path = path)
        file.rename(paste0(path,"Heatmap_top40_object_MCL_B_cells_X4clusters_Legend.jpeg"),
                    paste0(path,"Heatmap_top40_Pt10_LN2Pd_X4clusters.jpeg"))
}

#==== Figure ===========
path <- "Yang/Figure Sources/tSNE plots/Groups/"
if(!dir.exists(path)) dir.create(path, recursive = T)

# load data
(load(file="data/MCL_41_harmony_20191231.Rda"))
# reduce normal sample 
object$orig.ident %<>% gsub("N01|N02|N03","Normal",.)
all_normal = which(object$orig.ident %in% "Normal")
rm_normal = sample(all_normal,size = length(all_normal)/3*2, replace = F)
object = object[,-rm_normal]

# order the sample
df_samples <- readxl::read_excel("doc/191120_scRNAseq_info.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
table(object$orig.ident)

(df_samples = df_samples[df_samples$`sample name` %in% object$orig.ident,])
(samples = c("Normal",df_samples$`sample name`[df_samples$`sample name` %in% object$orig.ident]))
object$orig.ident %<>% factor(levels = samples)

# preprocess
Idents(object) = "Doublets"
object %<>% subset(idents = "Singlet")
Idents(object) = "cell.types"
object %<>% subset(idents = c("HSC/progenitors","Nonhematopoietic cells"), invert = TRUE)
table(Idents(object))

object$cell_types <- plyr::mapvalues(object@meta.data$cell.types,
                                     from = c("B_cells","MCL",
                                              "Myeloid cells",
                                              "NK_cells","T_cells:CD4+",
                                              "T_cells:CD8+"),
                                     to = c("B","MCL",
                                            "Monocytes",
                                            "NK","CD4 T",
                                            "CD8 T"))
object$cell_types %<>% as.character()
object$cell_types.colors = object$cell.types.colors


Idents(object) = "groups"
(groups = Idents(object) %>% unique %>% as.character %>% sort)
for(i in 3:length(groups)){
        sub_object <- subset(object, idents = groups[i])
        Idents(sub_object) = "cell_types"
        TSNEPlot.1(object = sub_object, label = F, label.repel = F, group.by = "cell_types",
                   cols = ExtractMetaColor(sub_object),no.legend = F,border = T,
                   pt.size = 1, do.print = T,do.return = F,legend.size = 25,
                   unique.name = "groups",
                   title.size = 20,title = paste("tSNE plot of",groups[i],"samples"),
                   units= "in",width=9, height=7,hjust =0.5, save.path = path)
        Progress(i-2,length(groups)-2)
}
path <- "Yang/Figure Sources/tSNE plots/Groups_split/"
if(!dir.exists(path)) dir.create(path, recursive = T)

Idents(object) = "groups"
for(i in 3:length(groups)){
        sub_object <- subset(object, idents = c("Normal",groups[i]))
        Idents(sub_object) = "cell_types"
        TSNEPlot.1(object = sub_object, label = F, label.repel = F, 
                   group.by = "cell_types",split.by = "orig.ident",
                   cols = ExtractMetaColor(sub_object),no.legend = T,border = T,
                   pt.size = 0.2, do.print = T,do.return = F,legend.size = 12,
                   unique.name = "groups",
                   title.size = 14,title = paste("tSNE plot of Normal and",groups[i],"samples"),
                   units= "in",
                   width=length(unique(sub_object$orig.ident))*2.5,
                   height = 3,hjust =0.5, save.path = path)
        Progress(i-2,length(groups)-2)
}
