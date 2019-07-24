########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(magrittr)
library(gplots)
library(cowplot)
source("../R/Seurat3_functions.R")
source("R/util.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#3.1  Compare DE across all major cell types==================
#We would need the data for all clusters, as well the subclusters.
#detect changes in gene expression between young and aged, 
#in the different cell types and subtypes. 
#It will also be interesting to check if there is some subtype enriched in young compared to aged or viceversa. 

# 3.1.1 load data
# Rename ident
(load(file="data/MCL_V3_Harmony_43_20190627.Rda"))
# B cells only ================
Idents(object) <-  "Doublets"
object %<>% subset(idents = "Singlet")
Idents(object) <-  "manual"
object <- sortIdent(object)
TSNEPlot.1(object,label = F,no.legend=F,cols = ExtractMetaColor(object),repel = T,
           do.print = F)
Idents(object) <-  "res.0.6"
table(Idents(object))
TSNEPlot.1(object,label = T,no.legend=F,repel = T, do.print = F)
B_cells_MCL <- subset(object, idents = c(0,1,5,6,9,13,14,18,19,20))
Idents(B_cells_MCL) <-  "singler1main"
B_cells_MCL <- subset(B_cells_MCL, idents = c("B_cells","MCL"))
table(B_cells_MCL@meta.data$singler1main) %>% as.data.frame %>%
        .[.[,"Freq"] >0,]
Idents(B_cells_MCL) <-  "manual"
table(Idents(B_cells_MCL))
B_cells_MCL <- subset(B_cells_MCL, idents = c("B_cells","MCL"))

B_cells_MCL@meta.data$manual %<>% droplevels
B_cells_MCL@meta.data$manual.colors %<>% droplevels
B_cells_MCL@meta.data$manual.colors %<>% factor(levels = c("#E6AB02","#2055da"))
object <- sortIdent(object)
TSNEPlot.1(object = B_cells_MCL, label = F, group.by = "ident",
           do.return = TRUE, no.legend = F, 
           cols = ExtractMetaColor(B_cells_MCL),
           pt.size = 0.2,label.size = 3,title = "Tsne plot of all B and MCL cells")

##############################
# re-scale, PCA, harmony and tsne
##############################
#B_cells_MCL <- NormalizeData(B_cells_MCL)
#B_cells_MCL <- FindVariableGenes(object = B_cells_MCL, mean.function = ExpMean, 
#                            dispersion.function = LogVMR, 
#                            x.low.cutoff = 0.1, x.high.cutoff = 8, y.cutoff = 0.5)
#length(B_cells_MCL@var.genes)

#B_cells_MCL %<>% ScaleData
#B_cells_MCL %<>% RunPCA(pc.genes = B_cells_MCL@var.genes, pcs.compute = 100, do.print = F)

#pcs = 1:75
#jpeg(paste0(path,"S1_RunHarmony_B.jpeg"), units="in", width=10, height=7,res=600)
#system.time(B_cells_MCL %<>% RunHarmony("orig.ident", dims.use = pcs,
#                                   theta = 2, plot_convergence = TRUE,
#                                   nclust = 50, max.iter.cluster = 100))
#dev.off()

#system.time(
#        B_cells_MCL <- RunTSNE(B_cells_MCL, reduction.use = "harmony", dims.use = 1:75,
#                                  perplexity = 30, do.fast = TRUE))
#npcs =75
#B_cells_MCL <- FindNeighbors(B_cells_MCL, reduction = "harmony",dims = 1:npcs) %>%
#        FindClusters(reduction = "harmony",resolution = 0.8,
#                     dims.use = 1:npcs,print.output = FALSE)
#Idents(B_cells_MCL) <- "RNA_snn_res.0.8"
#B_cells_MCL <- sortIdent(B_cells_MCL,numeric = T)
#TSNEPlot.1(B_cells_MCL, pt.size = 1,label = T,do.print = F,
#           label.size = 4, repel = T,title = "all clusters")
save(B_cells_MCL, file = "data/B_cells_MCL_43_20190713.Rda")

(load("data/B_cells_MCL_43_20190713.Rda"))
B_cells_MCL@reductions$tsne@cell.embeddings %<>% .[colnames(B_cells_MCL),]
meta.data = cbind.data.frame(B_cells_MCL@meta.data, B_cells_MCL@reductions$tsne@cell.embeddings)
meta.data[(meta.data$tsne_1 > -10) & (meta.data$res.0.6 == 9),"res.0.6"]=21
meta.data[(meta.data$tsne_1 < -17.5) & (meta.data$res.0.6 == 9),"res.0.6"]=21
B_cells_MCL@meta.data = meta.data
Idents(B_cells_MCL) <- "res.0.6"
B_cells_MCL %<>% subset(idents = 21, invert = TRUE)
B_cells_MCL %<>% sortIdent(numeric = T)
g1 <- TSNEPlot(B_cells_MCL, pt.size = 1,label = T,do.print = F,
                 label.size = 4, repel = T)

#B_cells_MCL@meta.data[(B_cells_MCL$tSNE_2 < -15 & B_cells_MCL$RNA_snn_res.0.8 %in% 5),"res.0.6"]=22
B_cells_MCL@meta.data$res.0.6 = as.numeric(as.character(B_cells_MCL@meta.data$res.0.6))
B_cells_MCL@meta.data[(B_cells_MCL$tsne_1 < -3 & B_cells_MCL$res.0.6 %in% 5),"res.0.6"]=6
B_cells_MCL@meta.data[(B_cells_MCL$tsne_2 > -10 & B_cells_MCL$res.0.6 %in% 6),"res.0.6"]=2
B_cells_MCL@meta.data$X5_clusters <- plyr::mapvalues(x = B_cells_MCL@meta.data$res.0.6,
                                                     from = c(0,1,5,6,9,13,14,18,19,20),
                                                     to =   c(1,2,3,4,5,2, 1, 1, 2, 1))
Idents(B_cells_MCL) <- "X5_clusters"
B_cells_MCL %<>% sortIdent
g2 <- TSNEPlot.1(B_cells_MCL, pt.size = 1,label = T,do.print = F,
                 label.size = 4, repel = T)
jpeg(paste0(path,"rename_tsne_B.jpeg"), units="in", width=10, height=7,res=600)
cowplot::plot_grid(g1+NoLegend(),g2+NoLegend()) +  ggtitle("Rename the clusters")+
        theme(text = element_text(size=15),
              plot.title = element_text(hjust = 0.5))
dev.off()

TSNEPlot.1(B_cells_MCL, pt.size = 0.2,label = F,do.print = T,no.legend = F,
           label.size = 4, repel = T, title = "5 clusters in B/MCL cells")
save(B_cells_MCL, file = "data/B_cells_MCL_43_20190713.Rda")
(load(file = "data/B_cells_MCL_43_20190713.Rda"))


##############
# DE genes between Clusters 5 cluster top 50
###############
Idents(B_cells_MCL) <- "X5_clusters"
B_cells_MCL <- sortIdent(B_cells_MCL,numeric = T)
table(Idents(B_cells_MCL))
X5_clusters_markers <- FindAllMarkers.UMI(B_cells_MCL,logfc.threshold = 0.1,only.pos = FALSE, 
                                          min.pct = 0.1,return.thresh = 0.05)
X5_clusters_markers <- FindAllMarkers.UMI(B_cells_MCL,logfc.threshold = 0.2,only.pos = T, 
                                          min.pct = 0.1,return.thresh = 0.01)
write.csv(X5_clusters_markers,paste0(path,"X5_clusters_FC0.2_markers.csv"))

X5_clusters_markers = read.csv("output/20190717/X5_clusters_FC0.2_markers.csv",row.names = 1)
#B_cells_MCL %<>% ScaleData()

markers <- FilterGenes(B_cells_MCL,c("CCND1","CDK4","RB1","CD19","SOX11","CD5","CD3D","CD8A"))
(MT_gene <- grep("^MT-",X5_clusters_markers$gene))
X5_clusters_markers = X5_clusters_markers[-MT_gene,]
Top_n = 40
top = X5_clusters_markers %>% group_by(cluster) %>% top_n(Top_n, avg_logFC)
B_cells_MCL %<>% ScaleData(features=unique(c(as.character(top$gene),markers)))
DoHeatmap.1(B_cells_MCL, marker_df = X5_clusters_markers, add.genes = markers, Top_n = Top_n, 
            do.print=T, angle = 0, group.bar = F, title.size = 20, no.legend = F,size=5,hjust = 0.5,
            group.bar.height = 0, label=F, cex.row= 2, legend.size = 0,width=10, height=6.5,
            pal_gsea = FALSE,
            title = paste("Top",Top_n,"differentially expressed genes in B and MCL clusters"))

# remove cluster with less than 3 cells======
# no scale down, keep the normal cells.
B_cells_MCL@assays$RNA@scale.data = matrix(0,0,0);GC()
B_cells_MCL@meta.data$old.ident = B_cells_MCL@meta.data$orig.ident
B_cells_MCL@meta.data$orig.ident = gsub("BH|DJ|MD|NZ","Normal",B_cells_MCL@meta.data$orig.ident)
B_cells_MCL@meta.data$X5_orig.ident = paste(B_cells_MCL@meta.data$orig.ident,
                                            B_cells_MCL@meta.data$X5_clusters, sep = "_")
B_cells_MCL@meta.data$X5_orig.ident = gsub('^Normal_.*', 'Normal', B_cells_MCL@meta.data$X5_orig.ident)

table(B_cells_MCL@meta.data$X5_orig.ident)
Idents(B_cells_MCL) <- "X5_orig.ident"

# add color to cluster
B_cells_MCL <- AddMetaColor(object = B_cells_MCL, label= "X5_clusters", 
                            colors = scales::hue_pal()(length(unique(B_cells_MCL$X5_clusters))))

Idents(B_cells_MCL)  = "orig.ident"
df_samples <- readxl::read_excel("doc/190626_scRNAseq_info.xlsx")
colnames(df_samples) <- tolower(colnames(df_samples))
(keep = df_samples$sample %in% unique(B_cells_MCL$orig.ident))
df_samples = df_samples[keep,]

tests <- paste0("test",c(2:12))
for(test in tests){
        sample_n = which(df_samples$tests %in% test)
        df <- as.data.frame(df_samples[sample_n,])
        samples <- unique(df$sample)
        rownames(df) = samples
        
        samples <- c(ifelse(length(samples)>5,NA,"Normal"),df$sample[order(df$tsne)])
        print(samples <- samples[!is.na(samples)])
        
        subset.MCL <- subset(B_cells_MCL,idents = samples)
        subset.MCL$orig.ident = factor(subset.MCL$orig.ident, levels = samples)
        Idents(subset.MCL) = "X5_clusters"
        subset.MCL %<>% sortIdent()
        TSNEPlot.1(subset.MCL, split.by = "orig.ident",pt.size = 1,label = F,
                   cols = ExtractMetaColor(subset.MCL),do.print = T, unique.name = T)
}

for(test in tests){
        sample_n = which(df_samples$tests %in% test)
        df <- as.data.frame(df_samples[sample_n,])
        samples <- unique(df$sample)
        rownames(df) = samples
        
        samples <- c(ifelse(length(samples)>5,NA,"Normal"),df$sample[order(df$tsne)])
        print(samples <- samples[!is.na(samples)])
        
        subset.MCL <- subset(B_cells_MCL,idents = samples)
        subset.MCL$orig.ident = factor(subset.MCL$orig.ident, levels = samples)
        Idents(subset.MCL) = "cell.type"
        subset.MCL <- sortIdent(subset.MCL)
        TSNEPlot.1(subset.MCL, split.by = "orig.ident",pt.size = 1,label = F,
                   cols = ExtractMetaColor(subset.MCL),do.print = T, unique.name = T)
}

table_B_cells_MCL <- table(B_cells_MCL@meta.data$X5_orig.ident) %>% as.data.frame
(remove <- table_B_cells_MCL[table_B_cells_MCL$Freq < 2,"Var1"] %>% as.character())
Idents(B_cells_MCL) = "X5_orig.ident"
B_cells_MCL <- subset(B_cells_MCL, idents = remove, invert = TRUE)
B_cells_MCL_exp <- AverageExpression(B_cells_MCL)
write.csv(B_cells_MCL_exp,paste0(path,"B_MCL_exp.csv"))

###############################
# Doheatmap for Normal / MCL
###############################
df_samples <- readxl::read_excel("doc/190626_scRNAseq_info.xlsx",sheet = "heatmap")
list_samples <- df2list(df_samples)
all(list_samples %>% unlist %>% as.vector %>% unique %in% 
              B_cells_MCL@meta.data$orig.ident)

Idents(B_cells_MCL) = "orig.ident"
markers <- FilterGenes(B_cells_MCL,c("CCND1","CDK4","RB1","CD19","SOX11","CD5","CD3D","CD8A"))
(block <- VariableFeatures(B_cells_MCL) %>% tail(1))
for(sample in list_samples$MCL){
        subset.MCL <- subset(B_cells_MCL, idents = c("Normal",sample))

        # SplitTSNEPlot======
        Idents(subset.MCL) = "X5_orig.ident"
        subset.MCL %<>% sortIdent()
        TSNEPlot.1(subset.MCL, split.by = "orig.ident",pt.size = 1,label = F,
                   do.print = F, unique.name = T)
        
        # remove cluster with less than 3 cells======
        table_subset.MCL <- table(subset.MCL$X5_orig.ident) %>% as.data.frame
        keep.MCL <- table_subset.MCL[table_subset.MCL$Freq > 2,"Var1"] %>% as.character()
        
        (X5_cluster <- keep.MCL %>% unique %>% 
                gsub('.*\\_',"",.) %>% sub("Normal","",.) %>% as.numeric %>% sort )
        
        subset.MCL <- subset(subset.MCL, idents = keep.MCL)
        print(ident.1 <- rep("Normal",length(X5_cluster)))
        print(ident.2 <- paste(sample,X5_cluster,sep="_"))

        # FindAllMarkers.UMI======

        #gde.markers <- FindPairMarkers(subset.MCL, ident.1 = c(ident.1,ident.2),
        #                               ident.2 = c(ident.2,ident.1), only.pos = T,
        #                               logfc.threshold = 0.1,min.cells.group =3,
        #                               min.pct = 0.1,return.thresh = 0.05,
        #                               save.files = FALSE)
        #write.csv(gde.markers, paste0(path,"Normal_vs_",sample,".csv"))
        gde.markers = read.csv(paste0(path,"Normal_vs_",sample,".csv"),row.names = 1)
        (mito.genes <- grep(pattern = "^MT-", x = gde.markers$gene))
        if(length(mito.genes)>0) gde.markers = gde.markers[-mito.genes,]
        GC()
        #DoHeatmap.1======
        Top_n = 40
        top <-  gde.markers %>% group_by(cluster1.vs.cluster2) %>% 
                top_n(Top_n, avg_logFC) %>% as.data.frame()
        add.genes = unique(c(as.character(top$gene),block,markers))
        subset.MCL %<>% ScaleData(features= add.genes)
        DoHeatmap.1(subset.MCL, add.genes = add.genes,
                    Top_n = Top_n, do.print=T, angle = 90,
                    group.bar = F, title.size = 20, no.legend = F,size=5,hjust = 0.5,
                    group.bar.height = 0,label=T, cex.row= 500/length(add.genes), legend.size = 0,
                    width=10, height=6.5,unique.name = T, pal_gsea = F,
                    title = paste("Top",Top_n,"DE genes in",sample,"/Normal B and MCL cells"))
        GC()
}


# Doheatmap for MCL.1 / MCL.2 ================
df_samples <- readxl::read_excel("doc/190626_scRNAseq_info.xlsx",sheet = "heatmap")
list_samples <- df2list(df_samples)
print(list_samples %>% unlist %>% as.vector %>% unique %in% 
              B_cells_MCL@meta.data$orig.ident)

Idents(B_cells_MCL) = "orig.ident"
for(i in 1:length(list_samples$MCL.1)){
        
        (samples1 = list_samples$MCL.1[i])
        (samples2 = list_samples$MCL.2[i])
        
        subset.MCL <- subset(B_cells_MCL, idents = c(samples1,samples2))
        subset.MCL@meta.data$orig.ident %<>% factor(levels = c(samples1,samples2))
        # remove cluster with less than 3 cells======
        
        table_subset.MCL <- table(subset.MCL@meta.data$X5_orig.ident) %>% as.data.frame
        (keep.MCL <- table_subset.MCL[table_subset.MCL$Freq > 2,"Var1"] %>% as.character())
        (X5_cluster <- keep.MCL %>% unique %>% 
                        gsub('.*\\_',"",.) %>% as.numeric %>% sort %>% .[duplicated(.)])
        
        print(ident.1 <- paste(samples1,X5_cluster,sep="_"))
        print(ident.2 <- paste(samples2,X5_cluster,sep="_"))
        
        #---SplitTSNEPlot----
        Idents(subset.MCL) = "X5_orig.ident"
        subset.MCL <- subset(subset.MCL, idents = keep.MCL)
        
        Idents(subset.MCL) %<>% factor(levels = c(ident.1,ident.2))
        TSNEPlot.1(subset.MCL, split.by = "orig.ident",pt.size = 1,label = F,do.return = F,
                   do.print = F, unique.name = T)

        #gde.markers <- FindPairMarkers(subset.MCL, ident.1 = c(ident.1,ident.2), 
        #                               ident.2 = c(ident.2,ident.1), only.pos = T,
        #                              logfc.threshold = 0.1,min.cells.group =3,
        #                               min.pct = 0.1,return.thresh = 0.05,
        #                               save.files = FALSE)
        #write.csv(gde.markers, paste0(path,samples1,"_vs_",samples2,".csv"))
        gde.markers = read.csv(paste0(path,samples1,"_vs_",samples2,".csv"),row.names = 1)
        print(table(gde.markers$cluster1.vs.cluster2))
        (mito.genes <- grep(pattern = "^MT-", x = gde.markers$gene))
        if(length(mito.genes)>0) gde.markers = gde.markers[-mito.genes,]
        GC()
        #DoHeatmap.1======
        Top_n = 40
        top <-  gde.markers %>% group_by(cluster1.vs.cluster2) %>% 
                top_n(Top_n, avg_logFC) %>% as.data.frame()
        add.genes = unique(c(as.character(top$gene),block,markers))
        subset.MCL %<>% ScaleData(features= add.genes)
        DoHeatmap.1(subset.MCL, add.genes = add.genes,
                    Top_n = Top_n, do.print=T, pal_gsea = F,
                    group.bar = F, title.size = 20, no.legend = F,size=5,hjust = 0.5,
                    group.bar.height = 0,label=T, cex.row=500/length(add.genes), legend.size = 0,
                    width=10, height=6.5,unique.name = T,
                    title = paste0("Top ",Top_n," DE genes in ",samples2," /",samples1, " B and MCL cells"))
        GC()
}

# Doheatmap for MCL longitudinal by different clusters ================

df_samples <- readxl::read_excel("doc/190626_scRNAseq_info.xlsx")
colnames(df_samples) =  tolower(colnames(df_samples))
df_samples[df_samples$sample %in% "MD","tsne"] = 0
df_samples[df_samples$sample %in% "MD","sample"] = "Normal"

Idents(B_cells_MCL) = "groups"
groups = c("Untreated","Pt-11","Pt-17","AFT-03","AFT-04","Pt-AA13","Pt-25","Pt-27")
clusters = list(c(1,2))#,3,c(4,6))
groups = c("Untreated","Pt-17","Pt-25")
for(i in 1:length(groups)){
        
        subset_MCL <- subset(B_cells_MCL, idents = c("Normal",groups[i]))
        Idents(subset_MCL) = "X5_clusters"
        for(k in 1:length(clusters)){
                subset.MCL <- subset(subset_MCL,idents = clusters[[k]])
                
                (samples = unique(subset.MCL$orig.ident))
                df = df_samples[df_samples$sample %in% samples,]
                subset.MCL@meta.data$orig.ident = factor(subset.MCL@meta.data$orig.ident, 
                                                         levels = df$sample[order(df$tsne)])
                Idents(subset.MCL) = "orig.ident"
                gde.markers <- FindAllMarkers.UMI(subset.MCL,logfc.threshold = 0.3, only.pos = T,
                                                  test.use = "MAST")
                write.csv(gde.markers,paste0(path,"B/B_MCL_DE/B_",groups[i],"_",
                                             paste(clusters[[k]],collapse = "_"),".csv"))
                (mito.genes <- grep(pattern = "^MT-", x = gde.markers$gene))
                if(length(mito.genes)>0) gde.markers = gde.markers[-mito.genes,]
                GC()
                #DoHeatmap.1======
                subset.MCL %<>% ScaleData(features= unique(gde.markers$gene))
                Top_n = 20
                DoHeatmap.1(subset.MCL, marker_df = gde.markers, Top_n = Top_n, do.print=T, angle = 0,
                            group.bar = T, title.size = 20, no.legend = F,size=5,hjust = 0.5,
                            label=F, cex.row=4, legend.size = NULL,width=10, height=7,unique.name = T,
                            title = paste("Top",Top_n,"DE genes in longitudinal",groups[i],
                                          "B/MCL cells cluster",paste(clusters[[k]],collapse = " and ")))
                # rename file
                v <- paste(unique(subset.MCL$orig.ident),collapse = "_")
                v = paste0(v,"_",FindIdentLabel(subset.MCL))
                old.name = paste0(path,"Doheatmap_top",Top_n,"_",v,".jpeg")
                file.rename(old.name, paste0(path,"Doheatmap_top",Top_n,"_",v,"cluster",
                                             paste(clusters[[k]],collapse = "_"),".jpeg"))
                remove(subset.MCL);GC()
        }
        remove(subset_MCL);GC()
}


# Doheatmap for MCL longitudinal ================
df_samples <- readxl::read_excel("doc/190626_scRNAseq_info.xlsx")
colnames(df_samples) =  tolower(colnames(df_samples))
df_samples[df_samples$sample %in% "MD","tsne"] = 0
df_samples[df_samples$sample %in% "MD","sample"] = "Normal"

Idents(B_cells_MCL) = "groups"
#groups = c("Untreated","Pt-11","Pt-17","AFT-03","AFT-04","Pt-AA13","Pt-25","Pt-27")
groups = c("Untreated","Pt-17","Pt-25")
for(i in 1:length(groups)){
        
        subset.MCL <- subset(B_cells_MCL, idents = c("Normal",groups[i]))
        
        (samples = unique(subset.MCL$orig.ident))
        df = df_samples[df_samples$sample %in% samples,]
        subset.MCL@meta.data$orig.ident = factor(subset.MCL@meta.data$orig.ident, 
                                                 levels = df$sample[order(df$tsne)])
        Idents(subset.MCL) %<>% factor()
        Idents(subset.MCL) = "orig.ident"
        samples = samples[-which(samples %in% "Normal")]
        gde.markers_list <- list()
        for(k in 1:length(samples)){
                gde.markers_list[[k]] <- FindMarkers.UMI(subset.MCL,
                                                         ident.1 = samples[k],
                                                         ident.2 = "Normal",
                                                         logfc.threshold = 0.1, only.pos = F,
                                                         test.use = "MAST")
                gde.markers_list[[k]]$cluster <- samples[k]
                gde.markers_list[[k]]$gene <- rownames(x = gde.markers_list[[k]])
        }
        gde.markers <- do.call(rbind, gde.markers_list)
        write.csv(gde.markers,paste0(path,"B/B_MCL_DE/B_",groups[i],".csv"))
        
        gde.markers = read.csv(file = paste0("output/20190622/B/B_MCL_DE/B_",groups[i],".csv"))
        (mito.genes <- grep(pattern = "^MT-", x = gde.markers$gene))
        if(length(mito.genes)>0) gde.markers = gde.markers[-mito.genes,]
        GC()
        #DoHeatmap.1======
        subset.MCL %<>% ScaleData(features= unique(gde.markers$gene))
        Top_n = 20
        DoHeatmap.1(subset.MCL, marker_df = gde.markers, Top_n = Top_n, do.print=T, angle = 0,
                    group.bar = T, title.size = 20, no.legend = F,size=5,hjust = 0.5,
                    label=T, cex.row=5, legend.size = NULL,width=10, height=7,unique.name = T,
                    title = paste("Top",Top_n,"DE genes of",groups[i],
                                  "MCL cells over Normal B cells"))
}

#######################
# heatmap.2
#######################
(load(file = "data/B_cells_MCL_43_20190615.Rda"))

markers <- c("BTK","CCND1","CCND2","CCND3","CD3D","CD5","CD8A","CDK4","IL2RA",
             "MS4A1","PDCD1","RB1","SOX11","CD19","CD79A","CD40")
markers <-  FilterGenes(B_cells_MCL,markers)
B_cells_MCL <- FindVariableFeatures(object = B_cells_MCL, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top20 <- head(x = VariableFeatures(object = B_cells_MCL), 20)

# plot variable features with labels
plot1 <- VariableFeaturePlot(object = B_cells_MCL)
plot1 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
jpeg(paste0(path,"VariableFeaturePlot_B_MCL.jpeg"), units="in", width=10, height=7,res=600)
print(plot1+  ggtitle("Variable genes in B and MCL cells")+
              theme(text = element_text(size=15),
                    plot.title = element_text(hjust = 0.5)))
dev.off()

## Column clustering (adjust here distance/linkage methods to what you need!)
B_cells_MCL %<>% ScaleData(features = VariableFeatures(B_cells_MCL))
y = B_cells_MCL@assays$RNA@scale.data[VariableFeatures(B_cells_MCL),]
#write.csv(exp,paste0(path,"exp1.csv"))
exp = read.csv(paste0(path,"exp.csv"),row.names = 1)
remove(B_cells_MCL);GC()
system.time(cor.mtr <- HiClimR::fastCor(y,nSplit = 3, optBLAS = TRUE))
system.time(distance <- as.dist(1-cor.mtr))
system.time(hc <- hclust(distance,method="ward.D2"))
save(hc, file = paste0("output/hc.Rda"))
load(file = paste0("output/hc.Rda"))

clusters = B_cells_MCL@meta.data[hc$labels,"orig.ident"]

cc <- plyr::mapvalues(x = clusters,
                      from = unique(clusters),
                      to = scales::hue_pal()(length(unique(clusters))))
clusters =factor(clusters, levels = sort(unique(clusters)))
clusters =factor(clusters, levels = 1:length(unique(clusters)))
jpeg(paste0(path,"Heatmap2_MCL.jpeg"), units="in", width=10, height=7,res=600)
heatmap.2(as.matrix(y),
          Colv = as.dendrogram(hc), 
          Rowv= FALSE,
          ColSideColors = cc, 
          trace ="none",labCol = FALSE,dendrogram = "column",
          key.xlab = "scale log nUMI",
          cexRow = 0.5,
          margins = c(2,5),
          #scale = "row",
          breaks = seq(-3,3,length.out = 101),
          col = bluered,
          main = "All B and MCL cells")
par(lend = 1)           # square line ends for the color legend
legend(-0.1, 0.8,       # location of the legend on the heatmap plot
       legend = unique(clusters), # category labels
       col = scales::hue_pal()(length(unique(clusters))),  # color key
       lty= 1,             # line style
       lwd = 10            # line width
)
dev.off()


# heatmap.2 for MCL longitudinal ================
df_samples <- readxl::read_excel("doc/190429_scRNAseq_info.xlsx")
colnames(df_samples) =  tolower(colnames(df_samples))

Idents(B_cells_MCL) = "groups"
groups = c("Pt-11","Pt-17","AFT-03","AFT-04","Pt-AA13","Pt-25","Pt-27")
for(i in 2:length(groups)){
        
        subset.MCL <- subset(B_cells_MCL, idents = groups[i])

        (samples = unique(subset.MCL$orig.ident))
        df = df_samples[df_samples$sample %in% samples,]
        
        subset.MCL <- FindVariableFeatures(object = subset.MCL, selection.method = "vst", 
                                            nfeatures = 2000)
        top20 <- head(x = VariableFeatures(object = subset.MCL), 20)
        
        # plot variable features with labels
        plot1 <- VariableFeaturePlot(object = subset.MCL)
        plot1 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
        jpeg(paste0(path,"VariableFeaturePlot_B_MCL",groups[i],".jpeg"), units="in", width=10, height=7,res=600)
        print(plot1+  ggtitle(paste("Variable genes in",groups[i],"B and MCL cells"))+
                      theme(text = element_text(size=15),
                            plot.title = element_text(hjust = 0.5)))
        dev.off()
        ## Column clustering (adjust here distance/linkage methods to what you need!)
        subset.MCL %<>% ScaleData(features = VariableFeatures(subset.MCL))
        y = subset.MCL@assays$RNA@scale.data[VariableFeatures(subset.MCL),]
        system.time(cor.mtr <- HiClimR::fastCor(y,nSplit = 2, optBLAS = TRUE))
        distance <- as.dist(1-cor.mtr)
        hc <- hclust(distance,method="ward.D2")
        clusters = subset.MCL@meta.data[hc$labels,c("orig.ident","X5_clusters.colors")]
        #clusters = subset.MCL@meta.data[hc$labels,c("orig.ident")]
        
        cc <- as.character(clusters$X5_clusters.colors)
        
        #cc <- plyr::mapvalues(x = clusters,
        #                      from = unique(clusters),
        #                      to = scales::hue_pal()(length(unique(clusters))))
        jpeg(paste0(path,"Heatmap2_MCL_",groups[i],".jpeg"), units="in", width=10, height=7,res=600)
        heatmap.2(as.matrix(y),
                  Colv = as.dendrogram(hc), 
                  Rowv= FALSE,
                  ColSideColors = cc, 
                  trace ="none",labCol = FALSE,dendrogram = "column",
                  key.xlab = "scale log nUMI",
                  cexRow = 0.5,
                  margins = c(2,5),
                  #scale = "row",
                  breaks = seq(-3,3,length.out = 101),
                  col = bluered,
                  main = paste("All B and MCL cells in",groups[i]))
        par(lend = 1)           # square line ends for the color legend
        legend(-0.01, 0.8,       # location of the legend on the heatmap plot
               legend = sort(unique(clusters)), # category labels
               col = scales::hue_pal()(length(unique(cc))),  # color key
               lty= 1,             # line style
               lwd = 10            # line width
        )
        dev.off()
}

        
# heatmap.2 for Normal / MCL ================

control <- "MD"
tests <- c("test3","test4")
for(test in tests){
        sample_n = which(df_samples$tests %in% test)
        samples <- df_samples$samples[sample_n]
        cell.use <- rownames(B_cells_MCL@meta.data)[B_cells_MCL@meta.data$orig.ident %in% c(control,sample)] #
        subset.MCL <- SubsetData(B_cells_MCL, cells.use = cell.use)
        test_markers <- FindAllMarkers.UMI(subset.MCL,logfc.threshold = 0.1, only.pos = T,
                                           test.use = "MAST")
        top <- test_markers %>% group_by(cluster) %>% top_n(50, avg_logFC)
        y = subset.MCL@scale.data[unique(c(markers,top$gene)),]
        ## Column clustering (adjust here distance/linkage methods to what you need!)
        hc <- hclust(as.dist(1-cor(as.matrix(y), method="spearman")), method="complete")
        cc = gsub("_.*","",hc$labels)
        cc = gsub(control,"#B3DE69",cc)
        cc = gsub(sample,"#195016",cc)
        
        jpeg(paste0(path,"/Heatmap2_",control,"_",sample,".jpeg"), units="in", width=10, height=7,res=600)
        heatmap.2(as.matrix(y),
                  Colv = as.dendrogram(hc), Rowv= FALSE,
                  ColSideColors = cc, trace ="none",labCol = FALSE,dendrogram = "column",#scale="row",
                  key.xlab = "scale log nUMI",
                  cexRow = 0.5,
                  margins = c(2,5),
                  #scale = "row",
                  breaks = seq(-3,3,length.out = 101),
                  col = bluered,
                  main = paste(control, "vs.",sample, "in B cells and MCL cells"))
        par(lend = 1)           # square line ends for the color legend
        legend(0, 0.8,       # location of the legend on the heatmap plot
               legend = c(control, sample), # category labels
               col = c("#B3DE69", "#195016"),  # color key
               lty= 1,             # line style
               lwd = 10            # line width
        )
        dev.off()
}


# heatmap.2 ================
markers <-  HumanGenes(B_cells_MCL,c("CD19","MS4A1","CD79A","CD5","CD40","CDK4"))
tests <- c("test3","test4")
for(test in tests){
        sample_n = which(df_samples$tests %in% test)
        df_samples[sample_n,] %>% kable() %>% kable_styling()
        samples <- df_samples$samples[sample_n]
        cell.use <- rownames(B_cells_MCL@meta.data)[B_cells_MCL@meta.data$orig.ident %in% samples]
        subset.MCL <- SubsetData(B_cells_MCL, cells.use = cell.use)
        subset.MCL <- SetAllIdent(subset.MCL,id = "orig.ident")
        test_markers <- FindAllMarkers.UMI(subset.MCL,logfc.threshold = 0.1, only.pos = T,
                                           test.use = "MAST")
        top <- test_markers %>% group_by(cluster) %>% top_n(30, avg_logFC)
        y = subset.MCL@scale.data[unique(c(markers,top$gene)),]
        ## Column clustering (adjust here distance/linkage methods to what you need!)
        hc <- hclust(as.dist(1-cor(as.matrix(y), method="spearman")), method="complete")
        cc = gsub("_.*","",hc$labels)
        cc = gsub(samples[1],"#B3DE69",cc)
        cc = gsub(samples[3],"#E31A1C",cc)
        cc = gsub(samples[2],"#195016",cc)
        
        jpeg(paste0(path,"/Heatmap2_B_MLC.jpeg"), units="in", width=10, height=7,res=600)
        heatmap.2(as.matrix(y),
                  Colv = as.dendrogram(hc), Rowv= FALSE,
                  ColSideColors = cc, trace ="none",labCol = FALSE,dendrogram = "column",
                  key.xlab = "scale log nUMI",
                  cexRow = 0.5,
                  margins = c(2,5),
                  #scale = "row",
                  breaks = seq(-3,3,length.out = 101),
                  col = bluered,
                  main = paste("all samples in B cells"))
        #par(lend = 1)           # square line ends for the color legend
        #legend(0, 0.8,       # location of the legend on the heatmap plot
        #       legend = c(samples[1], samples[2], samples[3]), # category labels
        #       col = c("#B3DE69", "#195016","#E31A1C"),  # color key
        #       lty= 1,             # line style
        #       lwd = 10            # line width
        #)
        dev.off()
}


# heatmap.2 for MCL bulk RNA ================
markers <-  c("CD19","MS4A1","CD79A","CD5","CD40","CDK4","SOX11","ITGA4","CCND1","CCND2")
markers <-  unique(MCL_markers$gene)[1:1000]
markers <- markers[markers %in% rownames(X181120_MCL_WTS)]

y = X181120_MCL_WTS[markers,]
## Column clustering (adjust here distance/linkage methods to what you need!)
hc <- hclust(dist(1-cor(as.matrix(y), method="spearman")), method="complete")

jpeg(paste0(path,"/Heatmap2_MCL_bulk_10.jpeg"), units="in", width=10, height=7,res=600)
heatmap.2(as.matrix(y),
          Colv = as.dendrogram(hc),
          trace ="none",dendrogram = "both",
          cexRow = 1,
          margins = c(5,5),
          scale = "row",
          col = bluered,
          main = paste("MCL bulk RNA-seq with marker genes"))
dev.off()

# heatmap.2 for MCL scRNA================
MCL_markers_list <- list()
cell.use <- rownames(B_cells_MCL@meta.data)[B_cells_MCL@meta.data$orig.ident %in% samples]
subset.MCL <- SubsetData(B_cells_MCL, cells.use = cell.use)
subset.MCL <- SetAllIdent(subset.MCL,id = "orig.ident")
TSNEPlot(subset.MCL)
MCL_exp <- AverageExpression(subset.MCL)

(samples <- c("DJ","MD","Pt-1294","Pt-RM","Pt-MS","Pt-LM","Pt-1475"))
markers <-  unique(MCL_markers$gene)[1:100]
markers <- markers[markers %in% rownames(MCL_exp)]
y = MCL_exp[markers,]

## Column clustering (adjust here distance/linkage methods to what you need!)
#hc <- hclust(dist(1-cor(as.matrix(y), method="spearman")), method="complete")

jpeg(paste0(path,"Heatmap2_MCL_sc_100.jpeg"), units="in", width=10, height=7,res=600)
heatmap.2(as.matrix(y),
          #Colv = as.dendrogram(hc),
          Colv=FALSE,
          trace ="none",dendrogram = "row",
          cexRow = 0.3,
          margins = c(5,5),
          scale = "row",
          col = bluered,
          main = paste("MCL scRNA RNA top 100 genes"))
dev.off()


# heatmap.2 for MCL bulk RNA + scRNA================
All_MCL <- merge(MCL_exp,log1p(X181120_MCL_WTS), by="row.names")
rownames(All_MCL) <-All_MCL$Row.names
All_MCL <- All_MCL[,-1]
column.sum <- colSums(All_MCL)

jpeg(paste0(path,"/All_MCL.jpeg"), units="in", width=10, height=7,res=600)
par(mfrow=c(2,1))
boxplot(All_MCL, ylab= "log UMI or FPKM")
title(main = "merge MCL bulk and RNA")
plot(x = column.sum)
dev.off()

#markers <-  c("CD19","MS4A1","CD79A","CD5","CD40","CDK4","SOX11","ITGA4","CCND1","CCND2")
markers <-  unique(MCL_markers$gene)[1:1000]
markers <- markers[markers %in% rownames(All_MCL)]

y = All_MCL[markers,]
## Column clustering (adjust here distance/linkage methods to what you need!)
hc <- hclust(dist(1-cor(as.matrix(y), method="spearman")), method="complete")

jpeg(paste0(path,"/Heatmap2_MCL_sc_bulk_1000.jpeg"), units="in", width=10, height=7,res=600)
heatmap.2(as.matrix(y),
          Colv = as.dendrogram(hc),
          trace ="none",dendrogram = "both",
          cexRow = 0.03,
          margins = c(5,5),
          scale = "row",
          col = bluered,
          main = paste("MCL bulk and scRNA top 1000 genes"))
dev.off()

# check HSC===========
B_cells_MCL <- SetAllIdent(B_cells_MCL,id='cell.type')
cell.type_markers <- FindAllMarkers(B_cells_MCL,logfc.threshold=0.5)
PROM1_markers <- FindMarkers.UMI(B_cells_MCL,logfc.threshold = 0, genes.use = "PROM1",
                                 min.pct = 0,ident.1 = "MCL/HSC",return.thresh = 1)

write.table(cell.type_markers,paste0(path,"cell.type_markers.txt"))


object@meta.data$cell.type = gsub("B_cells.*","B_cells",object@meta.data$singler1sub)
object@meta.data$cell.type = gsub("T_cells.*","T_cells",object@meta.data$cell.type)
object@meta.data$cell.type = gsub("MPP|MEP|CLP|HSC|CMP|GMP","MCL.HSC",object@meta.data$cell.type)

table_subset.MCL <- table(object@meta.data$cell.type) %>% as.data.frame 
(keep.MCL <- table_subset.MCL[table_subset.MCL$Freq > 2,"Var1"] %>% as.character())
object <- SetAllIdent(object, id="cell.type")
table(object@ident)
major_cell <- SubsetData(object,ident.use=keep.MCL)
major_cell_exp <- AverageExpression(major_cell)
write.csv(major_cell_raw_exp, paste0(path,"major_cell_raw_exp.csv"))



PROM1_markers <- FindMarkers.UMI(object,logfc.threshold = 0, genes.use = "PROM1",
                                 min.pct = 0,ident.1 = "MCL.HSC",return.thresh = 1)

