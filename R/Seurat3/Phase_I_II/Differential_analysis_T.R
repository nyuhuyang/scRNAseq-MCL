########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(SingleR)
library(dplyr)
library(plyr)
library(tidyr)
library(kableExtra)
library(magrittr)
library(harmony)
library(gplots)
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
Idents(object) <-  "Doublets"
object %<>% subset(idents = "Singlet")
# T cells only ================

Idents(object) <-  "manual"
object <- sortIdent(object)
table(Idents(object))
TSNEPlot.1(object,label = F, repel = T, no.legend = T,pt.size = 0.1,
           cols = ExtractMetaColor(object),do.return = T,do.print = F,
           title = "Cell type labeling by Blueprint + Encode + MCL")

Idents(object) <-  "res.0.6"
object <- sortIdent(object,numeric=T)
table(Idents(object))
TSNEPlot.1(object,label = T, repel = T, no.legend = F,pt.size = 0.1,
           do.return = T,do.print = F,
           title = "Unsupervised clustering")

T_NK_cells <- subset(object,  idents = c(2,3,4,10,12,15,16,17))
Idents(T_NK_cells) = "singler1main"
table(Idents(T_NK_cells))
T_NK_cells <- subset(T_NK_cells, idents = c("NK_cells","T_cells"))
table(T_NK_cells@meta.data$manual) %>% as.data.frame %>%
        .[.[,"Freq"] >0,]
Idents(T_NK_cells) = "manual"
T_NK_cells <- subset(T_NK_cells, idents = c("T_cells:CD4+","T_cells:CD8+","NK_cells"))
T_NK_cells$manual %<>% factor(levels = c("T_cells:CD4+","T_cells:CD8+","NK_cells"))
T_NK_cells$manual.colors %<>% factor(levels = c("#B3DE69","#F0027F","#A65628"))


T_NK_cells[['tSNE_1']] <- T_NK_cells@reductions$tsne@cell.embeddings[,"tsne_1"]

T_NK_cells <- subset(T_NK_cells, subset = tSNE_1 >5)

g0 <- TSNEPlot.1(T_NK_cells,do.label = F,no.legend=T,do.print = T,do.return =T, pt.size = 0.1,
           cols = ExtractMetaColor(T_NK_cells),
           title="T cells and NK cells")
##############################
# re-scale, PCA, harmony and tsne
##############################

#T_NK_cells <- NormalizeData(T_NK_cells)
#T_NK_cells <- FindVariableGenes(object = T_NK_cells, mean.function = ExpMean, 
#                            dispersion.function = LogVMR, 
#                            x.low.cutoff = 0.1, x.high.cutoff = 8, y.cutoff = 0.5)
#length(T_NK_cells@var.genes)

#T_NK_cells %<>% ScaleData
#T_NK_cells %<>% RunPCA(pc.genes = T_NK_cells@var.genes, pcs.compute = 100, do.print = F)

#pcs = 1:75
#jpeg(paste0(path,"S1_RunHarmony_B.jpeg"), units="in", width=10, height=7,res=600)
#system.time(T_NK_cells %<>% RunHarmony("orig.ident", dims.use = pcs,
#                                   theta = 2, plot_convergence = TRUE,
#                                   nclust = 50, max.iter.cluster = 100))
#dev.off()

#system.time(
#        T_NK_cells <- RunTSNE(T_NK_cells, reduction.use = "harmony", dims.use = 1:75,
#                                  perplexity = 30, do.fast = TRUE))

system.time(T_NK_cells %<>% FindNeighbors(reduction = "harmony",dims = 1:75,force.recalc = T))
system.time(T_NK_cells %<>% FindClusters(resolution = 0.3))
Idents(T_NK_cells) <- 'RNA_snn_res.0.3'
table(Idents(T_NK_cells))
g1 <- TSNEPlot.1(T_NK_cells,label = T,no.legend=T,do.print = T,pt.size = 0.2,
           title="unsupervised clustering of T and NK cells")
T_NK_cells@meta.data$X3_clusters <- plyr::mapvalues(x = T_NK_cells@meta.data$RNA_snn_res.0.3,
                                                     from = c(0,1,2,3,4,5,6,7,8),
                                                     to =   c(1,3,2,2,3,2,1,1,1))
T_NK_cells@meta.data$X3_clusters %<>% factor(levels = 1:3)

T_NK_cells@meta.data$X3_clusters.colors <- mapvalues(x = T_NK_cells@meta.data$X3_clusters,
                                                   from = 1:3,
                                                   to =   c("#B3DE69","#F0027F",
                                                            "#A65628"))
Idents(T_NK_cells) <- 'X3_clusters'
T_NK_cells <- sortIdent(T_NK_cells)
table(Idents(T_NK_cells))
# tsne plot
g2 <- TSNEPlot.1(T_NK_cells,label = F,no.legend=F,do.print = T,pt.size = 0.2,
                 cols = ExtractMetaColor(T_NK_cells),
                 title="unsupervised clustering")
jpeg(paste0(path,"rename_tsne_T_NK.jpeg"), units="in", width=10, height=7,res=600)
cowplot::plot_grid(g0,g2+NoLegend()) +  ggtitle("Rename the clusters")+
        theme(text = element_text(size=15),
              plot.title = element_text(hjust = 0.5,size=15))
dev.off()
TSNEPlot.1(T_NK_cells,label = F,no.legend=F,do.print = T,pt.size = 1,do.return = T,
           title="Merged clusters of T and NK cells")
T_NK_cells_exp <- AverageExpression(T_NK_cells)
write.csv(T_NK_cells_exp,paste0(path,"T_NK_cells_exp.csv"))

save(T_NK_cells, file = "data/T_NK_cells_43_20190717.Rda")
(load(file = "data/T_NK_cells_43_20190717.Rda"))

# Doheatmap for MCL longitudinal ================

df_samples <- readxl::read_excel("doc/190626_scRNAseq_info.xlsx")
colnames(df_samples) =  tolower(colnames(df_samples))
df_samples[df_samples$sample %in% "MD","tsne"] = 0
df_samples[df_samples$sample %in% "MD","sample"] = "Normal"

T_NK_cells@meta.data$old.ident = T_NK_cells@meta.data$orig.ident
T_NK_cells@meta.data$orig.ident = gsub("BH|DJ|MD|NZ","Normal",T_NK_cells@meta.data$orig.ident)

Idents(T_NK_cells) = "groups"

#groups = c("Untreated","Pt-11","Pt-17","AFT-03","AFT-04","Pt-AA13","Pt-25","Pt-27")
groups = c("Untreated","Pt-17","Pt-25")
for(i in 1:length(groups)){
        
        subset.MCL <- subset(T_NK_cells, idents = c("Normal",groups[i]))
        
        (samples = unique(subset.MCL$orig.ident))
        df = df_samples[df_samples$sample %in% samples,]
        subset.MCL@meta.data$orig.ident = factor(subset.MCL@meta.data$orig.ident, 
                                                 levels = df$sample[order(df$tsne)])
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
        write.csv(gde.markers,paste0(path,"/T/T_NK_DE/T_",groups[i],".csv"))
        
        gde.markers = read.csv(file = paste0("output/20190622/T/T_NK_DE/T_",groups[i],".csv"))
        (mito.genes <- grep(pattern = "^MT-", x = gde.markers$gene))
        if(length(mito.genes)>0) gde.markers = gde.markers[-mito.genes,]
        GC()
        #DoHeatmap.1======
        subset.MCL %<>% ScaleData(features= unique(gde.markers$gene))
        Top_n = 20
        DoHeatmap.1(subset.MCL, marker_df = gde.markers, Top_n = Top_n, do.print=T, angle = 0,
                    group.bar = T, title.size = 20, no.legend = F,size=5,hjust = 0.5,
                    label=F, cex.row=5, legend.size = NULL,width=10, height=7,unique.name = T,
                    title = paste("Top",Top_n,"DE genes of",groups[i],
                                  "Patient's T and NK cells over Normal cells"))
}


##############
# DE genes between Clusters 5 cluster top 50
###############
Idents(T_NK_cells) %<>% factor(levels = c("T_cells:CD4","T_cells:CD8",
                                          "NK_cells","HSC/progenitors"))
#T_NK_cells <- subset(T_NK_cells, idents = 1:2)
X3_clusters_markers <- FindAllMarkers.UMI(T_NK_cells,logfc.threshold = 0.1, only.pos = T, 
                                          min.pct = 0.1,return.thresh = 0.05)

X3_clusters_markers <- FindAllMarkers.UMI(T_NK_cells,logfc.threshold = -Inf,only.pos = FALSE, 
                                          min.pct = 0.01,return.thresh = 1)

write.csv(X3_clusters_markers,paste0(path,"T_NK_cells_X4clusters_FC0.1_markers.csv"))

X3_clusters_markers = read.csv("output/20190611/T_NK_cells_X4clusters_FC0.1_markers.csv",
                               row.names = 1)
X3_clusters_markers = X3_clusters_markers[(X3_clusters_markers$cluster %in% c("T_cells:CD4",
                                                                              "T_cells:CD8",
                                                                              "NK_cells")),]
X3_clusters_markers$cluster = factor(X3_clusters_markers$cluster,
                                     levels = c("T_cells:CD4","T_cells:CD8",
                                                "NK_cells"))
table(X3_clusters_markers$cluster)
T_NK_cells <- subset(T_NK_cells, idents = c("T_cells:CD4","T_cells:CD8","NK_cells"))
T_NK_cells %<>% ScaleData(features=rownames(T_NK_cells))

T_NK_cells$manual %<>% factor(levels = c("T_cells:CD4","T_cells:CD8","NK_cells"))

(MT_gene <- grep("^MT-",X3_clusters_markers$gene))
X3_clusters_markers = X3_clusters_markers[-MT_gene,]

DoHeatmap.1(T_NK_cells, marker_df = X3_clusters_markers, Top_n = 300, do.print=T, angle = 0,
            group.bar = T, title.size = 13, no.legend = F,size=5,hjust = 0.5,
            label=T, cex.row=0, legend.size = 15,width=7, height=20,
            title = "Top 300 differentially expressed genes in T and NK clusters")

MakeCorlorBar(T_NK_cells,  marker_df = X3_clusters_markers,Top_n = 300,do.return =F,do.print = T,
              legend.size = 15,width=7, height=20)
# remove cluster with less than 3 cells======
# no scale down, keep the normal cells.
T_NK_cells@meta.data$old.ident = T_NK_cells@meta.data$orig.ident
T_NK_cells@meta.data$orig.ident = gsub("BH|DJ|MD|NZ","Normal",T_NK_cells@meta.data$orig.ident)
T_NK_cells@meta.data$X4_orig.ident = paste(T_NK_cells@meta.data$orig.ident,
                                            T_NK_cells@meta.data$X3_clusters, sep = "_")
T_NK_cells@meta.data$X4_orig.ident = gsub('^Normal_.*', 'Normal', T_NK_cells@meta.data$X4_orig.ident)

# add color to cluster
T_NK_cells <- AddMetaColor(object = T_NK_cells, label= "X3_clusters", 
                           colors = scales::hue_pal()((length(unique(T_NK_cells$X3_clusters)))))

Idents(T_NK_cells) = "orig.ident"
df_samples <- readxl::read_excel("doc/190429_scRNAseq_info.xlsx")
colnames(df_samples) <- tolower(colnames(df_samples))
(keep = df_samples$sample %in% unique(T_NK_cells$orig.ident))
df_samples = df_samples[keep,]

tests <- paste0("test",c(2:12))

for(test in tests){
        sample_n = which(df_samples$tests %in% test)
        df <- as.data.frame(df_samples[sample_n,])
        samples <- unique(df$sample)
        rownames(df) = samples
        
        samples <- c(ifelse(length(samples)>5,NA,"Normal"),df$sample[order(df$tsne)])
        print(samples <- samples[!is.na(samples)])
        
        subset.MCL <- subset(T_NK_cells,idents = samples)
        subset.MCL$orig.ident = factor(subset.MCL$orig.ident, levels = samples)
        Idents(subset.MCL) = "X3_clusters"
        Idents(subset.MCL) %<>% factor(levels = c("T_cells:CD4","T_cells:CD8",
                                                  "NK_cells"))
        TSNEPlot.1(subset.MCL, split.by = "orig.ident",pt.size = 1,label = F,do.return = F,
                   cols = ExtractMetaColor(subset.MCL),do.print = T, unique.name = T)
}

for(test in tests){
        sample_n = which(df_samples$tests %in% test)
        df <- as.data.frame(df_samples[sample_n,])
        samples <- unique(df$sample)
        rownames(df) = samples
        
        samples <- c(ifelse(length(samples)>5,NA,"Normal"),df$sample[order(df$tsne)])
        print(samples <- samples[!is.na(samples)])
        
        subset.MCL <- subset(T_NK_cells,idents = samples)
        subset.MCL$orig.ident = factor(subset.MCL$orig.ident, levels = samples)
        Idents(subset.MCL) = "cell.type"
        subset.MCL <- sortIdent(subset.MCL)
        TSNEPlot.1(subset.MCL, split.by = "orig.ident",pt.size = 1,label = F,do.return = F,
                   cols = ExtractMetaColor(subset.MCL),do.print = T, unique.name = T)
}
table(T_NK_cells@meta.data$X4_orig.ident)
Idents(T_NK_cells) <- "X4_orig.ident"
table_T_NK_cells <- table(T_NK_cells@meta.data$X4_orig.ident) %>% as.data.frame
(keep.MCL <- table_T_NK_cells[table_T_NK_cells$Freq > 2,"Var1"] %>% as.character())
T_NK_cells <- subset(T_NK_cells, ident.use = keep.MCL)
T_NK_cells_exp <- AverageExpression(T_NK_cells)
write.csv(T_NK_cells_exp,paste0(path,"T_NK_cells_by_samples.csv"))


# Doheatmap for MCL longitudinal by samples ================
(load(file = "data/T_NK_cells_43_20190611.Rda"))
df_samples <- readxl::read_excel("doc/190429_scRNAseq_info.xlsx")
colnames(df_samples) =  tolower(colnames(df_samples))
df_samples[df_samples$sample %in% "MD","tsne"] = 0
df_samples[df_samples$sample %in% "MD","sample"] = "Normal"

T_NK_cells@meta.data$old.ident = T_NK_cells@meta.data$orig.ident
T_NK_cells@meta.data$orig.ident = gsub("BH|DJ|MD|NZ","Normal",T_NK_cells@meta.data$orig.ident)

Idents(T_NK_cells) = "groups"
table(Idents(T_NK_cells))
groups = c("Untreated","Pt-17","Pt-25")#,"Pt-11","AFT-03","AFT-04","Pt-AA13","Pt-27","Pt-13")
for(i in 1:length(groups)){
        
        subset.MCL <- subset(T_NK_cells, idents = c("Normal",groups[i]))
        (samples = unique(subset.MCL$orig.ident))
        df = df_samples[df_samples$sample %in% samples,]
        subset.MCL@meta.data$orig.ident = factor(subset.MCL@meta.data$orig.ident, 
                                                 levels = df$sample[order(df$tsne)])
        Idents(subset.MCL) = "orig.ident"
        samples = samples[-which(samples %in% "Normal")]
        gde.markers_list <- list()
        for(k in 1:(length(samples))){
                gde.markers_list[[k]] <- FindMarkers.UMI(subset.MCL,
                                                         ident.1 = samples[k],
                                                         ident.2 = "Normal",
                                                         logfc.threshold = 0.1, only.pos = F,
                                                         test.use = "MAST")
                gde.markers_list[[k]]$cluster <- samples[k]
                gde.markers_list[[k]]$gene <- rownames(x = gde.markers_list[[k]])
        }
        gde.markers <- do.call(rbind, gde.markers_list)
        write.csv(gde.markers,paste0(path,"T/T_NK_DE/T_",groups[i],".csv"))
}
        (mito.genes <- grep(pattern = "^MT-", x = gde.markers$gene))
        if(length(mito.genes)>0) gde.markers = gde.markers[-mito.genes,]
        GC()
        #DoHeatmap.1======
        subset.MCL %<>% ScaleData(features= unique(gde.markers$gene))
        Top_n =20
        DoHeatmap.1(subset.MCL, marker_df = gde.markers, Top_n = Top_n, do.print=T, angle = 0,
                    group.bar = T, title.size = 20, no.legend = F,size=5,hjust = 0.5,
                    label=F, cex.row=5, legend.size = NULL,width=10, height=7,unique.name = T,
                    title = paste("Top",Top_n,"DE genes of",groups[i],
                                  "over Normal T and NK cells"))
        remove(subset.MCL);GC()
}

# Doheatmap for MCL longitudinal by clusters ================

Idents(T_NK_cells) = "groups"
table(Idents(T_NK_cells))
groups = c("Untreated","Pt-11","Pt-17","AFT-03","AFT-04","Pt-AA13","Pt-25","Pt-27","Pt-13")
clusters = list("T_cells:CD4","T_cells:CD8","NK_cells")
for(i in 1:length(groups)){
        
        subset_MCL <- subset(T_NK_cells, idents = c("Normal",groups[i]))
        Idents(subset_MCL) = "X3_clusters"
        for(k in 1:length(clusters)){
                subset.MCL <- subset(subset_MCL,idents = clusters[[k]])
                
                (samples = unique(subset.MCL$orig.ident))
                df = df_samples[df_samples$sample %in% samples,]
                subset.MCL@meta.data$orig.ident = factor(subset.MCL@meta.data$orig.ident, 
                                                         levels = df$sample[order(df$tsne)])
                Idents(subset.MCL) = "orig.ident"
                gde.markers <- FindAllMarkers.UMI(subset.MCL,logfc.threshold = 0.3, only.pos = T,
                                                  test.use = "MAST")
                write.csv(gde.markers,paste0(path,"T_",groups[i],"_",clusters[[k]],".csv"))
                (mito.genes <- grep(pattern = "^MT-", x = gde.markers$gene))
                if(length(mito.genes)>0) gde.markers = gde.markers[-mito.genes,]
                GC()
                #DoHeatmap.1======
                subset.MCL %<>% ScaleData(features= unique(gde.markers$gene))
                Top_n =20
                DoHeatmap.1(subset.MCL, marker_df = gde.markers, Top_n = Top_n, do.print=T, angle = 0,
                            group.bar = T, title.size = 20, no.legend = F,size=5,hjust = 0.5,
                            label=F, cex.row=5, legend.size = NULL,width=10, height=7,unique.name = T,
                            title = paste("Top",Top_n,clusters[[k]],"DE genes in Normal and",groups[i]))
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

#============ Pt-17 =================
df_samples <- readxl::read_excel("doc/190626_scRNAseq_info.xlsx")
colnames(df_samples) <- tolower(colnames(df_samples))
tests <- paste0("test",4)
sample_n = which(df_samples$tests %in% tests)
df <- as.data.frame(df_samples[sample_n,])
(samples <- c("Normal",unique(df$sample)))
cell.type <- c("T_cells:CD4+","T_cells:CD8+","NK_cells")

T_NK_cells@meta.data$orig.ident = gsub("BH|DJ|MD|NZ","Normal",T_NK_cells@meta.data$orig.ident)
T_NK_cells@meta.data$orig.ident_cell.type = paste0(T_NK_cells@meta.data$orig.ident,".",
                                                   T_NK_cells@meta.data$manual)

Idents(T_NK_cells) = "groups"
Pt_17 <- subset(T_NK_cells, idents = c("Normal", "Pt-17"))
Pt_17@meta.data$orig.ident %<>% factor(levels = samples)

TSNEPlot.1(Pt_17,do.label = F,no.legend=T,do.print = T,do.return =T, pt.size = 0.1,
           group.by = "manual",cols =  ExtractMetaColor(Pt_17),
           split.by = 'orig.ident',ncol = length(samples), 
           title="T cells and NK cells in Patient 17",
           width=length(samples)*2+2, height=3)

Idents(Pt_17) = "orig.ident"
for (i in 2:length(samples)) {
        print(paste(samples[i],"vs. Normal"))
        subset.MCL <- subset(Pt_17, idents = c("Normal",samples[i]))
        
        Idents(subset.MCL) = "manual"
        TSNEPlot.1(subset.MCL, group.by = "manual",split.by = "orig.ident",pt.size = 0.3,
                   cols = ExtractMetaColor(subset.MCL),
                   label = F, do.print = F,do.return = T, unique.name = T)
        
        # remove cluster with less than 3 cells======
        table_subset.MCL <- table(subset.MCL$orig.ident_cell.type) %>% as.data.frame
        keep.MCL <- table_subset.MCL[table_subset.MCL$Freq > 2,"Var1"] %>% as.character()
        
        (cell.type1 <- keep.MCL %>% .[grep("Normal",.,invert = T)] %>%
                        gsub('.*\\.',"",.) %>% unique)

        print(ident.1 <- paste0(samples[i],".",cell.type1))
        print(ident.2 <- paste0("Normal.",cell.type1))
        
        Idents(subset.MCL) = "orig.ident_cell.type"
        gde.markers <- FindPairMarkers(subset.MCL, ident.1 = ident.1,
                                       ident.2 = ident.2, only.pos = F,
                                       logfc.threshold = 0.1,min.cells.group =3,
                                       min.pct = 0.1,return.thresh = 0.05,
                                       save.files = FALSE)
        write.csv(gde.markers, paste0(path,samples[i],"_vs_Normal",".csv"))
        }
