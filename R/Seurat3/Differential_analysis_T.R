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
(load(file="data/MCL_V3_Harmony_43_20190610.Rda"))

# T cells only ================

Idents(object) <-  "singler1sub"
object <- sortIdent(object)
table(Idents(object))
TSNEPlot.1(object,label = F, repel = T, no.legend = T,pt.size = 1,
           cols = ExtractMetaColor(object),do.return = T,do.print = F,
           title = "Cell type labeling by Blueprint + Encode + MCL")

Idents(object) <-  "res.0.6"
object <- sortIdent(object,numeric=T)
table(Idents(object))
TSNEPlot.1(object,label = T, repel = T, no.legend = F,pt.size = 1,
           cols = singler.colors,do.return = T,do.print = F,
           title = "Unsupervised clustering")

T_NK_cells <- subset(object,  idents = c(2,3,4,10,12,15,16,17))
Idents(T_NK_cells) = "singler1main"
table(Idents(T_NK_cells))
T_NK_cells <- subset(T_NK_cells, idents = c("HSC","NK_cells","T_cells"))
table(T_NK_cells@meta.data$singler1sub) %>% as.data.frame %>%
        .[.[,"Freq"] >0,]
Idents(T_NK_cells) = "singler1sub"
HSC <- c("CLP","CMP","GMP","HSC","MEP","MPP")
(T_cell_type <- grep("T_cells:.*$",unique(T_NK_cells$singler1sub),value= T))
T_NK_cells <- subset(T_NK_cells, idents = c(HSC,T_cell_type,"NK_cells"))

T_NK_cells[['tSNE_1']] <- T_NK_cells@reductions$tsne@cell.embeddings[colnames(T_NK_cells),'tSNE_1']

T_NK_cells <- subset(T_NK_cells, subset = tSNE_1 >5)

T_NK_cells <- sortIdent(T_NK_cells)
table(Idents(T_NK_cells))
TSNEPlot.1(object = T_NK_cells,label = F, repel = T, no.legend = F,pt.size = 1,
           cols = ExtractMetaColor(T_NK_cells),do.return = T,do.print =F,
           title = "Tsne plot of all minor T cells")

T_NK_cells@meta.data$cell.type = gsub("T_cells:CD4.*$","T_cells:CD4",T_NK_cells@meta.data$singler1sub)
T_NK_cells@meta.data$cell.type = gsub("T_cells:CD8.*$","T_cells:CD8",T_NK_cells@meta.data$cell.type)
T_NK_cells@meta.data$cell.type = gsub("CLP|CMP|GMP|HSC|MEP|MPP","HSC/progenitors",T_NK_cells@meta.data$cell.type)
table(T_NK_cells@meta.data$cell.type)

T_NK_cells@meta.data$cell.type.colors <- mapvalues(T_NK_cells@meta.data$cell.type,
                                            from = c("HSC/progenitors","NK_cells",
                                                     "T_cells:CD4","T_cells:CD8",
                                                     "T_cells:Tregs"),
                                            to =   c("#4038b0","#A65628",
                                                     "#FDB462","#F0027F",
                                                     "#7570B3"))
Idents(T_NK_cells) <- "cell.type"
T_NK_cells <- sortIdent(T_NK_cells)
table(Idents(T_NK_cells))

g0 <- TSNEPlot.1(T_NK_cells,do.label = F,no.legend=T,do.print = T,do.return =T, pt.size = 1,
           cols = ExtractMetaColor(T_NK_cells),
           title="Cell types")


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
g1 <- TSNEPlot.1(T_NK_cells,label = T,no.legend=T,do.print = F,pt.size = 1,
           title="unsupervised clustering of T and NK cells")

T_NK_cells@meta.data$X4_clusters <- plyr::mapvalues(x = T_NK_cells@meta.data$RNA_snn_res.0.3,
                                                     from = c(0,1,2,3,4,5,6,7,8),
                                                     to =   c(1,3,2,2,2,2,1,1,4)) %>% as.character()
T_NK_cells@meta.data$X4_clusters <- plyr::mapvalues(x = T_NK_cells@meta.data$X4_clusters,
                                               from = 1:4,
                                               to = c("T_cells:CD4","T_cells:CD8",
                                                      "NK_cells","HSC/progenitors"))  %>% as.character()
T_NK_cells@meta.data$X4_clusters.colors <- mapvalues(x = T_NK_cells@meta.data$X4_clusters,
                                                   from = c("T_cells:CD4","T_cells:CD8",
                                                            "NK_cells","HSC/progenitors"),
                                                   to =   c("#FDB462","#F0027F",
                                                            "#A65628","#4038b0")) %>% as.character()
Idents(T_NK_cells) <- 'X4_clusters'
T_NK_cells <- sortIdent(T_NK_cells)
table(Idents(T_NK_cells))
# tsne plot
g2 <- TSNEPlot.1(T_NK_cells,label = F,no.legend=F,do.print = T,pt.size = 1,
                 cols = ExtractMetaColor(T_NK_cells),
                 title="unsupervised clustering")
jpeg(paste0(path,"rename_tsne_T_NK.jpeg"), units="in", width=10, height=7,res=600)
cowplot::plot_grid(g0,g2) +  ggtitle("Rename the clusters")+
        theme(text = element_text(size=15),
              plot.title = element_text(hjust = 0.5,size=15))
dev.off()
TSNEPlot.1(T_NK_cells,label = F,no.legend=F,do.print = T,pt.size = 1,do.return = T,
           title="Merged clusters of T and NK cells")
T_NK_cells_exp <- AverageExpression(T_NK_cells)
write.csv(T_NK_cells_exp,paste0(path,"T_NK_cells_exp.csv"))

save(T_NK_cells, file = "data/T_NK_cells_43_20190611.Rda")
(load(file = "data/T_NK_cells_43_20190611.Rda"))

##############
# DE genes between Clusters 5 cluster top 50
###############
Idents(T_NK_cells) %<>% factor(levels = c("T_cells:CD4","T_cells:CD8",
                                          "NK_cells","HSC/progenitors"))
#T_NK_cells <- subset(T_NK_cells, idents = 1:2)
X4_clusters_markers <- FindAllMarkers.UMI(T_NK_cells,logfc.threshold = 0.1, only.pos = T, 
                                          min.pct = 0.1,return.thresh = 0.05)

X4_clusters_markers <- FindAllMarkers.UMI(T_NK_cells,logfc.threshold = -Inf,only.pos = FALSE, 
                                          min.pct = 0.01,return.thresh = 1)

write.csv(X4_clusters_markers,paste0(path,"T_NK_cells_X4clusters_FC0.1_markers.csv"))

X4_clusters_markers = read.csv("output/20190611/T_NK_cells_X4clusters_FC0.1_markers.csv",
                               row.names = 1)
X4_clusters_markers = X4_clusters_markers[(X4_clusters_markers$cluster %in% c("T_cells:CD4",
                                                                              "T_cells:CD8",
                                                                              "NK_cells")),]
X4_clusters_markers$cluster = factor(X4_clusters_markers$cluster,
                                     levels = c("T_cells:CD4","T_cells:CD8",
                                                "NK_cells"))
table(X4_clusters_markers$cluster)
T_NK_cells <- subset(T_NK_cells, idents = c("T_cells:CD4","T_cells:CD8","NK_cells"))
T_NK_cells %<>% ScaleData(features=rownames(T_NK_cells))

Idents(T_NK_cells) %<>% factor(levels = c("T_cells:CD4","T_cells:CD8","NK_cells"))

(MT_gene <- grep("^MT-",X4_clusters_markers$gene))
X4_clusters_markers = X4_clusters_markers[-MT_gene,]

DoHeatmap.1(T_NK_cells, marker_df = X4_clusters_markers, Top_n = 300, do.print=T, angle = 0,
            group.bar = T, title.size = 13, no.legend = F,size=5,hjust = 0.5,
            label=T, cex.row=0, legend.size = 15,width=7, height=20,
            title = "Top 300 differentially expressed genes in T and NK clusters")

MakeCorlorBar(T_NK_cells,  marker_df = X4_clusters_markers,Top_n = 300,do.return =F,do.print = T,
              legend.size = 15,width=7, height=20)
# remove cluster with less than 3 cells======
# no scale down, keep the normal cells.
T_NK_cells@meta.data$old.ident = T_NK_cells@meta.data$orig.ident
T_NK_cells@meta.data$orig.ident = gsub("BH|DJ|MD|NZ","Normal",T_NK_cells@meta.data$orig.ident)
T_NK_cells@meta.data$X4_orig.ident = paste(T_NK_cells@meta.data$orig.ident,
                                            T_NK_cells@meta.data$X4_clusters, sep = "_")
T_NK_cells@meta.data$X4_orig.ident = gsub('^Normal_.*', 'Normal', T_NK_cells@meta.data$X4_orig.ident)

table(T_NK_cells@meta.data$X4_orig.ident)
Idents(T_NK_cells) <- "X4_orig.ident"
table_T_NK_cells <- table(T_NK_cells@meta.data$X4_orig.ident) %>% as.data.frame
(keep.MCL <- table_T_NK_cells[table_T_NK_cells$Freq > 2,"Var1"] %>% as.character())
T_NK_cells <- subset(T_NK_cells, ident.use = keep.MCL)
T_NK_cells_exp <- AverageExpression(T_NK_cells)
write.csv(T_NK_cells_exp,paste0(path,"T_NK_cells_by_samples.csv"))

# add color to cluster
gg_color_hue <- function(n) {
        hues = seq(15, 375, length = n + 1)
        grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
}
gg_color_hue(8)
T_NK_cells <- AddMetaColor(object = T_NK_cells, label= "X4_clusters", colors = gg_color_hue(8))

Idents(T_NK_cells) = "orig.ident"
df_samples <- readxl::read_excel("doc/190429_scRNAseq_info.xlsx")
colnames(df_samples) <- tolower(colnames(df_samples))

#meta.data = T_NK_cells@meta.data
#(remove <- colnames(meta.data) %in% "singler1main" %>% which)
#meta.data = meta.data[,-remove]
#T_NK_cells@meta.data = meta.data
tests <- paste0("test",c(2:12))
tsnes <- c("singler1sub")#"cell.type")#,"X4_clusters")
for(tsne in tsnes){
        folder_path <- paste0(path,tsne,"/")
        if(!dir.exists(folder_path)) dir.create(folder_path, recursive = T)
        for(test in tests){
                sample_n = which(df_samples$tests %in% test)
                df <- as.data.frame(df_samples[sample_n,])
                samples <- unique(df$sample)
                rownames(df) = samples
                
                samples <- c(ifelse(length(samples)>5,NA,"Normal"),df$sample[order(df$tsne)])
                print(samples <- samples[!is.na(samples)])
                
                g <- lapply(samples,function(sample) {
                        temp <- subset(T_NK_cells, idents = sample)
                        Idents(temp) <- tsne
                        TSNEPlot.1(temp, no.legend = T,label = F, label.size=3,size=20,
                                   cols = ExtractMetaColor(temp),pt.size =1,
                                   do.return = T, repel = F, do.print = F)+
                        ggtitle(sample)+theme(text = element_text(size=15),
                                              plot.title = element_text(hjust = 0.5))
                })
                jpeg(paste0(folder_path,test,"_Plots.jpeg"), units="in", width=10, height=7,
                     res=600)
                print(do.call(cowplot::plot_grid, c(g, nrow = ifelse(length(samples)>2,2,1))))
                dev.off()
        }
}

###############################
# Doheatmap for Normal / MCL
###############################
df_samples <- readxl::read_excel("doc/190429_scRNAseq_info.xlsx",sheet = "heatmap")
list_samples <- df2list(df_samples)
print(list_samples %>% unlist %>% as.vector %>% unique %in% 
              T_NK_cells@meta.data$orig.ident)

T_NK_cells %<>% SetAllIdent(id = "orig.ident")
for(sample in list_samples$MCL){
        subset.MCL <- SubsetData(T_NK_cells, ident.use = c("Normal",sample))

        # SplitTSNEPlot======
        subset.MCL <- SetAllIdent(subset.MCL,id = "X3_orig.ident")
        SplitTSNEPlot(subset.MCL, select.plots = order(c("Normal",sample)),
                      do.return = FALSE,do.print = TRUE)
        
        # remove cluster with less than 3 cells======
        subset.MCL %<>% SetAllIdent(id = "X3_orig.ident")
        table_subset.MCL <- table(subset.MCL@meta.data$X5_orig.ident) %>% as.data.frame
        keep.MCL <- table_subset.MCL[table_subset.MCL$Freq > 2,"Var1"] %>% as.character()
        subset.MCL <- SubsetData(subset.MCL, ident.use = keep.MCL)
        
        X5_cluster <- subset.MCL@ident %>% unique
        X5_cluster = X5_cluster[-grep("^Normal",X5_cluster)] %>% as.character %>% 
                gsub('.*\\_',"",.) %>% as.numeric %>% sort
        print(ident.1 <- rep("Normal",length(X5_cluster)))
        print(ident.2 <- paste(sample,X5_cluster,sep="_"))
        # FindAllMarkers.UMI======
        subfolder <- paste0(path,"Normal_vs_",sample,"/")
        
        gde.markers <- FindPairMarkers(subset.MCL, ident.1 = c(ident.1,ident.2),
                                       ident.2 = c(ident.2,ident.1), only.pos = T,
                                       logfc.threshold = 1.005,min.cells.group =3,
                                       min.pct = 0.01,return.thresh = 0.05,
                                       save.path = subfolder)
        (mito.genes <- grep(pattern = "^MT-", x = gde.markers$gene))
        if(length(mito.genes)>0) gde.markers = gde.markers[-mito.genes,]
        GC()
        #DoHeatmap.1======
        subset.MCL <- ScaleData(subset.MCL)
        Top_n =25
        g1 <- DoHeatmap.1(subset.MCL, gde.markers, add.genes = markers, Top_n = 25,
                         use.scaled = T,group.order = c("Normal",ident.2),
                         title=paste("Top",Top_n,"DE genes","in Normal B/MCL cells vs.",sample, "B/MCL cells"),
                         group.label.rot = T,cex.row = 3,remove.key =F,title.size = 12)
        jpeg(paste0(subfolder,"DoHeatmap_Normal_",sample,".jpeg"), units="in", width=10, height=7,
             res=600)
        print(g1)
        dev.off()
        #---DoHeatmap vertical bar----
        g2 <- MakeCorlorBar(subset.MCL, gde.markers,Top_n = 25, do.print = F,
                            do.return = T)
        jpeg(paste0(subfolder,"DoHeatmap_Normal_",sample,"_legend.jpeg"),units="in", width=10, height=7,res=600)
        print(g2)
        dev.off()
}


# Doheatmap for MCL.1 / MCL.2 ================
df_samples <- readxl::read_excel("doc/190429_scRNAseq_info.xlsx",sheet = "heatmap")
list_samples <- df2list(df_samples)
print(list_samples %>% unlist %>% as.vector %>% unique %in% 
              T_NK_cells@meta.data$orig.ident)

T_NK_cells %<>% SetAllIdent(id = "orig.ident")
for(i in 1:length(list_samples$MCL.1)){
        
        samples1 = list_samples$MCL.1[i]
        samples2 = list_samples$MCL.2[i]
        
        subset.MCL <- SubsetData(T_NK_cells, ident.use = c(samples1,samples2))

        #---SplitTSNEPlot----
        subset.MCL %<>% SetAllIdent(id = "X5_orig.ident")
        select.plots =1:2
        if(all(c(samples1,samples1) != sort(c(samples1,samples2)))){
                select.plots =2:1
        }
        SplitTSNEPlot(subset.MCL, do.return = FALSE,do.print = TRUE,select.plots = select.plots)
        # remove cluster with less than 3 cells======

        table_subset.MCL <- table(subset.MCL@meta.data$X5_orig.ident) %>% as.data.frame
        keep.MCL <- table_subset.MCL[table_subset.MCL$Freq > 2,"Var1"] %>% as.character()
        
        X5_cluster <- subset.MCL@ident %>% unique %>% 
                gsub('.*\\_',"",.) %>% as.numeric %>% sort %>% .[duplicated(.)]
        
        print(ident.1 <- paste(samples1,X5_cluster,sep="_"))
        print(ident.2 <- paste(samples2,X5_cluster,sep="_"))
        
        subset.MCL <- SubsetData(subset.MCL, ident.use = c(ident.1,ident.2))
        
        subfolder <- paste0(path,samples1,"_vs_",samples2,"/")
        gde.markers <- FindPairMarkers(subset.MCL, ident.1 = c(ident.1,ident.2), 
                                       ident.2 = c(ident.2,ident.1), only.pos = T,
                                       logfc.threshold = 0.125,min.cells.group =3,
                                       min.pct = 0.01,
                                       return.thresh = 0.5,
                                       save.path = subfolder)
        (mito.genes <- grep(pattern = "^MT-", x = gde.markers$gene))
        if(length(mito.genes)>0) gde.markers = gde.markers[-mito.genes,]
        GC()

        #---DoHeatmap.1----
        subset.MCL <- ScaleData(subset.MCL)
        markers <- c("BTK","CCND1","CCND2","CCND3","CD3D","CD5","CD8A","CDK4","IL2RA",
                     "MS4A1","PDCD1","RB1","SOX11")
        Top_n = 25
        g <- DoHeatmap.1(subset.MCL, gde.markers, add.genes = markers, Top_n = 25,
                         group.order = c(ident.1,ident.2),
                         title=paste("Top",Top_n,"DE genes","in",samples1,"vs.",samples2, "in B/MCL cells"),
                         group.label.rot = T,cex.row = 3,remove.key =F,title.size = 12)
        jpeg(paste0(subfolder,"DoHeatmap_",paste(c(samples1,samples2),collapse = "_"),
                    ".jpeg"), units="in", width=10, height=7,res=600)
        print(g)
        dev.off()
        #---DoHeatmap vertical bar----
        g1 <- MakeCorlorBar(subset.MCL, gde.markers,Top_n = 25, add.genes = markers, do.print = F,
                            do.return = T)
        jpeg(paste0(subfolder,"DoHeatmap_",paste(c(samples1,samples2),collapse = "_"),
                    "_legend.jpeg"),units="in", width=10, height=7,res=600)
        print(g1)
        dev.off()
}

# heatmap.2 for Normal / MCL ================
markers <-  HumanGenes(T_NK_cells,c("CD19","MS4A1","CD79A","CD5","CD40","CDK4"))
control <- "MD"
tests <- c("test3","test4")
for(test in tests){
        sample_n = which(df_samples$tests %in% test)
        samples <- df_samples$samples[sample_n]
        cell.use <- rownames(T_NK_cells@meta.data)[T_NK_cells@meta.data$orig.ident %in% c(control,sample)] #
        subset.MCL <- SubsetData(T_NK_cells, cells.use = cell.use)
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
markers <-  HumanGenes(T_NK_cells,c("CD19","MS4A1","CD79A","CD5","CD40","CDK4"))
tests <- c("test3","test4")
for(test in tests){
        sample_n = which(df_samples$tests %in% test)
        df_samples[sample_n,] %>% kable() %>% kable_styling()
        samples <- df_samples$samples[sample_n]
        cell.use <- rownames(T_NK_cells@meta.data)[T_NK_cells@meta.data$orig.ident %in% samples]
        subset.MCL <- SubsetData(T_NK_cells, cells.use = cell.use)
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
        
        jpeg(paste0(path,"/Heatmap2_MD_",samples[2],"_",samples[3],".jpeg"), units="in", width=10, height=7,res=600)
        heatmap.2(as.matrix(y),
                  Colv = as.dendrogram(hc), Rowv= FALSE,
                  ColSideColors = cc, trace ="none",labCol = FALSE,dendrogram = "column",
                  key.xlab = "scale log nUMI",
                  cexRow = 0.5,
                  margins = c(2,5),
                  #scale = "row",
                  breaks = seq(-3,3,length.out = 101),
                  col = bluered,
                  main = paste(samples[1],"vs.",samples[2],"vs.",samples[3], "in B cells"))
        par(lend = 1)           # square line ends for the color legend
        legend(0, 0.8,       # location of the legend on the heatmap plot
               legend = c(samples[1], samples[2], samples[3]), # category labels
               col = c("#B3DE69", "#195016","#E31A1C"),  # color key
               lty= 1,             # line style
               lwd = 10            # line width
        )
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
cell.use <- rownames(T_NK_cells@meta.data)[T_NK_cells@meta.data$orig.ident %in% samples]
subset.MCL <- SubsetData(T_NK_cells, cells.use = cell.use)
subset.MCL <- SetAllIdent(subset.MCL,id = "orig.ident")
TSNEPlot(subset.MCL)
MCL_exp <- AverageExpression(subset.MCL)

(samples <- c("DJ","MD","Pt-1294","Pt-RM","Pt-MS","Pt-LM","Pt-1475"))
markers <-  unique(MCL_markers$gene)[1:100]
markers <- markers[markers %in% rownames(MCL_exp)]
y = MCL_exp[markers,]

## Column clustering (adjust here distance/linkage methods to what you need!)
#hc <- hclust(dist(1-cor(as.matrix(y), method="spearman")), method="complete")

jpeg(paste0(path,"/Heatmap2_MCL_sc_100.jpeg"), units="in", width=10, height=7,res=600)
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
T_NK_cells <- SetAllIdent(T_NK_cells,id='cell.type')
cell.type_markers <- FindAllMarkers(T_NK_cells,logfc.threshold=0.5)
PROM1_markers <- FindMarkers.UMI(T_NK_cells,logfc.threshold = 0, genes.use = "PROM1",
                                 min.pct = 0,ident.1 = "MCL/HSC",return.thresh = 1)

write.table(cell.type_markers,paste0(path,"cell.type_markers.txt"))


object@meta.data$cell.type = gsub("B_cells.*","B_cells",object@meta.data$singler1sub)
object@meta.data$cell.type = gsub("T_NK_cells.*","T_NK_cells",object@meta.data$cell.type)
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
