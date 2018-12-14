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
source("../R/Seurat_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#3.1  Compare DE across all major cell types==================
#We would need the data for all clusters, as well the subclusters.
#detect changes in gene expression between young and aged, 
#in the different cell types and subtypes. 
#It will also be interesting to check if there is some subtype enriched in young compared to aged or viceversa. 

# 3.1.1 load data
# Rename ident
(load(file="data/MCL_Harmony_12_20181121.Rda"))

# B cell and MCL cells only ================
MCL <- SetAllIdent(MCL, id="res.0.6")
table(MCL@ident)
TSNEPlot.1(MCL,do.label = T)
B_cells_MCL <- SubsetData(MCL, ident.use = c(0,2,3,7,8,10,12))
B_cells_MCL <- SetAllIdent(B_cells_MCL, id="singler2main")
table(B_cells_MCL@ident)

B_cells_MCL <- SubsetData(B_cells_MCL, ident.use = c("B_cells","HSC","MCL"))
table(B_cells_MCL@meta.data$singler1main)
B_cells_MCL <- SetAllIdent(B_cells_MCL, id="singler1main")
B_cells_MCL <- SubsetData(B_cells_MCL, ident.use = c("B_cells"))
B_cells_MCL <- SetAllIdent(B_cells_MCL, id="singler2sub")
p4 <- TSNEPlot.1(B_cells_MCL, do.label = F, do.return = T, pt.size = 0.5, 
                 colors.use = ExtractMetaColor(B_cells_MCL), no.legend =T)

B_cells_MCL <- SetAllIdent(B_cells_MCL, id="res.0.6")
idents <- as.data.frame(table(B_cells_MCL@ident))
old.ident.ids <- idents$Var1
new.cluster.ids <- c(0,2,1,1,3,4,5)
B_cells_MCL@ident <- plyr::mapvalues(x = B_cells_MCL@ident,
                              from = old.ident.ids, to = new.cluster.ids)
B_cells_MCL@ident <- factor(B_cells_MCL@ident, levels = 0:5)
B_cells_MCL <- StashIdent(object = B_cells_MCL, save.name = "5_clusters")

p3 <- TSNEPlot.1(B_cells_MCL, do.return = T, pt.size = 0.5, do.label = T, 
                 group.by = "ident",no.legend =T )

jpeg(paste0(path,"/S1_TSNEPlot.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p3, p4)
dev.off()


# Doheatmap for Normal / MCL ================
markers <-  HumanGenes(B_cells_MCL,c("CCND1","CCND2","CCND3","CD19","MS4A1","CD79A","CD5","CD40",
                                     "CDK4","CDK6","PCNA","CDK1","SOX11",
                                     "RB1","TP53","ATM","MYC","MTAP",
                                     "FOXO1","FOXO3"))
(samples <- df_samples$sample[5:7])
control <- "MD"
for(sample in samples){
        cell.use <- rownames(B_cells_MCL@meta.data)[B_cells_MCL@meta.data$orig.ident %in% c(control,sample)] #
        subset.MCL <- SubsetData(B_cells_MCL, cells.use = cell.use)
        #---plot_grid(p1, p2)----
        p1 <- TSNEPlot.1(subset.MCL, do.return = T, pt.size = 0.5, do.label = T, 
                         group.by = "ident",no.legend =T )
        subset.MCL <- SetAllIdent(subset.MCL, id= "singler2sub")
        p2 <- TSNEPlot.1(subset.MCL, do.label = F, do.return = T, pt.size = 0.5, 
                         colors.use = ExtractMetaColor(subset.MCL), no.legend =T)
        jpeg(paste0(path,control,"_",sample,"_","TSNEPlot.jpeg"), units="in", width=10, height=7,res=600)
        print(plot_grid(p1, p2)+
                      ggtitle(paste(control, "vs.",sample, "in B and MCL cells"))+
                      theme(text = element_text(size=15),							
                            plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) )
        dev.off()
        #---FindAllMarkers.UMI----
        subset.MCL@meta.data$X5_clusters <- paste(subset.MCL@meta.data$orig.ident, 
                                                  subset.MCL@meta.data$X5_clusters, sep = "_")
        subset.MCL <- SetAllIdent(subset.MCL,id = "X5_clusters")
        print(table(subset.MCL@ident))
        #---DoHeatmap.1----
        test_markers <- FindAllMarkers.UMI(subset.MCL,logfc.threshold = 0.4, only.pos = T,
                                           test.use = "MAST")

        g <- DoHeatmap.1(subset.MCL, test_markers, add.genes = markers, Top_n = 25,
                         ident.use = paste(control, "vs.",sample, "in B and MCL cells"),
                         group.label.rot = T,cex.row = 2,remove.key =F,title.size = 12)
        jpeg(paste0(path,"/DoHeatmap_",control,"_",sample,".jpeg"), units="in", width=10, height=7,
             res=600)
        print(g)
        dev.off()
        #---DoHeatmap vertical bar----
        g1 <- MakeCorlorBar(subset.MCL, test_markers,Top_n = 25, add.genes = markers, do.print = F,
                            do.return = T)
        jpeg(paste0(path,"/DoHeatmap_",control,"_",sample,"_legend.jpeg"),
             units="in", width=10, height=7,res=600)
        print(g1)
        dev.off()
}


# Doheatmap for Normal / MCL.1 / MCL.2 ================
control <- "MD"
tests <- c("test3","test4")
for(test in tests){
        sample_n = which(df_samples$tests %in% test)
        samples <- df_samples$sample[sample_n]
        cell.use <- rownames(B_cells_MCL@meta.data)[B_cells_MCL@meta.data$orig.ident %in%
                                                            c(control,samples)]
        subset.MCL <- SubsetData(B_cells_MCL, cells.use = cell.use)
        #---plot_grid(p1, p2)----
        p1 <- TSNEPlot.1(subset.MCL, do.return = T, pt.size = 0.5, do.label = T, 
                         group.by = "ident",no.legend =T )
        subset.MCL <- SetAllIdent(subset.MCL, id= "singler2sub")
        p2 <- TSNEPlot.1(subset.MCL, do.label = F, do.return = T, pt.size = 0.5, 
                         colors.use = ExtractMetaColor(subset.MCL), no.legend =T)
        jpeg(paste0(path,paste0(samples,collapse = "_"),"TSNEPlot.jpeg"), units="in", width=10, height=7,res=600)
        print(plot_grid(p1, p2)+
                      ggtitle(paste(control, "vs.",paste0(samples,collapse = " vs. "), "in B and MCL cells"))+
                      theme(text = element_text(size=15),							
                            plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) )
        dev.off()
        #---FindAllMarkers.UMI----
        subset.MCL@meta.data$X5_clusters <- paste(subset.MCL@meta.data$orig.ident, 
                                                  subset.MCL@meta.data$X5_clusters, sep = ".")
        subset.MCL <- SetAllIdent(subset.MCL,id = "X5_clusters")
        table(subset.MCL@ident)
        test_markers <- FindAllMarkers.UMI(subset.MCL,logfc.threshold = 0.4, only.pos = T,
                                           test.use = "MAST")
        subset.MCL <- SetAllIdent(subset.MCL, id="5_clusters")
        #---DoHeatmap.1----
        g <- DoHeatmap.1(subset.MCL, test_markers, add.genes = markers, Top_n = 15,
                         ident.use = paste(control, "vs.",paste0(samples,collapse = " & "), "in MCL"),
                         group.label.rot = T,cex.row = 3,remove.key =F,
                         title.size = 12)
        jpeg(paste0(path,"/DoHeatmap_",control,"_",paste0(samples,collapse = "_"),".jpeg"),
                    units="in", width=10, height=7,res=600)
        print(g)
        dev.off()
        #---DoHeatmap vertical bar----
        g1 <- MakeCorlorBar(subset.MCL, test_markers,Top_n = 15, add.genes = markers, do.print = F,
                            do.return = T)
        jpeg(paste0(path,"/DoHeatmap_",control,"_",paste0(samples,collapse = "_"),"legend.jpeg"),
             units="in", width=10, height=7,res=600)
        print(g1)
        dev.off()
}


# heatmap.2 for Normal / MCL ================
markers <-  HumanGenes(B_cells_MCL,c("CD19","MS4A1","CD79A","CD5","CD40","CDK4"))
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

