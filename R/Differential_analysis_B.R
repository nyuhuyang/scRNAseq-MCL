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
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#3.1  Compare DE across all major cell types==================
#We would need the data for all clusters, as well the subclusters.
#detect changes in gene expression between young and aged, 
#in the different cell types and subtypes. 
#It will also be interesting to check if there is some subtype enriched in young compared to aged or viceversa. 

# 3.1.1 load data
# Rename ident
(load(file="data/MCL_Harmony_30_20190320.Rda"))

# select 1/4 of cell from control
object <- ScaleDown(object = object)

# B cells only ================
object <- SetAllIdent(object, id="res.0.6")
table(object@ident)
TSNEPlot(object,do.label = T)
B_cells_MCL <- SubsetData(object, ident.use = c(0,1,4,5,10,12))
B_cells_MCL <- SetAllIdent(B_cells_MCL, id="singler1main")
table(B_cells_MCL@ident)
B_cells_MCL <- SubsetData(B_cells_MCL, ident.use = c("B_cells","HSC","MCL"))
table(B_cells_MCL@meta.data$singler1sub) %>% as.data.frame %>%
        .[.[,"Freq"] >0,]
B_cells_MCL <- SetAllIdent(B_cells_MCL, id="singler1sub") 
(ident.remove <- grep("T_cells.",B_cells_MCL@meta.data$singler1sub, value = T) %>% unique)
B_cells_MCL <- SubsetData(B_cells_MCL, ident.remove =  ident.remove)
remove(object);GC()
B_cells_MCL %<>% FindClusters(reduction.type = "harmony", resolution = 0.2, 
                              dims.use = 1:75,
                              save.SNN = TRUE, n.start = 10, nn.eps = 0.5,
                              force.recalc = TRUE, print.output = FALSE)
TSNEPlot(B_cells_MCL,do.label = T,label.size=5)
B_cells_MCL@ident <- plyr::mapvalues(x = B_cells_MCL@ident,
                                     from = c(0,1,2,3,4),
                                     to = c(1,2,3,4,5))
B_cells_MCL@ident %<>% factor(levels = 1:5)
B_cells_MCL <- StashIdent(object = B_cells_MCL, save.name = "X5_clusters")
B_cells_MCL@meta.data$X5_orig.ident = paste(B_cells_MCL@meta.data$orig.ident,
                                            B_cells_MCL@meta.data$X5_clusters, sep = "_")
B_cells_MCL@meta.data$X5_orig.ident = gsub('^Normal_.*', 'Normal', B_cells_MCL@meta.data$X5_orig.ident)

# tsne plot
B_cells_MCL %<>% SetAllIdent(id="X5_clusters")
TSNEPlot(B_cells_MCL, do.label = T)

###############################
# Doheatmap for Normal / MCL
###############################
df_samples <- readxl::read_excel("doc/190320_scRNAseq_info.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
sample_n = which(df_samples$tests %in% paste0("test",8))
(samples <- df_samples$sample[sample_n])
# remove samples with low B cells======
table_df <- table(B_cells_MCL@meta.data$orig.ident) %>% as.data.frame
keep <- table_df[table_df$Freq > 100,"Var1"] %>% as.character()
(samples <- samples[samples %in% keep])
B_cells_MCL %<>% SetAllIdent(id = "orig.ident")
for(sample in samples[1:length(samples)]){
        subset.MCL <- SubsetData(B_cells_MCL, ident.use = c(sample,"Normal"))

        # SplitTSNEPlot======
        subset.MCL <- SetAllIdent(subset.MCL,id = "X5_clusters")
        SplitTSNEPlot(subset.MCL, do.return = FALSE,do.print = TRUE)
        
        # remove cluster with less than 3 cells======
        subset.MCL %<>% SetAllIdent(id = "X5_orig.ident")
        table_subset.MCL <- table(subset.MCL@meta.data$X5_orig.ident) %>% as.data.frame
        keep.MCL <- table_subset.MCL[table_subset.MCL$Freq > 2,"Var1"] %>% as.character()
        subset.MCL <- SubsetData(subset.MCL, ident.use = keep.MCL)
        
        x6_cluster <- subset.MCL@ident %>% unique
        x6_cluster = x6_cluster[-grep("^Normal",x6_cluster)] %>% as.character %>% 
                gsub('.*\\_',"",.) %>% as.numeric %>% sort
        print(ident.1 <- rep("Normal",length(ident.1)))
        print(ident.2 <- paste(sample,x6_cluster,sep="_"))
        # FindAllMarkers.UMI======
        subfolder <- paste0(path,sample,"_vs_Normal/")
        
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
        markers <- c("BTK","CCND2","CCND3","CD3D","CD5","CD8A","CDK4","IL2RA",
                     "MS4A1","PDCD1","RB1","SOX11")
        g <- DoHeatmap.1(subset.MCL, gde.markers, add.genes = markers, Top_n = 25,
                         use.scaled = T,
                         ident.use = paste("Normal vs.",sample, "in B and MCL cells"),
                         group.label.rot = T,cex.row = 4,remove.key =F,title.size = 12)
        jpeg(paste0(path,"/DoHeatmap_Normal_",sample,".jpeg"), units="in", width=10, height=7,
             res=600)
        print(g)
        dev.off()
        #---DoHeatmap vertical bar----
        g1 <- MakeCorlorBar(subset.MCL, gde.markers,Top_n = 25, add.genes = markers, do.print = F,
                            do.return = T)
        jpeg(paste0(path,"/DoHeatmap_Normal_",sample,"_legend.jpeg"),
             units="in", width=10, height=7,res=600)
        print(g1)
        dev.off()
}


# Doheatmap for MCL.1 / MCL.2 ================
samples1 <- c("Pt-11-C14","Pt-25-C24","Pt-25-C24")
samples2 <- c("Pt-11-C28","Pt-25-C25","Pt-25-AMB-C25")
B_cells_MCL %<>% SetAllIdent(id="orig.ident")

for(i in 1:length(samples1)){
        subset.MCL <- SubsetData(B_cells_MCL, ident.use = c(samples1[i],samples2[i]))

        #---SplitTSNEPlot----
        subset.MCL %<>% SetAllIdent(id = "X5_orig.ident")
        if(all(c(samples1[i],samples2[i]) != sort(c(samples1[i],samples2[i])))){
                select.plots =2:1
        }
        SplitTSNEPlot(subset.MCL, do.return = FALSE,do.print = TRUE,select.plots = select.plots)
        # remove cluster with less than 3 cells======

        table_subset.MCL <- table(subset.MCL@meta.data$X5_orig.ident) %>% as.data.frame
        keep.MCL <- table_subset.MCL[table_subset.MCL$Freq > 2,"Var1"] %>% as.character()
        subset.MCL <- SubsetData(subset.MCL, ident.use = keep.MCL)
        
        x6_cluster <- subset.MCL@ident %>% unique %>% 
                gsub('.*\\_',"",.) %>% as.numeric %>% sort %>% .[duplicated(.)]
        
        print(ident.1 <- paste(samples1[i],x6_cluster,sep="_"))
        print(ident.2 <- paste(samples2[i],x6_cluster,sep="_"))
        
        subfolder <- paste0(path,samples1[i],"_vs_",samples2[i],"/")
        gde.markers <- FindPairMarkers(subset.MCL, ident.1 = c(ident.1,ident.2), 
                                       ident.2 = c(ident.2,ident.1), only.pos = T,
                                       logfc.threshold = 0.5,min.cells.group =3,
                                       min.pct = 0.01,
                                       return.thresh = 0.5,
                                       save.path = subfolder)
        (mito.genes <- grep(pattern = "^MT-", x = gde.markers$gene))
        if(length(mito.genes)>0) gde.markers = gde.markers[-mito.genes,]
        GC()

        #---DoHeatmap.1----
        subset.MCL <- ScaleData(subset.MCL)
        markers <- c("BTK","CCND2","CCND3","CD3D","CD5","CD8A","CDK4","IL2RA",
                     "MS4A1","PDCD1","RB1","SOX11")
        g <- DoHeatmap.1(subset.MCL, gde.markers, add.genes = markers, Top_n = 25,
                         group.order = c(ident.1,ident.2),
                         ident.use = paste(paste(c(samples1[i],samples2[i]),collapse = " vs. "),
                                           "in MCL cells"),
                         group.label.rot = T,cex.row = 3,remove.key =F,title.size = 12)
        jpeg(paste0(path,"/DoHeatmap_",paste(c(samples1[i],samples2[i]),collapse = "_"),
                    ".jpeg"), units="in", width=10, height=7,res=600)
        print(g)
        dev.off()
        #---DoHeatmap vertical bar----
        g1 <- MakeCorlorBar(subset.MCL, gde.markers,Top_n = 25, add.genes = markers, do.print = F,
                            do.return = T)
        jpeg(paste0(path,"/DoHeatmap_",paste(c(samples1[i],samples2[i]),collapse = "_"),
                    "_legend.jpeg"),units="in", width=10, height=7,res=600)
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

