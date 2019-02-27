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
(load(file="data/MCL_Harmony_24_20190128.Rda"))
df_samples <- readxl::read_excel("doc/190126_scRNAseq_info.xlsx")
colnames(df_samples) <- tolower(colnames(df_samples))
sample_n = which(df_samples$tests %in% c("control",paste0("test",2:7)))
df_samples[sample_n,] %>% kable() %>% kable_styling()
table(df_samples$tests);nrow(df_samples)
(samples <- df_samples$sample[sample_n])
(tests <- df_samples$tests[sample_n]) %>% unique
# cell population
table(object@meta.data$singler1main, object@meta.data$orig.ident) %>% 
        as.data.frame.matrix %>% .[,samples] %>% t %>% kable %>%
        kable_styling()

# B cell and MCL cells only ================
object <- SetAllIdent(object, id="res.0.6")
table(object@ident)
TSNEPlot.1(object,do.label = T)
B_cells_MCL <- SubsetData(object, ident.use = c(0,1,4,5,6,9))
B_cells_MCL <- SetAllIdent(B_cells_MCL, id="singler1main")
table(B_cells_MCL@ident)

B_cells_MCL <- SubsetData(B_cells_MCL, ident.use = c("B_cells","HSC","MCL"))
table(B_cells_MCL@meta.data$singler1sub)
B_cells_MCL <- SetAllIdent(B_cells_MCL, id="singler1sub")
(ident.remove <- grep("T_cells.",B_cells_MCL@meta.data$singler1sub, value = T) %>% unique)
B_cells_MCL <- SubsetData(B_cells_MCL, ident.remove =  ident.remove)

# cell population
table(B_cells_MCL@meta.data$singler1main, B_cells_MCL@meta.data$orig.ident) %>% 
        as.data.frame.matrix %>% .[,samples] %>% t %>% kable %>%
        kable_styling()
table(B_cells_MCL@meta.data$orig.ident) %>% .[samples] %>% kable %>% kable_styling()
# tsne plot
B_cells_MCL %<>% SetAllIdent(id="res.0.6")
p2 <- TSNEPlot.1(B_cells_MCL, do.label = T, do.return = T, pt.size = 0.5, no.legend =T)
                 #colors.use = ExtractMetaColor(B_cells_MCL))

B_cells_MCL %<>% FindClusters(reduction.type = "harmony", resolution = 0.2, dims.use = 1:50,
             save.SNN = TRUE, n.start = 10, nn.eps = 0.5,
             force.recalc = TRUE, print.output = FALSE)
idents <- as.data.frame(table(B_cells_MCL@ident))
(old.ident.ids <- idents$Var1)
new.cluster.ids <- c(1,3,4,5,2,6)
B_cells_MCL@ident <- plyr::mapvalues(x = B_cells_MCL@ident,
                              from = old.ident.ids, to = new.cluster.ids)
B_cells_MCL@ident <- factor(B_cells_MCL@ident, levels = 1:6)
B_cells_MCL <- StashIdent(object = B_cells_MCL, save.name = "X6_clusters")
B_cells_MCL %<>% SubsetData(ident.remove = "6")
p3 <- TSNEPlot.1(B_cells_MCL, do.return = T, pt.size = 0.5, do.label = T, 
                 group.by = "ident",no.legend =T )

jpeg(paste0(path,"/S1_MCL_TSNEPlot.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p2, p3)+ggtitle("B cell only") + 
        theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"))
dev.off()
remove(object);GC();GC();GC();GC();GC();GC();GC();GC();GC();

B_cells_MCL <- SetAllIdent(B_cells_MCL, id = "orig.ident")
Normal <- SubsetData(B_cells_MCL, ident.use = c("BH","DJ","MD","NZ"))
g1 <- TSNEPlot(Normal)
X1.4_normal <- SetAllIdent(X1.4_normal, id = "old.ident")
g2 <- TSNEPlot(X1.4_normal)

jpeg(paste0(path,"/X1.4_TSNEPlot.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(g1+xlim(-25, 20)+ylim(-30, 15), g2+xlim(-25, 20)+ylim(-30, 15))+
        ggtitle("Using 1/4  of B cells from each of the 4 normal samples") + 
        theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"))
dev.off()

# average UMI for each subpopulation in each patient
split_B_cells_MCL <- SplitSeurat(B_cells_MCL, split.by = "orig.ident")
ave.B_cells_MCL <- lapply(split_B_cells_MCL,AverageExpression)
for(i in 1:length(ave.B_cells_MCL)) write.csv(ave.B_cells_MCL[[i]],
                                             file = paste0(path,names(ave.B_cells_MCL[i]),"_exp.csv"))
# select 1/4 of cell from control
normal_cells = lapply(c("BH","DJ","MD","NZ"), function(x){
        rownames(B_cells_MCL@meta.data)[(B_cells_MCL@meta.data$orig.ident %in% x)]
})
remove_normal_cells = lapply(normal_cells, function(x) sample(x, size = length(x)*3/4)) %>%
        unlist
table(B_cells_MCL@cell.names %in% remove_normal_cells)
cell.use <- B_cells_MCL@cell.names[!(B_cells_MCL@cell.names %in% remove_normal_cells)]
B_cells_MCL <- SubsetData(B_cells_MCL, cells.use = cell.use)

X1.4_normal <- SubsetData(B_cells_MCL, ident.use = "normal")
X1.4_normal <- SetAllIdent(X1.4_normal, id = "old.ident")
TSNEPlot.1(X1.4_normal,do.label = T, do.print =T, do.return =F)

# Doheatmap for Normal / MCL ================
df_markers <- readxl::read_excel("doc/MCL-markers.xlsx")
(markers <- df_markers[,1] %>% as.matrix %>% as.character %>% HumanGenes(B_cells_MCL,marker.genes = .))
table(B_cells_MCL@meta.data$orig.ident)
table(B_cells_MCL@ident)

B_cells_MCL@meta.data$orig.ident = gsub("BH|DJ|MD|NZ","Normal",B_cells_MCL@meta.data$orig.ident)
Normal.cells <- rownames(B_cells_MCL@meta.data)[(B_cells_MCL@meta.data$orig.ident =="Normal")]
B_cells_MCL@meta.data[Normal.cells,"X6_clusters"] = "normal"
B_cells_MCL <- SetAllIdent(B_cells_MCL,id = "X6_clusters")
table(B_cells_MCL@ident)

sample_n = which(df_samples$tests %in% c("test2"))
(samples <- unique(df_samples$sample[sample_n])[2])

control <- "Normal"
for(sample in samples){
        cell.use <- rownames(B_cells_MCL@meta.data)[B_cells_MCL@meta.data$orig.ident %in% c(control,sample)] #
        subset.MCL <- SubsetData(B_cells_MCL, cells.use = cell.use)

        #---SplitTSNEPlot----
        subset.MCL <- SetAllIdent(subset.MCL, id= "singler1sub")
        SplitTSNEPlot(subset.MCL,do.return = FALSE,do.print = TRUE, do.label = F)
        subset.MCL <- SetAllIdent(subset.MCL,id = "X6_clusters")
        SplitTSNEPlot(subset.MCL,do.return = FALSE,do.print = TRUE)

        #---FindAllMarkers.UMI----
        MCL.dent.1 <- subset.MCL@ident %>% unique
        (MCL.dent.1 = MCL.dent.1[-which(MCL.dent.1 == "normal")] %>% droplevels %>% as.numeric %>% sort)
        subset.MCL@ident = factor(subset.MCL@ident, levels = c("normal",MCL.dent.1))
        print(table(subset.MCL@ident))
        gde.markers <- FindPairMarkers(subset.MCL, ident.1 = MCL.dent.1, 
                                       ident.2 = rep("normal",length(MCL.dent.1)),
                                       logfc.threshold = 0.2,min.cells.group =1,
                                       return.thresh = 0.05,save.files = FALSE)
        write.csv(gde.markers, paste0(path,control,"_",sample,".csv"))
        (mito.genes <- grep(pattern = "^MT-", x = gde.markers$gene))
        if(length(mito.genes)>0) gde.markers = gde.markers[-mito.genes,]
        GC();GC();GC();GC();GC();GC();GC();GC();
        subset.MCL <- ScaleData(subset.MCL)
        #---DoHeatmap.1----
        g <- DoHeatmap.1(subset.MCL, gde.markers, add.genes = markers, Top_n = 25,
                         use.scaled = FALSE,
                         ident.use = paste(control, "vs.",sample, "in B and MCL cells"),
                         group.label.rot = F,cex.row = 4,remove.key =F,title.size = 12)
        jpeg(paste0(path,"/DoHeatmap_",control,"_",sample,"~.jpeg"), units="in", width=10, height=7,
             res=600)
        print(g)
        dev.off()
        #---DoHeatmap vertical bar----
        g1 <- MakeCorlorBar(subset.MCL, gde.markers,Top_n = 25, add.genes = markers, do.print = F,
                            do.return = T)
        jpeg(paste0(path,"/DoHeatmap_",control,"_",sample,"_legend.jpeg"),
             units="in", width=10, height=7,res=600)
        print(g1)
        dev.off()
}

# Doheatmap for MCL.1 / MCL.2 ================
df_samples <- readxl::read_excel("doc/190126_scRNAseq_info.xlsx")
colnames(df_samples) <- tolower(colnames(df_samples))
Pairs <-list("11-C1-C14"=c("Pt-11-C1","Pt-11-C14"),
             "11-C14-C28"=c("Pt-11-C14","Pt-11-C28"),
             "11-C1-C28"=c("Pt-11-C1","Pt-11-C28"),
             "11-C1-C1-LN"=c("Pt-11-C1","Pt-11-LN-C1"),
             "11-C2-C7" = c("Pt-17-C2","Pt-17-C7"),
             "17-C2-C31" = c("Pt-17-C2","Pt-17-C31"),
             "17-C7-C31" = c("Pt-17-C7","Pt-17-C31"),
             "17-C1-LN-C2" = c("Pt-17-LN-C1","Pt-17-C2"))
B_cells_MCL %<>% SetAllIdent(id="orig.ident")

for(pair in Pairs){
        subset.MCL <- SubsetData(B_cells_MCL, ident.use = pair)
        select.plots = match(pair,unique(subset.MCL@meta.data$orig.ident %>% sort))
        #---SplitTSNEPlot----
        #subset.MCL <- SetAllIdent(subset.MCL, id= "singler1sub")
        #SplitTSNEPlot(subset.MCL,do.return = FALSE,do.print = TRUE, do.label = F)
        subset.MCL <- SetAllIdent(subset.MCL,id = "X6_clusters")
        SplitTSNEPlot(subset.MCL, select.plots = select.plots,
                      do.return = FALSE,do.print = TRUE)
        # remove low cell ident
        subset.MCL@meta.data$clusters_6 <- paste(subset.MCL@meta.data$orig.ident, 
                                                 subset.MCL@meta.data$X6_clusters, sep = ".")
        subset.MCL %<>% SetAllIdent(id='clusters_6')
        
        keep = table(subset.MCL@meta.data$clusters_6) %>% as.data.frame %>% .[.$Freq >2,"Var1"]
        subset.MCL %<>% SubsetData(ident.use = keep)
        #---FindAllMarkers.UMI----
        MCL.dent = list()
        for(i in 1:length(select.plots)){
                ident = pair[select.plots[i]]
                X6_clusters = subset.MCL@meta.data$X6_clusters[(subset.MCL@meta.data$orig.ident %in%
                                                                        ident)] %>% unique %>% sort
                MCL.dent[[i]] = paste(ident, X6_clusters, sep = ".")
        }

        
        if(length(MCL.dent[[1]]) != length(MCL.dent[[2]])){
                MCL.dent[[1]] = MCL.dent[[1]][(gsub('.*\\.', '', MCL.dent[[1]]) %in% 
                                                       gsub('.*\\.', '', MCL.dent[[2]]))]
                MCL.dent[[2]] = MCL.dent[[2]][(gsub('.*\\.', '', MCL.dent[[2]]) %in% 
                                                       gsub('.*\\.', '', MCL.dent[[1]]))]
        }

       subset.MCL %<>% SetAllIdent(id='clusters_6')
        print(table(subset.MCL@ident))
        gde.markers <- FindPairMarkers(subset.MCL, ident.1 = unlist(MCL.dent), 
                                       ident.2 = c(MCL.dent[[2]],MCL.dent[[1]]),
                                       logfc.threshold = 0.4,min.cells.group =1,
                                       return.thresh = 0.05)
        (mito.genes <- grep(pattern = "^MT-", x = gde.markers$gene))
        if(length(mito.genes)>0) gde.markers = gde.markers[-mito.genes,]
        write.csv(gde.markers, paste0(path,paste(pair,collapse = "_"),".csv"))
        GC();GC();GC();GC();GC();GC();GC();GC();
        subset.MCL <- ScaleData(subset.MCL)
        #---DoHeatmap.1----
        g <- DoHeatmap.1(subset.MCL, gde.markers, add.genes = markers, Top_n = 15,
                         ident.use = paste(paste(pair,collapse = " vs. "), "in MCL cells"),
                         group.label.rot = T,cex.row = 3,remove.key =F,title.size = 12)
        jpeg(paste0(path,"/DoHeatmap_",paste(pair,collapse = "_"),
                    ".jpeg"), units="in", width=10, height=7,res=600)
        print(g)
        dev.off()
        #---DoHeatmap vertical bar----
        g1 <- MakeCorlorBar(subset.MCL, gde.markers,Top_n = 15, add.genes = markers, do.print = F,
                            do.return = T)
        jpeg(paste0(path,"/DoHeatmap_",paste(pair,collapse = "_"),
                    "_legend.jpeg"),units="in", width=10, height=7,res=600)
        print(g1)
        dev.off()
} #20190129


control = "Pt-AA13-Ib-p"
(test <- paste0("test",7))
sample_n = which(df_samples$tests %in% test)
print(samples <- df_samples$sample[sample_n][2])

for(sample in samples){

        cell.use <- rownames(B_cells_MCL@meta.data)[B_cells_MCL@meta.data$orig.ident %in%
                                                            c(control, sample)]
        subset.MCL <- SubsetData(B_cells_MCL, cells.use = cell.use)
        
        #---SplitTSNEPlot----
        #subset.MCL <- SetAllIdent(subset.MCL, id= "singler1sub")
        SplitTSNEPlot(subset.MCL,do.return = FALSE,select.plots = 1:2,do.print = TRUE, do.label = F)
        subset.MCL <- SetAllIdent(subset.MCL,id = "clusters_6")
        SplitTSNEPlot(subset.MCL,do.return = FALSE,select.plots = 1:2,do.print = TRUE)
        
        #---FindAllMarkers.UMI----
        subset.MCL@meta.data$X6_clusters <- paste(subset.MCL@meta.data$orig.ident, 
                                                  subset.MCL@meta.data$X6_clusters, sep = ".")
        subset.MCL <- SetAllIdent(subset.MCL,id='X6_clusters')
        
        (MCL.dent.1 <- subset.MCL@ident %>% unique %>% as.character %>% sort)
        (len <- length(MCL.dent.1))
        subset.MCL@ident = factor(subset.MCL@ident, levels = MCL.dent.1)
        print(table(subset.MCL@ident))
        gde.markers <- FindPairMarkers(subset.MCL, ident.1 = MCL.dent.1[-6], 
                                       ident.2 = MCL.dent.1[c(7:11,1:5)],
                                       logfc.threshold = 0.4,min.cells.group =1,
                                       return.thresh = 0.05)
        write.csv(gde.markers, paste0(path,paste(c(control,sample),collapse = "_"),".csv"))
        (mito.genes <- grep(pattern = "^MT-", x = gde.markers$gene))
        if(length(mito.genes)>0) gde.markers = gde.markers[-mito.genes,]
        GC();GC();GC();GC();GC();GC();GC();GC();
        subset.MCL <- ScaleData(subset.MCL)
        #---DoHeatmap.1----
        g <- DoHeatmap.1(subset.MCL, gde.markers, add.genes = markers, Top_n = 15,
                         ident.use = paste(paste(c(control,sample),collapse = " vs. "), "in B and MCL cells"),
                         group.label.rot = T,cex.row = 3,remove.key =F,title.size = 12)
        jpeg(paste0(path,"/DoHeatmap_",paste(c(control,sample),collapse = "_"),
                    ".jpeg"), units="in", width=10, height=7,res=600)
        print(g)
        dev.off()
        #---DoHeatmap vertical bar----
        g1 <- MakeCorlorBar(subset.MCL, gde.markers,Top_n = 15, add.genes = markers, do.print = F,
                            do.return = T)
        jpeg(paste0(path,"/DoHeatmap_",paste(c(control,sample),collapse = "_"),
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
