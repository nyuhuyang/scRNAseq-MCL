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
source("../R/SingleR_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#3.1  Compare DE across all major cell types==================
#We would need the data for all clusters, as well the subclusters.
#detect changes in gene expression between young and aged, 
#in the different cell types and subtypes. 
#It will also be interesting to check if there is some subtype enriched in young compared to aged or viceversa. 

# 3.1.1 load data
# Rename ident
(load(file="data/MCL_Harmony_20181121.Rda"))
MCL <- SetAllIdent(MCL, id="res.0.6")
table(MCL@ident)
TSNEPlot.1(MCL,do.label = T)
B_cells_MCL <- SubsetData(MCL, ident.use = c(0,2,3,7,8,10,12))
TSNEPlot.1(B_cells_MCL,do.label = T)
B_cells_MCL <- SetAllIdent(B_cells_MCL, id="singler2main")
table(B_cells_MCL@ident)

B_cells_MCL <- SubsetData(B_cells_MCL, ident.use = c("B_cells","HSC","MCL"))
table(B_cells_MCL@meta.data$singler1main)
B_cells_MCL <- SetAllIdent(B_cells_MCL, id="singler1main")
B_cells_MCL <- SubsetData(B_cells_MCL, ident.use = c("B_cells"))

# Compare FACS data using GenePlot =======
df_samples <- readxl::read_excel("doc/181002_Single_cell_sample list.xlsx")
sample_n = which(df_samples$tests %in% paste0("test",2:4))
table(df_samples$tests)
df_samples[sample_n,] %>% kable() %>% kable_styling()

tests <- paste0("test",2:4)
for(test in tests){
        sample_n = which(df_samples$tests %in% test)
        samples <- df_samples$samples[sample_n]
        print(samples)
        
        cell.use <- rownames(MCL@meta.data)[MCL@meta.data$orig.ident %in% samples]
        subset.MCL <- SubsetData(MCL, cells.use = cell.use)
        
        cells_use_list <- split(row.names(subset.MCL@meta.data), subset.MCL@meta.data[,"orig.ident"])
        G1 <- lapply(cells_use_list, function(cells_use) {
                single.MCL <- SubsetData(subset.MCL, cells_use_list[[2]])
                GenePlot.1(single.MCL, gene1 = "CD5", gene2 = "CD19",use.raw = F, 
                                 title = unique(single.MCL@meta.data$orig.ident))
        })
        G2 <- lapply(cells_use_list, function(cells_use) {
                single.MCL <- SubsetData(subset.MCL, cells_use)
                GenePlot.1(single.MCL, gene1 = "CD3E", gene2 = "CD19",use.raw = F, 
                                 title = unique(single.MCL@meta.data$orig.ident))
        })
        
        jpeg(paste0(path,test,"_CD5_CD3E_CD19.jpeg"), units="in", width=8, height=7,
             res=600)
        #print(plot_grid(G1[[1]],G1[[5]],G1[[4]],G1[[3]],G1[[2]],
        #          G2[[1]],G2[[5]],G2[[4]],G2[[3]],G2[[2]],nrow = 2))
        print(plot_grid(G1[[1]],G1[[2]],G2[[1]],G2[[2]],nrow = 2))
        dev.off()
}
#  count cells in geneplot=============================
CountCells <- function(object = subset.MCL, gene1="CD5", gene2 = "CD19", split.by ="orig.ident"){
        cells_use_list <- split(row.names(subset.MCL@meta.data), subset.MCL@meta.data[,split.by])
        cell_counts <- lapply(cells_use_list, function(cells_use) {
                single.MCL <- SubsetData(subset.MCL, cells_use)
                c(length(which(single.MCL@data[gene1,] == 0 & single.MCL@data[gene2,] == 0)),
                  length(which(single.MCL@data[gene1,] > 0 & single.MCL@data[gene2,] == 0)),
                  length(which(single.MCL@data[gene1,] == 0 & single.MCL@data[gene2,] > 0)),
                  length(which(single.MCL@data[gene1,] > 0 & single.MCL@data[gene2,] > 0)))/
                        length(single.MCL@cell.names)*100
        })
        df_cell_counts <- list2df(cell_counts)
        rownames(df_cell_counts) <- c("l.l","h.l","l.h","h.h")
        return(df_cell_counts)
}

tests <- paste0("test",2:4)
for(test in tests){
        sample_n = which(df_samples$tests %in% test)
        samples <- df_samples$samples[sample_n]
        print(samples)
        
        cell.use <- rownames(MCL@meta.data)[MCL@meta.data$orig.ident %in% samples]
        subset.MCL <- SubsetData(MCL, cells.use = cell.use)
        CD5_CD19 <- CountCells(subset.MCL, "CD5","CD19")
        CD3E_CD19 <- CountCells(subset.MCL, "CD3E","CD19")
        df <- merge(CD5_CD19,CD3E_CD19, by="row.names",all.x=TRUE)
        write.csv(df, paste0(path,test,"_countcell.csv"),row.names = F)
}

#  demonstrate cell cycle in geneplot=============================
tests <- paste0("test",3:4)
for(test in tests){
        sample_n = which(df_samples$tests %in% test)
        samples <- df_samples$samples[sample_n]
        print(samples)
        
        cell.use <- rownames(MCL@meta.data)[MCL@meta.data$orig.ident %in% samples]
        subset.MCL <- SubsetData(MCL, cells.use = cell.use)
        
        cells_use_list <- split(row.names(subset.MCL@meta.data), subset.MCL@meta.data[,"orig.ident"])
        G1 <- lapply(cells_use_list, function(cells_use) {
                single.MCL <- SubsetData(subset.MCL, cells_use)
                GenePlot.1(single.MCL, gene1 = "G2M.Score", gene2 = "S.Score",use.raw = F, 
                           title = unique(single.MCL@meta.data$orig.ident))+
                        xlim(0, 1.25)+ylim(0, 0.8)
        })
        jpeg(paste0(path,test,"_cellcyle_geneplot.jpeg"), units="in", width=8, height=7,
             res=600)
        #print(plot_grid(G1[[1]],G1[[5]],G1[[4]],G1[[3]],G1[[2]]))
        print(do.call(plot_grid,G1))
        dev.off()
}

# geom_density  ===========
markers <-  HumanGenes(MCL,c("CCND1","CD5","CD19","CDK4","MS4A1","SOX11"))
MCL <- B_cells_MCL
tests <- paste0("test",c(1,4))
for(test in tests){
        sample_n = which(df_samples$tests %in% test)
        samples <- unique(df_samples$samples[sample_n])
        print(samples)
        
        cell.use <- rownames(MCL@meta.data)[MCL@meta.data$orig.ident %in% samples]
        subset.MCL <- SubsetData(MCL, cells.use = cell.use)
        
        g <- split(rownames(subset.MCL@meta.data), subset.MCL@meta.data[,"orig.ident"]) %>% lapply(function(cells_use) {
                single.MCL <- SubsetData(subset.MCL, cells.use = cells_use)
                sample <- unique(single.MCL@meta.data$orig.ident)
                data.use <- single.MCL@data[markers,] %>% as.matrix %>% t %>% as.data.frame %>%
                       gather(key = markers, value = ave.expr)
                ggplot(data.use, aes(x = ave.expr, color = markers)) + 
                        geom_density(size = 1) +
                        scale_y_sqrt() + ylim(0, 1)+
                        xlab("Average expression (log nUMI)")+
                        ggtitle(sample)+
                        theme(text = element_text(size=15),
                              #legend.position="none", 
                              legend.position=c(0.3,0.85) ,
                              plot.title = element_text(hjust = 0.5,size = 15, face = "bold"))
        })
        jpeg(paste0(path,"density_",test,".jpeg"), units="in", width=10, height=7,res=600)
        print(do.call(plot_grid,c(g,nrow = 1))+ #   plot_grid(g[[1]],g[[5]],g[[4]],g[[3]],g[[2]])
                ggtitle("Density plot for typical markers in MCL B cells")+
                theme(text = element_text(size=15),							
                      plot.title = element_text(hjust = 0.5,size = 15, face = "bold")))
        dev.off()
}

#======== compare threshold ===============
sample = "Pt-1475"
marker = "MS4A1"
thresholds = c(0.01,0.1,0.5,1,2,3,4,5)
cell.use <- rownames(MCL@meta.data)[MCL@meta.data$orig.ident == sample]
subset.MCL <- SubsetData(MCL, cells.use = cell.use)

for(threshold in thresholds){
        g <- list()
        g[[1]] <- SingleFeaturePlot.1(object = MD.MCL, threshold=threshold,
                                      feature = marker,title = "MD")
        g[[2]] <- SingleFeaturePlot.1(object = subset.MCL, threshold=threshold,
                                      feature = marker,title = sample)
        jpeg(paste0(path,"Splited_MD_",sample,"_",marker,"_",threshold,".jpeg"), units="in", width=10, height=7,
             res=600)
        print(do.call(plot_grid, g)+
                      ggtitle(paste("threshold =",threshold))+
                      theme(text = element_text(size=20),							
                            plot.title = element_text(hjust = 0.5,size = 15, face = "bold")))
        print(paste0(which(threshold == thresholds),":",length(thresholds)))
        dev.off()
}

# Identify MCL marker genes===================
(samples <- df_samples$samples[6:9])
MCL_markers_list <- list()
for(i in 1:length(samples)){
        cell.use <- rownames(B_cells_MCL@meta.data)[B_cells_MCL@meta.data$orig.ident %in% c("MD",samples[i])]
        subset.MCL <- SubsetData(B_cells_MCL, cells.use = cell.use)
        subset.MCL <- SetAllIdent(subset.MCL,id = "orig.ident")
        MCL_markers_list[[i]] <- FindAllMarkers.UMI(subset.MCL,logfc.threshold = 0.2, only.pos = T,
                                           test.use = "MAST")
}
MCL_markers_list <- lapply(MCL_markers_list, function(df) {df[!(df$cluster == "MD"),]})
MCL_markers <- do.call(rbind.data.frame , MCL_markers_list)
MCL_markers <- MCL_markers[,-(1:2)]
MCL_markers <- MCL_markers[order(MCL_markers$avg_logFC,decreasing = T),]
write.csv(MCL_markers, paste0(path,"MCL_markers.csv"))


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

