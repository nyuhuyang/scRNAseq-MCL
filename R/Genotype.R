library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(gplots)
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 2.1 identify phenotype for each cluster  ==========================================
(load(file="data/MCL_Harmony_12_20181121.Rda"))
markers.to.plot <-  HumanGenes(MCL,c("CD","LYZ","S100A9","CST3","CD68","FCER1A","FCGR3A","MS4A7","VMO1",
                     "CD2","CD3G","CD3D","CD8A","CD19","MS4A1","CD79A","CD40","CD22",
                     "FCER2"))
markers <-  HumanGenes(MCL,c("CD19","MS4A1","SOX11","ITGA4","CD79A","CCND1","CCND2","CD5","CD40"))

# SingleFeaturePlot.1 for Normal / MCL ================
markers <-  HumanGenes(MCL,c("CCND1","CCND2","CD19","CD3G","CDK4","CDK6","PCNA","CDK1"))
markers <-  HumanGenes(MCL,c("CD5","IL2RA","CD40","TNFRSF8","BTK"))
markers <-  HumanGenes(MCL,c("BTK"))
markers <-  HumanGenes(MCL,c("CCND3","RB1","TP53","ATM","MYC","MTAP","CDKN2a","BCL6",
                             "EZH2","EZH1","FOXO3","PRMT5"))
markers <-  HumanGenes(MCL,c("CDK4","CDK6","CCND1","CCND3","PCNA",
                             "CDK1","CD69","PRF1","SOX11","GZMB",
                             "CXCR4","CXCL10","PDCD1","CD274"))
df_samples <- readxl::read_excel("doc/181128_scRNAseq_info.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
tests <- paste0("test",c(3:4))
control = "MD"
for(test in tests){
    sample_n = which(df_samples$tests %in% test)
    samples <- unique(df_samples$sample[sample_n])
    print(paste(c(control,samples), collapse = " "))
    
    cell.use <- rownames(MCL@meta.data)[MCL@meta.data$orig.ident %in% c(control,samples)]
    subset.MCL <- SubsetData(MCL, cells.use = cell.use)
    SplitSingleFeaturePlot(subset.MCL, 
                           #select.plots = c(1,5,4,3,2),
                           group.by = "ident",split.by = "orig.ident",
                           no.legend = T,label.size=3,do.print =T,markers = markers,
                           threshold = 0.1)
}

# Identify MCL marker genes===================
(samples <- df_samples$samples[6:9])
MCL_markers_list <- list()
for(i in 1:length(samples)){
    cell.use <- rownames(MCL@meta.data)[MCL@meta.data$orig.ident %in% c("MD",samples[i])]
    subset.MCL <- SubsetData(MCL, cells.use = cell.use)
    subset.MCL <- SetAllIdent(subset.MCL,id = "orig.ident")
    MCL_markers_list[[i]] <- FindAllMarkers.UMI(subset.MCL,logfc.threshold = 0.2, only.pos = T,
                                                test.use = "MAST")
}
MCL_markers_list <- lapply(MCL_markers_list, function(df) {df[!(df$cluster == "MD"),]})
MCL_markers <- do.call(rbind.data.frame , MCL_markers_list)
MCL_markers <- MCL_markers[,-(1:2)]
MCL_markers <- MCL_markers[order(MCL_markers$avg_logFC,decreasing = T),]
write.csv(MCL_markers, paste0(path,"MCL_markers.csv"))

#==================================================================

for(i in 1:length(test.markers)) {
    jpeg(paste0(path,test.markers[i],".jpeg"), units="in", width=10, height=7,
        res=600)
    p1 <- SingleFeaturePlot.1(object = MCL, feature = test.markers[i])
    print(p1)
    print(paste0(i,":",length(test.markers)))
    dev.off()
}

# 2. SingleFeaturePlot.1 for Normal / MCL.1 / MCL.2 ================
tests <- c("test3","test4")

for(test in tests){
    sample_n = which(df_samples$tests %in% test)
    df_samples[sample_n,] %>% kable() %>% kable_styling()
    samples <- df_samples$samples[sample_n]
    cell.use <- rownames(MCL@meta.data)[MCL@meta.data$orig.ident %in% samples]
    subset.MCL <- SubsetData(MCL, cells.use = cell.use)
    subset.MCL <- SplitSeurat(subset.MCL, split.by = "orig.ident")
    lvl <- subset.MCL[[3]]
    
    for(marker in markers){
        g <- list()
        g[[1]] <- SingleFeaturePlot.1(object = subset.MCL[[1]], threshold=0.5,
                                      feature = marker,title = lvl[1])
        g[[2]] <- SingleFeaturePlot.1(object = subset.MCL[[2]],  threshold=0.5,
                                      feature = marker,title = lvl[2])
        jpeg(paste0(path,"_Splited_",test,"_",marker,".jpeg"), units="in", width=10, height=7,
             res=600)
        print(do.call(plot_grid, g))
        print(paste0(which(markers == marker),":",length(markers)))
        dev.off()
    }
}

#######################
# FACS
#######################
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
#######################
# geom_density
#######################
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



#######################
# FACS
#######################
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
#######################
# geom_density
#######################
markers <-  HumanGenes(MCL,c("CCND1","CD5","CD19","CDK4","MS4A1","SOX11"))
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
#============Selina's 10/31 email=====

markers <-  HumanGenes(MCL,c("CD19","MS4A1","ITGA4","CD79A","CCND1","CCND2","CD5","SOX11"))
all.markers <- sort(unique(c(markers.to.plot,test.markers,markers)))

Q3 <- HumanGenes(MCL,c("CD72","HLA-DQA1","LILRA4","STMN1"))

object.subsets <- SplitSeurat(object = MCL, split.by = "orig.ident")
levels <- object.subsets[[length(object.subsets)]]
levels


for(marker in markers){
    g <- list()
    for(i in 1:length(samples)){
        g[[i]] <- SingleFeaturePlot.1(object = object.subsets[[i]], 
                                      feature = marker,title =levels[i])
    }
    jpeg(paste0(path,"Split_",marker,".jpeg"), units="in", width=10, height=7,
         res=600)
    print(do.call(plot_grid, g))
    print(paste0(which(all.markers == marker),":",length(all.markers)))
    dev.off()
}


# two normals (DJ and MD)
sample_n = which(df_samples$Tests %in% "test2")
samples <- df_samples$Samples[sample_n]
print(samples)
cell.use <- rownames(MCL@meta.data)[MCL@meta.data$orig.ident %in% samples]
subset.MCL <- SubsetData(MCL, cells.use = cell.use)

object.subsets <- SplitSeurat(object = subset.MCL, split.by = "orig.ident")
levels <- object.subsets[[length(object.subsets)]]

for(marker in markers){
    g <- list()
    for(i in 1:length(levels)){
        g[[i]] <- SingleFeaturePlot.1(object = object.subsets[[i]], 
                                      feature = marker,title =samples[i])
    }
    jpeg(paste0(path,"Split_",marker,"_MCL.jpeg"), units="in", width=10, height=7,
         res=600)
    print(do.call(plot_grid, g))
    dev.off()
}

markers <-  HumanGenes(MCL,c("CCND1","SOX11","MS4A1","FCER2","CD72","HLA-DQA1","STMN1"))
object.subsets <- SplitSeurat(object = MCL, split.by = "orig.ident")
levels <- object.subsets[[length(object.subsets)]]

p <- FeaturePlot.1(object.subsets[[2]],x = markers)
p1 <- do.call(plot_grid, p)
p1 <- p1 + ggtitle(levels[2])+
    theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"))
jpeg(paste0(path,levels[2],".jpeg"),
     units="in", width=10, height=7,res=600)
print(p1)
dev.off()
#====== 2.2 marker gene analysis ==========================================
blueprint_encode_main = read.csv("../SingleR/output/blueprint_encode_main.csv",row.names =1,header = T,
                     stringsAsFactors = F)
marker.list <- df2list(blueprint_encode_main)
marker.list <- lapply(marker.list, function(x) HumanGenes(MCL,x[1:18]) %>% .[1:9])

marker.list %>% list2df %>% t %>% kable() %>% kable_styling()

FeaturePlot.1 <- function(object = MCL, x){
    p <- FeaturePlot(object = object, 
                     reduction.use = "tsne",
                     features.plot = x, min.cutoff = NA, do.return =T,
                     cols.use = c("lightgrey","blue"), pt.size = 0.5)
    return(p)
}

dev.off()
for(i in 1:length(marker.list)){
    p <- FeaturePlot.1(object = MCL, x = marker.list[[i]][1:9])
    p1 <- do.call(plot_grid, p)
    p1 <- p1 + ggtitle(paste(names(marker.list)[i],"markers"))+
        theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"))
    jpeg(paste0(path,names(marker.list)[i],".jpeg"),
         units="in", width=10, height=7,res=600)
    print(p1)
    print(paste0(i,":",length(marker.list)))
    dev.off()
}
