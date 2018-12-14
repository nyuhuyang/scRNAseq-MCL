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
markers <-  HumanGenes(MCL,c("CD28"))
markers <-  HumanGenes(MCL,c("CCND3","RB1","TP53","ATM","MYC","MTAP","CDKN2a","BCL6",
                             "EZH2","EZH1","FOXO3","PRMT5"))

df_samples <- readxl::read_excel("doc/181128_scRNAseq_info.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
tests <- paste0("test",c(3,4))
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
                           threshold = 0.5)
}


#  MD compare all others pairwisely in tsne
for(sample in samples){
    
    cell.use <- rownames(MCL@meta.data)[MCL@meta.data$orig.ident == sample]
    subset.MCL <- SubsetData(MCL, cells.use = cell.use)
    
    for(marker in markers){
        g <- list()
        g[[1]] <- SingleFeaturePlot.1(object = MD.MCL, threshold=0.1,
                                      feature = marker,title = "MD")
        g[[2]] <- SingleFeaturePlot.1(object = subset.MCL, threshold=0.1,
                                      feature = marker,title = sample)
        jpeg(paste0(path,"Splited_MD_",sample,"_",marker,".jpeg"), units="in", width=10, height=7,
             res=600)
        print(do.call(plot_grid, g))
        print(paste0(which(markers == marker),":",length(markers)))
        dev.off()
    }
}
#========= test markers for Normal / MCL ================
blueprint_encode_main <- read_csv("../SingleR/output/blueprint_encode_main.csv")
(samples <- df_samples$samples[6:9])
(B_markers <- HumanGenes(MCL, blueprint_encode_main$B_cells[1:15]))

for(sample in samples){
    
    cell.use <- rownames(MCL@meta.data)[MCL@meta.data$orig.ident == sample]
    subset.MCL <- SubsetData(MCL, cells.use = cell.use)
    
    for(marker in markers){
        g <- list()
        g[[1]] <- SingleFeaturePlot.1(object = MD.MCL, threshold=0.1,
                                      feature = marker,title = "MD")
        g[[2]] <- SingleFeaturePlot.1(object = subset.MCL, threshold=0.1,
                                      feature = marker,title = sample)
        jpeg(paste0(path,"Splited_MD_",sample,"_",marker,".jpeg"), units="in", width=10, height=7,
             res=600)
        print(do.call(plot_grid, g))
        print(paste0(which(markers == marker),":",length(markers)))
        dev.off()
    }
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
immgen_main = read.csv("../SingleR/output/immgen_main.csv",row.names =1,header = T,
                       stringsAsFactors = F)
marker.list <- df2list(immgen_main)
marker.list <- lapply(marker.list, function(x) HumanGenes(MCL,x))

marker.list %>% list2df %>% head(15) %>% kable() %>% kable_styling()

FeaturePlot.1 <- function(object = MCL, x){
    p <- FeaturePlot(object = object, 
                     reduction.use = "tsne",
                     features.plot = x, min.cutoff = NA, do.return =T,
                     cols.use = c("lightgrey","blue"), pt.size = 0.5)
    return(p)
}

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




#--Hematopoietic----
Hematopoietic <- HumanGenes(MCL,c("PTPRC","LAPTM5","SRGN"))
#------Myeloid----
megakaryocytes <-  HumanGenes(MCL,c("PPBP","GNG11"))
erythrocyte <-  HumanGenes(MCL,c("HBA2","HBB"))
MastCells <- HumanGenes(MCL,c("Cma1","Mcpt4","Tpsb2","Cpa3"))
Neutrophil <- HumanGenes(MCL,c("ADAM8","MSMO1","FUT4","FCGR3A","CEACAM8"))
CD14_Monocytes <-  HumanGenes(MCL,c("CD14","LYZ","S100A9","CCL2","CCR2"))
CD16_Monocytes <- HumanGenes(MCL,c("FCGR3A","MS4A7","VMO1","CCR2"))
Macrophages <- HumanGenes(MCL,c("LYZ","CD68","MARCO","EMR1"))
DendriticCells <- HumanGenes(MCL,c("Itgax","TGAM","GPR183","CST3","HLA-DQA1","FCER1A","TSPAN13",
                                     "IL3RA","IGJ","TLR7","ZBTB46","CD1C"))
interferon <- HumanGenes(MCL,c("IFNG","IFNA1","IFNGR1","IRF1","IRF7"))
Myeloid <-  HumanGenes(MCL,c(megakaryocytes,erythrocyte,MastCells,
                             CD14_Monocytes,CD16_Monocytes,Macrophages,DendriticCells),unique=T)
#------Lymphoid----
Lymphoid <- HumanGenes(MCL,c("Cd19","CD79A","MS4A1",
                               "GNLY","Ncr1","CCL5","KLRD1","NKG7"))
# T cell
T_Cell <- HumanGenes(MCL,c("CD2","CD3G","CD3D","CD4","CD8A","IL2RA","FOXP3",
                           "IL7R","SELL","IL2RG","GIMAP5"))

Treg <- HumanGenes(MCL,c("FOXP3","CD4","IL2RA","CTLA4","PDCD1","ENTPD1","CD38",
                         "ICOS","TNFSF9","TNFRSF9"))
CD4_Naive_T <- HumanGenes(MCL,c("CD4","IL7R","GIMAP5","SELL","IL2RG"))
NK <- HumanGenes(MCL,c("NKG7","CCL5","NCAM1","FCGR3A","Ncr1","KLRD1"))
# B cell
B_Cell <-HumanGenes(MCL,c("CD19","MS4A1","CD79A","CD40","CD22","FCER2","HLA-DRB1",
                          "CXCR4","SOX11","CD5","PAX5","CD27","IL4R"))

B_StemCell <- HumanGenes(MCL,c("SPN","MS4A1"))
Pre_Pro_B <- HumanGenes(MCL,c("CD34","MME","CD38"))
Pro_B <- HumanGenes(MCL,c("MME","CD19","SPN","CD38","CD24","IL7","IL3RA"))
Pre_B <- HumanGenes(MCL,c("MME","CD19","MS4A1","CD24","CD38","IL7","IL3RA","IL4R"))
Immature_B <- HumanGenes(MCL,c("MME","CD19","MS4A1","CR2","CD40","CD24","CD38","IL4R"))
Transitional_B <- HumanGenes(MCL,c("CD19","MS4A1","CD5","CR2","CD24","CD38"))
Marginal_zone_B <- HumanGenes(MCL,c("CD1C","CD19","MS4A1","CR2","CD27"))
Regulatory_B <- HumanGenes(MCL,c("CD1D","CD5","CD19","CR2","CD24"))
Follicular_B <- HumanGenes(MCL,c("CD19","MS4A1","CR2","CD22","FCER2","CD24",
                                   "HLA-DRB1","HLA-DQB1","HLA-DRA","HLA-DQA1"))
Activated_B <- HumanGenes(MCL,c("CD27","CD19","MS4A1","IL2RA","TNFRSF8","CD69","CD80","CD86","FLT3"))
Germinal_center_B <- HumanGenes(MCL,c("MME","CD19","MS4A1","FCER2","CD27","CD38","TNFRSF17"))
Plasma_blast <- HumanGenes(MCL,c("CD19","CD38","CD27","TNFRSF17","HLA-DRB1"))
Plasma_cell_long_lived <- HumanGenes(MCL,c("CXCR4","CD27","CD38","CD138","CD269"))
Memory_B <- HumanGenes(MCL,c("CD19","MS4A1","CD40","CD27","CXCR4","CXCR5","ACKR3"))


Melanocytes <- HumanGenes(MCL,c("Pmel","Mlana"))
Mesenchymal <- HumanGenes(MCL,c("Pdgfrb","Vim","Has2","Dcn"))
Myelinating_Schwann_cells <- HumanGenes(MCL,c("MBP","MPZ"))
Pericytes <- HumanGenes(MCL,c("Pdgfrb","Cspg4","Anpep","Rgs5",
                                "Myh11","Mylk","Des","Vtn","Ifitm1"))
Smooth_muscle_cells <- HumanGenes(MCL,c("Acta2","Myh11"))
Stem_cell <- HumanGenes(MCL,c("POU5F1","FUT4","CD34","PROM1","ABCG2","Runx1","ATXN1",
                                "Nes","NCAM","NGFR"))
#GP38,PODP
Stromal_fibroblasts <- HumanGenes(MCL,c("DCN","COL6A1","TIMP3","PDGFRA"))
Neurons <- HumanGenes(MCL,c("Ihh","Gli1", "Ptch1", "Hhip"))
cellcycle <- HumanGenes(MCL,c("CCND1","CCND2","CCND3","CDK4","CDK6","PCNA","SOX11",
                              "RB1","E2F1","TK1","CCNA2","MKI67","CDK1"),unique = T)



markers.to.plot <- c(CD14_Monocytes[c(1:3)],DendriticCells[c(3,5)],
                     Macrophages[2],Natural_killer_T[2],
                     CD16_Monocytes[1:3],T_Cell[1:4],B_Cell[1:6])
markers.to.plot <- c("CD14","LYZ","S100A9","CST3","CD68","FCER1A","FCGR3A","MS4A7","VMO1",
                     "CD2","CD3G","CD3D","CD8A","CD19","MS4A1","CD79A","CD40","CD22",
                     "FCER2")
                     