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
(load(file="data/MCL_Harmony_20181121.Rda"))
markers.to.plot <-  HumanGenes(MCL,c("CD","LYZ","S100A9","CST3","CD68","FCER1A","FCGR3A","MS4A7","VMO1",
                     "CD2","CD3G","CD3D","CD8A","CD19","MS4A1","CD79A","CD40","CD22",
                     "FCER2"))

markers <-  HumanGenes(MCL,c("CD19","MS4A1","SOX11","ITGA4","CD79A","CCND1","CCND2","CD5","CD40"))
MCL <- SetAllIdent(object = MCL, id = "singler1sub")

df_samples <- readxl::read_excel("doc/181002_Single_cell_sample list.xlsx")
(samples <- df_samples$samples[6:9])

cell.use <- rownames(MCL@meta.data)[MCL@meta.data$orig.ident == "MD"]
MD.MCL <- SubsetData(MCL, cells.use = cell.use)

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
#========= test B cell markers=======
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

# DoHeatmap ================
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

p3 <- TSNEPlot.1(B_cells_MCL, do.return = T, pt.size = 0.5, group.by = "orig.ident",no.legend =T )
p4 <- TSNEPlot.1(B_cells_MCL, do.label = F, do.return = T, pt.size = 0.5, 
                 colors.use = ExtractMetaColor(B_cells_MCL), no.legend =T)
jpeg(paste0(path,"/S1_TSNEPlot.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p3, p4)
dev.off()


(samples <- df_samples$samples[6:9])

for(sample in samples){
    cell.use <- rownames(B_cells_MCL@meta.data)[B_cells_MCL@meta.data$orig.ident %in% c("MD",sample)]
    subset.MCL <- SubsetData(B_cells_MCL, cells.use = cell.use)
    subset.MCL <- SetAllIdent(subset.MCL,id = "orig.ident")
    test_markers <- FindAllMarkers.UMI(subset.MCL,logfc.threshold = 0.1, only.pos = T,
                                       test.use = "MAST")

    g <- DoHeatmap.1(subset.MCL,test_markers,Top_n = 50,
                     ident.use = paste("MD vs.",sample, "in B cells"),
                     group.label.rot = F,cex.row = 6,remove.key =T)
    jpeg(paste0(path,"_heatmap_MD_",sample,".jpeg"), units="in", width=10, height=7,
         res=600)
    print(g)
    dev.off()
}

# heatmap.2 ================
markers <-  HumanGenes(B_cells_MCL,c("CD19","MS4A1","CD79A","CD5","CD40","CDK4"))
(samples <- df_samples$samples[6:9])
for(sample in samples){
    cell.use <- rownames(B_cells_MCL@meta.data)[B_cells_MCL@meta.data$orig.ident %in% c("DJ",sample)] #
    subset.MCL <- SubsetData(B_cells_MCL, cells.use = cell.use)
    subset.MCL <- SetAllIdent(subset.MCL,id = "orig.ident")
    test_markers <- FindAllMarkers.UMI(subset.MCL,logfc.threshold = 0.1, only.pos = T,
                                       test.use = "MAST")
    top <- test_markers %>% group_by(cluster) %>% top_n(50, avg_logFC)
    y = subset.MCL@scale.data[unique(c(markers,top$gene)),]
    ## Column clustering (adjust here distance/linkage methods to what you need!)
    hc <- hclust(as.dist(1-cor(as.matrix(y), method="spearman")), method="complete")
    cc = gsub("_.*","",hc$labels)
    cc = gsub("MD","#B3DE69",cc)
    cc = gsub(sample,"#195016",cc)
    
    jpeg(paste0(path,"/Heatmap2_MD_",sample,".jpeg"), units="in", width=10, height=7,res=600)
    heatmap.2(as.matrix(y),
              Colv = as.dendrogram(hc), Rowv= FALSE,
              ColSideColors = cc, trace ="none",labCol = FALSE,dendrogram = "column",#scale="row",
              key.xlab = "scale log nUMI",
              cexRow = 0.5,
              margins = c(2,5),
              #scale = "row",
              breaks = seq(-3,3,length.out = 101),
              col = bluered,
              main = paste("MD vs.",sample, "in B cells"))
    par(lend = 1)           # square line ends for the color legend
    legend(0, 0.8,       # location of the legend on the heatmap plot
           legend = c("MD", sample), # category labels
           col = c("#B3DE69", "#195016"),  # color key
           lty= 1,             # line style
           lwd = 10            # line width
    )
    dev.off()
}

# 2. Compare Pt 17-C7 with C2, Compare Pt 11-C14 with C1 in tsne
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

# heat map of the top 50 genes that are differentially regulated in MCL cells in Pt17, C2 vs C14, and Pt 11, C1 vs C14.

for(test in tests){
    sample_n = which(df_samples$tests %in% test)
    df_samples[sample_n,] %>% kable() %>% kable_styling()
    samples <- df_samples$samples[sample_n]
    cell.use <- rownames(B_cells_MCL@meta.data)[B_cells_MCL@meta.data$orig.ident %in% samples]
    subset.MCL <- SubsetData(B_cells_MCL, cells.use = cell.use)
    subset.MCL <- SetAllIdent(subset.MCL,id = "orig.ident")
    test_markers <- FindAllMarkers.UMI(subset.MCL,logfc.threshold = 0.1, 
                       test.use = "MAST")
    write.csv(test_markers, paste0(path,test,"_DE_genes.csv"))
    g <- DoHeatmap.1(subset.MCL,test_markers,Top_n = 50,
                     ident.use = paste(paste(samples,collapse=" vs. "), "in B cells"),
              group.label.rot = T,cex.row = 6,remove.key =T)
    jpeg(paste0(path,"_heatmap_",test,".jpeg"), units="in", width=10, height=7,
         res=600)
    print(g)
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

#==================================================================

for(i in 1:length(test.markers)) {
    jpeg(paste0(path,test.markers[i],".jpeg"), units="in", width=10, height=7,
        res=600)
    p1 <- SingleFeaturePlot.1(object = MCL, feature = test.markers[i])
    print(p1)
    print(paste0(i,":",length(test.markers)))
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
                     