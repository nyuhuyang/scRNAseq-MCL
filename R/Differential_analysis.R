########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
source("../R/Seurat_functions.R")

#3.1  Compare DE across all major cell types==================
#We would need the data for all clusters, as well the subclusters.
#detect changes in gene expression between young and aged, 
#in the different cell types and subtypes. 
#It will also be interesting to check if there is some subtype enriched in young compared to aged or viceversa. 

# 3.1.1 load data
# Rename ident
lnames = load(file = "./data/MCL_alignment.Rda")
lnames
table(MCL@ident)
idents <- as.data.frame(table(MCL@ident))
old.ident.ids <- idents$Var1
new.cluster.ids <- c("T cells",
                     "CD14 Monocytes",
                     "B cells",
                     "CD8 T cells",
                     "B cells",
                     "NK T cells",
                     "CD16 Monocytes",
                     "CD8 T cells",
                     "Dendritic Cells",
                     "Myeloid cells")

MCL@ident <- plyr::mapvalues(x = MCL@ident,
                             from = old.ident.ids,
                             to = new.cluster.ids)
MCL@meta.data$conditions <- sub("_",".",MCL@meta.data$conditions)

#===========answer to 6/1's email
MCL.subsets <- SplitCells(MCL)
MCL.subsets[[3]]
MCL.patient <- MCL.subsets[[1]]
MCL.normal <- MCL.subsets[[2]]
g1 <- featureplot(Treg[1:2],object=MCL.patient,do.return=T)
g2 <- featureplot(Treg[1:2],object=MCL.normal,do.return=T)
g <- list()
g[[1]] <- g1[[1]];g[[2]] <- g1[[2]];g[[3]] <- g2[[1]];g[[4]] <- g2[[2]];
do.call(plot_grid, g)

#===========answer to 6/1's email=====
Treg.markers.to.plot <- c("FOXP3","CTLA4","PDCD1","ENTPD1","CD38","ICOS",
                     "TNFSF9","TNFRSF9")
markers.to.plot <- HumanGenes(MCL,markers.to.plot,unique=T)
SplitDotPlotGG.1(object=MCL, genes.plot = rev(Treg.markers.to.plot),
                 cols.use = c("blue","red"), x.lab.rot = T, 
                 plot.legend = T, dot.scale = 8, do.return = T, 
                 grouping.var = "conditions")
#===========answer to 6/5's email=====
markers.to.plot <- c("CD79A","CD3G","CD8A","CD40","CXCR4")
markers.to.plot <- HumanGenes(MCL,markers.to.plot,unique=T)
SplitDotPlotGG.1(object=MCL, genes.plot = rev(markers.to.plot),
                 cols.use = c("blue","red"), x.lab.rot = T, 
                 plot.legend = T, dot.scale = 8, do.return = T, 
                 grouping.var = "conditions")
cellcycle.to.plot <- c("CCND1","CCND2","CCND3","CDK4","CDK6",
                       "PCNA","SOX11")
SplitDotPlotGG.1(object=MCL, genes.plot = rev(cellcycle.to.plot),
                 cols.use = c("blue","red"), x.lab.rot = T, 
                 plot.legend = T, dot.scale = 8, do.return = T, 
                 grouping.var = "conditions")

gde.all <- FindAllMarkersAcrossConditions(MCL,return.thresh = 0.01)
write.csv(gde.all,"./doc/MCL_patient_vs_normal.csv")
#===========answer to 6/17 & 6/21's email=====
Featureplot <- function(x,object = object,...){
        p <- FeaturePlot(object = object, 
                         reduction.use = "tsne",
                         features.plot = x, min.cutoff = NA, 
                         cols.use = c("lightgrey","blue"), ...)
        return(p)
}
markers.to.plot <- HumanGenes(MCL,c("CD79A","CD3G","CD8A","CD40","CXCR4","CEACAM1"))
MCL.subsets <- SplitCells(MCL)
MCL.subsets[[3]]
MCL.patient <- MCL.subsets[[1]]
MCL.normal <- MCL.subsets[[2]]
# Featureplot
Featureplot(cellcycle[-5],MCL.patient) # cellcycle
Featureplot(cellcycle[-5],MCL.normal) # cellcycle
Featureplot(T_Cell[1:6],MCL.patient) 
Featureplot(T_Cell[1:6],MCL.normal) 
Featureplot(Treg[1:6],MCL.patient) # Treg
Featureplot(Treg[1:6],MCL.normal) # Treg
Featureplot(Treg[5:10],MCL.patient) # Treg
Featureplot(Treg[5:10],MCL.normal) # Treg
Featureplot(c(CD4_Naive_T,"CD8A"),MCL.patient) 
Featureplot(c(CD4_Naive_T,"CD8A"),MCL.normal) 
Featureplot(NK,MCL.patient) # Treg
Featureplot(NK,MCL.normal) # Treg

Featureplot(B_Cell[1:6],MCL.patient) 
Featureplot(B_Cell[1:6],MCL.normal) 
Featureplot(c("CD1C","CD5","CD27","CD38","IL4R","CD40"),MCL.patient) 
Featureplot(c("CD1C","CD5","CD27","CD38","IL4R","CD40"),MCL.normal) 
Featureplot(c(Neutrophil,B_Cell[1:2]),MCL.patient) 
Featureplot(c(Neutrophil,B_Cell[1:2]),MCL.normal)

Featureplot(markers.to.plot,MCL.patient) 
Featureplot(markers.to.plot,MCL.normal)

Featureplot(c(CD14_Monocytes[-4],CD16_Monocytes),MCL.patient) 
Featureplot(c(CD14_Monocytes[-4],CD16_Monocytes),MCL.normal)


Featureplot(Stem_cell,MCL.patient) 
Featureplot(Stem_cell,MCL.normal)
#===========answer to 6/24's email=====
# Featureplot
Featureplot(c("CCND1","CDK4","CCND2","CDK6","CCND3","RB1"),MCL.patient,nCol = 2) # cellcycle
Featureplot(c("CCND1","CDK4","CCND2","CDK6","CCND3","RB1"),MCL.normal,nCol = 2) # cellcycle

Featureplot(c("E2F1","TK1","CCNA2","MKI67","PCNA","CDK1"),MCL.patient,nCol = 2) # cellcycle
Featureplot(c("E2F1","TK1","CCNA2","MKI67","PCNA","CDK1"),MCL.normal,nCol = 2) # cellcycle

Featureplot(c("SOX11","MS4A1","CD19","CD79A","CD5","CD40"),MCL.patient) # B cells
Featureplot(c("SOX11","MS4A1","CD19","CD79A","CD5","CD40"),MCL.normal) # B cells

Featureplot(c("CD22","PAX5","FCER2","CXCR4","CD27","IL4R"),MCL.patient) # B cells-2
Featureplot(c("CD22","PAX5","FCER2","CXCR4","CD27","IL4R"),MCL.normal) # B cells-2

Featureplot(c("CD3D","CD4","CD3G","CD8A","CD2","IL2RA"),MCL.patient) # T cells
Featureplot(c("CD3D","CD4","CD3G","CD8A","CD2","IL2RA"),MCL.normal) # T cells

Featureplot(c("IL7R","SELL","IL2RG","GIMAP5","GIMAP5","GIMAP5"),MCL.patient) # T cells -2
Featureplot(c("IL7R","SELL","IL2RG","GIMAP5","GIMAP5","GIMAP5"),MCL.normal) # T cells -2

Featureplot(c("CD14","LYZ","FCGR3A","S100A9","MS4A7","VMO1"),MCL.patient) # Monocytes
Featureplot(c("CD14","LYZ","FCGR3A","S100A9","MS4A7","VMO1"),MCL.normal) # Monocytes

Featureplot(c("ITGAM","CEACAM1","ITGAX","ITGAX","CD38","CD8B"),MCL.patient)
Featureplot(c("ITGAM","CEACAM1","ITGAX","ITGAX","CD38","CD8B"),MCL.normal)

Featureplot(DendriticCells[c(-6,-8)],MCL.patient)
Featureplot(DendriticCells[c(-6,-8)],MCL.normal)

Featureplot(c(Macrophages,Macrophages[1:2]),MCL.patient)
Featureplot(c(Macrophages,Macrophages[1:2]),MCL.normal)

#===========answer to 7/4's email=====
# Featureplot
markers.to.plot <- HumanGenes(MCL,c("EZH2","SETD2","EZH1","CTCF","NSD1","WHSC1L1"),unique=T)
Featureplot(markers.to.plot,MCL.patient,pt.size = 1)
Featureplot(markers.to.plot,MCL.normal,pt.size = 1)

markers.to.plot <- HumanGenes(MCL,c("TP53","BCL6","ATM","PIK3CA","Myc","PIK3CD"),unique=T)
Featureplot(markers.to.plot,MCL.patient,pt.size = 1)
Featureplot(markers.to.plot,MCL.normal,pt.size = 1)

#===========answer to 7/6's email=====
# Featureplot
markers.to.plot <- HumanGenes(MCL,c("PLK1","DNMT1","FOXO1","DNMT3A","FOXO3","DNMT3L"),unique=T)
Featureplot(markers.to.plot,MCL.patient,pt.size = 1)
Featureplot(markers.to.plot,MCL.normal,pt.size = 1)

grepl("FOXO",MCL@raw.data@Dimnames[1])

#==========7/10=======
Featureplot(interferon,MCL.patient,pt.size = 1)
Featureplot(interferon,MCL.normal,pt.size = 1)

Featureplot(DendriticCells,MCL.patient,pt.size = 1)
Featureplot(DendriticCells,MCL.normal,pt.size = 1)

Featureplot(c(CD14_Monocytes,CD16_Monocytes),MCL.patient,pt.size = 1)
Featureplot(c(CD14_Monocytes,CD16_Monocytes),MCL.normal,pt.size = 1)


#===========answer to 8/8 's email
Hpca_Blueprint_encode_main <- read.csv(file="../SingleR/output/Hpca_Blueprint_encode_main.csv")
Summary <- SearchAllMarkers(df = Hpca_Blueprint_encode_main, 
                            markers = c("CD3G","CD3D","CD8A"))
Summary

lnames = load(file = "./data/singler_MCL.RData")
lnames
# split singler
singler$seurat = MCL
singler.subsets <- SplitSingler(singler = singler)
singler.subsets[[length(singler.subsets)]] # levels of conditions
singler.patient = singler.subsets[[1]]
singler.normal = singler.subsets[[2]]
# main types-------
out = SingleR.PlotTsne.1(singler.patient$singler[[1]]$SingleR.single,
                         singler.patient$meta.data$xy,do.label=T,
                         do.letters = F,labels = singler.patient$singler[[1]]$SingleR.single$labels,
                         label.size = 5, dot.size = 2,do.legend = F,alpha = 1,
                         label.repel = T,force=2)
out$p+  ggtitle("Supervised labeling for MCL patient's major cell types")+#ggplot title
        theme(text = element_text(size=20),     #larger text including legend title
              plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) #title in middle

out = SingleR.PlotTsne.1(singler.normal$singler[[2]]$SingleR.single.main,
                         singler.normal$meta.data$xy,do.label=T,
                         do.letters = F,labels = singler.normal$singler[[2]]$SingleR.single.main$labels,
                         label.size = 5, dot.size = 2,do.legend = F,alpha = 1,
                         label.repel = T,force=2)
out$p+  ggtitle("Supervised labeling for normal samples's major cell types")+#ggplot title
        theme(text = element_text(size=20),     #larger text including legend title
              plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) #title in middle

# extract Singler_B_NK and MCL_B_NK========
# 1st run
table(singler$singler[[2]]$SingleR.single.main$labels)
B_NK_index1 <- singler$singler[[2]]$SingleR.single.main$labels[,1] %in% c('B-cells','NK cells')
table(B_NK_index1)

table(singler$singler[[1]]$SingleR.single$labels)

# 2nd run
B_NK_index2 = grepl("B_cell",singler$singler[[1]]$SingleR.single$labels[,1]) | 
        grepl("NK_cell",singler$singler[[1]]$SingleR.single$labels[,1])
table(B_NK_index2)

B_NK_index3 = B_NK_index1 & B_NK_index2
table(B_NK_index3)

# 3th run
MCL_B_NK = SubsetData(MCL, cells.use = MCL@cell.names[B_NK_index3])
remove1 <- FeaturePlot(MCL_B_NK, features.plot = "RPL4",do.identify = T) # manually select
write.csv(remove1,"./output/Cell_remove_20180808.csv",na = "")
B_NK_cells = B_NK_cells[!(B_NK_cells %in% remove1)]
length(B_NK_cells)

# get index
B_NK_index <- (MCL@cell.names %in% B_NK_cells)

# SubsetData
MCL_B_NK = SubsetData(MCL, cells.use = B_NK_cells)
Singler_B_NK <- SingleR.Subset.1(singler=singler,which(B_NK_index))

Singler_B_NK$seurat = MCL_B_NK
output <- SplitSingleR.PlotTsne(singler = Singler_B_NK, split.by = "conditions",main =FALSE,
                                return.plots=T,do.label=F,do.legend = F,legend.size = 15,
                                alpha = 1,label.repel = T, force=1)
output[[1]]$p+  ggtitle("B cells and NK cells in MCL patient")+#ggplot
        theme(text = element_text(size=20),     #larger text including legend title
              plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) #title in middle
output[[2]]$p+  ggtitle("B cells and NK cells in normal sample")+#ggplot title
        theme(text = element_text(size=20),     #larger text including legend title
              plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) #title in middle

# split Seurat
MCL.subsets = SplitSeurat(MCL_B_NK)
MCL.subsets[[3]]
MCL.patient = MCL.subsets[[1]]
MCL.normal = MCL.subsets[[2]]

# select cell manually
B_cells_normal <- FeaturePlot(MCL.normal, features.plot = "RPL4",pt.size = 3,
                              do.identify = T) # manually select
Cell_names_20180808 <- list("top_NK_cells_normal" = top_NK_cells_normal,
                            "bottom_NK_cells_normal" = bottom_NK_cells_normal,
                            "Normal B cells" = B_cells_normal,
                            "top_NK_cells_MCL" = top_NK_cells_MCL,
                            "bottom_NK_cells_MCL" = bottom_NK_cells_MCL,
                            "MCL B cells middle" = middle_B_cells,
                            "MCL B cells left" = left_B_cells,
                            "MCL B cells right" = right_B_cells)
Cell_names_20180808 <- list2df(Cell_names_20180808)
write.csv(Cell.names_20180808,"./output/Cell_names_20180808.csv",na = "")

Cell_names_20180808 <- df2list(Cell_names_20180808)
for(i in 1:length(Cell_names_20180808)){
        MCL_B_NK <- RenameIdent.1(MCL_B_NK, new.ident.name = names(Cell_names_20180808)[i],
                     cells.use = Cell_names_20180808[[i]])
}

TSNEPlot(object = MCL_B_NK, no.legend = TRUE, do.label = TRUE, pt.size = 2,
         do.return = TRUE, label.size = 5)+
        ggtitle("TSNE plot of four major cell types in normal sample")+
        theme(text = element_text(size=20),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5)) #title in middle

MCL_B_NK =  SubsetData(MCL_B_NK, ident.remove = "4")
top_NK <- SubsetData(MCL_B_NK,ident.use = c("top_NK_cells_normal","top_NK_cells_MCL"))
top_NK.markers <- FindAllMarkers.UMI(top_NK)

bottom_NK <- SubsetData(MCL_B_NK,ident.use = c("bottom_NK_cells_normal","bottom_NK_cells_MCL"))
bottom_NK.markers <- FindAllMarkers.UMI(bottom_NK)

MCL_B <- SubsetData(MCL_B_NK,ident.use = c("Normal B cells","MCL B cells middle",
                                           "MCL B cells left","MCL B cells right"))
MCL_B.markers <- FindAllMarkers.UMI(MCL_B)

top <- MCL_B.markers %>% group_by(cluster) %>% top_n(20, avg_logFC)
write.csv(MCL_B.markers,"./output/MCL_B_markers.csv")
# setting slim.col.label to TRUE will print just the cluster IDS instead of
# every cell name
DoHeatmap.1(object = MCL_B, genes.use = c(top$gene,"CCND1","CD5","CD19"), 
          slim.col.label = TRUE, remove.key = T,cex.row = 6,
          group.order = c("Normal B cells","MCL B cells left","MCL B cells middle",
                          "MCL B cells right"),
          rotate.key = T,group.label.rot = T)+
        ggtitle("Expression heatmap of top 15 DE genes")+
        theme(text = element_text(size=20),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5)) #title in middle
