########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(kableExtra)
source("../R/Seurat_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
########################################################################
#
#  1 Seurat Alignment 
# 
# ######################################################################
#======1.1 Setup the Seurat objects =========================
# Load the mouse.eyes dataset
# setup Seurat objects since both count matrices have already filtered
# cells, we do no additional filtering here
df_samples <- readxl::read_excel("doc/181002_Single_cell_sample list.xlsx")
#sample_n = which(df_samples$tests %in% c("test1", "test2", "test3","test4"))
sample_n = which(df_samples$samples %in% c("Pt-MD","Pt-RM"))
table(df_samples$tests)
df_samples[sample_n,] %>% kable() %>% kable_styling()
samples <- df_samples$samples[sample_n]
projects <- df_samples$projects[sample_n]
conditions <- df_samples$conditions[sample_n]

MCL_raw <- list()
MCL_Seurat <- list()

for(i in 1:length(samples)){
    MCL_raw[[i]] <- Read10X(data.dir = paste0("./data/",
                                                     samples[i],"/outs/filtered_gene_bc_matrices/hg19/"))
    colnames(MCL_raw[[i]]) <- paste0(conditions[i],
                                            "_",colnames(MCL_raw[[i]]))
    MCL_Seurat[[i]] <- CreateSeuratObject(MCL_raw[[i]],
                                                 min.cells = 3,
                                                 min.genes = 200,
                                                 project = projects[i],
                                                 names.delim = "_")
    MCL_Seurat[[i]]@meta.data$conditions <- conditions[i]
}
MCL_Seurat <- lapply(MCL_Seurat, FilterCells, 
                            subset.names = "nGene", 
                            low.thresholds = 200, 
                            high.thresholds = Inf) %>%
                lapply(NormalizeData) %>%
                lapply(ScaleData) %>%
                lapply(FindVariableGenes, do.plot = FALSE)

# we will take the union of the top 1k variable genes in each dataset for
# alignment note that we use 1k genes in the manuscript examples, you can
# try this here with negligible changes to the overall results
g <- lapply(MCL_Seurat, function(x) head(rownames(x@hvg.info), 1200))
genes.use <- unique(c(g[[1]],g[[2]]))
for(i in 1:length(conditions)){
    genes.use <- intersect(genes.use, rownames(MCL_Seurat[[i]]@scale.data))
}
length(genes.use)

#  Perform a canonical correlation analysis (CCA) =========================
# run a canonical correlation analysis to identify common sources
# of variation between the two datasets.
remove(MCL_raw)
GC()
MCL <- RunCCA(MCL_Seurat[[1]],MCL_Seurat[[2]],
                     genes.use = genes.use,
                     num.cc = 30)
save(MCL, file = "./data/MCL_alignment20181031.Rda")
load("./data/MCL_alignment20181031.Rda")
#======1.2 QC, pre-processing and normalizing the data=========================
QC_list <- lapply(MCL_Seurat, function(x) as.matrix(x = x@raw.data))
lapply(QC_list, function(x) median(colSums(x))) # Median nUMI
lapply(QC_list, function(x) median(apply(x,2,function(y) sum(length(y[y>0]))))) # Median nGene

lapply(QC_list, function(x) min(colSums(x))) # min nUMI
lapply(QC_list, function(x) min(apply(x,2,function(y) sum(length(y[y>0]))))) # min nGene


# CCA plot CC1 versus CC2 and look at a violin plot
p1 <- DimPlot(object = MCL, reduction.use = "cca", group.by = "conditions", 
              pt.size = 0.5, do.return = TRUE)
p2 <- VlnPlot(object = MCL, features.plot = "CC1", group.by = "conditions", 
              do.return = TRUE)
jpeg(paste0(path,"/S1_CCA_vplot.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p1, p2)
dev.off()

p3 <- MetageneBicorPlot(MCL, grouping.var = "conditions", dims.eval = 1:30, 
                        display.progress = FALSE) # run on cluster
jpeg(paste0(path,"/S1_CCA_number.jpeg"), units="in", width=10, height=7,res=600)
p3 + geom_smooth(method = 'loess')
dev.off()

PrintDim(object = MCL, reduction.type = "cca", dims.print = 1:2, genes.print = 10)

DimHeatmap(object = MCL, reduction.type = "cca", cells.use = 500, dim.use = c(1:3,11:20), 
           do.balanced = TRUE)

DimHeatmap(object = MCL, reduction.type = "cca", cells.use = 500, dim.use = 10:18, 
           do.balanced = TRUE)

#======1.4 filteration =========================
mito.genes <- grep(pattern = "^MT-", x = rownames(x = MCL@data), value = TRUE)
percent.mito <- Matrix::colSums(MCL@raw.data[mito.genes, ])/Matrix::colSums(MCL@raw.data)
MCL <- AddMetaData(object = MCL, metadata = percent.mito, col.name = "percent.mito")


g1 <- VlnPlot(object = MCL, features.plot = c("nGene", "nUMI", "percent.mito"), 
              nCol = 1,point.size.use = 0.2,
              x.lab.rot = T, do.return = T,return.plotlist =T)

MCL <- FilterCells(object = MCL, subset.names = c("nGene",  "nUMI", "percent.mito"), 
                          low.thresholds = c(500,1500, -Inf), 
                   high.thresholds = c(Inf,Inf, 0.15))

par(mfrow = c(1, 2))
GenePlot(object = MCL, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = MCL, gene1 = "nUMI", gene2 = "nGene")

g2 <- VlnPlot(object = MCL, features.plot = c("nGene", "nUMI", "percent.mito"), 
              nCol = 1,point.size.use = 0.2,
              x.lab.rot = T, do.return = T,return.plotlist =T)
jpeg(paste0(path,"/S1_nGene.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[1]]+ggtitle("nGene in raw data")+ 
                    scale_y_log10(limits = c(200,10000)),#+ylim(c(0,1000)),
                g2[[1]]+ggtitle("nGene after filteration")+ 
                    scale_y_log10(limits = c(200,10000))))
dev.off()
jpeg(paste0(path,"/S1_nUMI.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[2]]+ggtitle("nUMI in raw data")+ 
                    scale_y_log10(limits = c(1000,100000)),#+ylim(c(0,1000)),
                g2[[2]]+ggtitle("nUMI after filteration")+ 
                    scale_y_log10(limits = c(1000,100000))))
dev.off()
jpeg(paste0(path,"/S1_mito.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[3]]+ggtitle("mito % in raw data")+ 
                    ylim(c(0,0.4)),
                g2[[3]]+ggtitle("mito % after filteration")+ 
                    ylim(c(0,0.4))))
dev.off()

# filter var.ratio.pca

MCL <- CalcVarExpRatio(object = MCL, reduction.type = "pca",
                       grouping.var = "conditions", dims.use = 1:20)
MCL <- SubsetData(MCL, subset.name = "var.ratio.pca",accept.low = 0.5)

#======1.4 align seurat objects =========================
#Now we can run a single integrated analysis on all cells!
#Now we align the CCA subspaces, which returns a new dimensional reduction called cca.aligned
set.seed(42)
MCL <- AlignSubspace(object = MCL, reduction.type = "cca", grouping.var = "conditions", 
                     dims.align = 1:20)

MCL <- RunTSNE(object = MCL, reduction.use = "cca.aligned", dims.use = 1:20, 
                      do.fast = TRUE)
MCL <- FindClusters(object = MCL, reduction.type = "cca.aligned", dims.use = 1:20, 
                    resolution = 0.6, force.recalc = T, save.SNN = TRUE)

MCL <- RunPCA(object = MCL, pcs.compute = 30, do.print = TRUE, 
              pcs.print = 1:5, genes.print = 5)
p1 <- TSNEPlot(MCL, do.return = T, pt.size = 1, group.by = "conditions")
p2 <- TSNEPlot(MCL, do.label = F, do.return = T, pt.size = 1)

jpeg(paste0(path,"/TSNEplot_alignment.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p1, p2)
dev.off()

p3 <- TSNEPlot.1(object = MCL, do.label = T, group.by = "ident", 
           do.return = TRUE, no.legend = T, 
           pt.size = 1,label.size = 6 )+
    ggtitle("Tsne plot for all clusters")+
    theme(text = element_text(size=15),							
          plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 

jpeg(paste0(path,"/TSNEplot.jpeg"), units="in", width=10, height=7,res=600)
p3
dev.off()

SplitTSNEPlot(object = MCL)
#dev.off()
save(MCL, file = "./data/MCL_alignment20181031.Rda")
