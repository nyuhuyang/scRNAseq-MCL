########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
source("./R/Seurat_functions.R")
########################################################################
#
#  1 Seurat Alignment 
# 
# ######################################################################
#======1.1 Setup the Seurat objects =========================
# Load the mouse.eyes dataset

# setup Seurat objects since both count matrices have already filtered
# cells, we do no additional filtering here
MCL_raw <- list()
MCL_Seurat <- list()

samples <- c("Pt_DJ","Pt_1294")
projects <- c("EC-MDL-5092")
conditions <- c("normal", "MCL_patient")

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
                            low.thresholds = 500, 
                            high.thresholds = Inf)
MCL_Seurat <- lapply(MCL_Seurat, NormalizeData)
MCL_Seurat <- lapply(MCL_Seurat, ScaleData)
MCL_Seurat <- lapply(MCL_Seurat, FindVariableGenes, do.plot = FALSE)

# we will take the union of the top 1k variable genes in each dataset for
# alignment note that we use 1k genes in the manuscript examples, you can
# try this here with negligible changes to the overall results
g <- lapply(MCL_Seurat, function(x) head(rownames(x@hvg.info), 1000))
genes.use <- unique(c(g[[1]],g[[2]]))
for(i in 1:length(conditions)){
    genes.use <- intersect(genes.use, rownames(MCL_Seurat[[i]]@scale.data))
}
length(genes.use)

#======1.2 Perform a canonical correlation analysis (CCA) =========================
# run a canonical correlation analysis to identify common sources
# of variation between the two datasets.
remove(MCL_raw)
GC()
MCL <- RunCCA(MCL_Seurat[[1]],MCL_Seurat[[2]],
                     genes.use = genes.use,
                     num.cc = 30)
save(MCL, file = "./data/MCL_alignment.Rda")

# CCA plot CC1 versus CC2 and look at a violin plot
p1 <- DimPlot(object = MCL, reduction.use = "cca", group.by = "conditions", 
              pt.size = 0.5, do.return = TRUE)
p2 <- VlnPlot(object = MCL, features.plot = "CC1", group.by = "conditions", 
              do.return = TRUE)
plot_grid(p1, p2)

p3 <- MetageneBicorPlot(MCL, grouping.var = "conditions", dims.eval = 1:30, 
                        display.progress = FALSE) # run on cluster
p3 + geom_smooth(method = 'loess')

PrintDim(object = MCL, reduction.type = "cca", dims.print = 1:2, genes.print = 10)

DimHeatmap(object = MCL, reduction.type = "cca", cells.use = 500, dim.use = c(1:3,11:20), 
           do.balanced = TRUE)

DimHeatmap(object = MCL, reduction.type = "cca", cells.use = 500, dim.use = 10:18, 
           do.balanced = TRUE)

#======1.3 QC =========================
MCL <- CalcVarExpRatio(object = MCL, reduction.type = "pca",
                              grouping.var = "conditions", dims.use = 1:20)
MCL <- SubsetData(MCL, subset.name = "var.ratio.pca",accept.low = 0.5)

mito.genes <- grep(pattern = "^MT-", x = rownames(x = MCL@data), value = TRUE)
percent.mito <- Matrix::colSums(MCL@raw.data[mito.genes, ])/Matrix::colSums(MCL@raw.data)
MCL <- AddMetaData(object = MCL, metadata = percent.mito, col.name = "percent.mito")
#MCL <- ScaleData(object = MCL, genes.use = genes.use, display.progress = FALSE, 
#                         vars.to.regress = "percent.mito")
#Now we can run a single integrated analysis on all cells!
VlnPlot(object = MCL, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

MCL <- FilterCells(object = MCL, subset.names = c("nGene", "percent.mito"), 
                          low.thresholds = c(500, -Inf), high.thresholds = c(5000, 0.10))

par(mfrow = c(1, 2))
GenePlot(object = MCL, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = MCL, gene1 = "nUMI", gene2 = "nGene")

#Benjamin <- AddMetaData(object = Benjamin, metadata = percent.mito, col.name = "percent.mito")
#======1.4 align seurat objects =========================
#Now we align the CCA subspaces, which returns a new dimensional reduction called cca.aligned
set.seed(42)
MCL <- AlignSubspace(object = MCL, reduction.type = "cca", grouping.var = "conditions", 
                            dims.align = 1:20)
#Now we can run a single integrated analysis on all cells!

MCL <- FindClusters(object = MCL, reduction.type = "cca.aligned", dims.use = 1:20, 
                           resolution = 0.6, force.recalc = T, save.SNN = TRUE)

MCL <- RunTSNE(object = MCL, reduction.use = "cca.aligned", dims.use = 1:20, 
                      do.fast = TRUE)
MCL <- RunPCA(object = MCL, pcs.compute = 30, do.print = TRUE, 
              pcs.print = 1:5, genes.print = 5)
p1 <- TSNEPlot(MCL, do.return = T, pt.size = 1, group.by = "conditions")
p2 <- TSNEPlot(MCL, do.label = F, do.return = T, pt.size = 1)
#png('./output/TSNESplot_alignment.png')
plot_grid(p1, p2)
#dev.off()

TSNEPlot(object = MCL,do.label = TRUE, group.by = "ident", 
         do.return = TRUE, no.legend = TRUE,
         pt.size = 1,label.size = 8 )+
    ggtitle("TSNEplot for all cell clusters")+
    theme(text = element_text(size=20),     #larger text including legend title							
          plot.title = element_text(hjust = 0.5)) #title in middle

SplitTSNEPlot(object = MCL)
#dev.off()
save(MCL, file = "./data/MCL_alignment.Rda")
