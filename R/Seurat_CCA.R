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
#======1.1 read sample file =========================
# Load the mouse.eyes dataset
# setup Seurat objects since both count matrices have already filtered
# cells, we do no additional filtering here
df_samples <- readxl::read_excel("doc/181002_Single_cell_sample list.xlsx")
#sample_n = which(df_samples$tests %in% c("test1", "test2", "test3","test4"))
sample_n = which(df_samples$tests %in% c("test2","test3","test4"))
table(df_samples$tests)
df_samples[sample_n,] %>% kable() %>% kable_styling()
samples <- df_samples$samples[sample_n]
projects <- df_samples$projects[sample_n]
conditions <- df_samples$conditions[sample_n]
tests <- df_samples$tests[sample_n]
#======1.2 load  SingleCellExperiment =========================
(load(file = "./data/sce_list_5_20181107.Rda"))
names(sce_list)
MCL_Seurat <- lapply(sce_list, as.seurat) %>%
                lapply(NormalizeData) %>%
                lapply(ScaleData) %>%
                lapply(FindVariableGenes, do.plot = FALSE)

for(i in 1:length(samples)){
    MCL_Seurat[[i]]@meta.data$conditions <- conditions[i]
    MCL_Seurat[[i]]@meta.data$tests <- tests[i]
}
# we will take the union of the top 1k variable genes in each dataset for
# alignment note that we use 1k genes in the manuscript examples, you can
# try this here with negligible changes to the overall results
g <- lapply(MCL_Seurat, function(x) head(rownames(x@hvg.info), 800))
genes.use <- unique(unlist(g))
for(i in 1:length(conditions)){
    genes.use <- intersect(genes.use, rownames(MCL_Seurat[[i]]@data))
}
length(genes.use)

#========1.3 merge ===================================
MCL <- Reduce(function(x, y) MergeSeurat(x, y, do.normalize = F), MCL_Seurat)
remove(sce_list,MCL_Seurat);GC()


#======1.4 mito, QC, filteration =========================
mito.genes <- grep(pattern = "^MT-", x = rownames(x = MCL@data), value = TRUE)
percent.mito <- Matrix::colSums(MCL@raw.data[mito.genes, ])/Matrix::colSums(MCL@raw.data)
MCL <- AddMetaData(object = MCL, metadata = percent.mito, col.name = "percent.mito")

(load(file = "./data/MCL_g1.Rda"))

MCL <- FilterCells(object = MCL, subset.names = c("nGene","nUMI","percent.mito"),
                   low.thresholds = c(500,1000, -Inf), 
                   high.thresholds = c(Inf,Inf, 0.5))

par(mfrow = c(1, 2))
GenePlot(object = MCL, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = MCL, gene1 = "nUMI", gene2 = "nGene")

MCL@ident = factor(MCL@ident,levels = samples)

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
                    ylim(c(0,0.5)),
                g2[[3]]+ggtitle("mito % after filteration")+ 
                    ylim(c(0,0.5))))
dev.off()
# After removing unwanted cells from the dataset, the next step is to normalize the data.
MCL <- NormalizeData(object = MCL, normalization.method = "LogNormalize", 
                     scale.factor = 10000)
MCL <- FindVariableGenes(object = MCL, mean.function = ExpMean, 
                         dispersion.function = LogVMR, do.plot = FALSE, 
                         x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(MCL@var.genes)
#======1.5 1st run of pca-tsne  =========================
MCL <- ScaleData(object = MCL) %>%
    RunPCA() %>%
    FindClusters(dims.use = 1:20, force.recalc = T, print.output = FALSE) %>%
    RunTSNE()


#======1.5 Add Cell-cycle score =========================
# Read in a list of cell cycle markers, from Tirosh et al, 2015
cc.genes <- readLines(con = "./data/seurat_resources/regev_lab_cell_cycle_genes.txt")
s.genes <- HumanGenes(MCL,cc.genes[1:43])
g2m.genes <- HumanGenes(MCL,cc.genes[44:97])
MCL <- CellCycleScoring(object = MCL, s.genes = s.genes, g2m.genes = g2m.genes, 
                        set.ident = TRUE)
RidgePlot(object = MCL, features.plot = HumanGenes(MCL,c("CCND1","CDK4","CCND2","CDK6","CCND3","RB1")), 
          nCol = 2)
MCL@meta.data$CC.Difference <- MCL@meta.data$S.Score - MCL@meta.data$G2M.Score
head(x = MCL@meta.data)

#======1.6 check batch effect =========================

SingleFeaturePlot.1(MCL,"nUMI",threshold=10000)
SingleFeaturePlot.1(MCL,"nGene",threshold=2000)
TSNEPlot(MCL,group.by = "orig.ident")
SingleFeaturePlot.1(MCL,"percent.mito",threshold=0.05)
SingleFeaturePlot.1(MCL,"CC.Difference",threshold=0.05)
MCL <- ScaleData(object = MCL, 
                 model.use = "linear", do.par=T, do.center = T, do.scale = T,
                 vars.to.regress = c("orig.ident"),
                 display.progress = T)

# === 1.7 Perform a canonical correlation analysis (CCA) =========================
# run a canonical correlation analysis to identify common sources
# of variation between the two datasets.
MCL.subsets <- SplitSeurat(MCL, split.by = "orig.ident")
remove(MCL)
GC();GC();GC();GC();GC();GC();GC();GC();GC();GC()
MCL <- RunMultiCCA(object.list = MCL.subsets[1:length(samples)], 
                     genes.use = genes.use,
                     niter = 25, num.ccs = 30,
                     standardize =TRUE)
remove(MCL.subsets);GC()
#======1.2 QC, pre-processing and normalizing the data=========================
# CCA plot CC1 versus CC2 and look at a violin plot
# filter var.ratio.pca

MCL <- CalcVarExpRatio(object = MCL, reduction.type = "pca",
                       grouping.var = "orig.ident", dims.use = 1:20)
#MCL.copy <- MCL
MCL <- SubsetData(MCL, subset.name = "var.ratio.pca",accept.low = 0.5)

p1 <- DimPlot(object = MCL, reduction.use = "cca", group.by = "orig.ident", 
              pt.size = 0.5, do.return = TRUE)
p2 <- VlnPlot(object = MCL, features.plot = "CC1", group.by = "orig.ident", 
              do.return = TRUE)
jpeg(paste0(path,"/S1_CCA.jpeg"), units="in", width=10, height=7,res=600)
p1
dev.off()

jpeg(paste0(path,"/S1_CCA_vplot.jpeg"), units="in", width=10, height=7,res=600)
p2
dev.off()


PrintDim(object = MCL, reduction.type = "cca", dims.print = 1:2, genes.print = 10)
DimHeatmap(object = MCL, reduction.type = "cca", cells.use = 500, dim.use = c(1:3,11:20), 
           do.balanced = TRUE)
DimHeatmap(object = MCL, reduction.type = "cca", cells.use = 500, dim.use = 10:18, 
           do.balanced = TRUE)

#======1.7 align seurat objects =========================
#Now we can run a single integrated analysis on all cells!
#Now we align the CCA subspaces, which returns a new dimensional reduction called cca.aligned
set.seed(42)
MCL <- AlignSubspace(object = MCL, reduction.type = "cca", grouping.var = "orig.ident", 
                     dims.align = 1:20)

MCL <- RunTSNE(object = MCL, reduction.use = "cca.aligned", dims.use = 1:20, 
               do.fast = TRUE)
MCL <- FindClusters(object = MCL, reduction.type = "cca.aligned", dims.use = 1:20, 
                    resolution = 0.6, force.recalc = T, save.SNN = TRUE)

MCL <- RunPCA(object = MCL, pcs.compute = 30, do.print = TRUE, 
              pcs.print = 1:5, genes.print = 5)
p1 <- TSNEPlot(MCL, do.return = T, pt.size = 1, group.by = "orig.ident")
p2 <- TSNEPlot(MCL, do.label = F, do.return = T, pt.size = 1)

jpeg(paste0(path,"/TSNEplot_CCA_alignment.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p1, p2)
dev.off()
MCL <- SetAllIdent(MCL,id='res.0.6')
MCL@ident <- factor(MCL@ident, levels = 0:18)
p3 <- TSNEPlot.1(object = MCL, do.label = T, group.by = "ident", 
                 do.return = TRUE, no.legend = T, 
                 pt.size = 1,label.size = 6 )+
    ggtitle("Tsne plot for all clusters")+
    theme(text = element_text(size=15),							
          plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 

jpeg(paste0(path,"/TSNEplot.jpeg"), units="in", width=10, height=7,res=600)
p3
dev.off()

table(MCL@meta.data$orig.ident)

jpeg(paste0(path,"/SplitTSNEPlot.jpeg"), units="in", width=10, height=7,res=600)
SplitTSNEPlot(object = MCL,split.by = "orig.ident",do.return=F,
              select.plots = c(7,6,8,9,3,1,2,4,5))
dev.off()

save(MCL, file = "./data/MCL_CCA_20181110.Rda")
