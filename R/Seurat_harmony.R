########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

suppressPackageStartupMessages({
        library(Seurat)
        library(magrittr)
        library(harmony)
        library(dplyr)
        library(kableExtra)
        source("../R/Seurat_functions.R")
})
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
########################################################################
#
#  1 harmony Alignment 
# 
# ######################################################################
#======1.1 read sample file =========================
# Load the mouse.eyes dataset
# setup Seurat objects since both count matrices have already filtered
# cells, we do no additional filtering here
df_samples <- readxl::read_excel("doc/181002_Single_cell_sample list.xlsx")
#sample_n = which(df_samples$tests %in% c("test1", "test2", "test3","test4"))
sample_n = which(df_samples$tests %in% paste0("test",1:4))
table(df_samples$tests)
df_samples[sample_n,] %>% kable() %>% kable_styling()
samples <- df_samples$samples[sample_n]
projects <- df_samples$projects[sample_n]
conditions <- df_samples$conditions[sample_n]
tests <- df_samples$tests[sample_n]
#======1.2 load  SingleCellExperiment =========================
(load(file = "./data/sce_list_20181121.Rda"))
names(sce_list)
MCL_Seurat <- lapply(sce_list, as.seurat) %>%
        lapply(NormalizeData) %>%
        #lapply(ScaleData) %>%
        lapply(FindVariableGenes, do.plot = FALSE)

for(i in 1:length(samples)){
        MCL_Seurat[[i]]@meta.data$conditions <- conditions[i]
        MCL_Seurat[[i]]@meta.data$tests <- tests[i]
}
# we will take the union of the top 1k variable genes in each dataset for alignment
genes.use <- MCL_Seurat %>% lapply(function(object) head(rownames(object@hvg.info), 1000)) %>%
                unlist %>% unique
length(genes.use)

#========1.3 merge ===================================
MCL <- Reduce(function(x, y) MergeSeurat(x, y, do.normalize = F), MCL_Seurat)
MCL@var.genes = genes.use
remove(sce_list,MCL_Seurat);GC()

#MCL@meta.data$orig.ident <- gsub("Pt-MD", "MD", MCL@meta.data$orig.ident)
#MCL = SetAllIdent(MCL, id = "orig.ident")
#======1.4 mito, QC, filteration =========================
mito.genes <- grep(pattern = "^MT-", x = rownames(x = MCL@data), value = TRUE)
percent.mito <- Matrix::colSums(MCL@raw.data[mito.genes, ])/Matrix::colSums(MCL@raw.data)
MCL <- AddMetaData(object = MCL, metadata = percent.mito, col.name = "percent.mito")

(load(file = "./data/MCL_g1.Rda"))

MCL <- FilterCells(object = MCL, subset.names = c("nGene","nUMI","percent.mito"),
                   low.thresholds = c(500,1000, -Inf), 
                   high.thresholds = c(Inf,Inf, 0.5))
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
MCL@meta.data$S.Score = MCL@meta.data$S.Score - min(MCL@meta.data$S.Score)
MCL@meta.data$G2M.Score = MCL@meta.data$G2M.Score - min(MCL@meta.data$G2M.Score)
head(x = MCL@meta.data)

#======1.5 1st run of pca-tsne  =========================
MCL %<>% NormalizeData %>% ScaleData %>%
         RunPCA(pc.genes = MCL@var.genes, pcs.compute = 50, do.print = F)
#MCL %<>% RunICA(ic.genes = MCL@var.genes, ics.compute = 50, print.results = F)
DimElbowPlot(MCL, reduction.type = "pca", dims.plot = 50)
system.time(MCL %<>% RunHarmony("orig.ident", 
                                theta = 2, plot_convergence = TRUE,
                                nclust = 50, max.iter.cluster = 100))

MCL@ident %<>% factor(levels = samples)
p1 <- DimPlot(object = MCL, reduction.use = "harmony", pt.size = .1, group.by = "orig.ident", do.return = T)
p2 <- VlnPlot(object = MCL, features.plot = "Harmony1", group.by = "orig.ident", do.return = TRUE,
              x.lab.rot = T)
jpeg(paste0(path,"/S1_Harmony_vplot.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p1,p2)
dev.off()

jpeg(paste0(path,"/S1_Harmony_DimHeatmap.jpeg"), units="in", width=10, height=7,res=600)
DimHeatmap(object = MCL, reduction.type = "harmony", cells.use = 500, dim.use = 1:6, do.balanced = TRUE)
dev.off()

#========1.6 Seurat tSNE Functions for Integrated Analysis Using Harmony Results=======
system.time({
        MCL %<>% RunTSNE(reduction.use = "harmony", dims.use = 1:20, do.fast = T)
        MCL %<>% FindClusters(reduction.type = "harmony", resolution = 0.6, dims.use = 1:20, 
                              force.recalc = TRUE, print.output = FALSE)
})

p3 <- TSNEPlot(MCL, do.return = T, pt.size = 0.5, group.by = "orig.ident")
p4 <- TSNEPlot(MCL, do.label = T, do.return = T, pt.size = 0.5)
jpeg(paste0(path,"/S1_Harmony_TSNEPlot.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p3, p4)
dev.off()

g_Harmony <- TSNEPlot.1(object = MCL, do.label = F, group.by = "ident", 
                 do.return = TRUE, no.legend = T, 
                 colors.use = ExtractMetaColor(MCL),
                 pt.size = 1,label.size = 6 )+
        ggtitle("Tsne plot of all cell types")+
        theme(text = element_text(size=15),							
              plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 

jpeg(paste0(path,"/TSNEplot-Harmony.jpeg"), units="in", width=10, height=7,res=600)
print(g_Harmony)
dev.off()

save(MCL, file = "./data/MCL_Harmony_20181121.Rda")

jpeg(paste0(path,"/TSNEplot-alignment~.jpeg"), units="in", width=10, height=7,res=600)
do.call(plot_grid, list(g_CCA,g_MNN,g_Harmony))
dev.off()
