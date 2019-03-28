########################################################################
#
#  0 setup environment, install libraries if nLynchessary, load libraries
# 
# ######################################################################

library(Seurat)
library(magrittr)
library(harmony)
library(dplyr)
library(kableExtra)
source("../R/Seurat_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
if(!dir.exists("data/")) dir.create("data")
########################################################################
#
#  1 Data preprocessing
# 
# ######################################################################
#======1.1 Load the data files and Set up Seurat object =========================
# read sample summary list
df_samples <- readxl::read_excel("doc/190320_scRNAseq_info.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
sample_n = which(df_samples$tests %in% c("control",paste0("test",2:10)))
df_samples = df_samples[sample_n,]
attach(df_samples)
df_samples %>% kable() %>% kable_styling()
samples = sample

#======1.2 load  SingleCellExperiment =========================
(load(file = "data/sce_30_20190320.Rda"))
names(sce_list)
object_list <- lapply(sce_list, as.seurat) %>%
        lapply(NormalizeData) %>%
        #lapply(ScaleData) %>%
        lapply(FindVariableGenes, do.plot = FALSE)

for(i in 1:length(samples)){
        object_list[[i]]@meta.data$tests <- tests[i]
        object_list[[i]]@meta.data$conditions <- conditions[i]
        object_list[[i]]@meta.data$projects <- project[i]
        object_list[[i]]@meta.data$groups <- group[i]
        object_list[[i]]@meta.data$tissues <- tissue[i]
        object_list[[i]]@meta.data$tsne <- tsne[i]
        
}
# we will take the union of the top 1k variable genes in each dataset for alignment
intersect.gene = 300
genes.use = c()
while(length(genes.use) < 2500){
        genes.use <- object_list %>% 
                lapply(function(object) head(rownames(object@hvg.info), 
                                             intersect.gene)) %>%
                unlist %>% unique
        intersect.gene =intersect.gene+100
}
print(paste0("genes.use: ",length(genes.use)))

#========1.3 merge ===================================
object <- Reduce(function(x, y) MergeSeurat(x, y, do.normalize = F), object_list)
object@var.genes = genes.use
remove(sce_list,object_list);GC()

object = SetAllIdent(object, id = "orig.ident")
#======1.4 mito, QC, filteration =========================
object@meta.data$percent.mito = object@meta.data$pct_counts_Mito/100

(remove <- which(colnames(object@meta.data) %in%c("is_cell_control",
                                           "pct_counts_in_top_500_features_Mito")))
meta.data = object@meta.data[,-seq(remove[1], remove[2], by=1)]
object@meta.data = meta.data 

remove(meta.data)
(load(file = paste0(path,"g1_30_20190320.Rda")))

object <- FilterCells(object = object, subset.names = c("nGene","nUMI","percent.mito"),
                   low.thresholds = c(500,800, -Inf), 
                   high.thresholds = c(Inf,Inf, 0.5))

object@ident = factor(object@ident,levels = samples)
g2 <- lapply(c("nGene", "nUMI", "percent.mito"), function(features){
        VlnPlot(object = object, features.plot = features, nCol = 3, 
                point.size.use = 0.2,size.x.use = 10, group.by = "ident",
                x.lab.rot = T, do.return = T)
})
save(g2,file= paste0(path,"g2_30_20190320.Rda"))

jpeg(paste0(path,"/S1_nGene.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[1]]+ggtitle("nGene in raw data")+ 
                        scale_y_log10(limits = c(100,10000)),
                g2[[1]]+ggtitle("nGene after filteration")+ 
                        scale_y_log10(limits = c(100,10000))))
dev.off()
jpeg(paste0(path,"/S1_nUMI.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[2]]+ggtitle("nUMI in raw data")+ 
                        scale_y_log10(limits = c(500,100000)),
                g2[[2]]+ggtitle("nUMI after filteration")+ 
                        scale_y_log10(limits = c(500,100000))))
dev.off()
jpeg(paste0(path,"/S1_mito.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[3]]+ggtitle("mito % in raw data")+ 
                        ylim(c(0,1)),
                g2[[3]]+ggtitle("mito % after filteration")+ 
                        ylim(c(0,1))))
dev.off()
#======1.5 FindVariableGenes=======================
object <- NormalizeData(object = object)
jpeg(paste0(path,"/S1_dispersion.jpeg"), units="in", width=10, height=7,res=600)
object <- FindVariableGenes(object = object, mean.function = ExpMean, 
                            dispersion.function = LogVMR, do.plot = F, 
                            x.low.cutoff = 0.1, x.high.cutoff = 8, y.cutoff = 0.5)
dev.off()
length(object@var.genes)
table(object@var.genes %in% genes.use)

#======1.5 Add Cell-cycle score =========================
# Read in a list of cell cycle markers, from Tirosh et al, 2015
cc.genes <- readLines(con = "../R/seurat_resources/regev_lab_cell_cycle_genes.txt")
s.genes <- FilterGenes(object,cc.genes[1:43])
g2m.genes <- FilterGenes(object,cc.genes[44:97])
object <- CellCycleScoring(object = object, s.genes = s.genes, g2m.genes = g2m.genes, 
                        set.ident = FALSE)
RidgePlot(object = object, features.plot = FilterGenes(object,c("CCND1","CDK4","CCND2","CDK6","CCND3","RB1")), 
          nCol = 2)
object@meta.data$CC.Difference <- object@meta.data$S.Score - object@meta.data$G2M.Score
object@meta.data$S.Score = object@meta.data$S.Score - min(object@meta.data$S.Score)
object@meta.data$G2M.Score = object@meta.data$G2M.Score - min(object@meta.data$G2M.Score)
tail(x = object@meta.data)


#======1.6 PCA =========================
object %<>% ScaleData
object@scale.data = readRDS("data/MCL.scale.data_Harmony_30_20190320.rds")
object %<>% RunPCA(pc.genes = object@var.genes, pcs.compute = 75, do.print = F)

jpeg(paste0(path,"/S1_PCElbowPlot.jpeg"), units="in", width=10, height=7,res=600)
PCElbowPlot(object, num.pc = 100)
dev.off()

jpeg(paste0(path,"/S1_PCHeatmap.jpeg"), units="in", width=10, height=7,res=600)
PCHeatmap(object, pc.use = c(1:3, 48:50,73:75), cells.use = 500, do.balanced = TRUE)
dev.off()

GC()
pcs =1:75
system.time({
        object %<>% RunTSNE(reduction.use = "pca", dims.use = pcs, do.fast = TRUE) %>%
                FindClusters(reduction.type = "pca", resolution = 0.6, dims.use = pcs,
                             save.SNN = TRUE, n.start = 10, nn.eps = 0.5,
                             force.recalc = TRUE, print.output = FALSE)
})
p0 <- DimPlot(object = object, reduction.use = "tsne", pt.size = 0.3, group.by = "orig.ident", do.return = T)
#======1.6 RunHarmony=======================
jpeg(paste0(path,"/S1_RunHarmony.jpeg"), units="in", width=10, height=7,res=600)
system.time(object %<>% RunHarmony("orig.ident", dims.use = pcs,
                                theta = 2, plot_convergence = TRUE,
                                nclust = 50, max.iter.cluster = 100))
dev.off()

object@ident %<>% factor(levels = samples)
p1 <- DimPlot(object = object, reduction.use = "harmony", pt.size = 0.3, group.by = "orig.ident", do.return = T)
p2 <- VlnPlot(object = object, features.plot = "Harmony1", group.by = "orig.ident", do.return = TRUE,
              x.lab.rot = T)
jpeg(paste0(path,"/S1_Harmony_vplot.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p1,p2)
dev.off()

jpeg(paste0(path,"/S1_Harmony_DimHeatmap.jpeg"), units="in", width=10, height=7,res=600)
DimHeatmap(object = object, reduction.type = "harmony", cells.use = 500, 
           dim.use =c(1:3, 48:50,73:75), do.balanced = TRUE)
dev.off()

#========1.6 Seurat tSNE Functions for Integrated Analysis Using Harmony Results=======
system.time(
        object %<>% RunTSNE(reduction.use = "harmony", dims.use = pcs, do.fast = TRUE))
system.time(
        object %<>% FindClusters(reduction.type = "harmony", resolution = 0.6, dims.use = pcs,
                              save.SNN = TRUE, n.start = 10, nn.eps = 0.5,
                              force.recalc = TRUE, print.output = FALSE))


p3 <- TSNEPlot(object, do.return = T, pt.size = 0.3, group.by = "orig.ident")
p4 <- TSNEPlot(object, do.label = T, do.return = T, pt.size = 0.3)

jpeg(paste0(path,"/S1_pca_vs_Harmony_TSNEPlot.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p0+ggtitle("Raw data")+
                  theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold")),
          p3+ggtitle("After alignment")+
                  theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) )
dev.off()

jpeg(paste0(path,"/S1_Harmony_TSNEPlot.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p3+ggtitle("group by samples")+
                  theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold")),
          p4+ggtitle("group by clusters")+
                  theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold")))
dev.off()

g_Harmony <- TSNEPlot.1(object = object, do.label = T, group.by = "ident",
                        do.return = TRUE, no.legend = F, 
                        #colors.use = ExtractMetaColor(object),
                        pt.size = 1,label.size = 6 )+
        ggtitle("Tsne plot of all clusters")+
        theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 

jpeg(paste0(path,"/TSNEplot-Harmony.jpeg"), units="in", width=10, height=7,res=600)
print(g_Harmony)
dev.off()

saveRDS(object@scale.data, file = "data/MCL.scale.data_Harmony_30_20190320.rds")
object@scale.data = NULL; GC()
save(object, file = "data/MCL_Harmony_30_20190320.Rda")
