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
df_samples <- readxl::read_excel("doc/181128_scRNAseq_info.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
sample_n = which(df_samples$tests %in% c("test",paste0("test",1:6)))
df_samples[sample_n,] %>% kable() %>% kable_styling()
table(df_samples$tests);nrow(df_samples)
samples <- df_samples$sample[sample_n]
projects <- df_samples$project[sample_n]
tests <- df_samples$tests[sample_n]

current <- list.files("data")[!grepl(".Rda|RData",list.files("data"))]
samples[!(samples %in% current)]

MCL_raw <- list()
MCL_Seurat <- list()
for(i in 1:length(samples)){
    MCL_raw[[i]] <- Read10X(data.dir = paste0("./data/",
                                              samples[i],"/outs/filtered_gene_bc_matrices/hg19/"))
    colnames(MCL_raw[[i]]) <- paste0(samples[i],
                                     "_",colnames(MCL_raw[[i]]))
    MCL_Seurat[[i]] <- CreateSeuratObject(MCL_raw[[i]],
                                          min.cells = 3,
                                          min.genes = 0,
                                          project = projects[i],
                                          names.delim = "_")
}
#======1.1.2 QC before merge =========================
cell.number <- sapply(MCL_Seurat, function(x) length(x@cell.names))
QC_list <- lapply(MCL_Seurat, function(x) as.matrix(x = x@raw.data))
median.nUMI <- sapply(QC_list, function(x) median(colSums(x)))
median.nGene <- sapply(QC_list, function(x) median(apply(x,2,function(y) sum(length(y[y>0])))))

min.nUMI <- sapply(QC_list, function(x) min(colSums(x)))
min.nGene <- sapply(QC_list, function(x) min(apply(x,2,function(y) sum(length(y[y>0])))))

QC.list <- cbind(df_samples,cell.number, median.nUMI,median.nGene,min.nUMI,min.nGene,
                 row.names = samples)
write.csv(QC.list,paste0(path,"QC_list.csv"))
QC.list %>% kable() %>% kable_styling()
remove(QC_list,median.nUMI,median.nGene,min.nUMI,min.nGene,QC.list);GC()

#========1.1.3 merge ===================================
MCL <- Reduce(function(x, y) MergeSeurat(x, y, do.normalize = F), MCL_Seurat)
remove(MCL_raw,MCL_Seurat);GC()

mito.genes <- grep(pattern = "^MT-", x = rownames(x = MCL@data), value = TRUE)
percent.mito <- Matrix::colSums(MCL@raw.data[mito.genes, ])/Matrix::colSums(MCL@raw.data)
MCL <- AddMetaData(object = MCL, metadata = percent.mito, col.name = "percent.mito")

MCL@ident = factor(MCL@ident,levels = samples)

g1 <- VlnPlot(object = MCL, features.plot = c("nGene", "nUMI", "percent.mito"),
              nCol = 1,point.size.use = 0.2,,size.x.use = 10, group.by = "ident",
              x.lab.rot = T, do.return = T,return.plotlist =T)

save(g1,file="./data/g1_23_20181205.Rda")

#======1.2 load  SingleCellExperiment =========================
(load(file = "./data/sce_23_20181205.Rda"))
names(sce_list)
MCL_Seurat <- lapply(sce_list, as.seurat) %>%
        lapply(NormalizeData) %>%
        #lapply(ScaleData) %>%
        lapply(FindVariableGenes, do.plot = FALSE)

for(i in 1:length(samples)){
        MCL_Seurat[[i]]@meta.data$tests <- tests[i]
}
# we will take the union of the top 1k variable genes in each dataset for alignment
genes.use <- MCL_Seurat %>% lapply(function(object) head(rownames(object@hvg.info), 500)) %>%
                unlist %>% unique
length(genes.use)

#========1.3 merge ===================================
MCL <- Reduce(function(x, y) MergeSeurat(x, y, do.normalize = F), MCL_Seurat)
#MCL@var.genes = genes.use
remove(sce_list,MCL_Seurat);GC()

#MCL@meta.data$orig.ident <- gsub("Pt-MD", "MD", MCL@meta.data$orig.ident)
#MCL = SetAllIdent(MCL, id = "orig.ident")
#======1.4 mito, QC, filteration =========================
mito.genes <- grep(pattern = "^MT-", x = rownames(x = MCL@data), value = TRUE)
percent.mito <- Matrix::colSums(MCL@raw.data[mito.genes, ])/Matrix::colSums(MCL@raw.data)
MCL <- AddMetaData(object = MCL, metadata = percent.mito, col.name = "percent.mito")

(load(file = "./data/g1_23_20181205.Rda"))

MCL <- FilterCells(object = MCL, subset.names = c("nGene","nUMI","percent.mito"),
                   low.thresholds = c(50,100, -Inf), 
                   high.thresholds = c(Inf,Inf, 0.5))

MCL@ident = factor(MCL@ident,levels = samples)
g2 <- VlnPlot(object = MCL, features.plot = c("nGene", "nUMI", "percent.mito"), 
              nCol = 1,point.size.use = 0.2,size.x.use = 10, group.by = "ident",
              x.lab.rot = T, do.return = T,return.plotlist =T)
save(g2,file = "./data/g2_23_20181205.Rda")
jpeg(paste0(path,"/S1_nGene.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[1]]+ggtitle("nGene in raw data")+ 
                    scale_y_log10(limits = c(10,10000)),#+ylim(c(0,1000)),
                g2[[1]]+ggtitle("nGene after filteration")+ 
                    scale_y_log10(limits = c(10,10000))))
dev.off()
jpeg(paste0(path,"/S1_nUMI.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[2]]+ggtitle("nUMI in raw data")+ 
                    scale_y_log10(limits = c(50,100000)),#+ylim(c(0,1000)),
                g2[[2]]+ggtitle("nUMI after filteration")+ 
                    scale_y_log10(limits = c(50,100000))))
dev.off()
jpeg(paste0(path,"/S1_mito.jpeg"), units="in", width=10, height=7,res=600)
print(plot_grid(g1[[3]]+ggtitle("mito % in raw data")+ 
                    ylim(c(0,0.5)),
                g2[[3]]+ggtitle("mito % after filteration")+ 
                    ylim(c(0,0.5))))
dev.off()

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

#======1.5 loom pca=======================
MCL <- NormalizeData(object = MCL)
jpeg(paste0(path,"/S1_dispersion.jpeg"), units="in", width=10, height=7,res=600)
MCL <- FindVariableGenes(object = MCL, mean.function = ExpMean, 
                         dispersion.function = LogVMR, do.plot = T, 
                         x.low.cutoff = 0.0125, x.high.cutoff = 8, y.cutoff = 1)
dev.off()
length(MCL@var.genes)

# Convert from Seurat to loom Convert takes and object in 'from', a name of
# a class in 'to', and, for conversions to loom, a filename

MCL %<>% ScaleData %>%
         RunPCA(pc.genes = MCL@var.genes, pcs.compute = 50, do.print = F)

jpeg(paste0(path,"/S1_DimElbowPlot_pca.jpeg"), units="in", width=10, height=7,res=600)
DimElbowPlot(MCL, reduction.type = "pca", dims.plot = 50)
dev.off()

jpeg(paste0(path,"/S1_PCHeatmap.jpeg"), units="in", width=10, height=7,res=600)
PCHeatmap(MCL, pc.use = c(1:3, 45:50), cells.use = 500, do.balanced = TRUE)
dev.off()

saveRDS(MCL@scale.data, file = "./data/MCL.scale.data_23_20181205.Rda")

#MCL@scale.data = readRDS("./data/MCL.scale.data_23_20181205.Rda")
#MCL <- RunPCA(MCL, pc.genes = MCL@var.genes, pcs.compute = 50, do.print = F)
GC()

#======1.6 RunHarmony=======================
jpeg(paste0(path,"/S1_RunHarmony.jpeg"), units="in", width=10, height=7,res=600)
system.time(MCL %<>% RunHarmony("orig.ident", dims.use = 1:50,
                                theta = 2, plot_convergence = TRUE,
                                nclust = 50, max.iter.cluster = 100))
dev.off()

MCL@ident %<>% factor(levels = samples)
p1 <- DimPlot(object = MCL, reduction.use = "harmony", pt.size = .1, group.by = "orig.ident", do.return = T)
p2 <- VlnPlot(object = MCL, features.plot = "Harmony1", group.by = "orig.ident", do.return = TRUE,
              x.lab.rot = T)
jpeg(paste0(path,"/S1_Harmony_vplot.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p1,p2)
dev.off()

jpeg(paste0(path,"/S1_Harmony_DimHeatmap.jpeg"), units="in", width=10, height=7,res=600)
DimHeatmap(object = MCL, reduction.type = "harmony", cells.use = 500, 
           dim.use = c(1:3,48:50), do.balanced = TRUE)
dev.off()

#========1.6 Seurat tSNE Functions for Integrated Analysis Using Harmony Results=======
system.time({
        MCL %<>% RunTSNE(reduction.use = "harmony", dims.use = 1:50, do.fast = TRUE)
        MCL %<>% FindClusters(reduction.type = "harmony", resolution = 0.6, dims.use = 1:50,
                              save.SNN = TRUE, n.start = 10, nn.eps = 0.5,
                              force.recalc = TRUE, print.output = FALSE)
})

p3 <- TSNEPlot(MCL, do.return = T, pt.size = 0.5, group.by = "orig.ident")
 np4 <- TSNEPlot(MCL, do.label = T, do.return = T, pt.size = 0.5)
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

save(MCL, file = "./data/MCL_Harmony_23_20181205.Rda")
