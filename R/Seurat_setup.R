########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(SingleR)
library(scran)
source("../R/Seurat_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
dir.create(path, recursive = T)
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
sample_n = which(df_samples$Tests %in% c("test1", "test2", "test3","test4"))
table(df_samples$Tests)
df_samples[sample_n,]
samples <- df_samples$Samples[sample_n]
projects <- df_samples$Projects[sample_n]
conditions <- df_samples$Conditions[sample_n]

MCL_raw <- list()
MCL_Seurat <- list()
for(i in 1:length(samples)){
    MCL_raw[[i]] <- Read10X(data.dir = paste0("./data/",
                                             samples[i],"/outs/filtered_gene_bc_matrices/hg19/"))
    colnames(MCL_raw[[i]]) <- paste0(samples[i],
                                            "_",colnames(MCL_raw[[i]]))
    MCL_Seurat[[i]] <- CreateSeuratObject(MCL_raw[[i]],
                                                 min.cells = 3,
                                                 min.genes = 200,
                                                 project = projects[i],
                                                 names.delim = "_")
    MCL_Seurat[[i]]@meta.data$conditions <- conditions[i]
}
MCL <- Reduce(function(x, y) MergeSeurat(x, y, do.normalize = F), MCL_Seurat)
remove(MCL_raw,MCL_Seurat);GC()
MCL <- FilterCells(MCL, subset.names = "nGene",
                    low.thresholds = 200,
                    high.thresholds = Inf) %>%
    NormalizeData() %>%
    ScaleData(display.progress = FALSE) %>%
    FindVariableGenes(do.plot = FALSE, display.progress = FALSE)
save(MCL, file = "./data/MCL_20181019.Rda")

#======1.2 QC, pre-processing and normalizing the data=========================
# 1.2.1 Calculate median UMI per cell
Iname = load(file = "./data/MCL_20181019.Rda")
MCL_raw_data <- as.matrix(x = MCL@raw.data)
mean(colSums(MCL_raw_data))
median(colSums(MCL_raw_data))
min(colSums(MCL_raw_data))
remove(MCL_raw_data);GC()

# 1.2.3 calculate mitochondria percentage
mito.genes <- grep(pattern = "^MT-", x = rownames(x = MCL@data), value = TRUE)
percent.mito <- Matrix::colSums(MCL@raw.data[mito.genes, ])/Matrix::colSums(MCL@raw.data)
MCL <- AddMetaData(object = MCL, metadata = percent.mito, col.name = "percent.mito")

MCL@ident = factor(MCL@ident,levels = samples)

g1 <- VlnPlot(object = MCL, features.plot = c("nGene", "nUMI", "percent.mito"), 
              nCol = 1,point.size.use = 0.2,
              x.lab.rot = T, do.return = T,return.plotlist =T)

jpeg(paste0(path,"/S1_nUMI1.jpeg"), units="in", width=10, height=7,res=600)
print(g1[[2]]+ scale_y_log10() )
dev.off()


MCL <- FilterCells(object = MCL, subset.names = c("nGene","nUMI","percent.mito"),
                    low.thresholds = c(500,1000, -Inf), 
                    high.thresholds = c(Inf,Inf, 0.5))

g2 <- VlnPlot(object = MCL, features.plot = c("nGene", "nUMI", "percent.mito"), 
              nCol = 1,point.size.use = 0.2,
              x.lab.rot = T, do.return = T,return.plotlist =T)

jpeg(paste0(path,"/S1_nGene2.jpeg"), units="in", width=10, height=7,res=600)
print(g2[[1]]+ scale_y_log10() )
dev.off()
######################################

# After removing unwanted cells from the dataset, the next step is to normalize the data.
MCL <- NormalizeData(object = MCL, normalization.method = "LogNormalize", 
                      scale.factor = 10000)
MCL <- FindVariableGenes(object = MCL, mean.function = ExpMean, 
                          dispersion.function = LogVMR, do.plot = FALSE, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(MCL@var.genes)
#======1.3 1st run of pca-tsne  =========================
MCL <- ScaleData(object = MCL) %>%
    RunPCA() %>%
    FindClusters(dims.use = 1:20, force.recalc = T, print.output = FALSE) %>%
    RunTSNE()
#MCL@meta.data$orig.ident <- gsub("PND18pre","PND18",MCL@meta.data$orig.ident)
g1 <- TSNEPlot(object = MCL, do.label = F, group.by = "orig.ident", 
         do.return = TRUE, no.legend = F, #colors.use = singler.colors,
         pt.size = 1,label.size = 8 )+
    ggtitle("Oiginal")+
    theme(text = element_text(size=15),							
          plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 

save(MCL, file = "./data/MCL_20181019.Rda")
Iname = load("./data/MCL_20181019.Rda")

#======1.4 Add Cell-cycle score =========================
# Read in a list of cell cycle markers, from Tirosh et al, 2015
cc.genes <- readLines(con = "./data/seurat_resources/regev_lab_cell_cycle_genes.txt")
# We can segregate this list into markers of G2/M phase and markers of S phase
s.genes <- HumanGenes(MCL,cc.genes[1:43])
g2m.genes <- HumanGenes(MCL,cc.genes[44:97])
# Assign Cell-Cycle Scores
MCL <- CellCycleScoring(object = MCL, s.genes = s.genes, g2m.genes = g2m.genes, 
                         set.ident = TRUE)
# Visualize the distribution of cell cycle markers across
RidgePlot(object = MCL, features.plot = HumanGenes(MCL,c("CCND1","CDK4","CCND2","CDK6","CCND3","RB1")), 
          nCol = 2)
# regressing out the difference between the G2M and S phase scores
MCL@meta.data$CC.Difference <- MCL@meta.data$S.Score - MCL@meta.data$G2M.Score
# view cell cycle scores and phase assignments
head(x = MCL@meta.data)

#======1.5 Add project id =========================
batchname = MCL@meta.data$orig.ident
batch.effect = rep(NA,length(batchname))
batch.effect[batchname %in% c("Pt-DJ","Pt-1294")] = 1
batch.effect[!(batchname %in% c("Pt-DJ","Pt-1294"))] = 2
names(batch.effect) = rownames(MCL@meta.data)
MCL <- AddMetaData(object = MCL, metadata = batch.effect, col.name = "project.id")
table(MCL@meta.data$project.id)
#------
batchname = MCL@meta.data$orig.ident
batch.effect = as.numeric(factor(batchname,levels = samples))
names(batch.effect) = rownames(MCL@meta.data)
MCL <- AddMetaData(object = MCL, metadata = batch.effect, col.name = "batch.effect")
table(MCL@meta.data$batch.effect)
head(x = MCL@meta.data)

#======1.6 vars.to.regress ScaleData =========================
SingleFeaturePlot.1(MCL,"nUMI",threshold=10000)
SingleFeaturePlot.1(MCL,"nGene",threshold=2000)
SingleFeaturePlot.1(MCL,"batch.effect",threshold=3.0)
SingleFeaturePlot.1(MCL,"percent.mito",threshold=0.05)
SingleFeaturePlot.1(MCL,"CC.Difference",threshold=0.05)
TSNEPlot(object = MCL, do.label = F, group.by = "batch.effect")
MCL@scale.data = NULL
MCL <- ScaleData(object = MCL, #genes.use = MCL@var.genes,
                  model.use = "linear", do.par=T, do.center = T, do.scale = T,
                  vars.to.regress = c("nUMI","percent.mito","batch.effect"),#"CC.Difference","percent.mito"--nogood,"nUMI"--nogood
                  display.progress = T)
#save(MCL, file = "./data/MCL_20181001.Rda") #do.center = F, do.scale = T
#======1.7 Performing MNN-based correction =========================
#https://bioconductor.org/packages/3.8/workflows/vignettes/simpleSingleCell/inst/doc/work-5-mnn.html#4_performing_mnn-based_correction
set.seed(100)
original <- lapply(samples, function(x) MCL@scale.data[MCL@var.genes, 
                                                 (MCL@meta.data$orig.ident %in% x)])
mnn.out <- do.call(fastMNN, c(original, list(k=20, d=50, auto.order=T,
                                             approximate=TRUE)))
dim(mnn.out$corrected)
rownames(mnn.out$corrected) = MCL@cell.names
colnames(mnn.out$corrected) = paste0("MNN_",1:ncol(mnn.out$corrected))
#Storing a new MNN
MCL <- SetDimReduction(object = MCL, reduction.type = "MNN", slot = "cell.embeddings",
                       new.data = mnn.out$corrected)
MCL <- SetDimReduction(object = MCL, reduction.type = "MNN", slot = "key", 
                       new.data = "MNN_")
remove(original);GC()
MCL <- SetAllIdent(MCL,id = "orig.ident")
DimPlot(object = MCL, reduction.use = "MNN", pt.size = 0.5)

#======1.7 unsupervised clustering based on MNN =========================
MCL <- RunPCA(object = MCL, pc.genes = MCL@var.genes, pcs.compute = 100, 
               do.print = TRUE, pcs.print = 1:5, genes.print = 5)
PCAPlot(object = MCL)
PCElbowPlot(object = MCL, num.pc = 100)
PCHeatmap(MCL, pc.use = c(1:3, 25:30), cells.use = 500, do.balanced = TRUE)

DimElbowPlot.1(object = MCL, reduction.type = "MNN", 
             dims.plot = 50,slot = "cell.embeddings")

MCL <- RunTSNE(object = MCL, reduction.use = "MNN", dims.use = 1:50, 
                do.fast = TRUE, perplexity= 30)

MCL <- FindClusters(object = MCL, reduction.type = "MNN", 
                    dims.use = 1:50, resolution = 0.6, 
                     k.param = 30,force.recalc = T,
                     save.SNN = TRUE, n.start = 100, nn.eps = 0, print.output = FALSE)

g2 <- TSNEPlot.1(object = MCL, do.label = F, group.by = "orig.ident", 
           do.return = TRUE, no.legend = T, 
           #colors.use = singler.colors,
           pt.size = 1,label.size = 4 )+
    ggtitle("Corrected")+
    theme(text = element_text(size=15),							
          plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 

jpeg(paste0(path,"TSNEPlot.jpeg"), units="in", width=10, height=7,res=600)
g2
dev.off()

jpeg(paste0(path,"remove_batch~.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(g1 +ggtitle("Original")+
              theme(text = element_text(size=15),
                    legend.position="none",
                    plot.title = element_text(hjust = 0.5,size = 18, face = "bold")),g2)
dev.off()
MCL <- StashIdent(object = MCL, save.name = "MNN_scale")
MCL <- SetAllIdent(object = MCL, id = "MNN_ident")
save(MCL, file = "./data/MCL_20181019.Rda")