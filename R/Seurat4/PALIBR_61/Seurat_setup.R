########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
#
# ######################################################################
# conda activate r4.1.1
#devtools::install_github("immunogenomics/harmony", ref= "ee0877a",force = T)
invisible(lapply(c("Seurat","dplyr","ggplot2","cowplot","pbapply","sctransform","harmony"), function(x) {
    suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)


########################################################################
#
#  1 Seurat Alignment
#
# ######################################################################
#======1.1 Setup the Seurat objects =========================
# read sample summary list
df_samples <- readxl::read_excel("output/20220318/20220318_scRNAseq_info.xlsx")
df_samples = as.data.frame(df_samples)
colnames(df_samples) %<>% tolower()
df_samples %<>% filter(sequence %in% "GEX") %>% filter(phase %in% "PALIBR_I")
nrow(df_samples)
df_samples$date %<>% gsub(" UTC","",.) %>% as.character()
#======1.2 load  Seurat =========================
object = readRDS(file = "data/MCL_51_20210724.rds")
DefaultAssay(object) = "RNA"
object[["SCT"]] = NULL
object[["integrated"]] = NULL
GC()
object1 = readRDS(file = "data/MCL_10_20220318.rds")
object %<>% merge(object1)
remove(object1);GC()

meta.data = object@meta.data
for(i in 1:length(df_samples$sample)){
    cells <- meta.data$orig.ident %in% df_samples$sample[i]
    print(table(cells))
    meta.data[cells,"tissue"] = df_samples$tissue[i]
    meta.data[cells,"project"] = df_samples$project[i]
    meta.data[cells,"date"] = as.character(df_samples$date[i])
    meta.data[cells,"Mean.Reads.per.Cell"] = df_samples$mean.reads.per.cell[i]
    meta.data[cells,"Number.of.Reads"] = df_samples$number.of.reads[i]
    meta.data[cells,"Sequencing.Saturation"] = df_samples$sequencing.saturation[i]
}
meta.data$date %<>% as.character() %>% gsub(" .*","",.)
meta.data$orig.ident %<>% factor(levels = df_samples$sample)
table(rownames(object@meta.data) == rownames(meta.data))
table(colnames(object) == rownames(meta.data))
object@meta.data = meta.data
Idents(object) = "orig.ident"
#======1.6 Performing SCTransform and integration =========================

set.seed(100)
object_list <- SplitObject(object, split.by = "orig.ident")
remove(object);GC()

object_list %<>% pblapply(SCTransform,method = "glmGamPoi",vars.to.regress = "percent.mt")
features <- SelectIntegrationFeatures(object.list = object_list)

options(future.globals.maxSize= object.size(object_list)*1.5)
object_list %<>% pblapply(FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
    x
})
anchors <- FindIntegrationAnchors(object.list = object_list,reference = c(1, 40), reduction = "rpca",
                                  dims = 1:50)
#remove(object_list)
GC()
# this command creates an 'integrated' data assay
object <- IntegrateData(anchorset = anchors,normalization.method = "SCT", dims = 1:50)
remove(anchors);GC()
format(object.size(object),unit = "GB")
saveRDS(object, file = "data/MCL_61_20220318.rds")
# Perform an integrated analysis
# Now we can run a single integrated analysis on all cells!

# specify that we will perform downstream analysis on the corrected data note that the original
# unmodified data still resides in the 'RNA' assay
DefaultAssay(object) <- "integrated"

# Run the standard cca workflow for umap & tsne visualization
object %<>% ScaleData(verbose = FALSE)
object %<>% RunPCA(npcs = 100, verbose = FALSE)
jpeg(paste0(path,"ElbowPlot.jpeg"), units="in", width=10, height=7,res=600)
print(ElbowPlot(object,ndims = 100))
dev.off()
npcs = 100

object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)
system.time(object %<>% RunTSNE(reduction = "pca", dims = 1:npcs))

object[["cca.umap"]] <- CreateDimReducObject(embeddings = object@reductions[["umap"]]@cell.embeddings,
                                             key = "ccaUMAP_", assay = DefaultAssay(object))
colnames(object[["cca.umap"]]@cell.embeddings) %<>% paste0("cca-",.)

object[["cca.tsne"]] <- CreateDimReducObject(embeddings = object@reductions[["tsne"]]@cell.embeddings,
                                             key = "cca-SNE_", assay = DefaultAssay(object))
colnames(object[["cca.tsne"]]@cell.embeddings) %<>% paste0("cca-",.)

saveRDS(object, file = "data/MCL_61_20220318.rds")
object@reductions$pca = NULL
object@reductions$umap = NULL
object@reductions$tsne = NULL
#======1.7 UMAP from raw pca =========================
DefaultAssay(object)  = "SCT"
object %<>% SCTransform(method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = TRUE)

object <- FindVariableFeatures(object = object, selection.method = "vst",
                               num.bin = 20, nfeatures = 2000,
                               mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))
object %<>% ScaleData(verbose = FALSE)
object %<>% RunPCA(verbose = T,npcs = 100)

jpeg(paste0(path,"S1_ElbowPlot_SCT.jpeg"), units="in", width=10, height=7,res=600)
ElbowPlot(object, ndims = 100)
dev.off()

npcs = 100
object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)
system.time(object %<>% RunTSNE(reduction = "pca", dims = 1:npcs))

object[["raw.umap"]] <- CreateDimReducObject(embeddings = object@reductions[["umap"]]@cell.embeddings,
                                             key = "rawUMAP_", assay = DefaultAssay(object))
colnames(object[["raw.umap"]]@cell.embeddings) %<>% paste0("raw-",.)

object[["raw.tsne"]] <- CreateDimReducObject(embeddings = object@reductions[["tsne"]]@cell.embeddings,
                                             key = "rawtSNE_", assay = DefaultAssay(object))
colnames(object[["raw.tsne"]]@cell.embeddings) %<>% paste0("raw-",.)
object@reductions$umap = NULL
object@reductions$tsne = NULL
saveRDS(object, file = "data/MCL_61_20220318.rds")


#======1.8 UMAP from harmony =========================

npcs = 100
jpeg(paste0(path,"S1_RunHarmony.jpeg"), units="in", width=10, height=7,res=600)
system.time(object %<>% RunHarmony.1(group.by = "orig.ident", dims.use = 1:npcs,
                                     theta = 2, plot_convergence = TRUE,
                                     nclust = 50, max.iter.cluster = 100))
dev.off()

object %<>% FindNeighbors(reduction = "harmony",dims = 1:npcs)
object %<>% FindClusters(resolution = 0.8)
object %<>% RunUMAP(reduction = "harmony", dims = 1:npcs)
system.time(object %<>% RunTSNE(reduction = "harmony", dims = 1:npcs))

saveRDS(object, file = "data/MCL_61_20220318.rds")

# ======1.8.5 Unimodal UMAP Projection =========================
object = readRDS("data/MCL_61_20220318.rds")
DefaultAssay(object) = "SCT"
object.reference = readRDS("data/MCL_SCT_51_20210724.rds")
object.reference %<>% RunUMAP(reduction = "harmony", dims = 1:100,return.model = TRUE)

length(setdiff(colnames(object),colnames(object_reference)))
newCells <- setdiff(colnames(object),colnames(object_reference))
object.query = object[,newCells]

anchors <- FindTransferAnchors(reference = object.reference, query = object.query,
                                        dims = 1:30, reference.reduction = "pca")
object.query <- MapQuery(anchorset = anchors,
                           reference = object.reference, query = object.query,
                           refdata = list(cell.types = "cell.types"),
                         reference.reduction = "pca", reduction.model = "umap")

object[["umap"]]@cell.embeddings[newCells,] = object.query[["umap"]]@cell.embeddings[newCells,]
oldCells <- intersect(colnames(object),colnames(object.reference))
object[["umap"]]@cell.embeddings[oldCells,] = object.reference[["umap"]]@cell.embeddings[oldCells,]
saveRDS(object, file = "data/MCL_61_20220318.rds")

#=======1.9 save SCT only =======================================
format(object.size(object),unit = "GB")

format(object.size(object@assays$RNA),unit = "GB")
format(object.size(object@assays$integrated),unit = "GB")
object[['RNA']] <- NULL
object[['integrated']] <- NULL
format(object.size(object),unit = "GB")

object[["SCT"]]@scale.data = matrix(0,0,0)
object[["SCT"]]@counts = matrix(0,0,0)

saveRDS(object, file = "data/MCL_SCT_61_20220318.rds")


