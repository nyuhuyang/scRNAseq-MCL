########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
#
# ######################################################################
# conda activate r4.0
#devtools::install_github("immunogenomics/harmony", ref= "ee0877a",force = T)
invisible(lapply(c("Seurat","dplyr","kableExtra","ggplot2","cowplot","sctransform",
                   "harmony"), function(x) {
    suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
if(!dir.exists("data")) dir.create("data")
if(!dir.exists("doc")) dir.create("doc")

########################################################################
#
#  1 Seurat Alignment
#
# ######################################################################
#======1.1 Setup the Seurat objects =========================
# read sample summary list
df_samples <- readxl::read_excel("doc/20210311_scRNAseq_info.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
(attach(df_samples))
sample = df_samples$`sample name`

#======1.2 load  SingleCellExperiment =========================
(load(file = "data/MCL_AIM_74_20210311.Rda"))
meta.data = object@meta.data
for(i in 1:length(samples)){
        cells <- meta.data$orig.ident %in% samples[i]
        meta.data[cells,"patient"]    = df_samples$patient[i]
        meta.data[cells,"project"]    = df_samples$project[i]
        meta.data[cells,"phase"]     = df_samples$phase[i]
        meta.data[cells,"tissue"]     = df_samples$tissue[i]
        meta.data[cells,"conditions"] = df_samples$conditions[i]
        meta.data[cells,"disease"]    = df_samples$disease[i]
}
meta.data$orig.ident %<>% factor(levels = df_samples$`sample name`)


table(rownames(object@meta.data) == rownames(meta.data))
table(colnames(object) == rownames(meta.data))

object@meta.data = meta.data

#======1.6 Performing SCTransform and integration =========================
set.seed(100)
object_list <- SplitObject(object, split.by = "orig.ident")
remove(object)GC()
#save(object_list, file = "data/object_list_27_20201008.Rda")
#(load(file = "data/object_list_27_20201008.Rda"))
object_list %<>% lapply(SCTransform)
object.features <- SelectIntegrationFeatures(object_list, nfeatures = 3000)
options(future.globals.maxSize= object.size(object_list)*1.5)
npcs = 30
object_list %<>% lapply(function(x) {
    x %<>% RunPCA(features = object.features, verbose = FALSE)
})

object_list <- PrepSCTIntegration(object.list = object_list, anchor.features = object.features,
                                  verbose = FALSE)
anchors <- FindIntegrationAnchors(object_list, normalization.method = "SCT",
                                  anchor.features = object.features,
                                  reference = c(1, 2), reduction = "rpca",
                                  dims = 1:npcs)
remove(object_list);GC()

object <- IntegrateData(anchorset = anchors, normalization.method = "SCT",dims = 1:npcs)

remove(anchors);GC()
save(object, file = "data/MCL_AIM_74_20210311.Rda")

#======1.8  Harmony =========================
(load(file = "data/MCL_AIM_74_20210311.Rda"))

DefaultAssay(object)  = "SCT"
object <- FindVariableFeatures(object = object, selection.method = "vst",
                               num.bin = 20, nfeatures = 2000,
                               mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))
object %<>% ScaleData
object %<>% RunPCA(verbose = T,npcs = 150)

jpeg(paste0(path,"S1_ElbowPlot.jpeg"), units="in", width=10, height=7,res=600)
ElbowPlot(object, ndims = 150)
dev.off()

object %<>% JackStraw(num.replicate = 20,dims = 150)
object %<>% ScoreJackStraw(dims = 1:150)
a <- seq(1,150, by = 10)
b <- a+9
for(i in seq_along(a)){
    jpeg(paste0(path,"JackStrawPlot_",i,"_",a[i],"_",min(b[i],100),".jpeg"), units="in", width=10, height=7,res=600)
    print(JackStrawPlot(object, dims = a[i]:min(b[i],100)))
    Progress(i,length(a))
    dev.off()
}

npcs = 85
jpeg(paste0(path,"S1_RunHarmony.jpeg"), units="in", width=10, height=7,res=600)
system.time(object %<>% RunHarmony.1(group.by = "orig.ident", dims.use = 1:npcs,
                                   theta = 2, plot_convergence = TRUE,
                                   nclust = 50, max.iter.cluster = 100))
dev.off()

object %<>% FindNeighbors(reduction = "harmony",dims = 1:npcs)
object %<>% FindClusters(resolution = 0.8)
system.time(object %<>% RunTSNE(reduction = "harmony", dims = 1:npcs))
object %<>% RunUMAP(reduction = "harmony", dims = 1:npcs)
p4 <- TSNEPlot.1(object, group.by="orig.ident",pt.size = 1,label = F,legend.size = 15,
                 label.size = 4, repel = T,title = "Harmony intergrated tSNE plot")
p5 <- UMAPPlot.1(object, group.by="orig.ident",pt.size = 1,label = F,legend.size = 15,
                 label.size = 4, repel = T,title = "Harmony intergrated UMAP plot")

#=======1.9 summary =======================================
lapply(c(TSNEPlot.1, UMAPPlot.1), function(fun)
    fun(object, group.by="orig.ident",pt.size = 0.5,label = F,
        label.repel = T,alpha = 0.9,cols = c(Singler.colors,Singler.colors),
        no.legend = T,label.size = 4, repel = T, title = "Harmony Integration",
        do.print = T, do.return = F))

tmp = object@reductions$tsne@cell.embeddings[,"tSNE_2"]
object@reductions$tsne@cell.embeddings[,"tSNE_2"] = object@reductions$tsne@cell.embeddings[,"tSNE_1"]
object@reductions$tsne@cell.embeddings[,"tSNE_1"] = tmp

tmp = object@reductions$umap@cell.embeddings[,"UMAP_2"]
object@reductions$umap@cell.embeddings[,"UMAP_2"] = object@reductions$umap@cell.embeddings[,"UMAP_1"]
object@reductions$umap@cell.embeddings[,"UMAP_1"] = tmp
object@reductions$umap@cell.embeddings[,"UMAP_1"] = - object@reductions$umap@cell.embeddings[,"UMAP_1"]

object@assays$RNA@scale.data = matrix(0,0,0)
object@assays$SCT@scale.data = matrix(0,0,0)
object@assays$integrated = NULL


save(object, file = "data/MCL_AIM_74_20210311.Rda")

format(object.size(object@assays$RNA),unit = "GB")
object@assays$RNA = NULL
save(object, file = "data/MCL_AIM_74_20210311_SCT.Rda")
