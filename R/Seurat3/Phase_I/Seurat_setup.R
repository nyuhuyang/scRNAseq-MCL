########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
#devtools::install_github("immunogenomics/harmony", ref= "ee0877a",force = T)
invisible(lapply(c("Seurat","dplyr","kableExtra","ggplot2","cowplot","sctransform",
                   "harmony"), function(x) {
    suppressPackageStartupMessages(library(x,character.only = T))
}))
source("../R/Seurat3_functions.R")
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
df_samples <- readxl::read_excel("doc/191120_scRNAseq_info.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
(attach(df_samples))
samples = df_samples$sample

#======1.2 load  SingleCellExperiment =========================
(load(file = "data/sce_41_20191205.Rda"))
names(sce_list)
object_list <- lapply(sce_list, as.Seurat)

for(i in 1:length(samples)){
        object_list[[i]]@meta.data$tests <- df_samples$tests[i]
        object_list[[i]]@meta.data$conditions <- df_samples$notes[i]
        object_list[[i]]@meta.data$projects <- df_samples$project[i]
        object_list[[i]]@meta.data$groups <- df_samples$`sample title`[i]
        object_list[[i]]@meta.data$tissues <- df_samples$tissue[i]
}
#========1.3 merge ===================================
object <- Reduce(function(x, y) merge(x, y, do.normalize = F), object_list)
object@assays$RNA@data = object@assays$RNA@data *log(2)
remove(sce_list,object_list);GC()

(remove <- which(colnames(object@meta.data) %in% "ident"))
meta.data = object@meta.data[,-remove]
object@meta.data = meta.data 
remove(meta.data);GC()


#======== 1.4 FindVariableFeatures ===================================
# After removing unwanted cells from the dataset, the next step is to normalize the data.
#object <- NormalizeData(object = object, normalization.method = "LogNormalize", 
#                      scale.factor = 10000)
object <- FindVariableFeatures(object = object, selection.method = "vst",
                               num.bin = 20,
                               mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))

# Identify the 20 most highly variable genes
top20 <- head(VariableFeatures(object), 20)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(object)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
jpeg(paste0(path,"VariableFeaturePlot.jpeg"), units="in", width=10, height=7,res=600)
print(plot2)
dev.off()

#======1.6 Performing SCTransform and integration =========================
set.seed(100)
object_list <- SplitObject(object, split.by = "orig.ident")
remove(object);GC()
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

object %<>% RunPCA(npcs = 100, verbose = FALSE)
object <- JackStraw(object, num.replicate = 20,dims = 100)
object <- ScoreJackStraw(object, dims = 1:100)

jpeg(paste0(path,"JackStrawPlot.jpeg"), units="in", width=10, height=7,res=600)
JackStrawPlot(object, dims = 80:90)+
    ggtitle("JackStrawPlot")+
    theme(text = element_text(size=15),	
          plot.title = element_text(hjust = 0.5,size = 18))
dev.off()
npcs =50
object %<>% FindNeighbors(reduction = "pca",dims = 1:npcs)
object %<>% FindClusters(resolution = 0.8)
object %<>% RunTSNE(reduction = "pca", dims = 1:npcs)
object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)
p2 <- TSNEPlot.1(object, group.by="orig.ident",pt.size = 1,label = F,legend.size = 15,
                 label.size = 4, repel = T,title = "CCA Intergrated tSNE plot")
p3 <- UMAPPlot.1(object, group.by="orig.ident",pt.size = 1,label = F,legend.size = 15,
                 label.size = 4, repel = T,title = "CCA Intergrated UMAP plot")
save(object, file = "data/MCL_41_20191205.Rda")

#======1.8 old Harmony =========================
(load(file = "data/MCL_41_20191205.Rda"))
DefaultAssay(object)  = "SCT"
object <- FindVariableFeatures(object = object, selection.method = "vst",
                               num.bin = 20, nfeatures = 2000,
                               mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))
object %<>% ScaleData
object %<>% RunPCA(verbose = T,npcs = 100)

jpeg(paste0(path,"S1_ElbowPlot.jpeg"), units="in", width=10, height=7,res=600)
ElbowPlot(object, ndims = 100)
dev.off()

object %<>% JackStraw(num.replicate = 20,dims = 100)
object %<>% ScoreJackStraw(dims = 1:100)
a <- seq(1,100, by = 10)
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
        label.repel = T,alpha = 0.9,cols = Singler.colors,
        no.legend = T,label.size = 4, repel = T, title = "Harmony Integration",
        do.print = T, do.return = F))
object@reductions$tsne@cell.embeddings = - object@reductions$tsne@cell.embeddings
object@assays$RNA@scale.data = matrix(0,0,0)
object@assays$SCT@scale.data = matrix(0,0,0)
object@assays$integrated@scale.data = matrix(0,0,0)

object$sample = object$orig.ident
object$orig.ident %<>% plyr::mapvalues(from = unique(df_samples$sample),
                                       to = unique(df_samples$`sample name`))
table(object$orig.ident)
object$groups = gsub("_.*", "", object$orig.ident)
object$groups %<>% gsub("N01|N02|N03","Normal",.)
object$groups %<>% gsub("PtU01|PtU02|PtU03|PtU04","Untreated",.)
table(object$groups)
save(object, file = "data/MCL_41_harmony_20191231.Rda")

format(object.size(object[["SCT"]]),unit = "GB")
