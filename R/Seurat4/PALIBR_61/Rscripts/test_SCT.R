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
# Need 64GB
set.seed(101)
# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.integer(as.character(slurm_arrayid))
print(paste0("slurm_arrayid=",args))

save.path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(save.path)) dir.create(save.path, recursive = T)

#==================
test_df = data.frame(FindVariableFeatures = rep(c(TRUE, TRUE, FALSE),times = 2),
                     nfeatures = rep(c(3000,2000,3000), times = 2),
                     ScaleData = rep(c(TRUE, FALSE),each = 3))
print(findVariableFeatures <- test_df[args,"FindVariableFeatures"])
print(nfeatures <- test_df[args,"nfeatures"])
print(scaleData <- test_df[args,"ScaleData"])

file.name = paste0("FindVF.",as.character(findVariableFeatures),"_nfeatures.",nfeatures,"_scaleData.",as.character(scaleData))
#======1.2 load  Seurat =========================
object = readRDS(file = "data/MCL_61_20220331.rds")
object@meta.data = readRDS(file = "MCL_61_20220331_metadata.rds")

#======1.7 UMAP from raw pca =========================
if(findVariableFeatures){
    object <- FindVariableFeatures(object = object, selection.method = "vst",
                                   num.bin = 20, nfeatures = nfeatures,
                                   mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))
}
if(scaleData){
    object %<>% ScaleData(verbose = FALSE)
}

object %<>% RunPCA(verbose = T,npcs = 100)

#======1.8 UMAP from harmony =========================

npcs = 100
system.time(object %<>% RunHarmony.1(group.by = "orig.ident", dims.use = 1:npcs,
                                     theta = 2, plot_convergence = TRUE,
                                     nclust = 50, max.iter.cluster = 100))

object %<>% RunUMAP(reduction = "harmony", dims = 1:npcs)

UMAPPlot.1(object, group.by = "cell.types",do.print = T,raster=FALSE,no.legend = T,
           title = file.name,file.name = paste0(file.name,".jpeg"),
           save.path = save.path)


