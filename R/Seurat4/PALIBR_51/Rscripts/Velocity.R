# conda activate r4.0.3
invisible(lapply(c("Seurat","SeuratDisk","SeuratWrappers","dplyr","plyr","velocyto.R",
                   "magrittr"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
#SBATCH --mem=64G
# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
i <- as.numeric(slurm_arrayid)
print(paste0("slurm_arrayid=",i))
# loom
loom_files = list.files("data/velocyto/",pattern = "merged.loom")
(loom_file = loom_files[i])
(sample = gsub("\\.loom","",loom_file))
if(!dir.exists(paste0(path, sample))) dir.create(paste0(path, sample), recursive = T)
ldat <- ReadVelocity(file = paste0("data/velocyto/",loom_file))
RNA_velocyto <- as.Seurat(x = ldat)

df_samples <- readxl::read_excel("doc/20210715_scRNAseq_info.xlsx", sheet = "fastq")
df_samples = as.data.frame(df_samples)
colnames(df_samples) %<>% tolower()
df_samples %<>% filter(sequence %in% "GEX") %>% filter(phase %in% "PALIBR_I")
RNA_velocyto$orig.ident = gsub(":.*","",colnames(RNA_velocyto))
RNA_velocyto$orig.ident %<>% mapvalues(from = df_samples$sample.old,
                                       to = df_samples$sample)
RNA_velocyto$barcode = gsub(".*:","",colnames(RNA_velocyto)) %>% gsub("x$","",.)
new_names = paste0(RNA_velocyto$orig.ident, "-", RNA_velocyto$barcode)
RNA_velocyto %<>% RenameCells(new.names = new_names)

# MCL
object = readRDS(file = "data/MCL_SCT_51_20210724.rds")

Idents(object) = "X4cluster"
object %<>% subset(idents = c("1","2","3","4"))
object %<>% AddMetaColor(label = "X4cluster", colors = gg_color_hue(n=4))


# clean cells

common = intersect(new_names, colnames(object))
length(common)
# subset Seurat
RNA_velocyto %<>% subset(cells = common)
object %<>% subset(cells = common)
RNA_velocyto@meta.data %<>% cbind(object@meta.data)


inherit = T
if(inherit) {
        # inherit umap
        RNA_velocyto[["RNA"]] = RNA_velocyto[["spliced"]]
        RNA_velocyto %<>% SCTransform
        RNA_velocyto@reductions = object@reductions
} else {
        RNA_velocyto %<>% SCTransform(assay = "spliced")
        RNA_velocyto %<>% RunPCA(verbose = F)
        RNA_velocyto %<>% FindNeighbors(dims = 1:20)
        RNA_velocyto %<>% FindClusters()
        RNA_velocyto %<>% RunUMAP(dims = 1:20)
        RNA_velocyto %<>% RunTSNE(dims = 1:20)
}

idents <- c("X4cluster","cell.types","orig.ident")
table(idents %in% colnames(RNA_velocyto@meta.data))
lapply(idents, function(ident){
        Idents(object) = ident
        TSNEPlot.1(RNA_velocyto, group.by=ident,pt.size = 3,label = T,
                   label.repel = T,alpha = 0.9,
                   do.return = F,
                   no.legend = T,label.size = 4, repel = T,
                   save.path = paste0(path, sample,"/"),
                   do.print = T)
})


saveRDS(RNA_velocyto, paste0(path,sample,"/velocity",sample,".rds"))
RNA_velocyto = readRDS(paste0(path,sample,"/velocity",sample,".rds"))
RNA_velocyto@reductions$pca = RNA_velocyto@reductions$raw.tsne
DefaultAssay(RNA_velocyto) <- "RNA"
file.remove(paste0(path,sample,"/velocity",sample,".h5ad"))
SaveH5Seurat(RNA_velocyto, filename = paste0(path,sample,"/velocity",sample,".h5Seurat"))
Convert(paste0(path,sample,"/velocity",sample,".h5Seurat"), dest = "h5ad")
file.remove(paste0(path,sample,"/velocity",sample,".h5Seurat"))
