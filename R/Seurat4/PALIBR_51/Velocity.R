#conda activate r4.1.0
invisible(lapply(c("Seurat","SeuratDisk","SeuratWrappers","dplyr","plyr","velocyto.R",
                   "magrittr"), function(x) {
                       suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
#=============
df_samples <- readxl::read_excel("doc/20210715_scRNAseq_info.xlsx", sheet = "fastq")
df_samples = as.data.frame(df_samples)
colnames(df_samples) %<>% tolower()
df_samples %<>% filter(sequence %in% "GEX") %>% filter(phase %in% "PALIBR_I") %>%
    filter(sample != "Pt11_31")

ldat <- ReadVelocity(file = "../scRNAseq-AIM/data/velocyto/MCL51_merged.loom")
RNA_velocyto <- as.Seurat(x = ldat)
rm(ldat);GC()
RNA_velocyto$orig.ident = gsub(":.*","",colnames(RNA_velocyto))
RNA_velocyto$orig.ident %<>% plyr::mapvalues(from = df_samples$sample.id, df_samples$sample)
RNA_velocyto$barcode = gsub(".*:","",colnames(RNA_velocyto)) %>% gsub("x$","",.)
new_names = paste0(RNA_velocyto$orig.ident, "-", RNA_velocyto$barcode)
RNA_velocyto %<>% RenameCells(new.names = new_names)
RNA_velocyto[["RNA"]] <- RNA_velocyto[["spliced"]]
DefaultAssay(RNA_velocyto) = "RNA"

# MCL
object = readRDS(file = "data/MCL_SCT_51_20210724.rds")


#Idents(object) = "X4cluster"
#object %<>% subset(idents = c("1","2","3","4"))
#object %<>% AddMetaColor(label = "X4cluster", colors = gg_color_hue(n=4))

common = intersect(new_names,colnames(object))
length(common)
# subset Seurat
RNA_velocyto %<>% subset(cells = common)
object %<>% subset(cells = common)
meta.data = object@meta.data
meta.data = meta.data[common,]
if(all(rownames(RNA_velocyto@meta.data) == rownames(meta.data))) {
    RNA_velocyto@meta.data %<>% cbind(meta.data)
}
RNA_velocyto

RNA_velocyto[["tsne"]] <- CreateDimReducObject(embeddings = object@reductions[["tsne"]]@cell.embeddings[common,],
                                               key = "tSNE_", assay = "RNA")

RNA_velocyto[["umap"]] <- CreateDimReducObject(embeddings = object@reductions[["umap"]]@cell.embeddings[common,],
                                               key = "UMAP_", assay = "RNA")

DefaultAssay(RNA_velocyto) <- "RNA"
file.remove("data/velocyto/MCL51_merged.h5ad")
SaveH5Seurat(RNA_velocyto, filename = "data/velocyto/MCL51_merged.h5Seurat")
Convert("data/velocyto/MCL51_merged.h5Seurat", dest = "h5ad")
file.remove("data/velocyto/MCL51_merged.h5Seurat")
