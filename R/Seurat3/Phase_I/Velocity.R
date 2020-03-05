library(Seurat)
library(velocyto.R)
library(SeuratWrappers)
library(dplyr)
library(magrittr)
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

# load data
(load(file="data/MCL_41_harmony_20200225.Rda"))
Idents(object) = "groups"
object %<>% subset(idents = "Pt25")
ldat <- ReadVelocity(file = "data/velocyto/Pt-25_merged.loom")
RNA_velocyto <- as.Seurat(x = ldat)

# rename Seurat
Barcode = colnames(object)
Barcode %<>% gsub("_",":",.)
Barcode %<>% paste0("x")
object %<>% RenameCells(new.names = Barcode)
RNA_velocyto@assays$SCT = object@assays$SCT
rownames(object@reductions$pca@cell.embeddings) = colnames(object)
rownames(object@reductions$tsne@cell.embeddings) = colnames(object)
rownames(object@reductions$umap@cell.embeddings) = colnames(object)
rownames(object@reductions$harmony@cell.embeddings) = colnames(object)

# subset RNA_velocyto
RNA_velocyto %<>% subset(cells = colnames(object))
RNA_velocyto@reductions$harmony = object@reductions$harmony
RNA_velocyto@reductions$umap = object@reductions$umap
RNA_velocyto@meta.data %<>% cbind(object@meta.data)


# RunVelocity
RNA_velocyto %<>% RunVelocity(deltaT = 1, kCells = 25, fit.quantile = 0.5, reduction = "harmony")
Idents(RNA_velocyto) = "cell.types"
RNA_velocyto %<>% sortIdent()
ident.colors = ExtractMetaColor(RNA_velocyto)
names(x = ident.colors) <- levels(x = RNA_velocyto)
cell.colors <- ident.colors[Idents(object = RNA_velocyto)]
names(x = cell.colors) <- colnames(x = RNA_velocyto)

RNA_velocyto %<>% prepare.velocity.on.embedding.cor(n = 200, reduction = "tsne",scale = "sqrt")

jpeg(paste0(path,"velocity_Pt25.jpeg"), units="in", width=10, height=7,res=600)
show.velocity.on.embedding.cor(emb = Embeddings(object = RNA_velocyto, reduction = "tsne"), 
                               vel = Tool(object = RNA_velocyto, slot = "RunVelocity"), n = 200, scale = "sqrt", 
                               cell.colors = ac(x = cell.colors, alpha = 0.5), 
                               cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, 
                               grid.n = 40, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0.1)
dev.off()
