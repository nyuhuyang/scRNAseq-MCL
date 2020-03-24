# conda activate r3.6
invisible(lapply(c("Seurat","velocyto.R","SeuratWrappers","dplyr",
                   "magrittr"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
i <- as.numeric(slurm_arrayid)
print(paste0("slurm_arrayid=",i))

loom_files = list.files("data/velocyto/",pattern = "loom")
(loom_file = loom_files[i])
(sample = gsub("\\.loom","",loom_file))

# load data
object = readRDS(file = "data/MCL_41_B_20200225.rds")
Idents(object) = "orig.ident"
object %<>% subset(idents = sample)

ldat <- ReadVelocity(file = paste0("data/velocyto/",loom_file))
RNA_velocyto <- as.Seurat(x = ldat)

# rename and subset Seurat
Barcode = gsub(".*_","",colnames(object))
Barcode %<>% paste0(sample,":",.,"x")
object %<>% RenameCells(new.names = Barcode)
RNA_velocyto %<>% subset(cells = colnames(object))
RNA_velocyto

# inherit umap
RNA_velocyto@assays$SCT = object@assays$SCT
rownames(object@reductions$pca@cell.embeddings) = colnames(object)
rownames(object@reductions$tsne@cell.embeddings) = colnames(object)
rownames(object@reductions$harmony@cell.embeddings) = colnames(object)
RNA_velocyto@reductions$harmony = object@reductions$harmony
RNA_velocyto@reductions$tsne = object@reductions$tsne
RNA_velocyto@meta.data %<>% cbind(object@meta.data)

# RunVelocity

RNA_velocyto %<>% RunVelocity(deltaT = 1, kCells = 25, fit.quantile = 0.5, reduction = "harmony")
idents <- c("X4clusters","cell.types")
RNA_velocyto %<>% subset(cells = c("PtU04:CTAATGGCATCGGGTCx"), invert = T)

for(ident in idents){
        Idents(RNA_velocyto) = ident
        RNA_velocyto %<>% sortIdent()
        ident.colors = ExtractMetaColor(RNA_velocyto)
        names(x = ident.colors) <- levels(x = RNA_velocyto)
        cell.colors <- ident.colors[Idents(object = RNA_velocyto)]
        names(x = cell.colors) <- colnames(x = RNA_velocyto)
        TSNEPlot.1(RNA_velocyto, group.by = ident, cols = unique(cell.colors),
                   label = T, label.repel = T,
                   do.print = T,do.return = F, save.path = paste0(path,sample,"/"))
        RNA_velocyto %<>% prepare.velocity.on.embedding.cor(n = 200, reduction = "tsne",scale = "sqrt")
        jpeg(paste0(path, sample,"/velocity_",sample,"_",ident,".jpeg"), units="in", width=10, height=7,res=600)
        show.velocity.on.embedding.cor(emb = Embeddings(object = RNA_velocyto, reduction = "tsne"), 
                                       vel = Tool(object = RNA_velocyto, slot = "RunVelocity"), n = 200, scale = "sqrt", 
                                       cell.colors = ac(x = cell.colors, alpha = 0.8), 
                                       cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, 
                                       grid.n = 40, arrow.lwd = 1, 
                                       do.par = FALSE, cell.border.alpha = 0.1)
        dev.off()  
}
saveRDS(RNA_velocyto, paste0(path,sample,"/velocity",sample,".rds"))
