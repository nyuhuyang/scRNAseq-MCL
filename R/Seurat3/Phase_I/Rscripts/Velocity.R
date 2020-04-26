# conda activate r3.6.2
invisible(lapply(c("Seurat","velocyto.R","SeuratWrappers","dplyr",
                   "magrittr"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("../R/Seurat3_functions.R")
path <- "Yang/20200409_Velocity/"
if(!dir.exists(path))dir.create(path, recursive = T)
#SBATCH --mem=64G
# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
i <- as.numeric(slurm_arrayid)
print(paste0("slurm_arrayid=",i))

loom_files = list.files("data/velocyto/",pattern = "loom")
(loom_file = loom_files[i])
(sample = gsub("\\.loom","",loom_file))
if(!dir.exists(paste0(path, sample))) dir.create(paste0(path, sample), recursive = T)

object = readRDS(file = "data/MCL_41_B_20200225.rds")
object %<>% AddMetaColor(label = "X4clusters", colors = gg_color_hue(n=4))
object %<>% AddMetaColor(label = "X5clusters", colors = c(gg_color_hue(n=4),"#D89000"))

Idents(object) = "orig.ident"
object %<>% subset(cells = c("Pt-U04_CTAATGGCATCGGGTC",
                                   "Pt-28-LN-C1D1_GCTGCTTGTTGGACCC"), invert = T)
object %<>% subset(idents = sample)

ldat <- ReadVelocity(file = paste0("data/velocyto/",loom_file))
RNA_velocyto <- as.Seurat(x = ldat)

# rename and subset Seurat
Barcode = gsub(".*_","",colnames(object))
Barcode %<>% paste0(sample,":",.,"x")
object %<>% RenameCells(new.names = Barcode)
RNA_velocyto %<>% subset(cells = colnames(object))
RNA_velocyto
RNA_velocyto %<>% subset(cells = c("PtU04:CTAATGGCATCGGGTCx",
                                   "Pt28_LN1:GCTGCTTGTTGGACCCx"), invert = T)
inherit = F
if(inherit) {
        # inherit umap
        RNA_velocyto@assays$SCT = object@assays$SCT
        rownames(object@reductions$pca@cell.embeddings) = colnames(object)
        rownames(object@reductions$tsne@cell.embeddings) = colnames(object)
        rownames(object@reductions$harmony@cell.embeddings) = colnames(object)
        RNA_velocyto@reductions$harmony = object@reductions$harmony
        RNA_velocyto@reductions$tsne = object@reductions$tsne
        RNA_velocyto %<>% RunVelocity(deltaT = 1, kCells = 25, fit.quantile = 0.02, reduction = "harmony")
} else {
        RNA_velocyto %<>% SCTransform(assay = "spliced")
        RNA_velocyto %<>% RunPCA(verbose = F)
        RNA_velocyto %<>% FindNeighbors(dims = 1:20)
        RNA_velocyto %<>% FindClusters()
        RNA_velocyto %<>% RunUMAP(dims = 1:20)
        RNA_velocyto %<>% RunTSNE(dims = 1:20)
}

RNA_velocyto@meta.data %<>% cbind(object@meta.data)
RNA_velocyto[["singler1main"]] =NULL
idents <- c("X4clusters","X5clusters","cell.types")
reductions = c("umap","tsne")
for(reduction in reductions){
        if(!inherit) RNA_velocyto %<>% RunVelocity(deltaT = 1, kCells = 25, fit.quantile = 0.02, reduction = reduction)
        for(ident in idents){
                Idents(RNA_velocyto) = ident
                RNA_velocyto %<>% sortIdent()
                ident.colors = ExtractMetaColor(RNA_velocyto)
                cell.colors <- ident.colors[RNA_velocyto@meta.data[,g]]
                names(x = cell.colors) <- colnames(x = RNA_velocyto)
                RNA_velocyto %<>% prepare.velocity.on.embedding.cor(n = 200, reduction = reduction,scale = "sqrt")
                jpeg(paste0(path, sample,"/velocity_",sample,"_",reduction,"_", ident,".jpeg"), units="in", width=10, height=7,res=600)
                show.velocity.on.embedding.cor(emb = Embeddings(object = RNA_velocyto, reduction = reduction), 
                                               vel = Tool(object = RNA_velocyto, slot = "RunVelocity"), n = 200, scale = "sqrt", 
                                               cell.colors = ac(x = cell.colors, alpha = 0.8), 
                                               cex = 1.5, arrow.scale = 1, show.grid.flow = TRUE, min.grid.cell.mass = 1, 
                                               fixed.arrow.length = FALSE,
                                               grid.n = 40, arrow.lwd = 1.5, max.grid.arrow.length =0.10,
                                               do.par = FALSE, cell.border.alpha = 0.1)
                dev.off()  
        }
}
for(ident in idents){
        lapply(c(TSNEPlot.1, UMAPPlot.1), function(fun){
                Idents(object) = ident
                fun(RNA_velocyto, group.by=ident,pt.size = 3,label = T,
                    label.repel = T,alpha = 0.9,
                    do.return = F,
                    no.legend = T,label.size = 4, repel = T, 
                    #title = paste("res =",res[i]," based on harmony"),
                    save.path = paste0(path, sample,"/"),
                    do.print = T)
        })

}

#saveRDS(RNA_velocyto, paste0(path,sample,"/velocity",sample,".rds"))

