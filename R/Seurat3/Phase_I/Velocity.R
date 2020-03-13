library(Seurat)
library(R.utils)
library(velocyto.R)
library(SeuratWrappers)
library(dplyr)
library(magrittr)
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

# on linux, make folder, move files
folder_list <- list.files("data/bam", pattern = "\\.bam$")
folder_list = gsub("\\.bam","",folder_list)
df_samples <- readxl::read_excel("doc/191120_scRNAseq_info.xlsx")
df_samples$`Sample Name`[df_samples$Sample.ID %in% folder_list]

for(f in folder_list) {
        sample <- df_samples$`Sample Name`[df_samples$Sample.ID %in% f]
        dir.create(paste0("data/bam/",sample), recursive = T)
        file.rename(paste0("data/bam/",f,".bam"), 
                  paste0("data/bam/",sample,"/",sample,".bam"))
        file.rename(paste0("data/bam/",f,".bam.bai"), 
                  paste0("data/bam/",sample,"/",sample,".bam.bai"))
}

# on mac, make folder, copy barcodes.tsv
df_samples <- readxl::read_excel("doc/191120_scRNAseq_info.xlsx")
folder_list <-  list.files("data/scRNA-seq")
table(df_samples$Sample.ID %in% folder_list)
for(f in folder_list) {
        s <- df_samples[df_samples$Sample.ID %in% f, "Sample Name"] %>% as.character()
        if(!dir.exists(paste0("data/velocyto/",s)))  dir.create(paste0("data/velocyto/",s), recursive = T)
        file.copy(paste0("data/scRNA-seq/",f,"/outs/filtered_gene_bc_matrices/hg19/barcodes.tsv"), 
                  paste0("data/velocyto/",s,"/barcodes.tsv"))
        file.copy(paste0("data/scRNA-seq/",f,"/outs/filtered_gene_bc_matrices/hg19/barcodes.tsv.gz"), 
                  paste0("data/velocyto/",s,"/barcodes.tsv.gz"))
        if(file.exists(paste0("data/velocyto/",s,"/barcodes.tsv.gz"))) gunzip(paste0("data/velocyto/",s,"/barcodes.tsv.gz")) 
        Progress(which(folder_list %in% f), length(folder_list))
}
# load data
object = readRDS(file = "data/MCL_41_B_20200225.rds")
Idents(object) = "groups"
object %<>% subset(idents = "Pt25")
ldat <- ReadVelocity(file = "data/velocyto/Pt-25_merged.loom")
RNA_velocyto <- as.Seurat(x = ldat)

# rename and subset Seurat
Barcode = colnames(object)
Barcode %<>% gsub("_",":",.)
Barcode %<>% paste0("x")
object %<>% RenameCells(new.names = Barcode)
RNA_velocyto %<>% subset(cells = colnames(object))

choose = c("inherit","redo")[2]
# inherit umap
if(choose == "inherit"){
        RNA_velocyto@assays$SCT = object@assays$SCT
        rownames(object@reductions$pca@cell.embeddings) = colnames(object)
        rownames(object@reductions$tsne@cell.embeddings) = colnames(object)
        rownames(object@reductions$umap@cell.embeddings) = colnames(object)
        rownames(object@reductions$harmony@cell.embeddings) = colnames(object)
        
        RNA_velocyto@reductions$harmony = object@reductions$harmony
        RNA_velocyto@reductions$umap = object@reductions$umap
        RNA_velocyto@meta.data %<>% cbind(object@meta.data)
}
# redo umap
if(choose == "redo"){
        table(rownames(RNA_velocyto@meta.data) == colnames(object))
        RNA_velocyto@meta.data = RNA_velocyto@meta.data[,-1]
        RNA_velocyto@meta.data %<>% cbind(object@meta.data)
        RNA_velocyto %<>% SCTransform(assay = "spliced")
        RNA_velocyto %<>% RunPCA(verbose = FALSE)
        RNA_velocyto %<>% FindNeighbors(dims = 1:50)
        RNA_velocyto %<>% FindClusters()
        RNA_velocyto %<>% RunUMAP(dims = 1:50)
        RNA_velocyto %<>% RunTSNE(dims = 1:50)
}
groups <- c("orig.ident","X4clusters","tissues","cell.types")
lapply(groups, function(g) {
        Idents(RNA_velocyto) = g
        TSNEPlot.1(RNA_velocyto, group.by = g, cols = ExtractMetaColor(RNA_velocyto),
                   do.print = T, save.path = paste0(path,"Pt-25/"))
})


# RunVelocity

RNA_velocyto %<>% RunVelocity(deltaT = 1, kCells = 25, fit.quantile = 0.5, reduction = "pca")
reductions = c("tsne", "umap")
for(redu in reductions){
        for(g in groups){
                Idents(RNA_velocyto) = g
                RNA_velocyto %<>% sortIdent()
                ident.colors = ExtractMetaColor(RNA_velocyto)
                names(x = ident.colors) <- levels(x = RNA_velocyto)
                cell.colors <- ident.colors[Idents(object = RNA_velocyto)]
                names(x = cell.colors) <- colnames(x = RNA_velocyto)
                TSNEPlot.1(RNA_velocyto, group.by = g, cols = unique(cell.colors),
                           label = T, label.repel = T,
                           do.print = T,do.return = F, save.path = paste0(path,"Pt-25/"))
                UMAPPlot.1(RNA_velocyto,do.return = F,  group.by = g, cols = unique(cell.colors),
                           do.print = T, save.path = paste0(path,"Pt-25/"))
                RNA_velocyto %<>% prepare.velocity.on.embedding.cor(n = 200, reduction = redu,scale = "sqrt")
                jpeg(paste0(path,"Pt-25/","velocity_Pt25_",g,"_",redu,".jpeg"), units="in", width=10, height=7,res=600)
                show.velocity.on.embedding.cor(emb = Embeddings(object = RNA_velocyto, reduction = redu), 
                                               vel = Tool(object = RNA_velocyto, slot = "RunVelocity"), n = 200, scale = "sqrt", 
                                               cell.colors = ac(x = cell.colors, alpha = 0.5), 
                                               cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, 
                                               grid.n = 40, arrow.lwd = 1, 
                                               do.par = FALSE, cell.border.alpha = 0.1)
                dev.off()  
        }

}
saveRDS(RNA_velocyto, paste0(path,"Pt-25/","velocity.rds"))
