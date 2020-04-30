library(Seurat)
library(R.utils)
library(velocyto.R)
library(SeuratWrappers)
library(dplyr)
library(magrittr)
source("../R/Seurat3_functions.R")
source("R/util.R")
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
Idents(object) = "orig.ident"

arg =1
groups = c("Pt-17","Pt-25")

sample_pairs = list(c("Pt17_LN1","Pt17_2","Pt17_7"),#7
                    c("Pt25_SB1","Pt25_24","Pt25_25Pd","Pt25_AMB25Pd"))
(sample = sample_pairs[[arg]])
object %<>% subset(idents = sample)
ldat <- ReadVelocity(file = paste0("data/velocyto/",groups[arg],"_merged.loom"))

RNA_velocyto <- as.Seurat(x = ldat)
rm(ldat)
save.path = paste0(path, paste(sample, collapse = "-"),"/")
if(!dir.exists(save.path)) dir.create(save.path, recursive = T)
# rename and subset Seurat
Barcode = colnames(object)
Barcode %<>% gsub("_",":",.)
Barcode %<>% paste0("x")
object %<>% RenameCells(new.names = Barcode)
RNA_velocyto %<>% subset(cells = colnames(object))

choose = c("inherit","redo","monocle2")[3]
# inherit umap
if(choose == "inherit"){
        RNA_velocyto@assays$SCT = object@assays$SCT
        rownames(object@reductions$pca@cell.embeddings) = colnames(object)
        rownames(object@reductions$tsne@cell.embeddings) = colnames(object)
        rownames(object@reductions$umap@cell.embeddings) = colnames(object)
        rownames(object@reductions$harmony@cell.embeddings) = colnames(object)
        RNA_velocyto@reductions$harmony = object@reductions$tsne
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
# inherit ddtree from monocle
if(choose == "monocle2"){
        RNA_velocyto@assays$SCT = object@assays$SCT
        monocle.path = paste0("Yang/20200426_monocle2/",paste(sample, collapse = "-"),"/")
        cds = readRDS(paste0(monocle.path,"monocle2_",paste(sample, collapse = "-"),"_cds.rds"))
        object@reductions$tsne@cell.embeddings = t(cds@reducedDimS)
        colnames(object@reductions$tsne@cell.embeddings) = paste0("tSNE_",1:2)
        rownames(object@reductions$pca@cell.embeddings) = colnames(object)
        rownames(object@reductions$tsne@cell.embeddings) = colnames(object)
        RNA_velocyto@reductions$pca = object@reductions$pca
        RNA_velocyto@reductions$tsne = object@reductions$tsne
        RNA_velocyto[["orig.ident"]] = NULL
        RNA_velocyto@meta.data %<>% cbind(object@meta.data)
        RNA_velocyto[["Pseudotime"]] = cds$Pseudotime
}
if(length(unique(cds$orig.ident)) == 4) orig.ident_color <- c("#6A3D9A","#BEAED4","#B2DF8A","#BF5B17")
if(length(unique(cds$orig.ident)) == 3) orig.ident_color <- c("#6A3D9A","#B2DF8A","#BF5B17")
group_by <- c("orig.ident","X4clusters","cell.types","Pseudotime")[1:3]
RNA_velocyto %<>% AddMetaColor(label= "X4clusters", colors = gg_color_hue(length(unique(RNA_velocyto$X4clusters))))
RNA_velocyto %<>% AddMetaColor(label= "orig.ident", colors = orig.ident_color[1:(length(unique(RNA_velocyto$orig.ident)))])

for(g in group_by) {
        Idents(RNA_velocyto) = g
        TSNEPlot.1(RNA_velocyto, group.by = g, cols = ExtractMetaColor(RNA_velocyto),
                   do.print = T, save.path = save.path, unique.name = F)  
}

# RunVelocity
RNA_velocyto %<>% RunVelocity(deltaT = 1, kCells = 25, fit.quantile = 0.5, reduction = "tsne") #harmony
for(g in group_by){
        RNA_velocyto[[g]] = droplevels(object[[g]])
        Idents(RNA_velocyto) = g
        RNA_velocyto %<>% sortIdent()
        ident.colors = ExtractMetaColor(RNA_velocyto)
        cell.colors <- ident.colors[RNA_velocyto@meta.data[,g]]
        names(x = cell.colors) <- colnames(x = RNA_velocyto)
        RNA_velocyto %<>% prepare.velocity.on.embedding.cor(n = 200, reduction = "tsne",scale = "sqrt")
        jpeg(paste0(save.path, "/","velocity_",groups[arg],"_",g,".jpeg"), units="in", width=8, height=8,res=600)
        show.velocity.on.embedding.cor.1(emb = Embeddings(object = RNA_velocyto, reduction = "tsne"), 
                                       vel = Tool(object = RNA_velocyto, slot = "RunVelocity"), n = 200, scale = "sqrt", 
                                       cell.colors = ac(x = cell.colors, alpha = 0.8), 
                                       cex = 1.2, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, 
                                       grid.n = 40, arrow.lwd = 0.8, arrow.colors = "grey40",
                                       do.par = FALSE, cell.border.alpha = 0.5)
        dev.off()
        Progress(which(g %in% group_by), length(group_by))
}

saveRDS(RNA_velocyto, paste0(save.path,"/velocity",groups[arg],".rds"))
#RNA_velocyto = readRDS(paste0(save.path,"/velocity",groups[arg],".rds"))
#Trajectory step 2: generate RNA velocyto on Trajectory plot by each sample
save.path.sub = paste0(save.path, "subset/")
if(!dir.exists(save.path.sub)) dir.create(save.path.sub, recursive = T)

samples <- unique(RNA_velocyto$orig.ident) %>% as.character()
group_by <- c("X4clusters", "cell.types")
Idents(RNA_velocyto) = "orig.ident"
for(s in seq_along(samples)){
        subset_object <- subset(RNA_velocyto, idents = samples[s])
        for(k in seq_along(group_by)){
                Idents(subset_object) = group_by[k]
                ident.colors = ExtractMetaColor(subset_object)
                cell.colors <- ident.colors[subset_object@meta.data[,group_by[k]]]
                names(x = cell.colors) <- colnames(x = subset_object)
                RNA_velocyto %<>% prepare.velocity.on.embedding.cor(n = 200, reduction = "tsne",scale = "sqrt")
                jpeg(paste0(save.path.sub,"velocity_",samples[s],"_",group_by[k],".jpeg"), units="in", width=8, height=8,res=600)
                show.velocity.on.embedding.cor.1(emb = Embeddings(object = subset_object, reduction = "tsne"), 
                                               vel = Tool(object = subset_object, slot = "RunVelocity"), n = 200, scale = "sqrt", 
                                               cell.colors = ac(x = cell.colors, alpha = 0.8), 
                                               cex = 1.2, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, 
                                               grid.n = 40, arrow.lwd = 0.8, arrow.colors = "grey40",
                                               do.par = FALSE, cell.border.alpha = 0.5)
                dev.off()
                Progress((s-1)*length(samples)+k, length(samples)*length(group_by))
        }
}


#Trajectory step 3: generate RNA velocyto on Trajectory plot by each cluster
X4clusters <- sort(unique(cds$X4clusters))
file_name = paste(c(samples[1],gsub(".*_","",samples[2:length(samples)])),collapse = "_")
sample_name = paste(c(samples[1],gsub(".*_","",samples[2:length(samples)])),collapse = ", ")
Idents(RNA_velocyto) = "X4clusters"
for(c in seq_along(X4clusters)){
        subset_object <- subset(RNA_velocyto, idents = X4clusters[c])
        Idents(subset_object) = "orig.ident"
        subset_object %<>% sortIdent()
        ident.colors = ExtractMetaColor(subset_object)
        cell.colors <- ident.colors[droplevels(subset_object@meta.data[,"orig.ident"])]
        names(x = cell.colors) <- colnames(x = subset_object)
        RNA_velocyto %<>% prepare.velocity.on.embedding.cor(n = 200, reduction = "tsne",scale = "sqrt")
        jpeg(paste0(save.path.sub,"velocity_",X4clusters[c],"_",file_name,".jpeg"), units="in", width=8, height=8,res=600)
        show.velocity.on.embedding.cor.1(emb = Embeddings(object = subset_object, reduction = "tsne"), 
                                       vel = Tool(object = subset_object, slot = "RunVelocity"), n = 200, scale = "sqrt", 
                                       cell.colors = ac(x = cell.colors, alpha = 0.8), 
                                       cex = 1.2, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, 
                                       grid.n = 40, arrow.lwd = 0.8, arrow.colors = "grey40",
                                       do.par = FALSE, cell.border.alpha = 0.5)
        dev.off()
        Progress(c, length(X4clusters))
}


path <- "Yang/20200325_Velocity/"
(samples = list.files(path))
for(s in seq_along(samples)){
        file.remove(paste0(path,samples[s],"/velocity",samples[s],".rds"))
        file.remove(paste0(path,samples[s],"/TSNEPlot_RNA_velocyto_cell.types_Legend.jpeg"))
        file.remove(paste0(path,samples[s],"/TSNEPlot_RNA_velocyto_X4clusters_Legend.jpeg"))
}

path <- "Yang/20200325_monocle2/"
(samples = list.files(path))
for(s in seq_along(samples)){
        file.remove(paste0(path,samples[s],"/monocle2_",samples[s],"_DE.rds"))
}

# The colors assigned to the clusters look great.  
# But using the same colors for different samples becomes confusing.   
# Can you  please use a different set of 4 colors for SB1, C24, C25, C25-AMB for Pt25?   
(n <- length(unique(object$orig.ident)))
object = readRDS(paste0("Yang/20200421_Velocity/",groups[arg],"-original/velocity.rds"))
Idents(object) = "orig.ident"
object %<>% sortIdent()
for(i in 1:100){
        colors = sample(Singler.colors[1:18], n)
        object %<>% AddMetaColor(label= "orig.ident", colors = colors)
        TSNEPlot.1(object, title = paste(c("Colors:",colors), collapse = ", "),
                   do.print = T, save.path = paste0(path,i,"-"))
        Progress(i, 100)
}
