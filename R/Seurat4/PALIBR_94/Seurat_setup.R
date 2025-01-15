
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
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)


########################################################################
#
#  1 Seurat Alignment
#
# ######################################################################
#======1.1 Setup the Seurat objects =========================
# read sample summary list
df_samples <- readxl::read_excel("output/20240614/20240525_scRNAseq_info_MD_SCK.xlsx")
df_samples = as.data.frame(df_samples)
colnames(df_samples) %<>% tolower()
nrow(df_samples)
#======1.2 load  Seurat =========================
object = readRDS(file = "data/MCL_94_20240614.rds")

meta.data = object@meta.data
for(i in 1:length(df_samples$sample)){
    cells <- meta.data$orig.ident %in% df_samples$sample_id[i]
    print(table(cells))
    meta.data[cells,"tissue"] = df_samples$tissue[i]
    meta.data[cells,"patient"] = df_samples$patient[i]
    meta.data[cells,"treatment"] = df_samples$treatment[i]
    meta.data[cells,"response"] = df_samples$`response at sampling`[i]
    meta.data[cells,"phase"] = df_samples$phase[i]
    meta.data[cells,"additional_groups"] = df_samples$additional_groups[i]
    meta.data[cells,"Mean.Reads.per.Cell"] = df_samples$mean.reads.per.cell[i]
    meta.data[cells,"Number.of.Reads"] = df_samples$number.of.reads[i]
    meta.data[cells,"Sequencing.Saturation"] = df_samples$sequencing.saturation[i]
}

meta.data$date %<>% as.character() %>% gsub(" .*","",.)
meta.data$orig.ident %<>% factor(levels = df_samples$sample_id)
table(rownames(object@meta.data) == rownames(meta.data))
table(colnames(object) == rownames(meta.data))
object@meta.data = meta.data
Idents(object) = "orig.ident"

# Determine the ‘dimensionality’ of the dataset  =========
npcs = 150

DefaultAssay(object) <- "RNA"
object %<>% NormalizeData()
object <- FindVariableFeatures(object = object, selection.method = "vst",
                               num.bin = 20, nfeatures = 3000,
                               mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))
object %<>% ScaleData(verbose = FALSE)
object %<>% RunPCA(npcs = npcs, verbose = FALSE)
jpeg(paste0(path,"ElbowPlot_RNA.jpeg"), units="in", width=10, height=7,res=600)
print(ElbowPlot(object,ndims = npcs))
dev.off()
object <- JackStraw(object, num.replicate = 20,dims = npcs)
object <- ScoreJackStraw(object, dims = 1:npcs)

for(i in 0:14){
    a = i*10+1; b = (i+1)*10
    jpeg(paste0(path,"JackStrawPlot_",a,"_",b,".jpeg"), units="in", width=10, height=7,res=600)
    print(JackStrawPlot(object, dims = a:b))
    dev.off()
    Progress(i, 14)
}
p.values = object[["pca"]]@jackstraw@overall.p.values
print(npcs <- max(which(p.values[,"Score"] <=0.05)))

npcs <- 110


#======1.6 Performing SCTransform and integration =========================

format(object.size(object),unit = "GB")
options(future.globals.maxSize= object.size(object)*15)
object %<>% SCTransform(method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = TRUE)
object %<>% ScaleData(verbose = FALSE)
object %<>% RunPCA(verbose = T,npcs = npcs)

jpeg(paste0(path,"S1_ElbowPlot_SCT.jpeg"), units="in", width=10, height=7,res=600)
ElbowPlot(object, ndims = npcs)
dev.off()


object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)
object[["raw.umap"]] <- CreateDimReducObject(embeddings = object@reductions[["umap"]]@cell.embeddings,
                                             key = "rawUMAP_", assay = DefaultAssay(object))
colnames(object[["raw.umap"]]@cell.embeddings) %<>% paste0("raw-",.)

#======1.8 UMAP from harmony =========================

jpeg(paste0(path,"S1_RunHarmony.jpeg"), units="in", width=10, height=7,res=600)
system.time(object %<>% RunHarmony(group.by = "orig.ident", dims.use = 1:npcs,
                                     theta = 2, plot_convergence = TRUE,#do_pca = FALSE,
                                     nclust = 50, max.iter.cluster = 100))
dev.off()

object %<>% FindNeighbors(reduction = "harmony",dims = 1:npcs)
object %<>% FindClusters(resolution = 0.8)
object %<>% RunUMAP(reduction = "harmony", dims = 1:npcs)
object[['RNA']]@scale.data = matrix(0,0,0)
object[["SCT"]]@scale.data = matrix(0,0,0)
format(object.size(object),unit = "GB")

saveRDS(object, file = "data/MCL_94_20240614.rds")

#======1.5 Add Cell-cycle score =========================
# Read in a list of cell cycle markers, from Tirosh et al, 2015
cc.genes <- readLines(con = "../seurat_resources/regev_lab_cell_cycle_genes.txt")
s.genes <- FilterGenes(object,cc.genes[1:43])
g2m.genes <- FilterGenes(object,cc.genes[44:97])
object <- CellCycleScoring(object = object, s.features = s.genes, g2m.features = g2m.genes,
                           set.ident = FALSE)
colnames(object@meta.data) %<>% sub("Phase","cell cycle phase",.)
saveRDS(object@meta.data, file = "output/MCL_94_20240614.rds")

plots <- UMAPPlot(object, group.by="label1.blue_encode",raster=FALSE,cols = cols)
jpeg(paste0(path, "/", "umap.jpeg"), units="in", width=10, height=7,res=600)
print(plots)
dev.off()

#=======1.9 save SCT only =======================================
format(object.size(object),unit = "GB")

format(object.size(object@assays$RNA),unit = "GB")
object@assays[["RNA"]] = NULL
format(object.size(object),unit = "GB")

object[["SCT"]]@scale.data = matrix(0,0,0)
object[["SCT"]]@counts = matrix(0,0,0)

saveRDS(object, file = "data/MCL_94_SCT_20240615.rds")


