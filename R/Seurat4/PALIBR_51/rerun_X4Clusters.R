invisible(lapply(c("Seurat","dplyr","ggplot2","cowplot","pbapply","sctransform","harmony","magrittr"), function(x) {
    suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)


#=========================
object <- readRDS(file = "data/MCL_51_20210724.rds")
meta.data <- readRDS(file = "output/MCL_51_20210724_meta.data_v2.rds")

if(all(colnames(object) == rownames(meta.data))){
    print("all cellID match!")
    object@meta.data = meta.data
}

object %<>% subset(subset =  discard == FALSE)
# CNV
CNV_meta = read.csv("shinyApp/PALIBR_I_51/cnv_meta_data.csv",row.names = 1)
if(all(colnames(object) == rownames(CNV_meta))){
    print("all cellID match!")
    object[["CNV_burden"]] = round(CNV_meta$cnv_score,digits = 6)
    lvl = sort(unique(object$`CNV_burden`))
    object$`CNV_burden` %<>% factor(levels = lvl)
}
object[["harmony.umap"]] <- CreateDimReducObject(embeddings = object@reductions[["umap"]]@cell.embeddings,
                                                    key = "harmonyUMAP_", assay = DefaultAssay(object))
object %<>% subset(subset =  UMAP_1 > 0
              & discard == FALSE
              & tSNE_2 < tSNE_1 +1
              & UMAP_2 <10
              #& orig.ident != "Pt2_30"
              & Doublets == "Singlet"
              #& X4cluster %in% c("1","2","3","4")
              & cell.types  %in% c("B_cells","MCL")
)
npcs <- 100
object %<>% FindNeighbors(reduction = "harmony",dims = 1:npcs,k.param = 20)
resolutions = c(seq(0.01,0.09, by = 0.01),seq(0.1,0.9, by = 0.1))
for(i in 1:length(resolutions)){
    object %<>% FindClusters(resolution = resolutions[i], algorithm = 1)
    UMAPPlot.1(object, group.by=paste0("SCT_snn_res.",resolutions[i]),
               #split.by = paste0("SCT_snn_res.",resolutions[i]),
               pt.size = 0.3,label = T,
               label.repel = T,alpha = 0.9,
               do.return = F,
               no.legend = T,label.size = 4, repel = T,
               title = paste("res =",resolutions[i]),
               do.print = T, save.path = path)
    Progress(i,length(resolutions))
}
saveRDS(object@meta.data, file = "output/MCL_B_51_20230113_meta.data.rds")
saveRDS(object@reductions[["umap"]], file = "output/MCL_B_51_20210724_umap.rds")

# Determine the ‘dimensionality’ of the dataset  =========
npcs <- 100

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

for(i in 0:9){
    a = i*10+1; b = (i+1)*10
    jpeg(paste0(path,"JackStrawPlot_",a,"_",b,".jpeg"), units="in", width=10, height=7,res=600)
    print(JackStrawPlot(object, dims = a:b))
    dev.off()
    Progress(i, 9)
}
p.values = object[["pca"]]@jackstraw@overall.p.values
print(npcs <- max(which(p.values[,"Score"] <=0.05)))
npcs <- 80
object[['RNA']]@scale.data = matrix(0,0,0)

#======= SCT ============================================

mito <- "^MT-"
message("mito.genes:")

(mito.features <- grep(pattern = mito, x = rownames(object), value = TRUE))
object[["percent.mt"]] <- PercentageFeatureSet(object = object, pattern = mito)
format(object.size(object),unit = "GB")
options(future.globals.maxSize= object.size(object)*50)
object[["SCT"]] <- NULL
object@commands$SCTransform.RNA <- NULL
object %<>% SCTransform(method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = TRUE)

#object <- FindVariableFeatures(object = object, selection.method = "vst",
#                               num.bin = 20, nfeatures = 3000,
#                               mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))
#object %<>% ScaleData(verbose = FALSE)
object %<>% RunPCA(verbose = T,npcs = npcs)

s.genes <- cc.genes$s.genes %>% gsub("MLF1IP","CENPU",.)
g2m.genes <- cc.genes$g2m.genes %>% plyr::mapvalues(from = c("FAM64A", "HN1"),
                                                    to = c("PIMREG","JPT1"))
object <- CellCycleScoring(object, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
colnames(object@meta.data) %<>% sub("Phase","cell.cycle.phase",.)

#================================================
object@reductions$umap <- NULL
object %<>% RunUMAP(reduction = "pca", dims = 1:npcs,return.model = TRUE)
object[["raw.umap.v3"]] <- CreateDimReducObject(embeddings = object@reductions[["umap"]]@cell.embeddings,
                                             key = "rawUMAP_", assay = DefaultAssay(object))
object@reductions$harmony <- NULL
jpeg(paste0(path,"S1_RunHarmony.jpeg"), units="in", width=10, height=7,res=600)
system.time(object %<>% RunHarmony.1(group.by = "orig.ident", dims.use = 1:npcs,
                                     theta = 2, plot_convergence = TRUE,
                                     nclust = 50, max.iter.cluster = 100))
dev.off()
object@reductions$umap <- NULL
object %<>% RunUMAP(reduction = "harmony", dims = 1:npcs,return.model = TRUE)
object %<>% FindNeighbors(reduction = "harmony",dims = 1:npcs,k.param = 20)
resolutions = c(seq(0.01,0.09, by = 0.01),seq(0.1,0.9, by = 0.1))
for(i in 1:length(resolutions)){
    object %<>% FindClusters(resolution = resolutions[i], algorithm = 1)
    UMAPPlot.1(object, group.by=paste0("SCT_snn_res.",resolutions[i]),
               #split.by = paste0("SCT_snn_res.",resolutions[i]),
               pt.size = 0.3,label = T,
               label.repel = T,alpha = 0.9,
               do.return = F,
               no.legend = T,label.size = 4, repel = T,
               title = paste("res =",resolutions[i]),
               do.print = T, save.path = path)
    Progress(i,length(resolutions))
}

UMAPPlot.1(object, group.by="patient",#paste0("SCT_snn_res.",resolutions[i]),
           #split.by = paste0("SCT_snn_res.",resolutions[i]),
           pt.size = 0.3,label = T,
           label.repel = T,alpha = 0.9,
           do.return = F,
           no.legend = T,label.size = 4, repel = T,
           title = paste("res =",resolutions[i]),
           do.print = T, save.path = path)
object$reRun_SCT_snn_res.0.05 <- object$SCT_snn_res.0.05
cl <- grep("^SCT_snn_res.",colnames(object@meta.data))
object@meta.data <- object@meta.data[,-cl]

object[['RNA']]@scale.data = matrix(0,0,0)
object[['SCT']]@scale.data = matrix(0,0,0)
format(object.size(object),unit = "GB")

saveRDS(object, file = "data/MCL_B_51_20230113.rds")



MCL$X4cluster = as.integer(MCL$SCT_snn_res.0.006)
MCL = subset(MCL, subset =  X4cluster != 5)

MCL$X4cluster %<>% gsub("6|7|8|9|10","1",.)
MCL$SCT_snn_res.0.006 %<>% as.integer() %>% as.factor()

object$X4cluster = object$cell.types
object@meta.data[rownames(MCL@meta.data),"X4cluster"] = MCL$X4cluster


