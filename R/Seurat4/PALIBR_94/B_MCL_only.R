invisible(lapply(c("Seurat","dplyr","ggplot2","cowplot","pbapply","sctransform","harmony"), function(x) {
    suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

object <- readRDS(file = "data/MCL_94_20240615.rds")
meta.data <- readRDS(file = "output/MCL_SCT_94_20240615_meta.data_v1.rds")
if(all(colnames(object) == rownames(meta.data))){
    print("all cellID match!")
    object@meta.data = meta.data
}

object@meta.data %<>% cbind(object[["umap"]]@cell.embeddings)
seurat_obj  <- subset(object, UMAP_2 >0  & cell.types %in% c("B cells","MCL"))

seurat_obj@meta.data[,"SCT_snn_res.0.8"] = NULL

seurat_obj %<>% SCTransform(method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = TRUE)
seurat_obj %<>% ScaleData(verbose = FALSE)
seurat_obj %<>% RunPCA(verbose = T,npcs = 50)

system.time(seurat_obj %<>% RunHarmony(group.by = "orig.ident", dims.use = 1:50,
                                   theta = 2, plot_convergence = TRUE,#do_pca = FALSE,
                                   nclust = 50, max.iter.cluster = 100))
seurat_obj %<>% FindNeighbors(reduction = "umap",dims = 1:2)

resolutions = c(seq(0.02,0.08, by = 0.02),seq(0.2,0.8, by = 0.2))
for(i in 1:length(resolutions)){
    seurat_obj %<>% FindClusters(resolution = resolutions[i], algorithm = 1,verbose = F)
    Progress(i, length(resolutions))
}

library(clustree)
library(pheatmap)

# Calculate average expression for each cluster
for (res in paste0("SCT_snn_res.",resolutions)){
    Idents(seurat_obj) = res
    cluster_averages <- AverageExpression(seurat_obj, return.seurat = TRUE)

    # Compute pairwise correlations between clusters
    cluster_correlation <- cor(cluster_averages@assays$RNA@data)

    # Visualize the correlation matrix
    jpeg(paste0(path,"cluster_correlation_",res,".jpeg"), units="in", width=7, height=7,res=600)
    pheatmap(cluster_correlation)
    dev.off()
}
for (res in paste0("SCT_snn_res.",resolutions)){
    Idents(seurat_obj) = res
    jpeg(paste0(path,"UMAPPlot_",res,".jpeg"), units="in", width=7, height=7,res=600)
    print(UMAPPlot(seurat_obj))
    dev.off()
    }
seurat_obj$X7cluster <- "1"
seurat_obj$X7cluster[seurat_obj$SCT_snn_res.0.04 %in% c("5","13")] = "2"
seurat_obj$X7cluster[seurat_obj$SCT_snn_res.0.04 %in% c("3")] = "3"
seurat_obj$X7cluster[seurat_obj$SCT_snn_res.0.04 %in% "6"] = "4"
seurat_obj$X7cluster[seurat_obj$SCT_snn_res.0.04 %in% "8"] = "5"
seurat_obj$X7cluster[seurat_obj$SCT_snn_res.0.04 %in% "9"] = "6"
seurat_obj$X7cluster[seurat_obj$SCT_snn_res.0.04 %in% "7"] = "7"
table(seurat_obj$X7cluster)
seurat_obj <- readRDS(file = "data/MCL_B_only_94_SCT_20240615.rds")
saveRDS(seurat_obj, file = "data/MCL_B_only_94_SCT_20240615.rds")
saveRDS(seurat_obj@meta.data, file = "output/MCL_B_only_94_SCT_20240615_meta.data_v1.rds")
