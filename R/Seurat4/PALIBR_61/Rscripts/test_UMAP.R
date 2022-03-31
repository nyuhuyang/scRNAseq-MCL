# test
invisible(lapply(c("Seurat","dplyr","ggplot2","cowplot","pbapply","sctransform","harmony"), function(x) {
    suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
save.path <- paste0("output/",gsub("-","",Sys.Date()))
if(!dir.exists(save.path)) dir.create(save.path, recursive = T)
# Need 64GB ?
set.seed(101)
# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.integer(as.character(slurm_arrayid))
print(paste0("slurm_arrayid=",args))

#=========
test_df = data.frame(min_dist = rep(c(0.2,0.4,0.6),each = 3),
                     spread = rep(c(0.6,1.0,1.4),times = 3),
                     npcs = rep(c(50,60,70,80,90,100),each = 9))

print(spread <- test_df[args,"spread"])
print(min.dist <- test_df[args,"min_dist"])
print(npcs <- test_df[args,"npcs"])

file.name = paste0("cs",npcs,"_dist.",min.dist,"_spread.",spread)

object = readRDS(file = "data/MCL_61_20220318.rds")

#object %<>% SCTransform(method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = TRUE)
DefaultAssay(object) = "SCT"
print(length(VariableFeatures(object)))
#object %<>% ScaleData(verbose = FALSE)
object %<>% RunPCA(verbose = T,npcs = npcs)
jpeg(paste0(save.path,"S1_RunHarmony.jpeg"), units="in", width=10, height=7,res=600)
system.time(object %<>% RunHarmony.1(group.by = "orig.ident", dims.use = 1:npcs,
                                     theta = 2, plot_convergence = TRUE,
                                     nclust = 50, max.iter.cluster = 100))
dev.off()
object %<>% RunUMAP(reduction = "harmony", dims = 1:npcs,return.model = TRUE)

UMAPPlot.1(object, group.by = "orig.ident", do.print = T,title = file.name,file.name = paste0(file.name,".jpeg"))
object %<>% FindNeighbors(reduction = "umap",dims = 1:2)
object[[paste0("umap_",file.name)]] <- CreateDimReducObject(embeddings = object@reductions[["umap"]]@cell.embeddings,
                                                            key = paste0(args,"UMAP_"), assay = DefaultAssay(object))
umap = object@reductions[paste0("umap_",file.name)]
saveRDS(umap, file = paste0(save.path, "/umap_",file.name,".rds"))
