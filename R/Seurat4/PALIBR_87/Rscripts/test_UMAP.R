# r4.1.1
invisible(lapply(c("Seurat","dplyr","ggplot2","cowplot","pbapply","magrittr","sctransform","harmony"), function(x) {
    suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")

# Need 64GB
set.seed(101)
# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.integer(as.character(slurm_arrayid))
print(paste0("slurm_arrayid=",args))

#=========
opts = c("test_npcs_dist_spread","densmap")[1]
if(opts == "test_npcs_dist_spread"){
    test_df = data.frame(min_dist = rep(c(0.2,0.4,0.6),each = 3),
                         spread = rep(c(0.6,1.0,1.4),times = 3),
                         npcs = rep(c(50,60,70,80,90,100),each = 9))
    test_df = bind_rows(list(test_df,test_df))
    test_df$nfeatures = rep(c(2000,3000),each =54)

    print(nfeatures <- test_df[args,"nfeatures"])
    print(spread <- test_df[args,"spread"])
    print(min.dist <- test_df[args,"min_dist"])
    print(npcs <- test_df[args,"npcs"])

    file.name = paste0("npcs",npcs,"_dist.",min.dist,"_spread.",spread)

    save.path <- paste0("output/",gsub("-","",Sys.Date()),"/",nfeatures)
    if(!dir.exists(save.path)) dir.create(save.path, recursive = T)


    object <- readRDS(file = "data/MCL_SCT_87_20220901.rds")
    DefaultAssay(object) = "SCT"

    object <- FindVariableFeatures(object = object, selection.method = "vst",
                                   num.bin = 20, nfeatures = nfeatures,
                                   mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))
    print(length(VariableFeatures(object)))
    object %<>% ScaleData(verbose = FALSE)
    object %<>% RunPCA(verbose = T,npcs = npcs)
    system.time(object %<>% RunHarmony.1(group.by = "orig.ident", dims.use = 1:npcs,
                                         theta = 2, plot_convergence = TRUE,
                                         nclust = 50, max.iter.cluster = 100))
    object %<>% RunUMAP(reduction = "harmony", dims = 1:npcs,min.dist = min.dist,spread = spread,
                        return.model = TRUE)

    saveRDS(object@reductions, file = paste0(save.path, "/reductions_",file.name,".rds"))
}



if(opts == "densmap"){
    test_df = data.frame(min_dist = rep(c(0.2,0.4,0.6),each = 3),
                         spread = rep(c(0.6,1.0,1.4),times = 3),
                         dens_lambda = rep(c(1,2,5),each = 9),
                         npcs = rep(c(70,80),each = 27))

    print(min.dist <- test_df[args,"min_dist"])
    print(spread <- test_df[args,"spread"])
    print(dens_lambda <- test_df[args,"dens_lambda"])
    print(npcs <- test_df[args,"npcs"])

    file.name = paste0("npcs",npcs,"_dist.",min.dist,"_spread.",spread,"_dens_lambda.",dens_lambda)

    save.path <- paste0("output/",gsub("-","",Sys.Date()),"/")
    if(!dir.exists(save.path)) dir.create(save.path, recursive = T)


    object = readRDS(file = "data/MCL_AIM_93_20220519.rds")
    object@meta.data = readRDS(file = "output/MCL_AIM_93_20220519_metadata_v2.rds")
    DefaultAssay(object) = "SCT"
    object@reductions = list()

    object <- FindVariableFeatures(object = object, selection.method = "vst",
                                   num.bin = 20, nfeatures = 3000,
                                   mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))
    print(length(VariableFeatures(object)))
    object %<>% ScaleData(verbose = FALSE)
    object %<>% RunPCA(verbose = T,npcs = npcs)
    system.time(object %<>% RunHarmony.1(group.by = "orig.ident", dims.use = 1:npcs,
                                         theta = 2, plot_convergence = TRUE,
                                         nclust = 50, max.iter.cluster = 100))
    object %<>% RunUMAP(reduction = "harmony", dims = 1:npcs,min.dist = min.dist,spread = spread,
                        return.model = TRUE,
                        umap.method = "umap-learn",densmap = TRUE,dens.lambda = dens_lambda, dens.frac = 0.3,dens_var_shift=0.1)

    object[[paste0("umap_",file.name)]] <- CreateDimReducObject(embeddings = object@reductions[["umap"]]@cell.embeddings,
                                                                key = paste0(args,"UMAP_"), assay = DefaultAssay(object))
    saveRDS(object@reductions, file = paste0(save.path, "umap_",file.name,".rds"))
}
