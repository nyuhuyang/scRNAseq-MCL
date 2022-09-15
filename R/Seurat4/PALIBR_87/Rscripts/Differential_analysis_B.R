########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
#
# ######################################################################
####################################
invisible(lapply(c("Seurat","dplyr","magrittr","tidyr",
                   "future","gplots"), function(x) {
        suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_differential_expression.R")
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)


# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.integer(as.character(slurm_arrayid))
print(paste0("slurm_arrayid=",args))


object = readRDS("data/MCL_AIM_93_20220519.rds")
object@meta.data = readRDS("output/MCL_AIM_93_20220519_metadata_v2.rds")

#object %<>% subset(subset = Doublets == "Singlet")
object[["umap"]] = NULL
file.name = paste0("output/20220526/3000/umap_npcs70_dist.0.4_spread.0.6.rds")
umap  = readRDS(file.name)[[1]]
umap@key = "UMAP_"
colnames(umap@cell.embeddings) = c("UMAP_1","UMAP_2")
object[["umap"]] <- umap

DefaultAssay(object) = "SCT"
object@meta.data %<>% cbind(umap@cell.embeddings)


step = c("SCT_snn_res.0.8")[1]
# choose == "orig.ident_X6_vs_orig.ident_X6"
if(step == "SCT_snn_res.0.8"){ # need 32 GB
        object = subset(object, subset =  UMAP_1 > -2.5 & UMAP_2 < 1.5 & 
                                label1.blue_encode %in% c("B cells","MCL")
)
        object$SCT_snn_res.0.8 %<>% droplevels()
        df <- table(object$SCT_snn_res.0.8) %>% as.data.frame.table() %>%
                filter(Freq >= 10) %>% arrange(desc(Freq))
        
        print(opt <- df[args,])

        markers = FindMarkers_UMI(object,
                                  ident.1 = opt$Var1,
                                  group.by = "SCT_snn_res.0.8",
                                  assay = "SCT",
                                  min.pct = 0.1,
                                  logfc.threshold = 0.5,
                                  only.pos = T)
        markers$gene = rownames(markers)
        markers$cluster = opt$Var1

        arg = args
        if(args < 10) arg = paste0("0",arg)
        if(args < 100) arg = paste0("0",arg)

        write.csv(markers,paste0(path,arg,"SCT_snn_res.0.8","_Cluster",opt$Var1, ".csv"))
}
