########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
####################################
invisible(lapply(c("Seurat","dplyr","magrittr","tidyr",
                   "MAST","future","gplots"), function(x) {
        suppressPackageStartupMessages(library(x,character.only = T))
}))
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

# change the current plan to access parallelization

# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
i <- as.numeric(slurm_arrayid)
print(paste0("slurm_arrayid=",i))

step = 3 
# choose == "X4clusters"
if(step == 1){ # need 32 GB
        args = data.frame(only.pos = rep(c(T,  T,   T,   F),  each = 4),
                          logfc =  rep(c(0.25, 0.1, 0.05, 0), each = 4),
                          ident.1 = rep(paste0("C",1:4),      time = 4))
        
        (arg = args[i,])
        object = readRDS(file = "data/MCL_41_B_20200225.rds")
        DefaultAssay(object)  = "SCT"
        Idents(object) = "orig.ident"
        object %<>% subset(idents = "Pt2_30Pd",invert = T)
        Idents(object) = "cell.types"
        object <- subset(object, idents= "MCL") 
        Idents(object) = "X4clusters"
        system.time(MCL_markers <- FindMarkers.UMI(object, 
                                                   ident.1 = as.character(arg$ident.1),
                                                   ident.2 = NULL,
                                                   logfc.threshold = arg$logfc, 
                                                   only.pos = arg$only.pos,
                                                   test.use = "MAST",
                                                   latent.vars = "nCount_SCT"))
        write.csv(MCL_markers,paste0(path,"MCL_only_41-FC",arg$logfc,"_",arg$ident.1,".csv"))
}
# choose == "X4clusters_vs_Normal"
if(step == 2){ # need 32 GB
        args = data.frame(only.pos = rep(c(T,  T,   T,   F),  each = 5),
                          logfc =  rep(c(0.25, 0.1, 0.05, 0), each = 5),
                          ident.1 = rep(c("B_cells",paste0("C",1:4)),      time = 4))
        
        (arg = args[i,])
        object = readRDS(file = "data/MCL_41_B_20200225.rds")
        DefaultAssay(object)  = "SCT"
        Idents(object) = "orig.ident"
        object %<>% subset(idents = "Pt2_30Pd",invert = T)
        object$X4clusters_normal = as.character(object$X4clusters)
        object$X4clusters_normal %<>% paste(object$cell.types, sep = "_")
        object$X4clusters_normal %<>% gsub(".*_B_cells","B_cells",.)
        object$X4clusters_normal %<>% gsub("_MCL","",.)
        normal <- grepl("N01|N02|N03",object$orig.ident)
        
        object@meta.data[normal,"X4clusters_normal"] = "Normal"
        Idents(object) = "X4clusters_normal"
        object %<>% sortIdent()
        table(Idents(object))
        system.time(MCL_markers <- FindMarkers.UMI(object, 
                                                   ident.1 = as.character(arg$ident.1),
                                                   ident.2 = "Normal",
                                                   logfc.threshold = arg$logfc, 
                                                   only.pos = arg$only.pos,
                                                   test.use = "MAST",
                                                   latent.vars = "nCount_SCT"))
        write.csv(MCL_markers,paste0(path,"MCL_Normal_41-FC",arg$logfc,"_",arg$ident.1,".csv"))
}
# choose == "X4clusters_vs_B_cells"
if(step == 3){ # need 32 GB
        args = data.frame(only.pos = rep(c(T,  T,   T,   F),  each = 4),
                          logfc =  rep(c(0.25, 0.1, 0.05, 0), each = 4),
                          ident.1 = rep(paste0("C",1:4),      time = 4))
        
        (arg = args[i,])
        object = readRDS(file = "data/MCL_41_B_20200225.rds")
        DefaultAssay(object)  = "SCT"
        Idents(object) = "orig.ident"
        object %<>% subset(idents = c("Pt2_30Pd","N01","N02","N03"),invert = T)
        object$X4clusters_B = as.character(object$X4clusters)
        object$X4clusters_B %<>% paste(object$cell.types, sep = "_")
        object$X4clusters_B %<>% gsub(".*_B_cells","B_cells",.)
        object$X4clusters_B %<>% gsub("_MCL","",.)
        Idents(object) = "X4clusters_B"
        object %<>% sortIdent()
        table(Idents(object))
        system.time(MCL_markers <- FindMarkers.UMI(object, 
                                                   ident.1 = as.character(arg$ident.1),
                                                   ident.2 = "B_cells",
                                                   logfc.threshold = arg$logfc, 
                                                   only.pos = arg$only.pos,
                                                   test.use = "MAST",
                                                   latent.vars = "nCount_SCT"))
        write.csv(MCL_markers,paste0(path,"MCL_B_41-FC",arg$logfc,"_",arg$ident.1,".csv"))
}