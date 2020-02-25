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

args = data.frame(only.pos = rep(c(T,  T,   T,   F),  each = 4),
                  logfc =  rep(c(0.25, 0.1, 0.05, 0), each = 4),
                  ident.1 = rep(paste0("C",1:4),      time = 4))

(arg = args[i,])
step = 1 # choose == "X4clusters"
if(step == 1){ # need 128 GB
        object = readRDS(file = "data/MCL_41_B_20200225.rds")
        DefaultAssay(object)  = "SCT"
        Idents(object) = "orig.ident"
        object %<>% subset(idents = "Pt2_30Pd",invert = T)
        Idents(object) = "X4clusters"
        system.time(MCL_markers <- FindMarkers.UMI(object, ident.1 = as.character(arg$ident.1),
                                                   ident.2 = NULL,
                                                   logfc.threshold = arg$logfc, 
                                                   only.pos = arg$only.pos,
                                                   test.use = "MAST",
                                                   return.thresh = 1))
        write.csv(MCL_markers,paste0(path,"MCL_41-FC",arg$logfc,"_",arg$ident.1,".csv"))
}
