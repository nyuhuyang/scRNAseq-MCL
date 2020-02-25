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
args <- as.numeric(slurm_arrayid)
print(paste0("slurm_arrayid=",args))
logfc = list(c(T,0.25),
             c(T,0.1),
             c(T,0.05),
             c(F,0))
(r <- logfc[[args]])
step = 1
if(step == 1){ # need 128 GB
        object = readRDS(file = "data/MCL_41_B_20200207.rds")
        DefaultAssay(object)  = "SCT"
        Idents(object) = "orig.ident"
        object %<>% subset(idents = "Pt2_30Pd",invert = T)
        Idents(object) = "X4clusters"
        system.time(MCL_markers <- FindAllMarkers.UMI(object, 
                                                       logfc.threshold = r[2], 
                                                      only.pos = ifelse(r[1],T,F),
                                                      test.use = "MAST",
                                                      return.thresh = 1))
        write.csv(MCL_markers,paste0(path,"MCL_41-FC",r[2],".csv"))
}
