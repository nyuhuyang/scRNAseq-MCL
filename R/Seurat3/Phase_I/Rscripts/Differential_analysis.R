########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
#
# ######################################################################
invisible(lapply(c("Seurat","dplyr","cowplot","fgsea",
                   "magrittr","data.table","future","ggplot2","tidyr"), function(x) {
                       suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

# change the current plan to access parallelization

# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.integer(as.character(slurm_arrayid))
print(paste0("slurm_arrayid=",args))

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
# Need 64GB
# load files
(load(file = "data/MCL_41_harmony_20210426.Rda"))
Idents(object) = "Doublets"
object %<>% subset(idents = "Singlet")
DefaultAssay(object)  = "SCT"
idents.1 = sort(unique(object$clusters))

Idents(object) = "clusters"


cluster_markers = FindMarkers.UMI(object = object,ident.1 = idents.1[args], #0-30
                                  group.by = "clusters",
                                  logfc.threshold = 0.75,
                                  only.pos = F,
                                  test.use = "MAST",
                                  latent.vars = "nFeature_SCT")
id = sub("_.*",idents.1[args]) %>% as.integer()
id = idents.1[args]
if(id < 10) id =  paste0("0",idents.1[args])
write.csv(cluster_markers,file = paste0(path,id,"-markers_FC0_",idents.1[args],".csv"))
