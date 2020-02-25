########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(cowplot)
library(magrittr)
library(DoubletFinder)
library(kableExtra)
source("../R/Seurat3_functions.R")
source("R/util.R")

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
i <- as.numeric(slurm_arrayid)
print(paste0("slurm_arrayid=",i))

########################################################################
#
#  2. DoubletFinder 
# 
# ######################################################################

# samples

(load(file = "data/MCL_41_harmony_20200225.Rda"))
(samples = unique(object$orig.ident))
object_list <- SplitObject(object,split.by = "orig.ident")
rm(object);GC()
(load(file = "output/MCL_41_harmony_20200203_sweep.res_list.Rda"))
sweep_list <- lapply(sweep.res_list, function(x) summarizeSweep(x, GT = FALSE))
bcmvn_list <- lapply(sweep_list,find.pK)


(maximal_pk <- sapply(bcmvn_list,function(x) {
    as.numeric(as.character(x[find.localMaxima(x$BCmetric),"pK"]))
    }))
maximal_pk

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
print(paste("processing",unique(object_list[[i]]$orig.ident)))
homotypic.prop <- modelHomotypic(object_list[[i]]@meta.data$cell.types)
nExp_poi <- round(Multiplet_Rate(object_list[[i]])*length(colnames(object_list[[i]])))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
object_list[[i]] <- doubletFinder_v3(object_list[[i]], PCs = 1:50, 
                                     pN = 0.25, pK = maximal_pk[i], 
                                     nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
object_list[[i]] <- doubletFinder_v3(object_list[[i]], PCs = 1:50, 
                                     pN = 0.25, pK = maximal_pk[i],
                                     nExp = nExp_poi.adj, 
                                     reuse.pANN = grep("pANN",colnames(object_list[[i]]@meta.data),value = T),
                                     sct = TRUE)
colName = colnames(object_list[[i]]@meta.data)
colName[grep("DF.classifications",colName)] = c("Low_confident_doublets",
                                                "High_confident_doublets")
colnames(object_list[[i]]@meta.data) = colName

saveRDS(object_list[[i]]@meta.data,file=paste0(path,"object_list_",i,"_meta.data.rds"))
