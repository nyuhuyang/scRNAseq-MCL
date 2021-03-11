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
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
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

# need 16GB

(load(file = "data/MCL_41_harmony_20200225.Rda"))
(samples = unique(object$orig.ident))
Idents(object) = "orig.ident"
sub_object <- subset(object, idents = samples[i])
rm(object);GC()
(load(file = "output/MCL_41_harmony_20200203_sweep.res_list.Rda"))
sweep_list <- lapply(sweep.res_list, function(x) summarizeSweep(x, GT = FALSE))
bcmvn_list <- lapply(sweep_list,find.pK)


(maximal_pk <- sapply(bcmvn_list,function(x) {
    as.numeric(as.character(x[find.localMaxima(x$BCmetric),"pK"]))
    }))

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
print(paste("processing",unique(sub_object$orig.ident)))
homotypic.prop <- modelHomotypic(sub_object@meta.data$cell.types)
nExp_poi <- round(Multiplet_Rate(sub_object)*length(colnames(sub_object)))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
sub_object <- doubletFinder_v3(sub_object, PCs = 1:50, 
                                     pN = 0.25, pK = maximal_pk[i], 
                                     nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
sub_object <- doubletFinder_v3(sub_object, PCs = 1:50, 
                                     pN = 0.25, pK = maximal_pk[i],
                                     nExp = nExp_poi.adj, 
                                     reuse.pANN = grep("pANN",colnames(sub_object@meta.data),value = T),
                                     sct = TRUE)
colName = colnames(sub_object@meta.data)
colName[grep("DF.classifications",colName)] = c("Low_confident_doublets",
                                                "High_confident_doublets")
colnames(sub_object@meta.data) = colName

saveRDS(sub_object@meta.data,file=paste0(path,"object_list_",i,"_meta.data.rds"))
