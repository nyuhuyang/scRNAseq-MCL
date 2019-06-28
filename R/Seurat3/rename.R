########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(magrittr)
library(gplots)
source("../R/Seurat3_functions.R")
source("R/util.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)


# 3.1.1 load data
# Rename ident
(load(file="data/MCL_V3_Harmony_43_20190610.Rda"))

df_samples <- readxl::read_excel("doc/190626_scRNAseq_info.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
sample_n = which(df_samples$tests %in% c("control",paste0("test",2:12)))
df_samples = df_samples[sample_n,]


object@meta.data$orig.ident %<>% plyr::mapvalues(from = unique(df_samples$sample),
                                               to = unique(df_samples$publication.id))
object@meta.data$groups %<>% plyr::mapvalues(from = unique(df_samples$group),
                                               to = unique(df_samples$group.id))
NewNames = paste0(object@meta.data$orig.ident,"_",object@meta.data$Barcode)
object %<>% RenameCells(new.names = NewNames)
rownames(object@reductions$tsne@cell.embeddings) = colnames(object)
gsub("_.*","",rownames(object@reductions$tsne@cell.embeddings)) %>% table

Idents(object) = "res.0.6"
table(Idents(object))
object %<>% sortIdent(numeric = T)

TSNEPlot.1(object)
save(object, file = "data/MCL_V3_Harmony_43_20190627.Rda")

