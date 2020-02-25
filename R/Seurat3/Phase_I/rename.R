########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
library(Seurat)
library(dplyr)
library(kableExtra)
library(magrittr)
library(ggplot2)
library(cowplot)
library(fgsea)
library(tibble)
library(ggsci)
source("../R/Seurat3_functions.R")

# load data
choose = c("All_samples","B")[2]
if(choose == "All_samples") (load(file="data/MCL_41_harmony_20191231.Rda"))
if(choose == "B") object = readRDS(file = "data/MCL_41_B_20200207.rds")

head(object@meta.data)
object$orig.ident %<>% gsub("Pt28_25", "Pt28_28",.)
object$sample %<>% gsub("Pt-28-PB-C25D1", "Pt-28-PB-C28D1",.)

Barcode = colnames(object)
Barcode %<>% gsub(".*_","",.)
NewNames = paste0(object$sample,"_",Barcode)
object %<>% RenameCells(new.names = NewNames)
rownames(object@reductions$pca@cell.embeddings) = colnames(object)
rownames(object@reductions$tsne@cell.embeddings) = colnames(object)
rownames(object@reductions$umap@cell.embeddings) = colnames(object)
rownames(object@reductions$harmony@cell.embeddings) = colnames(object)

df_samples <- readxl::read_excel("doc/191120_scRNAseq_info.xlsx")
colnames(df_samples) %<>% tolower
df_samples$`sample name`

(df_samples = df_samples[df_samples$`sample name` %in% object$orig.ident,])
(samples = df_samples$`sample name`[df_samples$`sample name` %in% object$orig.ident])
object$orig.ident %<>% factor(levels = samples)
if(choose == "All_samples") save(object, file = "data/MCL_41_harmony_20191231.Rda")
if(choose == "B") saveRDS(object, file = "data/MCL_41_B_20200207.rds")
