########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(SingleR)
library(kableExtra)
source("../R/Seurat_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
########################################################################
#
#  1 Seurat Alignment 
# 
# ######################################################################
#======1.1 Setup the Seurat objects =========================
# Load the mouse.eyes dataset

# setup Seurat objects since both count matrices have already filtered
# cells, we do no additional filtering here
df_samples <- readxl::read_excel("doc/181002_Single_cell_sample list.xlsx")
sample_n = which(df_samples$tests %in% c("test2", "test3","test4"))
table(df_samples$tests)
df_samples[sample_n,]
samples <- df_samples$samples[sample_n]
projects <- df_samples$projects[sample_n]
conditions <- df_samples$conditions[sample_n]

MCL_raw <- list()
MCL_Seurat <- list()
for(i in 1:length(samples)){
    MCL_raw[[i]] <- Read10X(data.dir = paste0("./data/",
                                             samples[i],"/outs/filtered_gene_bc_matrices/hg19/"))
    colnames(MCL_raw[[i]]) <- paste0(samples[i],
                                            "_",colnames(MCL_raw[[i]]))
    MCL_Seurat[[i]] <- CreateSeuratObject(MCL_raw[[i]],
                                                 min.cells = 3,
                                                 min.genes = 200,
                                                 project = projects[i],
                                                 names.delim = "_")
    MCL_Seurat[[i]]@meta.data$conditions <- conditions[i]
}
#======1.1.2 QC before merge =========================
cell.number <- sapply(MCL_Seurat, function(x) length(x@cell.names))
QC_list <- lapply(MCL_Seurat, function(x) as.matrix(x = x@raw.data))
median.nUMI <- sapply(QC_list, function(x) median(colSums(x)))
median.nGene <- sapply(QC_list, function(x) median(apply(x,2,function(y) sum(length(y[y>0])))))

min.nUMI <- sapply(QC_list, function(x) min(colSums(x)))
min.nGene <- sapply(QC_list, function(x) min(apply(x,2,function(y) sum(length(y[y>0])))))

QC.list <- cbind(df_samples,cell.number, median.nUMI,median.nGene,min.nUMI,min.nGene,
                 row.names = samples)
write.csv(QC.list,paste0(path,"QC_list.csv"))
QC.list %>% kable() %>% kable_styling()
remove(QC_list);GC()

#========1.1.3 merge ===================================
MCL <- Reduce(function(x, y) MergeSeurat(x, y, do.normalize = F), MCL_Seurat)
remove(MCL_raw,MCL_Seurat);GC()
#MCL <- FilterCells(MCL, subset.names = "nGene",low.thresholds = 200,high.thresholds = Inf) %>%
#    NormalizeData() #%>%
    #ScaleData(display.progress = FALSE) %>%
    #FindVariableGenes(do.plot = FALSE, display.progress = FALSE)
mito.genes <- grep(pattern = "^MT-", x = rownames(x = MCL@data), value = TRUE)
percent.mito <- Matrix::colSums(MCL@raw.data[mito.genes, ])/Matrix::colSums(MCL@raw.data)
MCL <- AddMetaData(object = MCL, metadata = percent.mito, col.name = "percent.mito")

MCL@ident = factor(MCL@ident,levels = samples)

g1 <- VlnPlot(object = MCL, features.plot = c("nGene", "nUMI", "percent.mito"), 
              nCol = 1,point.size.use = 0.2,
              x.lab.rot = T, do.return = T,return.plotlist =T)

save(g1, file = "./data/MCL_g1.Rda")
