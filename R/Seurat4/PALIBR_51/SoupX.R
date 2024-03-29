library(SoupX)
library(magrittr)
library(Seurat) # Seurat 4
library(sctransform)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
#https://cellgeni.github.io/notebooks/html/new-10kPBMC-SoupX.html
#======1.1 Load the data files and Set up Seurat object =========================
# read sample summary list
df_samples <- readxl::read_excel("doc/20210715_scRNAseq_info.xlsx",sheet = "fastq")
df_samples = as.data.frame(df_samples)
colnames(df_samples) %<>% tolower()
df_samples %<>% filter(sequence %in% "GEX") %>% filter(phase %in% "PALIBR_I") %>%
    filter(sample != "Pt11_31")
# check missing data
read.path = "data/scRNAseq/counts"
current <- list.files(read.path)
(current <- current[!grepl(".Rda|RData",current)])
(missing_data <- df_samples$sample.id[!(df_samples$sample.id %in% current)])


#============== filtered counts & re-clustering ====================
adj.matrix_list <- pbapply::pblapply(df_samples$sample.id, function(s){
    readDir <- file.path(read.path,as.character(s),"outs")
    filt.matrix <- Seurat::Read10X(file.path(readDir, "filtered_feature_bc_matrix"),strip.suffix = TRUE)
    raw.matrix <- Seurat::Read10X(file.path(readDir, "raw_feature_bc_matrix"),strip.suffix = TRUE)
    soup.channel  <- SoupChannel(raw.matrix, filt.matrix)

    obj <- CreateSeuratObject(filt.matrix,min.cells = 0,names.delim = "-",min.features = 0) %>%
        SCTransform(method = "glmGamPoi",verbose = FALSE) %>%
        #FindVariableFeatures(verbose = FALSE) %>%
        RunPCA(verbose = FALSE) %>%
        FindNeighbors(reduction = "pca",dims = 1:30) %>%
        FindClusters(resolution = 0.8, algorithm = 1,verbose = F)


    soup.channel  <- setClusters(soup.channel, setNames(obj$SCT_snn_res.0.8, colnames(obj)))
    soup.channel  <- autoEstCont(soup.channel, priorRhoStdDev = 0.3)
    adj.matrix  <- adjustCounts(soup.channel, roundToInt = T)
    return(adj.matrix)
})

names(adj.matrix_list) = df_samples$sample

for (s in df_samples$sample) {
    colnames(adj.matrix_list[[s]]) = paste0(s,"-", colnames(adj.matrix_list[[s]]))
}

adj.matrix <- do.call(cbind, adj.matrix_list)
meta.data <- readRDS(file = "output/MCL_51_20210724_meta.data_v1.rds")
table(colnames(adj.matrix) %in% rownames(meta.data))
table(rownames(meta.data) %in% colnames(adj.matrix))

object <- CreateSeuratObject(adj.matrix[,rownames(meta.data)],
                             min.cells = 0,names.delim = ".",min.features = 0)
mito <- "^MT-"
message("mito.genes:")

(mito.features <- grep(pattern = mito, x = rownames(object), value = TRUE))
object[["percent.mt"]] <- PercentageFeatureSet(object = object, pattern = mito)

format(object.size(object)*60,unit = "GB")
options(future.globals.maxSize= object.size(object)*60)
object %<>% SCTransform(method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = TRUE)
object[["SCT"]]@scale.data = matrix(0,0,0)

object_orig = readRDS(file = "data/MCL_SCT_51_20210724.rds")
table(colnames(object) == colnames(object_orig))
object@graphs = object_orig@graphs
object@reductions = object_orig@reductions
object@tools = object_orig@tools
object_orig@commands$SCTransform.RNA = object@commands
object@commands = object_orig@commands
object$orig.ident <- gsub("-.*","",colnames(object))
meta.data <- meta.data[,-which(colnames(meta.data) %in% colnames(object@meta.data))]

if(all(colnames(object) == rownames(meta.data))){
    print("all cellID match!")
    object@meta.data %<>% cbind(meta.data)
}

format(object.size(object),unit = "GB")

s.genes <- cc.genes$s.genes %>% gsub("MLF1IP","CENPU",.)
g2m.genes <- cc.genes$g2m.genes %>% plyr::mapvalues(from = c("FAM64A", "HN1"),
                                                    to = c("PIMREG","JPT1"))
object %<>% CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

saveRDS(object, file = "data/MCL_51_20210724.rds")
saveRDS(object@meta.data, file = "output/MCL_51_20210724_meta.data_v2.rds")
