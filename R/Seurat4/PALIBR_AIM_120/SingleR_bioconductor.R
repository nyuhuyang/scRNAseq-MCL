#====== 3.1 Create Singler Object  ==========================================
# conda activate r4.1.1 linux
invisible(lapply(c("Seurat","SingleR","SingleCellExperiment","magrittr","data.table","Matrix"), function(x) {
    suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/SingleR_functions.R")

# Need 64GB ?
set.seed(101)
# ====== load reference =============
reference = c("MCL+blue_encode","MCL+azimuth_PBMC")[2]

# ====== load single cell =============
object = readRDS(file = "data/MCL_120_SCT_20230106.rds")
sce <- SingleCellExperiment(list(logcounts=object[["SCT"]]@data),
                            colData=DataFrame(object@meta.data))
rm(object);GC()

# 2. check and prepare MCL data==============================
#remove duplicate rownames with lower rowsumns
#' @param mat input as data.frame with gene name
#' @export mat matrix with gene as rownames, no duplicated genes
RemoveDup <- function(mat){
    gene_id <- as.matrix(mat[,1])
    mat <- mat[,-1]
    if(!is.matrix(mat)) mat <- sapply(mat,as.numeric)
    rownames(mat) <- 1:nrow(mat)
    mat[is.na(mat)] = 0
    mat <- cbind(mat, "rowSums" = rowSums(mat))
    mat <- mat[order(mat[,"rowSums"],decreasing = T),]
    gene_id <- gene_id[as.numeric(rownames(mat))]
    remove_index <- duplicated(gene_id)
    mat <- mat[!remove_index,]
    rownames(mat) <- gene_id[!remove_index]
    return(mat[,-ncol(mat)])
}
if(reference == "MCL+blue_encode"){
    #MCL_bulk
    MCL_bulk <- read.csv(file="../scRNAseq-MCL/data/RNAseq/MCL_bulk_191006.csv") %>%
        RemoveDup() %>%
        log1p()
    meta.data = data.frame("label.main"= rep("MCL",ncol(MCL_bulk)),
                           "label.fine" = rep("MCL",ncol(MCL_bulk)),
                           "label.ont" = colnames(MCL_bulk))
    MCL_bulk_sce <- SummarizedExperiment(list(logcounts=MCL_bulk),
                                         colData=DataFrame(meta.data))
    #blue_encode
    blue_encode <- BlueprintEncodeData()
    common <- Reduce(intersect, list(rownames(sce),
                                     rownames(MCL_bulk_sce),
                                     rownames(blue_encode)
    ))
    length(common)
    combine_ref = do.call("cbind", list(blue_encode[common,],
                                        MCL_bulk_sce[common,]))
    table(combine_ref$label.fine)
    system.time(trained <- trainSingleR(ref = combine_ref,
                                        labels=combine_ref$label.fine))
    system.time(pred <- classifySingleR(sce[common,], trained))
    saveRDS(object = pred, file = "output/MCL_120_20230106_blueEncode_singleR_pred.rds")
}

if(reference == "MCL+azimuth_PBMC"){
    # MCL_bulk
    MCL_bulk <- read.csv(file="../scRNAseq-MCL/data/RNAseq/MCL_bulk_191006.csv") %>%
        RemoveDup() %>%
        log1p()
    meta.data1 = data.frame("celltype.l3" = rep("MCL",ncol(MCL_bulk)))
    MCL_bulk_sce <- SummarizedExperiment(list(logcounts=MCL_bulk),
                                         colData=DataFrame(meta.data1))
    # azimuth/PBMC
    path = "../seurat_resources/azimuth/PBMC/"
    counts <- Read10X(paste0(path, "GSE164378/GSM5008740_RNA_5P"))
    meta.data = read.csv(paste0(path,"GSE164378/GSE164378_sc.meta.data_5P.csv"),row.names =1)
    table(rownames(meta.data) == colnames(counts))

    object = CreateSeuratObject(counts,min.cells = 3,names.delim = "-",min.features = 3,
                                meta.data = meta.data)
    object %<>% NormalizeData()
    object$orig.ident <- object$donor
    object$celltype.l3_sample <- paste0(object$celltype.l3,"@",object$orig.ident)
    exp = AverageExpression(object,group.by = "celltype.l3_sample",assays = "RNA")
    exp = log1p(exp$RNA)
    meta.data <- object@meta.data
    meta.data = meta.data[!duplicated(meta.data$celltype.l3_sample),]
    rownames(meta.data) = meta.data$celltype.l3_sample
    PBMC <- SingleCellExperiment(list(logcounts=exp),
                                 colData=DataFrame("celltype.l3"=gsub("@.*","",colnames(exp)),
                                                   row.names = colnames(exp)))

    rm(counts,meta.data,object,exp);GC()

    common <- Reduce(intersect, list(rownames(sce),
                                     rownames(MCL_bulk_sce),
                                     rownames(PBMC)
    ))
    length(common)
    combine_ref = do.call("cbind", list(MCL_bulk_sce[common,],
                                        PBMC[common,]))
    table(combine_ref$celltype.l3)
    system.time(trained <- trainSingleR(ref = combine_ref,
                                        labels=combine_ref$celltype.l3))
    system.time(pred <- classifySingleR(sce[common,], trained))
    saveRDS(object = pred, file = "output/MCL_120_20230106_azimuth_MCL_singleR_pred.rds")
}

