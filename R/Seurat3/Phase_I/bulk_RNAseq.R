library(pheatmap)
library(genefilter)
library(dplyr)
library(magrittr)
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")

# 2. check and prepare MCL data==============================
MCL_bulk <- read.csv(file="data/RNAseq/MCL_bulk_191006.csv", stringsAsFactors = F)
head(MCL_bulk);dim(MCL_bulk)
B_bulk <- read.csv(file="data/RNAseq/B_bulk_191006.csv", stringsAsFactors = F)
head(B_bulk);dim(B_bulk)
meta.data <- read.csv(file="data/RNAseq/191006 PALIBR WTS_meta.data.csv",
                      stringsAsFactors = F,row.names = 1)
head(meta.data);dim(meta.data)

####functions===========
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
MCL_B_bulk <- inner_join(B_bulk, MCL_bulk, by = "X")
MCL_B_bulk <- RemoveDup(MCL_B_bulk);dim(MCL_B_bulk)
colnames(MCL_B_bulk) %<>% gsub("X","Pt",.)
#======= clean "1-Sep"  "1-Mar"  "10-Sep" "10-Mar" "1-Dec"  ==========
grep("^[0-9]{1,2}-",rownames(MCL_B_bulk),value = T) %>% sort
Excel_gene = data.frame("Excel_gene" = grep("^[0-9]{1,2}-",rownames(MCL_B_bulk),value = T) %>% sort)
Excel_gene$hg38 = c("BHLHE40",
                "MARCHF1",
                "SEPTIN1",
                "MARCHF10",
                "SEPTIN10",
                "MARCHF11",
                "SEPTIN11",
                "SEPTIN12",
                "SEPTIN14",
                "SELENOF",
                "MARCHF2",
                "SEPTIN2",
                "MARCHF3",
                "SEPTIN3",
                "MARCHF4",
                "SEPTIN4",
                "MARCHF5",
                "SEPTIN5",
                "MARCHF6",
                "SEPTIN6",
                "MARCHF7",
                "SEPTIN7",
                "MARCHF8",
                "SEPTIN8",
                "MARCHF9",
                "SEPTIN9")
rownames(MCL_B_bulk) %<>% plyr::mapvalues(from = Excel_gene$Excel_gene,
                                          to = Excel_gene$hg38)

colnames(meta.data) %<>% gsub("X","Pt",.)
meta.data["patient",] = gsub("\\..*","", colnames(meta.data))
meta.data %<>% t %>% as.data.frame()
genes.sd = genefilter::rowSds(MCL_B_bulk)
variable.genes = head(genes.sd[order(genes.sd,decreasing = T)], 4000)

labels = c("specimens","progress","patient")
for(i in seq_along(labels)){
        y = MCL_B_bulk[names(variable.genes),]
        colnames(y) %<>% plyr::mapvalues(from = colnames(y),
                                         to = as.character(meta.data[labels[i],]))
        system.time(cor.mtr <- cor(y, method = "spearman"))
        hc_c <- hclust(as.dist(1-cor(cor.mtr, method="spearman")), method="ward.D2")
        hc_r <- hclust(as.dist(1-cor(t(cor.mtr), method="spearman")), method="ward.D2")
        diag(cor.mtr) <-NA
        jpeg(paste0(path,"bulk_RNAseq_heatmap_",labels[i],".jpeg"), units="in", width=10, height=10,res=600)
        print(pheatmap(cor.mtr,cex=.9,
                       cluster_rows= hc_r,
                       cluster_cols = hc_c,
                       fontsize_row = 12,
                       fontsize_col = 12,
                       fontsize = 20,
                       main = ""))
        dev.off()
        Progress(i, length(labels))
}

bulk <- CreateSeuratObject(counts = MCL_B_bulk, meta.data = meta.data)
bulk <- NormalizeData(bulk, normalization.method = "LogNormalize", scale.factor = 10000)
bulk <- FindVariableFeatures(bulk, selection.method = "vst", nfeatures = 2000)
bulk <- ScaleData(bulk)
bulk <- RunPCA(bulk, features = VariableFeatures(object = bulk))

# These are now standard steps in the Seurat workflow for visualization and clustering
bulk <- RunTSNE(bulk, dims = 1:10, verbose = FALSE,perplexity = 20)
bulk <- RunUMAP(bulk, dims = 1:10, verbose = FALSE)
bulk <- FindNeighbors(bulk, dims = 1:10, verbose = FALSE)
bulk <- FindClusters(bulk, verbose = FALSE,  resolution = 1.1)
DimPlot(bulk, reduction = "tsne", group.by = "RNA_snn_res.1.1")
DimPlot(bulk, reduction = "umap", group.by = "RNA_snn_res.1.1")
FeaturePlot.1(bulk,features = c("CCND1","SOX11"))
bulk$orig.ident = rownames(bulk@meta.data)
saveRDS(bulk,"data/bulk_RNAseq_20190307.rds")
length(unique(bulk$patient))

sample_list <- list("bulk-seq" = sort(as.character(meta.data[2,])),
                    "scRNA-seq" = sort(as.character(unique(object$orig.ident))))
df <- list2df(sample_list)
write.xlsx(df, asTable = F,
           file = paste0(path,"scRNA_bulk_samples.xlsx"),
           borders = "surrounding")
