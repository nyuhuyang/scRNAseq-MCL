library(pheatmap)
library(genefilter)
library(dplyr)
library(magrittr)
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
source("../R/Seurat3_functions.R")

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

pt_names <- as.character(meta.data["sample",7:ncol(meta.data)])
meta.data["sample",7:ncol(meta.data)] = paste0("Pt",pt_names)
meta.data["patient",] = gsub(" .*","",meta.data["sample",])

genes.sd = genefilter::rowSds(MCL_B_bulk)
variable.genes = head(genes.sd[order(genes.sd,decreasing = T)], 4000)

y = MCL_B_bulk[names(variable.genes),];dim(y)
for(label in c("sample","cell.type","patient")){
        colnames(y) %<>% plyr::mapvalues(from = colnames(y), 
                                         to = as.character(meta.data[label,]))
        system.time(cor.mtr <- cor(y, method = "spearman"))
        hc_c <- hclust(as.dist(1-cor(cor.mtr, method="spearman")), method="ward.D2")
        hc_r <- hclust(as.dist(1-cor(t(cor.mtr), method="spearman")), method="ward.D2")
        diag(cor.mtr) <-NA
        jpeg(paste0(path,"bulk_RNAseq_heatmap_",label,".jpeg"), units="in", width=10, height=10,res=600)
        print(pheatmap(cor.mtr,cex=.9,
                       cluster_rows= hc_r,
                       cluster_cols = hc_c,
                       fontsize_row = 12,
                       fontsize_col = 12,
                       fontsize = 20,
                       main = ""))
        dev.off()
}