############################################
# combine hpca and blueprint_encode
############################################
#devtools::install_github('dviraran/SingleR')
library(SingleR)
library(genefilter)
library(dplyr)
library(magrittr)
source("../R/Seurat3_functions.R")
source("../R/SingleR_functions.R")
source("R/util.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
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


# 1. check and filter blueprint_encode data==============================
data("blueprint_encode")
names(blueprint_encode)
unique(blueprint_encode$types)
dim(blueprint_encode$data)
head(blueprint_encode$types)
length(blueprint_encode$types)
length(unique(blueprint_encode$types))
length(unique(blueprint_encode$main_types))
dim(blueprint_encode$data)
head(blueprint_encode$data[,1:5])
table(is.na(blueprint_encode$data))
blueprint_encode$data[is.na(blueprint_encode$data)] = 0
head(colSums(blueprint_encode$data))
#testMMM(blueprint_encode$data)
#boxplot(blueprint_encode$data, main="blueprint_encode")#slow!
#boxplot(blueprint_encode$data[,1:100])#slow!
# remove low quanlity blueprint_encode data
par(mfrow=c(2,1))
hist(colMeans(blueprint_encode$data),breaks=ncol(blueprint_encode$data))
quantile_75 <- apply(blueprint_encode$data,2,function(x) quantile(x,probs =0.75))
hist(quantile_75, breaks=ncol(blueprint_encode$data))
rm_samples <- names(quantile_75)[quantile_75<1]
(rm_index <- which(colnames(blueprint_encode$data) %in% rm_samples))
blueprint_encode_rm <- blueprint_encode$data[,-rm_index]
quantile_75_new <- apply(blueprint_encode_rm,2,function(x) quantile(x,probs =0.75))
hist(quantile_75, breaks=ncol(blueprint_encode$data))
hist(quantile_75_new, breaks=ncol(blueprint_encode_rm),xlim = c(0,4.1))
par(mfrow=c(1,1))
#boxplot(blueprint_encode_rm)#slow!
title(main="blueprint_encode_rm")


sort(unique(blueprint_encode$main_types))
types = FineTune(as.character(blueprint_encode$types[-rm_index]),
                       main.type = FALSE)
main_types = FineTune(as.character(blueprint_encode$main_types[-rm_index]),
                            main.type = TRUE)
sort(unique(main_types))
Blueprint_encode = CreateSinglerReference(name = 'Blueprint_encode',
                                          expr = blueprint_encode_rm,
                                          types = types, 
                                          main_types = main_types)
save(Blueprint_encode,file='../SingleR/data/Blueprint_encode.RData')
# 2. check and prepare MCL data==============================
MCL_bulk <- read.csv(file="data/RNAseq/MCL_bulk_191006.csv", stringsAsFactors = F)
head(MCL_bulk);dim(MCL_bulk)
B_bulk <- read.csv(file="data/RNAseq/B_bulk_191006.csv", stringsAsFactors = F)
head(B_bulk);dim(B_bulk)

B_MCL_bulk <- inner_join(B_bulk, MCL_bulk, by = "X")
B_MCL_bulk <- RemoveDup(B_MCL_bulk)
head(B_MCL_bulk);dim(B_MCL_bulk)

B_MCL_bulk <- log1p(B_MCL_bulk)

# 3. merge MCL and Blueprint_encode =====================
(load(file='../SingleR/data/Blueprint_encode.RData'))
dim(Blueprint_encode$data)
MCL_blue_encode <- merge(B_MCL_bulk, Blueprint_encode$data,
                         by="row.names",all=FALSE)
MCL_blue_encode <- RemoveDup(MCL_blue_encode)
testMMM(MCL_blue_encode)

colsum <- colSums(MCL_blue_encode)
#scale_factor = median(colsum)
#MCL_blue_encode = MCL_blue_encode/colsum * scale_factor
#testMMM(MCL_blue_encode)

jpeg(paste0(path,"boxplot_MCL_blue_encode.jpeg"), units="in", width=10, height=7,res=600)
par(mfrow=c(1,1))
boxplot(MCL_blue_encode) #slow
title(main = "boxplot for Blueprint + Encode + MCL")
dev.off()

# Create Singler Reference
ref = CreateSinglerReference(name = 'MCL_blue_encode',
                             expr = as.matrix(MCL_blue_encode), # the expression matrix
                             types = c(rep("B_cells:PB",ncol(B_bulk)-1),
                                       rep("MCL",ncol(MCL_bulk)-1),
                                       Blueprint_encode$types), 
                             main_types = c(rep("B_cells",ncol(B_bulk)-1),
                                            rep("MCL",ncol(MCL_bulk)-1),
                                            Blueprint_encode$main_types))

save(ref,file='data/ref_MCL_blue_encode_20200225.RData')