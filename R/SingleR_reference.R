############################################
# combine hpca and blueprint_encode
############################################
library(SingleR)
library(genefilter)
library(dplyr)
library(magrittr)
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
####functions===========

# check blueprint_encode data==============================
data("blueprint_encode")
names(blueprint_encode)
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
testMMM(blueprint_encode$data)
boxplot(blueprint_encode$data, main="blueprint_encode")#slow!
boxplot(blueprint_encode$data[,1:100])#slow!
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
boxplot(blueprint_encode_rm)#slow!
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

# check MCL data==============================
X181120_MCL_WTS <- readxl::read_excel("doc/181120 MCL WTS.xlsx", col_names = FALSE)

# remove NA columns
(remove_columns <- (X181120_MCL_WTS$X__1 == "gene") %>% which %>% .[1] %>% 
        X181120_MCL_WTS[.,] %>% is.na %>% as.vector)
X181120_MCL_WTS <- X181120_MCL_WTS[,!remove_columns]
head(X181120_MCL_WTS)
# remove NA rows
X181120_MCL_WTS <- X181120_MCL_WTS[!apply(X181120_MCL_WTS,1, function(x) all(is.na(x))),]

#rename column and remove gene row
(colnames(X181120_MCL_WTS) <- (X181120_MCL_WTS$X__1 == "gene") %>% which %>% .[1] %>% 
        X181120_MCL_WTS[.,] %>% as.vector)
head(X181120_MCL_WTS)

X181120_MCL_WTS <- X181120_MCL_WTS[-((X181120_MCL_WTS$gene == "gene") %>% which %>% .[1]),]
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
X181120_MCL_WTS <- RemoveDup(X181120_MCL_WTS)
dim(X181120_MCL_WTS)

testMMM(X181120_MCL_WTS)

# merge MCL and Blueprint_encode

(load(file='../SingleR/data/Blueprint_encode.RData'))
dim(Blueprint_encode$data)
MCL_blue_encode <- merge(log1p(X181120_MCL_WTS),Blueprint_encode$data,
                         by="row.names",all=FALSE)
rownames(MCL_blue_encode) = MCL_blue_encode$Row.names
MCL_blue_encode <- MCL_blue_encode[-which(colnames(MCL_blue_encode)=="Row.names")]
testMMM(MCL_blue_encode)

colsum <- colSums(MCL_blue_encode)
scale_factor = median(colsum)
MCL_blue_encode = MCL_blue_encode/colsum * scale_factor
testMMM(MCL_blue_encode)


jpeg(paste0(path,"boxplot_MCL_blue_encode.jpeg"), units="in", width=10, height=7,res=600)
boxplot(MCL_blue_encode) #slow
dev.off()

# Create Singler Reference
ref = CreateSinglerReference(name = 'MCL_blue_encode',
                             expr = as.matrix(MCL_blue_encode), # the expression matrix
                             types = c(paste0("MCL:",colnames(X181120_MCL_WTS)),
                                       Blueprint_encode$types), 
                             main_types = c(rep("MCL",ncol(X181120_MCL_WTS)),
                                            Blueprint_encode$main_types))

save(ref,file='data/ref_MCL_blue_encode.RData') # it is best to name the object and the file with the same name.