library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(magrittr)
library(gplots)
library(readxl)
library(pheatmap)
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
source("R/util.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#3.1  Compare DE across all major cell types==================
# 3.1.1 load bulk RNA-seq data
df_samples <- readxl::read_excel("doc/190126_scRNAseq_info.xlsx",
                                 sheet = "bulk")
colnames(df_samples) = tolower(colnames(df_samples))
df_samples$bulk.rnaseq.id = gsub(" ","_",df_samples$bulk.rnaseq.id)

MDL_6669 <- read_excel("data/GRCF-Elemento-MDL-6669-expression.genes.xls")
MDL_6392 <- read_excel("data/GRCF-Elemento-MDL-6392-expression.genes.xls")

(sample_6669 <- colnames(MDL_6669)[colnames(MDL_6669) %in% df_samples$bulk.rnaseq.id])
MDL_6669 = MDL_6669[,c("Gene_name",sample_6669)]

(sample_6392 <- colnames(MDL_6392)[colnames(MDL_6392) %in% df_samples$bulk.rnaseq.id])
MDL_6392 = MDL_6392[,c("Gene_name",sample_6392)]

colnames(MDL_6392)[-1] = paste0(colnames(MDL_6392)[-1],"_bulk")

MCL_bulk <- merge(MDL_6669, MDL_6392, by = "Gene_name")
MCL_bulk <- RemoveDup(MCL_bulk)
log_MCL_bulk <- log1p(MCL_bulk)
testMMM(log_MCL_bulk)
# 3.1.2 load sc RNA-seq data
(load(file="data/MCL_Harmony_24_20190128.Rda"))
object %<>% SetAllIdent(id="orig.ident")
MCL_scRNA <- AverageExpression(object)
testMMM(MCL_scRNA)

# 3.1.3 merge bulk and sc
MCL_bulk_sc <- merge(log_MCL_bulk, MCL_scRNA, by="row.names")
table(duplicated(MCL_bulk_sc$Row.names))
rownames(MCL_bulk_sc) = MCL_bulk_sc$Row.names
MCL_bulk_sc = MCL_bulk_sc[,-which(colnames(MCL_bulk_sc) %in% "Row.names")]
testMMM(MCL_bulk_sc)
dim(MCL_bulk_sc)

var.genes <- rownames(MCL_bulk_sc)[rownames(MCL_bulk_sc) %in% object@var.genes]
length(var.genes)
MCL_bulk_sc = MCL_bulk_sc[var.genes,]

testMMM(MCL_bulk_sc)
MCL_bulk_sc <- apply(MCL_bulk_sc, 2, function(x) (x- mean(x))/sd(x))

boxplot(MCL_bulk_sc)

# 3.1.4 Spearman correlation  ==================
y <- cor(MCL_bulk_sc, method="spearman")
diag(c) <-NA
ident_num <- ncol(MCL_bulk)
object_c_gendata <- c[(ident_num+1):nrow(c),1:ident_num]

jpeg(paste0(path,"cluster heatmap.jpeg"), units="in", width=10, height=7,res=600)
print(pheatmap(object_c_gendata,cex=.9,
         cluster_rows= F,
         cluster_cols = T,
         fontsize_row = 15,
         fontsize_col = 15,
         fontsize = 20,
         main = "cluster heatmap"))
dev.off()
