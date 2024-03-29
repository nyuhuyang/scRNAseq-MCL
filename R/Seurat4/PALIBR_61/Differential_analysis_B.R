########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
#
# ######################################################################
####################################
invisible(lapply(c("Seurat","dplyr","magrittr","tidyr","data.table",
                   "future","gplots"), function(x) {
        suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_differential_expression.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

object = readRDS(file = "data/MCL_61_20220331.rds")

meta.data1 = readRDS(file = "data/MCL_61_20220331_metadata.rds")
meta.data2 = readRDS(file = "data/MCL_65_20220411_metadata.rds")
meta.data1 %<>% tibble::rownames_to_column("barcode")
meta.data2 %<>% tibble::rownames_to_column("barcode")

meta.data = left_join(meta.data1, meta.data2[,c("barcode","Doublets")],by = "barcode")
meta.data$barcode.1 = NULL
meta.data %<>% tibble::column_to_rownames("barcode")
meta.data$Doublets[is.na(meta.data$Doublets)] = "unknown"
saveRDS(meta.data, file = "data/MCL_61_20220331_metadata.rds")

if(all(colnames(object) == rownames(meta.data))) object@meta.data = meta.data


object = subset(object, subset =  Doublets == "Singlet"
                                & X6cluster  %in% c("1","2","3","4","5")
)
object$X6cluster %<>% paste0("C",.)
object$orig.ident_X6 = paste0(object$orig.ident,"_",object$X6cluster)
df <- as.data.frame.matrix(table(object$orig.ident,object$X6cluster))
exp <- AverageExpression(object,assays = "SCT",group.by = "orig.ident_X6")
exp = (log2(exp$SCT + 1))
exp %<>% as.data.table()
rownames(exp) = rownames(object)
fwrite(exp,paste0(path,"B_MCL_log2UMI.csv"),row.names = TRUE)
fwrite(df,paste0(path,"B_MCL_number.csv"),row.names = TRUE)


#==============
csv_names = sapply(opts,function(x) paste)
csv_index = list.files("output/20220420",pattern = ".csv") %>% gsub("-.*","",.) %>% as.integer()
table(1:13 %in% csv_index)
csv_names = list.files("output/20220420",pattern = ".csv")
deg_list <- pbapply::pblapply(csv_names, function(csv){
    tmp <- read.csv(paste0("output/20220420/",csv),row.names = 1)
    tmp = tmp[tmp$p_val_adj < 0.05,]
    tmp$gene = rownames(tmp)
    tmp %<>% group_by(cluster) %>% arrange(desc(avg_log2FC), .by_group = TRUE)
    tmp
})

deg = bind_rows(deg_list)
deg %<>% filter(p_val_adj < 0.05)
deg_list = split(deg, f = deg$cluster)

openxlsx::write.xlsx(deg_list, file = paste0(path,"AFT12_DEG.xlsx"),
           colNames = TRUE, borders = "surrounding")

