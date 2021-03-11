########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
#
# ######################################################################

library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(magrittr)
library(gplots)
library(cowplot)
library(eulerr)
library(openxlsx)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")

save.path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(save.path)) dir.create(save.path, recursive = T)

#=====
B_cells_MCL = readRDS(file = "data/MCL_41_B_20200225.rds")
Idents(B_cells_MCL) = "orig.ident"
object <- subset(B_cells_MCL,idents = c("N01","Pt25_SB1","Pt25_24","Pt25_25Pd","Pt25_AMB25Pd"))
object <- subset(B_cells_MCL,idents = c("PtB13_Ibp","PtB13_Ib1","PtB13_IbR"))

Idents(object) = "X4clusters"
object <- subset(object, idents = "C4")
gene_pairs = list(c("PCNA","TXNIP"),
                 c("SLC2A1","TXNIP"),
                 c("PCNA","SLC2A1"),
                 c("MYC","TXNIP"),
                 c("EZH2","HLA-A"))
for(gene_pair in gene_pairs) {
    exp <- FetchData(object,slot = "data", var = c(gene_pair,"orig.ident"))
    #cut_exp <- cut(exp[,1],breaks = 2)
    #print(gene_pair[1])
    #print(levels(cut_exp))
    colnames(exp) %<>% gsub("-","_",.)
    gene_pair  %<>% gsub("-","_",.)
    exp[,1] = plyr::mapvalues(x =exp[,1] > switch(EXPR = gene_pair[1],
                                                  "PCNA" = 1,
                                                  "SLC2A1" = 0,
                                                  "MYC" = 1,
                                                  "EZH2" = 0), #cut_exp,
                    from = c(FALSE,TRUE),#levels(cut_exp),
                    to = c("Low", "High"))
    p <- ggboxplot(exp, x = gene_pair[1], y = gene_pair[2],
                   color = gene_pair[1], palette = "jco",
                   facet.by = "orig.ident",
                   line.color = "gray", line.size = 0.4,
                   add = "jitter",nrow =1)
    #  Add p-value

    jpeg(paste0(save.path, gene_pair[1],"_",gene_pair[2],".jpeg"),units="in", width=10, height=7,res=300)
    print(p + stat_compare_means(method = "t.test"))
    dev.off()
}

