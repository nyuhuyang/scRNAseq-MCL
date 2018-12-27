########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
library(Seurat)
library(SingleR)
library(dplyr)
library(tidyr)
library(kableExtra)
library(magrittr)
library(gplots)
source("../R/Seurat_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

# 5.1 load data ==============
(load(file="data/MCL_Harmony_12_20181121.Rda"))
##################################
# select B cells only
##################################
MCL <- SetAllIdent(MCL, id="res.0.6")
table(MCL@ident)
TSNEPlot.1(MCL,do.label = T)
B_cells_MCL <- SubsetData(MCL, ident.use = c(0,2,3,7,8,10,12))
remove(MCL);GC()
B_cells_MCL <- SetAllIdent(B_cells_MCL, id="singler2main")
table(B_cells_MCL@ident)

B_cells_MCL <- SubsetData(B_cells_MCL, ident.use = c("B_cells","HSC","MCL"))
table(B_cells_MCL@meta.data$singler1main)
B_cells_MCL <- SetAllIdent(B_cells_MCL, id="singler1main")
B_cells_MCL <- SubsetData(B_cells_MCL, ident.use = c("B_cells"))
B_cells_MCL <- SetAllIdent(B_cells_MCL, id="singler2sub")
p4 <- TSNEPlot.1(B_cells_MCL, do.label = F, do.return = T, pt.size = 0.5, 
                 colors.use = ExtractMetaColor(B_cells_MCL), no.legend =T)

B_cells_MCL <- SetAllIdent(B_cells_MCL, id="res.0.6")
idents <- as.data.frame(table(B_cells_MCL@ident))
old.ident.ids <- idents$Var1
new.cluster.ids <- c(0,2,1,1,3,4,5)
B_cells_MCL@ident <- plyr::mapvalues(x = B_cells_MCL@ident,
                                     from = old.ident.ids, to = new.cluster.ids)
B_cells_MCL@ident <- factor(B_cells_MCL@ident, levels = 0:5)
B_cells_MCL <- StashIdent(object = B_cells_MCL, save.name = "5_clusters")

p3 <- TSNEPlot.1(B_cells_MCL, do.return = T, pt.size = 0.5, do.label = T, 
                 group.by = "ident",no.legend =T )

jpeg(paste0(path,"/S1_TSNEPlot.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p3, p4)
dev.off()

##############################
# Split PrepareGSEA
###############################

df_samples <- readxl::read_excel("doc/181128_scRNAseq_info.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
# subset Seurat module in SingleR 
B_cells_MCL <- SetAllIdent(B_cells_MCL, id="orig.ident")
tests <- paste0("test",c(3:4))
control = ""
for(test in tests){
        sample_n = which(df_samples$tests %in% test)
        samples <- unique(df_samples$sample[sample_n])
        print(c(control,samples))
        
        cell.use <- rownames(B_cells_MCL@meta.data)[B_cells_MCL@meta.data$orig.ident %in% 
                                                            c(control,samples)]
        subset.B_cells_MCL <- SubsetData(B_cells_MCL, cells.use = cell.use)
        #print(continuous.label <- unique(subset.B_cells_MCL@meta.data$orig.ident))#[c(1,4,3,2,5)])
        PrepareGSEA(subset.B_cells_MCL, k <- 1, continuous.label = NULL)
        file.rename(paste0(path,list.files(path,pattern=paste0("orig.ident",k,"*"))), 
                    paste0(path,test,"_B_",k,c(".cls",".txt")))
}

# Run GSEA and generate reports
(gsea_path <- paste0("~/gsea_home/output/",tolower(format(Sys.Date(), "%b%d")), 
                     "/test4_B_1.c6.all.Gsea.1545890081484"))
#(pos.xls.path <- list.files(gsea_path,pattern="gsea_report_for_.*pos.*xls"))
(pos.xls.path <- list.files(gsea_path,pattern="gsea_report_for_.*xls")[2])
GSEA_output <- readr::read_delim(paste0(gsea_path,"/",pos.xls.path),"\t", 
                                 escape_double = FALSE, trim_ws = TRUE)
GSEA_output %>% .[-c(2,3,12)] %>% head(50) %>% kable() %>% kable_styling()

(GSEA.plots <- sapply(GSEA_output$NAME[1:9], function(name) {
        paste0("enplot_",name, "_([0-9]+)*\\.png$")}) %>%
                sapply(function(x) list.files(path = gsea_path, pattern =x)) %>%
                .[sapply(.,length)>0] %>% #Remove empty elements from list with character(0)
                paste(gsea_path, ., sep = "/")) 
CombPngs(GSEA.plots, ncol = 3)
file.rename(paste0(path,list.files(path,pattern=paste0("GSEA.plots_CombPngs"))), 
            paste0(path,"test4_B_1.c6.all.Gsea.1545890081484",".jpeg"))
