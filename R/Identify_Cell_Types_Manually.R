library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(gplots)
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 2.1 identify phenotype for each cluster  ==========================================
(load(file="data/MCL_Harmony_20_20181231.Rda"))

#blueprint_encode_main = read.csv("../SingleR/output/blueprint_encode_main.csv",row.names =1,header = T,
#                                 stringsAsFactors = F)
df_markers <- readxl::read_excel("../SingleR/output/bio-rad-markers.xlsx")
df_markers <- readxl::read_excel("../SingleR/output/bio-rad-markers.xlsx",sheet = "Human.sub")
colnames(df_markers) = gsub(" ","_",colnames(df_markers))
colnames(df_markers) = gsub(":|\\/","_",colnames(df_markers))
colnames(df_markers) = gsub("\\+","",colnames(df_markers))
markers = df_markers[,-grep("Alias",colnames(df_markers))]
marker.list <- df2list(markers)

marker.list %<>% lapply(function(x) HumanGenes(object,x))  %>%
        lapply(function(x) x[1:12]) %>% lapply(function(x) x[!is.na(x)])

marker.list %>% list2df %>% kable() %>% kable_styling()

for(i in 1:length(marker.list)){
        p <- lapply(marker.list[[i]], function(marker) {
                SingleFeaturePlot.1(object = object, feature = marker,pt.size = 0.5,
                                    gradient.use = c("lightblue", "blue3"),threshold=0.1)+
                        ggtitle(paste0(marker,Alias(df = df_markers,gene = marker)))+
                        theme(plot.title = element_text(hjust = 0.5,size = 15, face = "bold"))
                })

        jpeg(paste0(path,names(marker.list)[i],".jpeg"),
             units="in", width=10, height=7,res=600)
        print(do.call(plot_grid, p)+ ggtitle(paste(names(marker.list)[i],"markers"))+
                      theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold")))
        dev.off()
        print(paste0(i,":",length(marker.list)))
}
#====== 2.2 focus on MCL heterogeneity, and T cell subsets and NK cells in the tSNE plots =====
df_markers <- readxl::read_excel("doc/MCL-markers.xlsx")
(markers <- df_markers[,1] %>% as.matrix %>% as.character %>% HumanGenes(object,marker.genes = .))
markers <- HumanGenes(object,marker.genes = c("CD19","SOX11","CD5"))
table(object@meta.data$orig.ident)
table(object@ident)
object@meta.data$orig.ident = gsub("BH|DJ|MD|NZ","Normal",object@meta.data$orig.ident)

df_samples <- readxl::read_excel("doc/181227_scRNAseq_info.xlsx")
colnames(df_samples) <- tolower(colnames(df_samples))
tests <- paste0("test",c(2))
for(test in tests){
        sample_n = which(df_samples$tests %in% c("control",test))
        samples <- unique(df_samples$sample[sample_n])
        
        cell.use <- rownames(object@meta.data)[object@meta.data$orig.ident %in% 
                                                       c("Normal",samples)]
        subset.object <- SubsetData(object, cells.use = cell.use)
        subset.object@meta.data$orig.ident %>% unique %>% sort %>% print
        
        SplitSingleFeaturePlot(subset.object, 
                               select.plots = c(1:4,8:5),#
                               #alias = df_markers, 
                               group.by = "ident",split.by = "orig.ident",
                               no.legend = T,label.size=3,do.print =T,nrow = 2,
                               markers = "SOX11", threshold = 1.02)
}
#' Find threshold using maximal UMI from normal control
#' @param object Seurat object
#' @param marker single gene name
#' @param control control sample name
#' @param celltype focus on specific celltype
#' @example Normal_control(subset.object, "CD5", celltype =c("B_cells","MCL"))
Normal_control <- function(object, marker, control, celltype = NULL){
        if(!is.null(celltype)){
                cell.use <- grep(paste(celltype,collapse = "|"),object@ident, value = T)
                object <- SubsetData(object, cells.use = names(cell.use))
        } 
        object@meta.data$orig.ident = gsub(paste(control,collapse = "|"),
                                                  "Normal",object@meta.data$orig.ident)
        object <- SetAllIdent(object, id = "orig.ident")
        subset.object <- SubsetData(object, ident.use = "Normal")
        
        normal.max <- subset.object@data[marker,] %>% max 
        
        return(normal.max+0.01)
}
Normal_control(subset.object, marker = "SOX11", control = c("BH","DJ","MD","NZ"), 
               celltype =c("B_cells","MCL"))



#object@meta.data$orig.ident =gsub("Pt-11-LN-C14","Pt-11-C14",object@meta.data$orig.ident)
