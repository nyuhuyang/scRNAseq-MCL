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
(load(file="data/MCL_Harmony_24_20190128.Rda"))

#blueprint_encode_main = read.csv("../SingleR/output/blueprint_encode_main.csv",row.names =1,header = T,
#                                 stringsAsFactors = F)
df_markers <- readxl::read_excel("../seurat_resources/bio-rad-markers.xlsx")
df_markers <- readxl::read_excel("../seurat_resources/bio-rad-markers.xlsx",sheet = "Human.sub")
colnames(df_markers) = gsub(" ","_",colnames(df_markers))
colnames(df_markers) = gsub(":|\\/","_",colnames(df_markers))
colnames(df_markers) = gsub("\\+","",colnames(df_markers))
markers = df_markers[,-grep("Alias",colnames(df_markers))]
marker.list <- df2list(markers)

marker.list %<>% lapply(function(x) x[1:16]) %>% 
        lapply(function(x) HumanGenes(object,x)) %>% 
        lapply(function(x) x[!is.na(x)])
marker.list %>% list2df %>% t %>% kable() %>% kable_styling()

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
# select 1/4 of cell from control
normal_cells = lapply(c("BH","DJ","MD","NZ"), function(x){
        rownames(object@meta.data)[(object@meta.data$orig.ident %in% x)]
})
remove_normal_cells = lapply(normal_cells, function(x) sample(x, size = length(x)*3/4)) %>%
        unlist
table(object@cell.names %in% remove_normal_cells)
cell.use <- object@cell.names[!(object@cell.names %in% remove_normal_cells)]
object <- SubsetData(object, cells.use = cell.use)

object@meta.data$orig.ident = gsub("BH|DJ|MD|NZ","Normal",object@meta.data$orig.ident)

# ===== markers======
df_markers <- readxl::read_excel("doc/MCL-markers.xlsx", sheet = "20190128")
(markers <- df_markers[,1] %>% as.matrix %>% as.character %>% HumanGenes(object,marker.genes = .))
table(object@meta.data$orig.ident)
table(object@ident)
# ===== sample list ======
df_samples <- readxl::read_excel("doc/190126_scRNAseq_info.xlsx")
colnames(df_samples) <- tolower(colnames(df_samples))
tests <- paste0("test",c(3))
for(test in tests){
        sample_n = which(df_samples$tests %in% c("control",test))
        samples <- unique(df_samples$sample[sample_n])
        
        cell.use <- rownames(object@meta.data)[object@meta.data$orig.ident %in% 
                                                       c("Normal",samples)]
        subset.object <- SubsetData(object, cells.use = cell.use)
        subset.object@meta.data$orig.ident %>% unique %>% sort %>% print
        SplitTSNEPlot(subset.object,do.label = F,select.plots = c(1,2,5,3,4), do.print = T, do.return=F)
        SplitSingleFeaturePlot(subset.object, 
                               select.plots = c(1,2,5,3,4),#
                               alias = df_markers,
                               group.by = "ident",split.by = "orig.ident",
                               no.legend = T,label.size=3,do.print =T,nrow = 2,
                               markers = markers, threshold = NULL)   
        
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

#######################
# RidgePlot
########################
cell.use <- grep(paste(c("B_cells","MCL"),collapse = "|"),subset.object@ident, value = T)
subset.subset.object <- SubsetData(subset.object, cells.use = names(cell.use))
subset.subset.object <- SplitSeurat(subset.subset.object, split.by = "CD5")

jpeg(paste0(path,"/RidgePlot~~.jpeg"), units="in", width=10, height=7,res=600)
RidgePlot(object = subset.subset.object, features.plot = HumanGenes(object,c("CD5")),
          group.by = "orig.ident",nCol = 1, y.log = T)
dev.off()

#object@meta.data$orig.ident =gsub("Pt-11-LN-C14","Pt-11-C14",object@meta.data$orig.ident)

for(i in seq(0.8,2.2,length.out = 15)){
        SplitSingleFeaturePlot(subset.object, 
                               select.plots = c(1:4,8:5),#
                               #alias = df_markers, 
                               group.by = "ident",split.by = "orig.ident",
                               no.legend = T,label.size=3,do.print =T,nrow = 2,
                               markers = "CCND1", threshold = i)
}

#######################
# RidgePlot
########################
object <- SetAllIdent(object, id = "orig.ident")
sample_n = which(df_samples$tests %in% c("control","test2"))
(samples <- unique(df_samples$sample[sample_n]))
subset.object <- SubsetData(object , ident.use = samples )

object.CCND1 <- SplitSeurat(object,split.by = "CCND1")
jpeg(paste0(path,"RidgePlot_CD5.jpeg"), units="in", width=10, height=7,res=600)
RidgePlot(object,features.plot = "CD5",ident.include = samples)
dev.off()

SplitSingleFeaturePlot(subset.object, 
                       select.plots = c(1:4,8:5),#
                       #alias = df_markers,
                       group.by = "ident",split.by = "orig.ident",
                       no.legend = T,label.size=3,do.print =T,nrow = 2,
                       markers = "CD5", threshold = 0.001)
