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
markers = df_markers[,-grep("Alias",colnames(df_markers))]

marker.list <- df2list(markers)
marker.list %<>% lapply(function(x) HumanGenes(object,x))  %>%
        lapply(function(x) x[1:12]) %>% lapply(function(x) x[!is.na(x)])

marker.list %>% list2df %>% kable() %>% kable_styling()

for(i in 1:length(marker.list)){
        p <- lapply(marker.list[[i]], function(marker) {
                SingleFeaturePlot.1(object = object, feature = marker,pt.size = 0.5,
                                    gradient.use = c("lightblue", "blue3"))+
                        ggtitle(paste0(marker,Alias(gene = marker)))+
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
(markers = df_markers[,1] %>% as.matrix %>% as.character %>% HumanGenes(object,marker.genes = .))

table(object@meta.data$orig.ident)
table(object@ident)
object@meta.data$orig.ident = gsub("BH|DJ|MD|NZ","Normal",object@meta.data$orig.ident)

df_samples <- readxl::read_excel("doc/181227_scRNAseq_info.xlsx")
colnames(df_samples) <- tolower(colnames(df_samples))
tests <- paste0("test",c(6))
for(test in tests){
        sample_n = which(df_samples$tests %in% c("control",test))
        samples <- unique(df_samples$sample[sample_n])
        
        cell.use <- rownames(object@meta.data)[object@meta.data$orig.ident %in% 
                                                       c("Normal",samples)]
        subset.object <- SubsetData(object, cells.use = cell.use)
        subset.object@meta.data$orig.ident %>% unique %>% sort %>% print
        
        SplitSingleFeaturePlot(subset.object, 
                               select.plots = c(1,3,2,4),#c(6:8,1:5)
                               alias = df_markers, 
                               group.by = "ident",split.by = "orig.ident",
                               no.legend = T,label.size=3,do.print =T,
                               markers = markers, threshold = 0.1)
}
