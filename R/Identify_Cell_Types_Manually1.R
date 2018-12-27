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
(load(file="data/MCL_Harmony_12_20181121.Rda"))

blueprint_encode_main = read.csv("../SingleR/output/blueprint_encode_main.csv",row.names =1,header = T,
                                 stringsAsFactors = F)
blueprint_encode_sub = read.csv("../SingleR/output/blueprint_encode_sub.csv",row.names =1,header = T,
                                 stringsAsFactors = F)

MCL_markers = read.csv("doc/MCL_markers.csv",row.names =1,header = T,
                                stringsAsFactors = F)

#marker.list <- df2list(blueprint_encode_main)
marker.list <- df2list(blueprint_encode_sub)
marker.list <- df2list(MCL_markers)
marker.list %<>% lapply(function(x) HumanGenes(MCL,x[1:80]))  %>%
        lapply(function(x) x[1:9]) %>% lapply(function(x) x[!is.na(x)])

marker.list <- grep("B_cells", names(marker.list),
                    value = T) %>% marker.list[.]
marker.list %>% list2df %>% t %>% kable() %>% kable_styling()

blueprint_encode_sub = read.csv("../SingleR/output/blueprint_encode_sub.csv",row.names =1,header = T,
                                stringsAsFactors = F)

FeaturePlot.1 <- function(object = MCL, x){
        p <- FeaturePlot(object = object, 
                         reduction.use = "tsne",
                         features.plot = x, min.cutoff = NA, do.return =T,
                         cols.use = c("lightgrey","blue"), pt.size = 0.5)
        return(p)
}

dev.off()

marker.list <- list("PD1" = c("PDCD1","CD274","PDCD1LG2","CTLA4"))
for(i in 1:length(marker.list)){
        p <- FeaturePlot.1(object = MCL, x = marker.list[[i]])
        p1 <- do.call(plot_grid, p)
        p1 <- p1 + ggtitle(paste(names(marker.list)[i],"markers"))+
                theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"))
        jpeg(paste0(path,names(marker.list)[i],".jpeg"),
             units="in", width=10, height=7,res=600)
        print(p1)
        print(paste0(i,":",length(marker.list)))
        dev.off()
}
