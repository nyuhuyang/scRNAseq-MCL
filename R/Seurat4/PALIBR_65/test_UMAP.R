invisible(lapply(c("Seurat","dplyr","ggplot2","cowplot","pbapply","sctransform","harmony"), function(x) {
    suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
save.path <- paste0("output/",gsub("-","",Sys.Date()))
if(!dir.exists(save.path)) dir.create(save.path, recursive = T)

#======1.1 Setup the Seurat objects =========================
# read sample summary list
df_samples <- readxl::read_excel("output/20220411/20220411_scRNAseq_info.xlsx")
df_samples = as.data.frame(df_samples)
colnames(df_samples) %<>% tolower()
df_samples %<>% filter(sequence %in% "GEX") %>% filter(phase %in% c("PALIBR_I","Cell_line")) %>%
    filter(sample != "Pt11_31")
nrow(df_samples)
#======1.2 load  Seurat =========================

object = readRDS(file = "data/MCL_65_20220411.rds")

object@meta.data = readRDS(file = "data/MCL_65_20220411_metadata.rds")
npcs = 83
object %<>% FindNeighbors(reduction = "pca",dims = 1:npcs)
resolutions = c( 0.01, 0.1, 0.2, 0.5,0.8)
for(i in 1:length(resolutions)){
    object %<>% FindClusters(resolution = resolutions[i], algorithm = 1)
    Progress(i,length(resolutions))
}

DefaultAssay(object) = "SCT"
s.genes <- cc.genes$s.genes %>% gsub("MLF1IP","CENPU",.)
g2m.genes <- cc.genes$g2m.genes %>% plyr::mapvalues(from = c("FAM64A", "HN1"),
                                                    to = c("PIMREG","JPT1"))
object %<>% CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
colnames(object@meta.data) %<>% sub("Phase","cell cycle phase",.)
meta.data = object@meta.data


meta.data$Mean.Reads.per.Cell %<>% gsub(",","",.) %>% as.integer()
meta.data$Number.of.Reads %<>% gsub(",","",.) %>% as.integer()
meta.data$Sequencing.Saturation %<>% gsub("%","",.) %>% as.numeric()
meta.data$barcode  = rownames(meta.data)
# find X6cluster
meta_data = readRDS(file = "data/MCL_61_20220331_metadata.rds")
meta_data$barcode  = rownames(meta_data)

meta.data %<>% left_join(meta_data[,c("X6cluster","X4cluster","barcode")], by = "barcode")
meta.data %<>% tibble::column_to_rownames(var = "barcode")

table(colnames(object) == rownames(meta.data))

meta.data[meta.data$phase %in% "Cell_line","X6cluster"] = "Cell_line"
saveRDS(meta.data, file = "data/MCL_65_20220411_metadata.rds")



