#devtools::install_github('dviraran/SingleR')
library(SingleR)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/SingleR_functions.R")
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.1 Create Singler Object  ==========================================
(load(file = "data/MCL_41_harmony_20191231.Rda"))
object_data <- object@assays$SCT@data
remove(object); GC()
#(load(file = 'data/ref_MCL_blue_encode_20200225.RData'))
ref = readRDS(file='data/ref_MCL_blue_encode_GSE107011_20200605.rds')
singler <- CreateBigSingleRObject.1(object_data, annot = NULL, project.name="EC-MDL",
                                    N = 5000, min.genes = 200, technology = "10X",
                                    species = "Human", citation = "", ref.list = list(ref),
                                    normalize.gene.length = F, variable.genes = "de", fine.tune = T,
                                    reduce.file.size = F, do.signatures = F, do.main.types = T,
                                    temp.dir = getwd(), numCores = SingleR.numCores)
saveRDS(singler,file="output/singlerT_MCL_41_20200605.rds")
