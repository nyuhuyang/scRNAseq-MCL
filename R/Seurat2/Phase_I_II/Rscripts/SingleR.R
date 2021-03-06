library(SingleR)
library(Seurat)
source("../R/SingleR_functions.R")
path <- paste0(getwd(),"/output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.1 Create Singler Object  ==========================================
(load(file="data/MCL_Harmony_43_20190430.Rda"))
object_data <- object@data
remove(object);GC()
(load(file = 'data/ref_MCL_blue_encode_20190315.RData'))
singler <- CreateBigSingleRObject.1(object_data, annot = NULL, project.name="EC-MDL",
                       N = 2500, min.genes = 200, technology = "10X",
                       species = "Human", citation = "", ref.list = list(ref),
                       normalize.gene.length = F, variable.genes = "de", fine.tune = T,
                       reduce.file.size = F, do.signatures = F, do.main.types = T,
                       temp.dir = getwd(), numCores = SingleR.numCores/2)
save(singler,file="output/singlerT_MCL_43_20190430.Rda")
