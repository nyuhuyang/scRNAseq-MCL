#devtools::install_github('dviraran/SingleR')
library(SingleR)
source("../R/SingleR_functions.R")
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.1 Create Singler Object  ==========================================
(load(file = "data/MCL_41_harmony_20191231.Rda"))
object_data <- object@assays$SCT@data
remove(object); GC()                                                             
(load(file = 'data/ref_MCL_blue_encode_20190916.RData'))
singler <- CreateBigSingleRObject.1(object_data, annot = NULL, project.name="EC-MDL",
                                    N = 5000, min.genes = 200, technology = "10X",
                                    species = "Human", citation = "", ref.list = list(ref),
                                    normalize.gene.length = F, variable.genes = "de", fine.tune = T,
                                    reduce.file.size = F, do.signatures = F, do.main.types = T,
                                    temp.dir = getwd(), numCores = SingleR.numCores)
save(singler,file="output/singlerT_MCL_41_20200203.Rda")
