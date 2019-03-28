library(SingleR)
source("../R/SingleR_functions.R")
path <- paste0(getwd(),"/output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.1 Create Singler Object  ==========================================
(load(file="data/MCL.raw.data_Harmony_30_20190320.Rda"))
(load(file = 'data/ref_MCL_blue_encode_20190315.RData'))
singler <- CreateBigSingleRObject.1(raw_data, annot = NULL, project.name="EC-MDL",
                       N = 10000, min.genes = 200, technology = "10X",
                       species = "Human", citation = "", ref.list = list(ref),
                       normalize.gene.length = F, variable.genes = "de", fine.tune = T,
                       reduce.file.size = F, do.signatures = F, do.main.types = T,
                       temp.dir = getwd(), numCores = SingleR.numCores)
save(singler,file="./output/singlerT_MCL_30_20190320.Rda")