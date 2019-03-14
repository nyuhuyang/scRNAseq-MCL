library(SingleR)
path <- paste0(getwd(),"/output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.1 Create Singler Object  ==========================================
(load(file="data/MCL.data_Harmony_24_20190128.Rda"))
(load(file = 'data/ref_MCL_blue_encode.RData'))
singler = CreateSinglerObject(object.data, annot = NULL, project.name="EC-MDL",
                              min.genes = 500,technology = "10X", species = "Human", citation = "",
                              ref.list = list(ref),normalize.gene.length = F, variable.genes = "de",
                              fine.tune = F, do.signatures = F, clusters = NULL,
                              numCores = 1)
save(singler,file="./output/singlerT_MCL_24_20190128.Rda")
