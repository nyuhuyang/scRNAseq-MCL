# test
invisible(lapply(c("Seurat","dplyr","ggplot2","cowplot","sctransform","magrittr",
                   "harmony"), function(x) {
                       suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
save.path <- paste0("output/",gsub("-","",Sys.Date()))
if(!dir.exists(save.path)) dir.create(save.path, recursive = T)
# Need 64GB ?
set.seed(101)
# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.integer(as.character(slurm_arrayid))
print(paste0("slurm_arrayid=",args))

test_df = data.frame(npcs = rep(c(60,85,100),time = 5),
                     perplexity = rep(c(30,50,100,175,300),each =3))

print(paste("npcs =",npcs <- test_df[args,"npcs"]))
print(paste("perplexity =",perplexity <- test_df[args,"perplexity"]))

object = readRDS(file = "data/MCL_52_20210715.rds")
DefaultAssay(object) = "SCT"
object@reductions$tsne = NULL

system.time(object %<>% RunTSNE(reduction = "harmony",dims = 1:npcs, seed.use = args, perplexity = perplexity))

Idents(object)  = "cell.types"
g <- TSNEPlot.1(object,group.by = "cell.types", title = seed.use,
            label = T,cols = ExtractMetaColor(object),raster=FALSE,
           do.print = T)

file.name = paste0("npcs=",npcs,"_perplexity=",perplexity)
g1 <- TSNEPlot.1(object,group.by = "cell.types", title = paste0("npcs = ",npcs,", perplexity = ",perplexity),
                 label = F,cols = ExtractMetaColor(object),raster=FALSE,
                 do.print = F,do.return = T)

jpeg(paste0(save.path, "/tSNE_", file.name,".jpeg"), units= "in",width=10, height=7,res=600)
print(g1)
dev.off()
reductions = object@reductions
saveRDS(reductions, file = paste0(save.path, "reductions_",file.name,"rds"))
