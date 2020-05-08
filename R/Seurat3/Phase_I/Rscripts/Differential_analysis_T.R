########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
####################################
invisible(lapply(c("Seurat","dplyr","magrittr","tidyr",
                   "MAST","future","gplots"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
i <- as.numeric(slurm_arrayid)
print(paste0("slurm_arrayid=",i))

opts = data.frame(cell.types = rep(c("T_cells:CD8+","T_cells:CD4+"),  each = 11),
                  ident.1 = rep(c("Pt17_7","Pt17_12","Pt17_31","Pt28_4","Pt28_28","Pt25_1",
                              "Pt25_1_8","Pt25_24","Pt25_25Pd","Pt25_25Pd","Pt11_28"),2),
                  ident.2 = rep(c("Pt17_2","Pt17_2","Pt17_2","Pt28_1","Pt28_1","N01",
                              "Pt25_1","Pt25_1","Pt25_1","Pt25_24","Pt11_14"),2),
                  stringsAsFactors = F)
(opt = opts[i,])
# load data
(load(file="data/MCL_41_harmony_20200225.Rda"))
DefaultAssay(object) = "SCT"
Idents(object) = "Doublets"
object %<>% subset(idents = "Singlet")
Idents(object) = "cell.types"
object %<>% subset(idents = opt$cell.types)
Idents(object) = "orig.ident"
object %<>% subset(idents = c(opt$ident.1,opt$ident.2))

T_markers <- FindMarkers.UMI(object, 
                            logfc.threshold = 0, 
                            ident.1 = opt$ident.1,
                            ident.2 = opt$ident.2,
                            only.pos = F,
                            test.use = "MAST",
                            return.thresh = 1, 
                            latent.vars = "nCount_SCT")
T_markers$cluster = paste0(opt$ident.1, " \\ ",opt$ident.2)
write.csv(T_markers,paste0(path,sub(":","_",opt$cell.types),opt$ident.1, "\\",opt$ident.2,".csv"))
