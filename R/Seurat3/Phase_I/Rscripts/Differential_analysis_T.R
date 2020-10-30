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
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
i <- as.numeric(slurm_arrayid)
print(paste0("slurm_arrayid=",i))
step <- switch(as.character(which(sapply(list(1:46, 47:198,199:200), function(x) args %in% x))),
               "1" = "1. pairwise T cell comparisions",
               "2" = "2. pairwise immune cell comparisions",
               "3" = "3. CD4 T vs CD8 T")
print(step)

# load data
(load(file="data/MCL_41_harmony_20200225.Rda"))
DefaultAssay(object) = "SCT"
Idents(object) = "Doublets"
object %<>% subset(idents = "Singlet")
object@meta.data[object$SCT_snn_res.0.8 %in% c(5,16), "cell.types"] = "Monocytes"
Idents(object) = "cell.types"

if(step == "1. pairwise T cell comparisions"){
        opts = data.frame(cell.types = rep(c("T_cells:CD8+","T_cells:CD4+"),  each = 23),
                          ident.1 = rep(c("Pt17_7","Pt17_12","Pt17_31","Pt28_4","Pt28_28",
                                          "Pt25_1_8","Pt25_24","Pt25_25Pd","Pt25_25Pd","Pt11_28",
                                          "Pt25_1","PtU01","PtU02","PtU03","PtU04",
                                          "Pt17_LN1","Pt17_2","Pt17_7","Pt17_12","Pt17_31",
                                          "PtB13_Ibp","PtB13_Ib1","PtB13_IbR"),2),
                          ident.2 = rep(c("Pt17_2","Pt17_2","Pt17_2","Pt28_1","Pt28_1",
                                          "Pt25_1","Pt25_1","Pt25_1","Pt25_24","Pt11_14",
                                          rep("N01",13)),2),
                          stringsAsFactors = F)
        (opt = opts[i,])  # need 32 GB

        object %<>% subset(idents = opt$cell.types)
        Idents(object) = "orig.ident"
        object %<>% subset(idents = c(opt$ident.1,opt$ident.2))

        markers <- FindAllMarkers.UMI(object,
                                      logfc.threshold = 0,
                                      only.pos = F,
                                      return.thresh = 1,
                                      test.use = "MAST",
                                      latent.vars = "nFeature_SCT")
        save.path = paste0(path,step,"/")
        if(!dir.exists(save.path))dir.create(save.path, recursive = T)
        write.csv(markers,paste0(save.path,sub(":","_",opt$cell.types),opt$ident.1, "-",opt$ident.2,".csv"))
}

if(step == "2. pairwise immune cell comparisions"){
        samples <- c("N04","PtU01","PtU02","PtU03","PtU04",
                     "Pt2_30Pd","Pt10_LN2Pd","Pt11_LN1","Pt11_1",
                     "Pt11_14","Pt11_28","Pt13_BMA1","Pt13_1a",
                     "Pt13_1b","Pt16_3Pd","Pt17_LN1","Pt17_2",
                     "Pt17_7","Pt17_12","Pt17_31","Pt19_BM2Pd",
                     "Pt20_1","Pt25_SB1","Pt25_1","Pt25_1_8",
                     "Pt25_24","Pt25_25Pd","Pt25_AMB25Pd","Pt27_1",
                     "Pt27_1_8","Pt27_12","PtB13_Ibp","PtB13_Ib1",
                     "PtB13_IbR","Pt28_LN1","Pt28_1","Pt28_4",
                     "Pt28_28")
        opts = data.frame(cell.types = rep(c("T_cells:CD8+","T_cells:CD4+","Monocytes","NK_cells"),  each = 38),
                          ident.1 = rep(samples, time = 4),
                          ident.2 = rep("N01", 38*4),
                          stringsAsFactors = F)
        i = i-46
        (opt = opts[i,])  # need 32 GB

        object %<>% subset(idents = opt$cell.types)
        Idents(object) = "orig.ident"
        object %<>% subset(idents = c(opt$ident.1,opt$ident.2))

        markers <- FindAllMarkers.UMI(object,
                                      logfc.threshold = 0,
                                      only.pos = F,
                                      return.thresh = 1,
                                      test.use = "MAST",
                                      latent.vars = "nFeature_SCT")
        save.path = paste0(path,step,"/")
        if(!dir.exists(save.path))dir.create(save.path, recursive = T)
        write.csv(markers,paste0(save.path,sub(":","_",opt$cell.types),opt$ident.1, "-",opt$ident.2,".csv"))
}

if(step == "3. CD4 T vs CD8 T"){
    opts = data.frame(ident.1 = c("T_cells:CD8+","T_cells:CD4+"),
                      ident.2 = c("T_cells:CD4+","T_cells:CD8+"),
                      stringsAsFactors = F)
    i = i-198
    (opt = opts[i,])  # need 32 GB

    object %<>% subset(idents = opt)
    markers <- FindMarkers.UMI(object,
                               ident.1 = opt$ident.1, ident.2 = opt$ident.2,
                               logfc.threshold = 0,
                               only.pos = F,
                               test.use = "MAST",
                               latent.vars = "nFeature_SCT")
    save.path = paste0(path,step,"/")
    if(!dir.exists(save.path))dir.create(save.path, recursive = T)
    write.csv(markers,paste0(save.path,opt$ident.1, "-",opt$ident.2,".csv"))
}
