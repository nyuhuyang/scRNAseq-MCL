########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
#
# ######################################################################
####################################
invisible(lapply(c("Seurat","dplyr","magrittr","tidyr",
                   "future","gplots"), function(x) {
        suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_differential_expression.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)


# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.integer(as.character(slurm_arrayid))
print(paste0("slurm_arrayid=",args))

object = readRDS(file = "data/MCL_61_20220331.rds")
meta.data = readRDS(file = "data/MCL_61_20220331_metadata.rds")
if(all(colnames(object) == rownames(meta.data))) object@meta.data = meta.data

DefaultAssay(object) = "SCT"

object = subset(object, subset =  Doublets == "Singlet")


step = 1
# choose == "orig.ident_X6_vs_orig.ident_X6"
if(step == 1){ # need 32 GB
        opts = list(
                    c("AFT12_1_C5", "AFT12_1_C2"),
                    c("AFT12_1_C5", "AFT12_1_C3"),
                    c("AFT12_1_C2", "AFT12_1_C4"),
                    c("AFT12_1_C5", "AFT12_3_21_C5"),
                    c("AFT12_3_21_C5", "AFT12_18_1_C5"),
                    c("AFT12_18_1_C5", "AFT12_18_8_C5"),
                    c("AFT12_18_8_C5", "AFT12_week3_C5"),
                    c("AFT12_1_C2", "AFT12_3_21_C2"),
                    c("AFT12_3_21_C2", "AFT12_18_1_C2"),
                    c("AFT12_18_1_C2", "AFT12_18_8_C2"),
                    c("AFT12_18_8_C2", "AFT12_week3_C2"),
                    c("AFT12_18_8_C2", "Pt25_25_C2"),
                    c("AFT12_18_8_C3", "Pt25_25_C3")
                    )#13
        print(opt <- opts[[args]])
        object$orig.ident_X6 = paste0(object$orig.ident,"_C",object$X6cluster)

        object %<>% subset(subset = orig.ident_X6 %in% opt)

        markers = FindMarkers_UMI(object,
                                  ident.1 = opt[1],
                                  ident.2 = opt[2],
                                  group.by = "orig.ident_X6",
                                  assay = "SCT",
                                  min.pct = 0.1,
                                  logfc.threshold = 0.1,
                                  only.pos = F)
        markers$gene = rownames(markers)
        markers$cluster = paste(opt,collapse = "_vs_")

        arg = args
        if(args < 10) arg = paste0("0",arg)
        if(args < 100) arg = paste0("0",arg)

        # remove MT-
        MT <- grepl("^MT-",markers$gene)
        if(any(MT)) markers = markers[!MT,]

        write.csv(markers,paste0(save_path,arg,"-",opt$ident,"-",num,".",opt$type, ".csv"))


        g <- VolcanoPlots(markers, cut_off_value = 0.05, cut_off = "p_val", cut_off_logFC = 0.1,top = 20,
                          alpha=1, size=4,
                          legend.size = 12)+ theme(legend.position="bottom")
        g = g + ggtitle(paste(opt,collapse = "_vs_"))
        g = g + TitleCenter()#+theme_bw()

        jpeg(paste0(path,"VolcanoPlots_",opt[1],"_vs_",opt[2],".jpeg"),
             units="in", width=10, height=10,res=600)
        print(g)
        dev.off()

}
