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

object = readRDS(file = "data/MCL_41_B_20200225.rds")
DefaultAssay(object) = "SCT"
Idents(object) = "orig.ident"

step = 5
# choose == "MCL_vs_B_cells"
if(step == 0){  # need 32 GB
        # load data
        samples = as.character(unique(object$orig.ident))
        opts = data.frame(ident.1 = samples[2:length(samples)],
                          ident.2 = rep("N01", length(samples)-1),
                          stringsAsFactors = F)
        (opt = opts[i,])
        object %<>% subset(idents = c(opt$ident.1,"N01"))
        
        MCL_markers <- FindAllMarkers.UMI(object, 
                                        logfc.threshold = 0, 
                                        only.pos = F,
                                        return.thresh = 1,
                                        test.use = "MAST",
                                        latent.vars = "nCount_SCT")
        
        write.csv(MCL_markers,paste0(path,"MCL_B_",opt$ident.1, "-N01.csv"))
}
# choose == "X4clusters"
if(step == 1){ # need 32 GB
        opts = data.frame(only.pos = rep(c(T,  T,   T,   F),  each = 4),
                          logfc =  rep(c(0.25, 0.1, 0.05, 0), each = 4),
                          ident.1 = rep(paste0("C",1:4),      time = 4))
        
        (opt = opts[i,])
        object %<>% subset(idents = "Pt2_30Pd",invert = T)
        Idents(object) = "cell.types"
        object <- subset(object, idents= "MCL") 
        Idents(object) = "X4clusters"
        system.time(MCL_markers <- FindMarkers.UMI(object, 
                                                   ident.1 = as.character(opt$ident.1),
                                                   ident.2 = NULL,
                                                   logfc.threshold = opt$logfc, 
                                                   only.pos = opt$only.pos,
                                                   test.use = "MAST",
                                                   latent.vars = "nCount_SCT"))
        write.csv(MCL_markers,paste0(path,"MCL_only_41-FC",opt$logfc,"_",opt$ident.1,".csv"))
}
# choose == "X4clusters_vs_Normal"
if(step == 2){ # need 32 GB
        opts = data.frame(only.pos = rep(c(T,  T,   T,   F),  each = 5),
                          logfc =  rep(c(0.25, 0.1, 0.05, 0), each = 5),
                          ident.1 = rep(c("B_cells",paste0("C",1:4)),      time = 4))
        
        (opt = opts[i,])
        object %<>% subset(idents = "Pt2_30Pd",invert = T)
        object$X4clusters_normal = as.character(object$X4clusters)
        object$X4clusters_normal %<>% paste(object$cell.types, sep = "_")
        object$X4clusters_normal %<>% gsub(".*_B_cells","B_cells",.)
        object$X4clusters_normal %<>% gsub("_MCL","",.)
        normal <- grepl("N01|N02|N03",object$orig.ident)
        
        object@meta.data[normal,"X4clusters_normal"] = "Normal"
        Idents(object) = "X4clusters_normal"
        object %<>% sortIdent()
        table(Idents(object))
        system.time(MCL_markers <- FindMarkers.UMI(object, 
                                                   ident.1 = as.character(opt$ident.1),
                                                   ident.2 = "Normal",
                                                   logfc.threshold = opt$logfc, 
                                                   only.pos = opt$only.pos,
                                                   test.use = "MAST",
                                                   latent.vars = "nCount_SCT"))
        write.csv(MCL_markers,paste0(path,"MCL_Normal_41-FC",opt$logfc,"_",opt$ident.1,".csv"))
}
# choose == "X4clusters_vs_B_cells"
if(step == 3){ # need 32 GB
        opts = data.frame(only.pos = rep(c(T,  T,   T,   F),  each = 4),
                          logfc =  rep(c(0.25, 0.1, 0.05, 0), each = 4),
                          ident.1 = rep(paste0("C",1:4),      time = 4))
        
        (opt = opts[i,])
        object %<>% subset(idents = c("Pt2_30Pd","N01","N02","N03"),invert = T)
        object$X4clusters_B = as.character(object$X4clusters)
        object$X4clusters_B %<>% paste(object$cell.types, sep = "_")
        object$X4clusters_B %<>% gsub(".*_B_cells","B_cells",.)
        object$X4clusters_B %<>% gsub("_MCL","",.)
        Idents(object) = "X4clusters_B"
        object %<>% sortIdent()
        table(Idents(object))
        system.time(MCL_markers <- FindMarkers.UMI(object, 
                                                   ident.1 = as.character(opt$ident.1),
                                                   ident.2 = "B_cells",
                                                   logfc.threshold = opt$logfc, 
                                                   only.pos = opt$only.pos,
                                                   test.use = "MAST",
                                                   latent.vars = "nCount_SCT"))
        write.csv(MCL_markers,paste0(path,"MCL_B_41-FC",opt$logfc,"_",opt$ident.1,".csv"))
}

# choose == "B_cells_vs_B_cells"
if(step == 4){ # need 32 GB
        opts = list(list(FALSE, 0, "N01"),
                   list(FALSE, 0, "N02"),
                   list(FALSE, 0, "N03"),
                   list(FALSE, 0, "N04"),
                   list(FALSE, 0, c("N01","N02","N03","N04")), 
                   list(FALSE, 0, "Pt25_1"),
                   list(FALSE, 0, "Pt25_24"),
                   list(FALSE, 0, "Pt25_25Pd"),
                   list(FALSE, 0, c("Pt25_24", "Pt25_1")))
        
        (opt = opts[[i]])
        names(opt) = c("only.pos","logfc","specimens")
        object %<>% subset(idents = opt$specimens)
        Idents(object) = "cell.types"
        object %<>% subset(idents = "B_cells")
        Idents(object) = "X4clusters"
        table(Idents(object))
        if(!identical(opt$specimens, c("Pt25_24", "Pt25_1"))) {
                ident.1 = "C1"
                ident.2 = "C2"
                object %<>% subset(idents = c(ident.1,ident.2))
        } else {
                object$X4clusters_orig.ident = paste0(object$X4clusters,"_",
                                                      object$orig.ident)
                Idents(object) = "X4clusters_orig.ident"
                ident.1 ="C2_Pt25_24"
                ident.2 = "C2_Pt25_1"
                object %<>% subset(idents = c(ident.1,ident.2))
        }
        system.time(B_markers <- FindAllMarkers.UMI(object, 
                                                   logfc.threshold = opt$logfc, 
                                                   only.pos = opt$only.pos,
                                                   test.use = "MAST",
                                                   return.thresh = 1, 
                                                   latent.vars = "nCount_SCT"))
        write.csv(B_markers,paste0(path,"B_41-FC",opt$logfc,"_",
                                   paste(opt$specimens,collapse = "-"),".csv"))
}
# choose == "orig.ident_X4clusters_vs_Normal"
if(step == 5){ # need 32 GB
        object$orig.ident %<>% gsub("N02|N01|N03","Normal",.)
        object %<>% subset(idents = "N04", invert = T)
        object$orig.ident_X4clusters = paste0(object$orig.ident,"_",object$X4clusters)
        object$orig.ident_X4clusters %<>% gsub("Normal_.*","Normal",.)
        df = table(object$orig.ident_X4clusters) %>% as.data.frame.table 
        keep = df$Var1[df$Freq >= 3] %>% as.character()
        Idents(object) = "orig.ident_X4clusters"
        object %<>% subset(idents = keep)
        orig.ident_X4clusters = keep[-grepl("Normal",keep)]
        s = orig.ident_X4clusters[i]
        object %<>% subset(idents = c("Normal",s))
        
        save.path = paste0(path,"X4clusters_vs_Normal/",gsub("_.*","",s),"/")
        if(!dir.exists(save.path))dir.create(save.path, recursive = T)
        system.time(B_markers <- FindAllMarkers.UMI(object, 
                                                    logfc.threshold = 0, 
                                                    only.pos = T,
                                                    test.use = "MAST",
                                                    return.thresh = 1, 
                                                    latent.vars = "nCount_SCT"))
        write.csv(B_markers, paste0(save.path,"DE_FC0_",s,"-Normal",".csv"))
        normal = B_markers$cluster %in% "Normal"
        B_markers$avg_logFC[normal] = -1*B_markers$avg_logFC[normal]
        g <- VolcanoPlots(B_markers, cut_off_value = 0.05, cut_off = "p_val", cut_off_logFC = 0.1,top = 20,
                          cols = c("#2a52be","#d2dae2","#d9321f"),alpha=1, size=2,
                          legend.size = 12)+ theme(legend.position="bottom")
        g = g + ggtitle(paste(s, "/ Normal in B and MCL"))
        g = g + TitleCenter()#+theme_bw()
        
        jpeg(paste0(save.path,"VolcanoPlots_",s,"-Normal",".jpeg"), units="in", width=10, height=10,res=600)
        print(g)
        dev.off()
}