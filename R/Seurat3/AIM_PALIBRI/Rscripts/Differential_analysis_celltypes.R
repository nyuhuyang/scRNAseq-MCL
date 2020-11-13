########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
#
# ######################################################################
invisible(lapply(c("Seurat","dplyr","cowplot","fgsea",
                   "magrittr","data.table","future","ggplot2","tidyr"), function(x) {
                       suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
# change the current plan to access parallelization
plan("multiprocess", workers = 4)
plan()

# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.integer(as.character(slurm_arrayid))
print(paste0("slurm_arrayid=",args))

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
# Need 64GB
# load files
(load(file = "data/MCL_AIM_58_20201009.Rda"))

# Need 64GB
DefaultAssay(object) = "SCT"
orig.ident = levels(object$orig.ident)
object$orig.ident %<>% gsub("N01|N02|N03|N04","Normal",.)
orig.ident %<>% gsub("N01|N02|N03|N04","Normal",.)
object$orig.ident %<>% as.factor %>% factor(levels = unique(orig.ident))

object$label = object$label.fine
object$label %<>% gsub("T cells, CD4\\+.*","CD4 T",.)
object$label %<>% gsub("T cells, CD8\\+.*","CD8 T",.)
object$label %<>% gsub("Monocytes, CD14\\+.*","CD14 Monocytes",.)
object$label %<>% gsub("Monocytes, CD16\\+.*","CD16 Monocytes",.)
object$label %<>% gsub("B cells,.*","B cells",.)

object$patient = gsub("_.*","",object$orig.ident)
#Idents(object) = "Doublets"
#object <- subset(object, idents = "Singlet")

opts = data.frame(label = c(rep("CD4 T",10),
                            rep("CD8 T",10),
                            rep("CD14 Monocytes",10),
                            rep("CD16 Monocytes",10),
                            rep("NK cells",10),
                            rep("B cells", 10),
                            rep("MCL", 10)),
                  patient = rep(c("Pt11","Pt13","Pt17","Pt25","Pt27","Pt28","PtB13",
                              "AIM13","AIM17","AIM24"),times = 7),
                  stringsAsFactors = F)
set.seed(101)
print(opt <- opts[args,])
label = opt$label
patient = opt$patient

save.path = paste0(path,label,"/")
if(!dir.exists(save.path))dir.create(save.path, recursive = T)

# subset
Normal_cells =  colnames(object)[object$orig.ident %in% "Normal"]
remove_normal_cells = sample(Normal_cells, size = 3*round(length(Normal_cells)/4))
object %<>% subset(cells = remove_normal_cells, invert = T)
Idents(object) = "patient"
object %<>% subset(idents = c("Normal",patient))
Idents(object) = "label"
object %<>% subset(idents = label)
object$orig.ident %<>% droplevels()

Idents(object) = "orig.ident"
orig.ident.num = as.data.frame.vector(table(object$orig.ident))
if(any(orig.ident.num[,1] < 4)) {
    keep = rownames(orig.ident.num)[orig.ident.num[,1] >= 4]
    object %<>% subset(idents = keep)
}
GC()
run_DE = TRUE

if(!dir.exists(paste0(save.path,"DEGs")))dir.create(paste0(save.path,"DEGs"), recursive = T)
if(run_DE) {
    gde.markers <- FindAllMarkers.UMI(object,
                                      logfc.threshold = 0,
                                      only.pos = F,
                                      return.thresh = 1,
                                      test.use = "MAST",
                                      latent.vars = "nFeature_SCT")
    write.csv(gde.markers,paste0(save.path,"DEGs/",patient,"_longitudinal_F0_DEGs.csv"))
} else gde.markers = read.csv(paste0(save.path,"DEGs/",patient,"_longitudinal_F0_DEGs.csv"),
                              row.names = 1, stringsAsFactors = F)


#DoHeatmap.1======
(mito.genes <- grep(pattern = "^MT-", x = gde.markers$gene))
if(length(mito.genes)>0) gde.markers = gde.markers[-mito.genes,]
Top_n = 25
top = gde.markers %>% group_by(cluster) %>% top_n(Top_n, avg_logFC)

markers <- FilterGenes(object,c("CD19","CCND1","CD3D","CD3E","CD4","CD8A",
                                "NCAM1","NKG7","S100A8","S100A9"))

features = c(as.character(top$gene),
             tail(VariableFeatures(object = object), 2),
             markers)
object %<>% ScaleData(features=features)
featuresNum <- make.unique(features, sep = ".")
object %<>% MakeUniqueGenes(features = features)
DoHeatmap.1(object =object, features = featuresNum, Top_n = Top_n,
            do.print=T, angle = 0, group.bar = F, title.size = 20, no.legend = F,size=5,hjust = 0.5,
            group.bar.height = 0, label=F, cex.row= 4, legend.size = 0,width=10, height=6.5,
            pal_gsea = FALSE,
            title = paste("Top",Top_n,"DE genes in longitudinal",patient,label),
            save.path = paste0(save.path,"Heatmaps"),
            file.name = paste0(patient,"_Heatmap_nolabel.jpeg"))

DoHeatmap.1(object =object, features = featuresNum, Top_n = Top_n,
            do.print=T, angle = 45, group.bar = T, title.size = 20, no.legend = F,size=5,hjust = 0.5,
            group.bar.height = 0.05, label=T, cex.row= 4, legend.size = 0,width=10, height=6.5,
            pal_gsea = FALSE,
            title = paste("Top",Top_n,"DE genes in longitudinal",patient,label),
            save.path = paste0(save.path,"Heatmaps"),
            file.name = paste0(patient,"_Heatmap_label.jpeg"))
# fgsea
res = gde.markers
res = res[order(res["p_val_adj"]),]
head(res, 20)
(clusters <- unique(res$cluster))
hallmark <- fgsea::gmtPathways("../seurat_resources/msigdb/h.all.v6.2.symbols.gmt")
names(hallmark) = gsub("HALLMARK_","",names(hallmark))
names(hallmark) = gsub("\\_"," ",names(hallmark))

# Now, run the fgsea algorithm with 1000 permutations:
set.seed(100)
fgseaRes = FgseaDotPlot(stats=res, pathways=hallmark,
                        padj = 0.25,pval = 0.05,
                        order.yaxis.by = c("Normal","NES"),
                        decreasing = F,
                        order.xaxis = rownames(orig.ident.num),
                        Rowv = F,Colv = F,
                        size = " -log10(pval)", fill = "NES",
                        pathway.name = "Hallmark",rotate.x.text = T,
                        title = paste(label, "in", patient),
                        font.xtickslab=12, font.main=12, font.ytickslab = 10,
                        font.legend = list(size = 12),font.label = list(size = 12),
                        do.return = T,save.path = paste0(save.path,"DEGs"),
                        do.print = T,
                        width = 7,height = 6,hjust = 0.75)
write.csv(fgseaRes, file = paste0(save.path,"GSEA/",patient,"_FgseaDotPlot_FDR0.25_pval0.05.csv"))
