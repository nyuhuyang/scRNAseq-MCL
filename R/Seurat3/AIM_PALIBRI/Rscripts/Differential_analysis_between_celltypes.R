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

object$label = object$label.fine
object$label %<>% gsub("T cells, CD4\\+.*","CD4 T",.)
object$label %<>% gsub("T cells, CD8\\+.*","CD8 T",.)
object$label %<>% gsub("Monocytes, CD14\\+.*","CD14 Monocytes",.)
object$label %<>% gsub("Monocytes, CD16\\+.*","CD16 Monocytes",.)
object$label %<>% gsub("B cells,.*","B cells",.)

opts = data.frame(label = c("CD14 Monocytes","CD16 Monocytes"),
                  stringsAsFactors = F)
set.seed(101)
print(label <- opts$label)
file.name = "CD16_vs_CD14_monocytes"

save.path = paste0(path,file.name,"/")
if(!dir.exists(save.path))dir.create(save.path, recursive = T)

# subset

Idents(object) = "label"
object %<>% subset(idents = label)
GC()
run_DE = TRUE

if(run_DE) {
    gde.markers <- FindAllMarkers.UMI(object,
                                      logfc.threshold = 0,
                                      only.pos = F,
                                      return.thresh = 1,
                                      test.use = "MAST",
                                      latent.vars = "nFeature_SCT")
    file.path = paste0(save.path,file.name,"F0_DEGs.csv")
    write.csv(gde.markers,file = file.path)
} else gde.markers = read.csv(file = file.path,
                              row.names = 1, stringsAsFactors = F)


#DoHeatmap.1======
(mito.genes <- grep(pattern = "^MT-", x = gde.markers$gene))
if(length(mito.genes)>0) gde.markers = gde.markers[-mito.genes,]
Top_n = 40
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
            title = paste("Top",Top_n,"DE genes between CD14 monocytes and CD16 monocytes "),
            save.path = save.path,
            file.name = paste0("Heatmap_",file.name,"_noLabel.jpeg"))

DoHeatmap.1(object =object, features = featuresNum, Top_n = Top_n,
            do.print=T, angle = 45, group.bar = T, title.size = 20, no.legend = F,size=5,hjust = 0.5,
            group.bar.height = 0.05, label=T, cex.row= 4, legend.size = 0,width=10, height=6.5,
            pal_gsea = FALSE,
            title = paste("Top",Top_n,"DE genes between CD14 monocytes and CD16 monocytes "),
            save.path = save.path,
            file.name = paste0("Heatmap_",file.name,"_Label.jpeg"))
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
                        padj = 1,pval = 1,
                        order.yaxis.by = c("CD16 Monocytes","NES"),
                        decreasing = F,
                        order.xaxis = c("CD14 Monocytes","CD16 Monocytes"),
                        Rowv = F,
                        Colv = F,
                        size = " -log10(pval)", fill = "NES",
                        pathway.name = "Hallmark",rotate.x.text = T,
                        title = file.name,
                        font.xtickslab=12, font.main=12, font.ytickslab = 10,
                        font.legend = list(size = 12),font.label = list(size = 12),
                        do.return = T,
                        save.path = save.path,
                        do.print = T,
                        width = 7,height = 6,hjust = 0.75)
write.csv(fgseaRes, file = paste0(save.path,"FgseaDotPlot_FDR0.25_pval0.05.csv"))
# VolcanoPlots
gde.markers = gde.markers[gde.markers$cluster %in% "CD16 Monocytes",]
g <- VolcanoPlots(gde.markers, cut_off_value = 0.05, cut_off = "p_val", cut_off_logFC = 0,top = 20,
                  cols = c("#2a52be","#d2dae2","#d9321f"),alpha=1, size=2,
                  legend.size = 12)+ theme(legend.position="bottom")
g = g + ggtitle("CD16 Monocytes vs CD14 Monocytes")
g = g + TitleCenter()+theme_bw()

jpeg(paste0(save.path,"VolcanoPlots_",file.name,".jpeg"), units="in", width=10, height=10,res=600)
print(g)
dev.off()
