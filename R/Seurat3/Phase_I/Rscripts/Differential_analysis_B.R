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
object %<>% subset(idents = "Pt2_30Pd", invert = T)
Idents(object) = "Doublets"
object %<>% subset(idents = "Singlet")
Idents(object) = "cell.types"
table(Idents(object))

object$X4clusters %<>% factor(levels=paste0("C",1:4))
markers <- FilterGenes(object,c("CCND1","CD19","CD5","CDK4","RB1","BTK","SOX11"))


step = 2
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
                                        latent.vars = "nFeature_SCT")

        write.csv(MCL_markers,paste0(path,"MCL_B_",opt$ident.1, "-N01.csv"))
}
# choose == "X4clusters"
if(step == 1){ # need 32 GB
        opts = data.frame(only.pos = rep(c(T,  T,   T,   F),  each = 4),
                          logfc =  rep(c(0.25, 0.1, 0.05, 0), each = 4),
                          ident.1 = rep(paste0("C",1:4),      time = 4))

        (opt = opts[i,])
        Idents(object) = "cell.types"
        object <- subset(object, idents= "MCL")
        Idents(object) = "X4clusters"
        system.time(MCL_markers <- FindMarkers.UMI(object,
                                                   ident.1 = as.character(opt$ident.1),
                                                   ident.2 = NULL,
                                                   logfc.threshold = opt$logfc,
                                                   only.pos = opt$only.pos,
                                                   test.use = "MAST",
                                                   latent.vars = "nFeature_SCT"))
        write.csv(MCL_markers,paste0(path,"MCL_only_41-FC",opt$logfc,"_",opt$ident.1,".csv"))
}
# choose == "X4clusters_vs_Normal"
if(step == 2){ # need 32 GB
        opts = data.frame(only.pos = rep(c(T,  T,   T,   F),  each = 5),
                          logfc =  rep(c(0.25, 0.1, 0.05, 0), each = 5),
                          ident.1 = rep(c("B_cells",paste0("C",1:4)),      time = 4))

        (opt = opts[i,])
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
                                                   latent.vars = "nFeature_SCT"))
        write.csv(MCL_markers,paste0(path,"MCL_Normal_41-FC",opt$logfc,"_",opt$ident.1,".csv"))
}
# choose == "X4clusters_vs_B_cells"
if(step == 3){ # need 32 GB
        opts = data.frame(only.pos = rep(c(T,  T,   T,   F),  each = 4),
                          logfc =  rep(c(0.25, 0.1, 0.05, 0), each = 4),
                          ident.1 = rep(paste0("C",1:4),      time = 4))

        (opt = opts[i,])
        object %<>% subset(idents = c("N01","N02","N03"),invert = T)
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
                                                   latent.vars = "nFeature_SCT"))
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
                                                   latent.vars = "nFeature_SCT"))
        write.csv(B_markers,paste0(path,"B_41-FC",opt$logfc,"_",
                                   paste(opt$specimens,collapse = "-"),".csv"))
}
# choose == "orig.ident_X4clusters_vs_Normal"
if(step == 5){ # need 32 GB
        object$orig.ident %<>% gsub("N01|N02|N03","Normal",.)
        object %<>% subset(idents = "N04",invert = T)
        object$orig.ident_X4clusters = paste0(object$orig.ident, "_", object$X4clusters)
        object$orig.ident_X4clusters %<>% gsub("Normal_.*", "Normal", .)
        df = table(object$orig.ident_X4clusters) %>% as.data.frame.table
        keep = df$Var1[df$Freq >= 3] %>% as.character()
        Idents(object) = "orig.ident_X4clusters"
        object %<>% subset(idents = keep)
        orig.ident_X4clusters = keep[-grepl("Normal", keep)]
        print(s <- orig.ident_X4clusters[i]) #127
        object %<>% subset(idents = c("Normal",s))

        save.path = paste0(path,"X4cluster_vs_Normal/",gsub("_.*","",s),"/")
        if(!dir.exists(save.path))dir.create(save.path, recursive = T)
        system.time(B_markers <- FindAllMarkers.UMI(object,
                                                    logfc.threshold = 0,
                                                    only.pos = T,
                                                    test.use = "MAST",
                                                    return.thresh = 1,
                                                    latent.vars = "nFeature_SCT"))
        write.csv(B_markers, paste0(save.path,"DE_FC0_",s,"-Normal",".csv"))
        normal = B_markers$cluster %in% "Normal"
        B_markers$avg_logFC[normal] = -1*B_markers$avg_logFC[normal]
        # remove MT-
        MT <- grepl("^MT-",B_markers$gene)
        if(any(MT)) B_markers = B_markers[!MT,]

        g <- VolcanoPlots(B_markers, cut_off_value = 0.05, cut_off = "p_val", cut_off_logFC = 0.1,top = 20,
                          cols = c("#2a52be","#d2dae2","#d9321f"),alpha=1, size=2,
                          legend.size = 12)+ theme(legend.position="bottom")
        g = g + ggtitle(paste(s, "/ Normal in B and MCL"))
        g = g + TitleCenter()#+theme_bw()

        jpeg(paste0(save.path,"VolcanoPlots_",s,"-Normal",".jpeg"), units="in", width=10, height=10,res=600)
        print(g)
        dev.off()
}
# choose == "orig.ident_X4clusters_vs_orig.ident_X4clusters"
if(step == 6){ # need 32 GB
        opts = list(c("Pt11_LN1","Pt11_1"),
                    c("Pt11_LN1","Pt11_14"),
                    c("Pt11_LN1","Pt11_28"),
                    c("Pt11_1","Pt11_14"),
                    c("Pt11_1","Pt11_28"),
                    c("Pt11_14","Pt11_28"),
                    c("Pt17_LN1","Pt17_2"),
                    c("Pt17_LN1","Pt17_7"),
                    c("Pt17_LN1","Pt17_12"),
                    c("Pt17_2","Pt17_7"),
                    c("Pt17_2","Pt17_12"),
                    c("Pt17_7","Pt17_12"),
                    c("Pt25_SB1","Pt25_AMB25Pd"),
                    c("Pt25_1","Pt25_1_8"),
                    c("Pt25_1","Pt25_24"),
                    c("Pt25_1","Pt25_25Pd"),
                    c("Pt25_1","Pt25_AMB25Pd"),
                    c("Pt25_1_8","Pt25_24"),
                    c("Pt25_1_8","Pt25_25Pd"),
                    c("Pt25_1_8","Pt25_AMB25Pd"),
                    c("Pt25_24","Pt25_25Pd"),
                    c("Pt25_24","Pt25_AMB25Pd"),
                    c("Pt25_25Pd","Pt25_AMB25Pd"),
                    c("Pt27_1","Pt27_1_8"),
                    c("Pt27_1","Pt27_12"),
                    c("Pt27_1_8","Pt27_12"),
                    c("PtB13_Ibp","PtB13_Ib1"),
                    c("PtB13_Ibp","PtB13_IbR"),
                    c("PtB13_Ib1","PtB13_IbR"),
                    c("Pt28_LN1","Pt28_1"),
                    c("Pt28_1","Pt28_4"),
                    c("Pt28_1","Pt28_28"),
                    c("Pt28_4","Pt28_28")
                    ) #33
        print(opt <- opts[[i]])
        object %<>% subset(idents = opt)
        object$orig.ident_X4clusters = paste0(object$orig.ident,"_",object$X4clusters)
        # remove cluster with less than 3 cells======
        table_subset.MCL <- table(object$orig.ident_X4clusters) %>% as.data.frame.table
        (keep.MCL <- table_subset.MCL[table_subset.MCL$Freq > 2,"Var1"] %>% as.character())
        (X4_cluster <- keep.MCL %>% unique %>%
                        gsub('.*_C',"",.) %>% as.numeric %>% sort %>% .[duplicated(.)])

        print(ident.2 <- paste(opt[1], X4_cluster, sep="_C"))
        print(ident.1 <- paste(opt[2], X4_cluster, sep="_C"))

        Idents(object) = "orig.ident_X4clusters"
        object %<>% subset(idents = c(ident.1, ident.2))

        save.path = paste0(path,"X4clusters_vs_X4clusters/",gsub("_.*","",opt[1]),"/")
        if(!dir.exists(save.path))dir.create(save.path, recursive = T)
        system.time(B_markers <- FindPairMarkers(object,
                                                 ident.1 = ident.1,
                                                 ident.2 = ident.2,
                                                 logfc.threshold = 0,
                                                 only.pos = F,
                                                 test.use = "MAST",
                                                 return.thresh = 1,
                                                 latent.vars = "nFeature_SCT"))
        write.csv(B_markers, paste0(save.path,"DE_FC0_",opt[2],"_",opt[1],".csv"))
        # remove MT-
        MT <- grepl("^MT-",B_markers$gene)
        if(any(MT)) B_markers = B_markers[!MT,]

        clusters = unique(B_markers$cluster1.vs.cluster2)
        for(c in clusters){
                cluster_markers  <- B_markers[B_markers$cluster1.vs.cluster2 %in% c,]
                g <- VolcanoPlots(cluster_markers, cut_off_value = 0.05, cut_off = "p_val", cut_off_logFC = 0.1,top = 20,
                                  cols = c("#2a52be","#d2dae2","#d9321f"),alpha=1, size=2,
                                  legend.size = 12)+ theme(legend.position="bottom")
                g = g + ggtitle(paste(clusters,"in B and MCL"))
                g = g + TitleCenter()#+theme_bw()

                jpeg(paste0(save.path,"VolcanoPlots_",opt[2],"_",opt[1],"_",gsub(".*_C","C",c),".jpeg"),
                     units="in", width=10, height=10,res=600)
                print(g)
                dev.off()
        }
}

# Doheatmap for MCL longitudinal X4 clusters ================
if(step == 7){
        opts = data.frame(group = rep(c("Untreated","Pt11","Pt17","Pt25","Pt27","PtB13"),  each = 4),
                          cluster = rep(paste0("C",1:4),6),
                          stringsAsFactors = F)
        print(opt <- opts[i,])

        Idents(object) = "groups"
        subset.MCL <- subset(object, idents = opt$group)

        Idents(subset.MCL) = "X4clusters"
        subset.MCL %<>% subset(idents = opt$cluster)

        # remove low cell sample
        subset.MCL$orig.ident %<>% droplevels
        df <- as.data.frame(table(subset.MCL$orig.ident) )
        if(any(df$Freq < 5)) subset.MCL %<>% subset(idents = df[df$Freq >= 5,"Var1"])
        Idents(subset.MCL) = "orig.ident"

        gde.markers <- FindAllMarkers.UMI(subset.MCL,
                                       logfc.threshold = 0.05,
                                       only.pos = T,
                                       test.use = "MAST")
        write.csv(gde.markers,paste0(path,opt$group,"_",opt$cluster,
                                     "_FC0.05_markers.csv"))
        (mito.genes <- grep(pattern = "^MT-", x = gde.markers$gene))
        if(length(mito.genes)>0) gde.markers = gde.markers[-mito.genes,]
        GC()
        #DoHeatmap.1======
        Top_n = 40
        top = gde.markers %>% group_by(cluster) %>% top_n(Top_n, avg_logFC)

        features = c(as.character(top$gene),
                     tail(VariableFeatures(object = subset.MCL), 2),
                     markers)
        subset.MCL %<>% ScaleData(features=features)
        featuresNum <- make.unique(features, sep = ".")
        subset.MCL %<>% MakeUniqueGenes(features = features)

        DoHeatmap.1(object =subset.MCL, features = featuresNum, Top_n = Top_n,
                    do.print=T, angle = 0, group.bar = F, title.size = 20, no.legend = F,size=5,hjust = 0.5,
                    group.bar.height = 0, label=F, cex.row= ifelse(i==2,4,2), legend.size = 0,width=10, height=6.5,
                    pal_gsea = FALSE,
                    title = paste("Top",Top_n,"DE genes in longitudinal",opt$group,
                                  "B/MCL cells cluster",opt$cluster))

        # rename file
        v <- UniqueName(object = subset.MCL, fileName = "subset.MCL",unique.name = T)
        v = paste0(v,"_",FindIdentLabel(object))
        old.name = paste0(path,"Heatmap_top",Top_n,"_",v,"_Legend.jpeg")
        file.rename(old.name, paste0(path,"Heatmap_top",Top_n,"_",opt$group,
                                     "_Cluster",opt$cluster, ".jpeg"))
        DoHeatmap.1(object =subset.MCL, features = featuresNum, Top_n = Top_n,
                    do.print=T, angle = 45, group.bar = T, title.size = 20, no.legend = F,size=5,hjust = 0.5,
                    group.bar.height = 0.05, label=T, cex.row= ifelse(i==2,4,2), legend.size = 0,width=10, height=6.5,
                    pal_gsea = FALSE,
                    title = paste("Top",Top_n,"DE genes in longitudinal",opt$group,
                                  "B/MCL cells cluster",opt$cluster))

}

# Dec 4, 2020
# 2. Volcano plots comparing a) cluster 1 MCL cells with cluster 1 normal B cells,  and b) cluster 2 MCL cells with cluster 2 B cells.

if(step == 8){
        path <- "Yang/Figure Sources/MCL_vs_B_inX4Clusters/"
        if(!dir.exists(path))dir.create(path, recursive = T)

        Idents(object) = "X4clusters"
        X4Clusters = c("C1","C2","C3","C4")
        X4Cluster = X4Clusters[i]
        save.path <- paste0(path,X4Cluster,"/")
        if(!dir.exists(save.path))dir.create(save.path, recursive = T)
        sub_object <- subset(object,idents = X4Cluster)
        Idents(sub_object) = "cell.types"

        MCL_markers <- FindAllMarkers.UMI(sub_object,
                                          logfc.threshold = 0,
                                          return.thresh = 1,
                                          only.pos = F,
                                          test.use = "MAST",
                                          latent.vars = "nFeature_SCT")
        write.csv(MCL_markers,paste0(save.path,"MCL_vs_B_within_",X4Cluster,
                                     "_FC0_markers.csv"))
        table(MCL_markers$cluster)
        markers <- FilterGenes(object,c("CCND1","CD19","CD5","CDK4","RB1","BTK","SOX11"))
        (MT_gene <- grep("^MT-",MCL_markers$gene))
        if(length(MT_gene) >0 ) MCL_markers = MCL_markers[-MT_gene,]
        Top_n = 40

        top = MCL_markers %>% group_by(cluster) %>%
                top_n(40, avg_logFC)
        unique(top$cluster)
        #top = top[order(top$cluster),]
        write.csv(top,paste0(save.path,"Top40_","MCL_vs_B_within_",
                             X4Cluster,".csv"))
        features = c(as.character(top$gene),
                     tail(VariableFeatures(object = sub_object), 2),
                     markers)
        #DoHeatmap.1======
        # raw heatmap
        featuresNum <- make.unique(features, sep = ".")
        exp = AverageExpression(sub_object[features,],
                                assays = "SCT") %>% .$SCT
        exp %<>% MakeUniqueGenes(features = features)
        exp[tail(VariableFeatures(object = sub_object), 2),] =0

        (group.by = c("B_cells","MCL"))
        DoHeatmap.matrix(exp, features = featuresNum,
                         group.by = group.by,
                         size = 6,angle = 0,label =F,
                         draw.lines =F, raster = FALSE,
                         pal_gsea = FALSE,
                         width=2.5, height=10,res=600,no.legend = F,
                         cex.row=5,
                         group.colors = gg_color_hue(2),
                         do.print = T,
                         save.path = save.path,
                         file.name = paste0("Heatmap_top40_MCL_vs_B_within_",
                                            X4Cluster,"_raw.jpeg")
        )
        # scale heatmap
        sub_object %<>% ScaleData(features = features)
        scale_exp = AverageExpression(sub_object[features,],
                                      assays = "SCT", slot = "scale.data") %>% .$SCT
        scale_exp %<>% MakeUniqueGenes(features = features)
        scale_exp[tail(VariableFeatures(object = sub_object), 2),] =0
        DoHeatmap.matrix(scale_exp, features = featuresNum,
                         group.by = group.by,
                         size = 6,angle = 0,label =F,
                         draw.lines =F, raster = FALSE,
                         pal_gsea = FALSE,
                         width=2.5, height=10,res=600,no.legend = F,
                         cex.row=5,
                         group.colors = gg_color_hue(2),
                         do.print = T,
                         save.path = save.path,
                         file.name = paste0("Heatmap_top40_MCL_vs_B_within_",
                                            X4Cluster,"_scaled.jpeg")
        )
        MCL_markers = MCL_markers[MCL_markers$cluster %in% "MCL",]
        #avg_logFC = MCL_markers[MCL_markers$cluster %in% opt$ident.2,"avg_logFC"]
        #MCL_markers[MCL_markers$cluster %in% opt$ident.2,"avg_logFC"] = avg_logFC * -1
        p <- VolcanoPlots(data = MCL_markers, cut_off_value = 0.05, cut_off = "p_val", cut_off_logFC = 0.1,
                          top = 20, cols = c("#2a52be","#d2dae2","#d9321f"),alpha=1, size=2,
                          legend.size = 12)+
                ggtitle(paste("MCL vs B cells within",X4Cluster,"Cluster"))+
                theme(plot.title = element_text(hjust = 0.5,size=15,face = "plain"),
                      legend.position="bottom")
        jpeg(paste0(save.path,"VolcanoPlots_MCL_vs_B_within_",
                    X4Cluster,".jpeg"), units="in", width=10, height=7,res=600)
        print(p)
        dev.off()
}
