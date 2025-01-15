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
i <- as.numeric(slurm_arrayid)
print(paste0("slurm_arrayid=",i))

object = readRDS(file = "data/MCL_51_20210724.rds")
DefaultAssay(object) = "SCT"
meta.data <- readRDS(file = "output/MCL_B_51_20230113_meta.data_v2.rds")
if(all(rownames(meta.data) %in% colnames(object))){
    print("all cellID within!")
    object %<>% subset(cells = rownames(meta.data))
    object@meta.data = meta.data
}

object = subset(object, subset =  Doublets == "Singlet"
                                & X4cluster  %in% c("1","2","3","4")
)


step = c("MCL_vs_B_cells","X4clusters","X4clusters_vs_Normal","X4clusters_vs_B_cells",
         "orig.ident_X4clusters_vs_orig.ident_X4clusters","longitudinal X4 clusters","MCL_vs_B_inX4Clusters",
         "response_X4clusters_vs_rest in CR and PD","X4clusters_v1,2,3_vs_Normal")[3]
if(step == "MCL_vs_B_cells"){  # need 32 GB
        # load data
        samples = as.character(unique(object$orig.ident))
        samples = samples[!samples %in% c("N01","N02","N03","N04")]
        print(ident.1 <- samples[i])

        Idents(object) = "orig.ident"
        object %<>% subset(idents = c(ident.1,"N01"))

        MCL_markers <- FindAllMarkers_UMI(object,
                                        logfc.threshold = 0,
                                        only.pos = F,
                                        return.thresh = 1,
                                        test.use = "MAST",
                                        latent.vars = "nFeature_SCT")

        write.csv(MCL_markers,paste0(path,"MCL_B_",ident.1, "-vs-N01.csv"))
}
if(step == "X4clusters"){ # need 32 GB

        #object = subset(object, subset =  orig.ident %in% c("Pt3_BMA72_6","Pt3_72_6"), invert = T)
        object = subset(object, subset =  orig.ident %in% c("N01","N02","N03","N04"), invert = T)

        Idents(object) = "X4cluster"
        table(Idents(object))
        object$X4cluster %<>% factor(levels=c("1","2","3","4"))

        system.time(MCL_markers <- FindMarkers_UMI(object,
                                                   ident.1 = as.character(i),
                                                   ident.2 = NULL,
                                                   return.thresh = 1,
                                                   logfc.threshold = 0,
                                                   only.pos = FALSE,
                                                   test.use = "MAST"))
        write.csv(MCL_markers,paste0(path,"MCL_only_51-FC0_cluster",i,".csv"))
}

if(step == "X4clusters_vs_Normal"){ # need 32 GB
        #object = subset(object, subset =  orig.ident %in% c("Pt3_BMA72_6","Pt3_72_6"), invert = T)
        object$X4clusters_normal = as.character(object$X4cluster)
        object@meta.data[object$orig.ident %in% c("N01","N02","N03"),
                         "X4clusters_normal"] = "Normal"
        object@meta.data[(!object$X4clusters_normal %in% "Normal") & object$cell.types == "B_cells",
                         "X4clusters_normal"] = "B cells"
        object %<>% subset(subset =  orig.ident == "N04" | X4clusters_normal == "B cells", invert = T)
        Idents(object) = "X4clusters_normal"
        table(Idents(object))
        opts = as.character(1:4)
        print(paste0(opts[i]," vs. Normal"))
        system.time(markers <- FindMarkers_UMI(object,
                                                   ident.1 = opts[i],
                                                   ident.2 = "Normal",
                                                   logfc.threshold = 0.1,
                                                   only.pos = FALSE,
                                                   test.use = "MAST",
                                                   latent.vars = "nFeature_SCT"))
        markers$gene <- rownames(markers)
        markers$cluster = paste0(opts[i]," vs. Normal")

        write.csv(markers,paste0(path,i,"_MCL_Normal_51-FC0.1_",opt$X4cluster,"_vs_Normal.csv"))
}

if(step == "X4clusters_vs_B_cells"){ # need 32 GB
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
        system.time(MCL_markers <- FindMarkers_UMI(object,
                                                   ident.1 = as.character(opt$ident.1),
                                                   ident.2 = "B_cells",
                                                   logfc.threshold = opt$logfc,
                                                   only.pos = opt$only.pos,
                                                   test.use = "MAST",
                                                   latent.vars = "nFeature_SCT"))
        write.csv(MCL_markers,paste0(path,"MCL_B_41-FC",opt$logfc,"_",opt$ident.1,".csv"))
}

if(step == "B_cells_vs_B_cells"){ # need 32 GB
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
        system.time(B_markers <- FindMarkers_UMI(object,
                                                   logfc.threshold = opt$logfc,
                                                   only.pos = opt$only.pos,
                                                   test.use = "MAST",
                                                   return.thresh = 1,
                                                   latent.vars = "nFeature_SCT"))
        write.csv(B_markers,paste0(path,"B_41-FC",opt$logfc,"_",
                                   paste(opt$specimens,collapse = "-"),".csv"))
}

if(step == "orig.ident_X4clusters_vs_Normal"){ # need 32 GB
        object$orig.ident %<>% gsub("N01|N02|N03","Normal",.)
        Idents(object) = "orig.ident"
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
        system.time(B_markers <- FindMarkers_UMI(object,
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

if(step == "orig.ident_X4clusters_vs_orig.ident_X4clusters"){ # need 32 GB
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
if(step == "longitudinal X4 clusters"){
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

        gde.markers <- FindMarkers_UMI(subset.MCL,
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

if(step == 'MCL_vs_B_inX4Clusters'){
        path <- "Yang/Figure Sources/MCL_vs_B_inX4Clusters/"
        if(!dir.exists(path))dir.create(path, recursive = T)

        Idents(object) = "X4clusters"
        X4Clusters = c("C1","C2","C3","C4")
        X4Cluster = X4Clusters[i]
        save.path <- paste0(path,X4Cluster,"/")
        if(!dir.exists(save.path))dir.create(save.path, recursive = T)
        sub_object <- subset(object,idents = X4Cluster)
        Idents(sub_object) = "cell.types"

        MCL_markers <- FindMarkers_UMI(sub_object,
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

# choose == "response_X4clusters_vs_rest in CR and PD"
if(step == "response_X4clusters_vs_rest in CR and PD"){ # need 32 GB
        object = subset(object, subset =  response %in% c("CR","PD"))
        object$response %<>% droplevels()
        object$response_X4clusters = paste0(object$response,"_",object$X4cluster)

        # remove cluster with less than 3 cells======
        opts = as.data.frame(table(object$response_X4clusters)) %>% filter(Freq >= 3) #8
        print(opt <- as.character(opts[i,1]))

        save.path = paste0(path,"response_X4clusters_vs_rest/")
        if(!dir.exists(save.path))dir.create(save.path, recursive = T)
        Idents(object) = "response_X4clusters"
        system.time(B_markers <- FindMarkers_UMI(object,
                                                 ident.1 = opt,
                                                 logfc.threshold = 0.1,
                                                 only.pos = F,
                                                 test.use = "MAST",
                                                 latent.vars = "nFeature_SCT"))
        B_markers$cluster = opt
        B_markers$gene = rownames(B_markers)
        write.csv(B_markers, paste0(save.path,opt,"_FC0.1_.csv"))
}

if(step == "X4clusters_v1,2,3_vs_Normal"){ # need 32 GB
    object %<>% subset(subset =  orig.ident %in% c("N04"), invert = T)
    object %<>% subset(subset =  cell.types == "MCL")
    object$X4cluster_normal = as.character(object$X4cluster)
    object$X4cluster_v2_normal = as.character(object$X4cluster_v2)
    object$X4cluster_v3_normal = as.character(object$X4cluster_v3)
    normal <- object$orig.ident %in% c("N01","N02","N03")
    object@meta.data[normal,"X4cluster_normal"] = "Normal"
    object@meta.data[normal,"X4cluster_v2_normal"] = "Normal"
    object@meta.data[normal,"X4cluster_v3_normal"] = "Normal"

    opts <- data.frame("X4cluster" = rep(as.character(1:4),3),
                       "ident" = paste0(rep(c("X4cluster","X4cluster_v2","X4cluster_v3"),each = 4),"_normal")) #1~12
    opt <- opts[i,]
    print(opt)
    Idents(object) = opt$ident
    table(Idents(object))

    print(paste(opt$X4cluster," vs. Normal in",opt$ident))
    system.time(markers <- FindMarkers_UMI(object,
                                           ident.1 = opt$X4cluster,
                                           ident.2 = "Normal",
                                           logfc.threshold = 0.1,
                                           only.pos = TRUE,
                                           test.use = "MAST",
                                           latent.vars = "nFeature_SCT"))
    markers$gene <- rownames(markers)
    markers$cluster = opt$X4cluster
    markers$ident = opt$ident

    num = i
    if(i < 10) num = paste0("0",num)

    write.csv(markers,paste0(path,num,"_MCL_Normal_51-FC0.1_",opt$X4cluster,"_vs_Normal_in_",opt$ident,".csv"))

}
