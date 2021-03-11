library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(magrittr)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
save.path <- "Yang/Figure 4 Xiangao/"
if(!dir.exists(save.path))dir.create(save.path, recursive = T)
#=============== multiple color in the single Featureplot===================================
B_cells_MCL = readRDS(file = "data/MCL_41_B_20200225.rds")
B_cells_MCL$orig.ident %<>% gsub("N01|N02|N03","Normal",.)
Idents(B_cells_MCL) = "orig.ident"
samples = c("PtB13_Ibp","PtB13_Ib1","PtB13_IbR",
            "Pt25_SB1","Pt25_1","Pt25_1_8","Pt25_24","Pt25_25Pd","Pt25_AMB25Pd",
            "Pt10_LN2Pd","All_samples","Normal",
            "Pt11_LN1", "Pt17_LN1","PtU01","PtU02","PtU03","PtU04")
features.list = lapply(list(c("EZH2", "SUZ12"),
                            c("EZH2", "EED"),
                            c("EZH1", "SUZ12"),
                            c("EZH1", "EED"),
                            c("IRF4", "CYTIP"),
                            c("IRF4", "CD40"),
                            c("IRF4", "MAP3K8"),
                            c("MAP3K8", "RELB"),
                            c("CD40", "MAP3K8"),
                            c("PIK3IP1", "EZH1"),
                            c("HLA-DRA", "E2F1"),
                            c("HLA-DRA", "EZH1"),
                            c("HLA-DRA", "EZH2"),
                            c("EZH2","EZH1"),
                            c("EZH1","CCND1"),
                            c("EZH2","CCND1"),
                            c("EZH1","E2F1"),
                            c("EZH2","E2F1"),
                            c("EZH1","PCNA"),
                            c("EZH2","PCNA"),
                            c("IRF4","PIK3IP1"),
                            c("HLA-DPA1", "MCM7"),
                            c("HLA-DPB1", "EZH2"),
                            c("HLA-DPB1", "EZH1"),
                            c("HLA-A", "MCM7"),
                            c("HLA-B", "MCM7"),
                            c("HLA-DPA1", "CCND1"),
                            c("HLA-A", "CCND1"),
                            c("HLA-B", "CCND1"),
                            c("IRF4", "BCL6"),
                            c("IRF4", "CCND1"),
                            c("BCL6", "CCND1"),
                            c("CCND1","MCM7"),
                            c("BCL6","IRF4"),
                            c("MYC", "BCL6"),
                            c("MEF2B", "BCL6"),
                            c("CDKN1A", "BCL6"),
                            c("NFKB2","MAP3K8"),
                            c("FOXM1","IRF4"),
                            c("RELB","IRF4"),
                            c("NFKB2","RELB"),
                            c("MAP3K8","RELB"),
                            c("NFKB2","IRF4"),
                            c("MYC", "IRF4"),
                            c("POLR2M", "CRBN"),
                            c("POLR2M", "IKZF1"),
                            c("POLR2M", "IKZF3"),
                            c("IKZF1", "CRBN"),
                            c("IKZF3", "CRBN"),
                            c("MBOAT7", "CRBN"),
                            c("MBOAT7", "IKZF1"),
                            c("MBOAT7", "IKZF3")),
                       function(x) FilterGenes(B_cells_MCL,x,unique = F))
Idents(B_cells_MCL) ="orig.ident"
cluster=F
ScatterPlot = F
#features_list <- unlist(features.list) %>% unique()
#max_UMI = B_cells_MCL[["SCT"]]@data %>% .[features_list,] %>% apply(1,max)

breaks = 0

for(s in 1:length(samples)){ #
        sample = samples[s]
        s_path <- paste0(save.path,sample,"/")
        if(!dir.exists(s_path)) dir.create(s_path, recursive = T)
        if(sample == "All_samples") {
                subset_object = B_cells_MCL
        } else subset_object = subset(B_cells_MCL, idents = sample)
        for(i in 1:length(features.list)){ #5:
                # FeaturePlot.2
                g <- FeaturePlot.2(object = subset_object, features = features.list[[i]],
                                   do.return = T,
                                   overlay = T,cols.use = c("#d8d8d8","#b88801","#2c568c", "#E31A1C"),
                                   pt.size = 3, alpha = 1, breaks = breaks)
                jpeg(paste0(s_path,sample,"_",paste(features.list[[i]],collapse = "_"),".jpeg"),
                     units="in", width=7, height=7,res=600)
                g = g+theme(plot.title = element_text(hjust = 0.5,size = 20),
                            legend.position="bottom",
                            legend.text=element_text(size=15))+
                        guides(colour = guide_legend(override.aes = list(size=5)),
                               shape = guide_legend(override.aes = list(size=5)))
                print(g)
                dev.off()
                # X4clusters
                Idents(subset_object) = "X4clusters"
                if(!cluster) {
                        features_var <- FetchData(subset_object,features.list[[i]])
                        df1 <- table(features_var[,1]==0,features_var[,2]==0) %>% #prop.table(margin = NULL) %>%
                                as.data.frame() %>% spread(Var2,Freq)
                        df = df1[,-1]
                        if(class(df) == "integer") next # if no expression
                        if(nrow(df) == 1) next # if only one expression
                        FISH <- fisher.test(df,workspace = 2000000)
                        df$p_value = FISH$p.value
                        df$p_val_adj = p.adjust(p = df$p_value, method = "bonferroni",
                                                n = nrow(df))
                        df[2,c("p_value","p_val_adj")] = ""
                        rownames(df)= paste(features.list[[i]][2],c("> 0","== 0"))
                        colnames(df)[1:2] = paste(features.list[[i]][2],c("> 0","== 0"))
                        write.csv(df,file = paste0(s_path,sample,"_",paste(features.list[[i]],
                                                                      collapse = "_"),".csv"))
                        if(ScatterPlot){
                                jpeg(paste0(s_path,sample,"_ScatterPlot_",paste(features.list[[i]],collapse = "_"),".jpeg"),
                                     units="in", width=7, height=7,res=600)
                                g <- FeatureScatter(subset_object, feature1 = features.list[[i]][1],
                                                    group.by = "orig.ident",
                                                    pt.size = 2,#cols = "black",
                                                    feature2 = features.list[[i]][2],slot = "data")+ NoLegend()
                                print(g)
                                dev.off()
                        }
                }
                if(cluster == T) for(k in unique(subset_object$X4clusters)){
                        subset_object_n <- subset(subset_object, idents = k)
                        features_var <- FetchData(subset_object_n,features.list[[i]])
                        df1 <- table(features_var[,1]>0,features_var[,2]>0) %>% #prop.table(margin = NULL) %>%
                                as.data.frame() %>% spread(Var2,Freq)
                        df = df1[,-1]
                        if(class(df) == "integer") next # if no expression
                        if(nrow(df) == 1) next # if only one expression
                        FISH <- fisher.test(df,workspace = 2000000)
                        df$p_value = FISH$p.value
                        df$p_val_adj = p.adjust(p = df$p_value, method = "bonferroni",
                                                n = nrow(df))
                        df[2,c("p_value","p_val_adj")] = ""
                        rownames(df)= plyr::mapvalues(df1$Var1, c(FALSE, TRUE),
                                                      paste(features.list[[i]][1],c("== 0","> 0")))
                        colnames(df)[1:2] = paste(features.list[[i]][2],c("== 0","> 0"))
                        write.csv(df,file = paste0(s_path,sample,"_",paste(features.list[[i]],
                                                                      collapse = "_"),"_cluster_",k,".csv"))

                        jpeg(paste0(s_path,sample,"_ScatterPlot_",paste(features.list[[i]],collapse = "_"),"_cluster_",k,".jpeg"),
                             units="in", width=7, height=7,res=600)
                        g <- FeatureScatter(subset_object_n, feature1 = features.list[[i]][1],
                                            pt.size = 3,
                                            feature2 = features.list[[i]][2],slot = "data")+ NoLegend()
                        print(g)
                        dev.off()
                }
        Progress(s,length(samples))}
}

Idents(B_cells_MCL) = "orig.ident"
for(N in c("B_cells", "MCL")){
        if(N == "B_cells") subset_object <- subset(B_cells_MCL,
                                                   idents = c("N02","N01","N03"))
        if(N == "MCL") subset_object <- subset(B_cells_MCL,
                                                   idents = c("PtU01","PtU02","PtU03","PtU04"))
                s_path <- paste0(path,N,"/")
        if(!dir.exists(s_path)) dir.create(s_path, recursive = T)
        for(i in 5:length(features.list)){
                # FeaturePlot.2
                g <- FeaturePlot.2(object = subset_object, features = features.list[[i]],do.return = T,
                                   overlay = T,cols = c("#d8d8d8","#b88801","#2c568c", "#E31A1C"),
                                   pt.size = 4, alpha = 0.75, breaks =8)
                jpeg(paste0(s_path,N,"_",paste(features.list[[i]],collapse = "_"),".jpeg"),
                     units="in", width=7, height=7,res=600)
                g = g+theme(plot.title = element_text(hjust = 0.5,size = 20),
                            legend.position="bottom",
                            legend.text=element_text(size=15))+
                        guides(colour = guide_legend(override.aes = list(size=5)),
                               shape = guide_legend(override.aes = list(size=5)))
                print(g)
                dev.off()
                # percentage
                Idents(subset_object) = "X4clusters"
                for(k in unique(subset_object$X4clusters)){
                        subset_object_n <- subset(subset_object, idents = k)
                        features_var <- FetchData(subset_object_n,features.list[[i]])
                        df1 <- table(features_var[,1]>0,features_var[,2]>0) %>% #prop.table(margin = NULL) %>%
                                as.data.frame() %>% spread(Var2,Freq)
                        df = df1[,-1]
                        if(class(df) == "integer") next
                        p_value = c()
                        ColSum <- colSums(df)
                        for(m in 1:nrow(df)){
                                (conting <- rbind(df[m,],ColSum-df[m,]))
                                FISH <- fisher.test(conting,workspace = 2000000)
                                (p_value[m] = FISH$p.value)
                                #CHI = chisq.test(conting, correct = T)
                                #chisq_p_value[i] = CHI$p.value
                        }
                        df$p_value = p_value
                        df$p_val_adj = p.adjust(p = df$p_value, method = "bonferroni",
                                                n = nrow(df))
                        rownames(df)= plyr::mapvalues(df1$Var1, c(FALSE, TRUE),
                                                      paste(features.list[[i]][1],c("== 0","> 0")))
                        colnames(df)[1:2] = paste(features.list[[i]][2],c("== 0","> 0"))
                        write.csv(df,file = paste0(s_path,N,"_",paste(features.list[[i]],
                                                                      collapse = "_"),"_cluster_",k,".csv"))

                        jpeg(paste0(s_path,N,"_ScatterPlot_",paste(features.list[[i]],collapse = "_"),"_cluster_",k,".jpeg"),
                             units="in", width=7, height=7,res=600)
                        g <- FeatureScatter(subset_object_n, feature1 = features.list[[i]][1],
                                            pt.size = 4,
                                            feature2 = features.list[[i]][2],slot = "data")
                        print(g)
                        dev.off()
                }
        }
}

Idents(subset_object) = "orig.ident"
RidgePlot(subset_object, features = c("EZH1", "EZH2"), ncol = 1)
FeatureScatter(B_cells_MCL, feature1 = "EZH2",feature2 = "CDK1",slot = "data")+NoLegend()
