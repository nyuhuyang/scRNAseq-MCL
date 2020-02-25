library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(magrittr)
source("../R/Seurat3_functions.R")
path <- paste0("Yang/Figure 4 Xiangao/")
if(!dir.exists(path))dir.create(path, recursive = T)

#=============== multiple color in the single Featureplot===================================
B_cells_MCL = readRDS(file = "data/MCL_41_B_20200207.rds")

Idents(B_cells_MCL) = "orig.ident"
samples = c("All_samples","Pt10_LN2Pd","Pt11_LN1", "Pt17_LN1","PtU01","PtU02","PtU03","PtU04")
features.list = lapply(list(c("EZH2","E2F1"),
                            c("EZH2","PCNA"),
                            c("EZH2","CDK1"),
                            c("EZH1","EZH2"),
                            c("POLR2M", "CRBN"),
                            c("POLR2M", "IKZF1"),
                            c("POLR2M", "IKZF3"),
                            c("IKZF1", "CRBN"),
                            c("IKZF3", "CRBN"),
                            c("MBOAT7", "CRBN"),
                            c("MBOAT7", "IKZF1"),
                            c("MBOAT7", "IKZF3")),
                       function(x) FilterGenes(B_cells_MCL,x,unique = F))
(cols.use.list = rep(list(c("#b88801","#2c568c", "#E31A1C")), length(features.list)))
Idents(B_cells_MCL) ="orig.ident"
cluster=F
for(s in samples[1:length(samples)]){
        s_path <- paste0(path,s,"/")
        if(!dir.exists(s_path)) dir.create(s_path, recursive = T)
        if(s == "All_samples") {
                subset_object = B_cells_MCL
        } else subset_object = subset(B_cells_MCL, idents = s)
        for(i in 1:length(features.list)){ #5:
                # FeaturePlot.2
                g <- FeaturePlot.2(object = subset_object, features = features.list[[i]],do.return = T,
                                   overlay = T,cols = c("#d8d8d8",cols.use.list[[i]]),
                                   pt.size = 2, alpha = 0.75, breaks =8)
                jpeg(paste0(s_path,s,"_",paste(features.list[[i]],collapse = "_"),".jpeg"), 
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
                        write.csv(df,file = paste0(s_path,s,"_",paste(features.list[[i]],
                                                                      collapse = "_"),".csv"))
                        
                        jpeg(paste0(s_path,s,"_ScatterPlot_",paste(features.list[[i]],collapse = "_"),".jpeg"), 
                             units="in", width=7, height=7,res=600)
                        g <- FeatureScatter(subset_object, feature1 = features.list[[i]][1],
                                            pt.size = 2,
                                            feature2 = features.list[[i]][2],slot = "data")
                        print(g)
                        dev.off()
                }
                if(cluster == T) for(k in unique(subset_object$X4clusters)){
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
                        write.csv(df,file = paste0(s_path,s,"_",paste(features.list[[i]],
                                                                      collapse = "_"),"_cluster_",k,".csv"))
                        
                        jpeg(paste0(s_path,s,"_ScatterPlot_",paste(features.list[[i]],collapse = "_"),"_cluster_",k,".jpeg"), 
                             units="in", width=7, height=7,res=600)
                        g <- FeatureScatter(subset_object_n, feature1 = features.list[[i]][1],
                                            pt.size = 4,
                                            feature2 = features.list[[i]][2],,slot = "data")
                        print(g)
                        dev.off()
                }
                
        
                Progress(i,length(features.list))}
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
                                   overlay = T,cols = c("#d8d8d8",cols.use.list[[i]]),
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