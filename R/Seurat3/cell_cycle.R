########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(plyr)
library(tidyr)
library(kableExtra)
library(magrittr)
library(ggpubr)
source("../R/Seurat3_functions.R")
source("R/util.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#3.1  Compare DE across all major cell types==================
#We would need the data for all clusters, as well the subclusters.
#detect changes in gene expression between young and aged, 
#in the different cell types and subtypes. 
#It will also be interesting to check if there is some subtype enriched in young compared to aged or viceversa. 

# 3.1.1 load data
# Rename ident
(load(file = "data/B_cells_MCL_43_20190713.Rda"))
(load(file = "data/T_NK_cells_43_20190717.Rda"))
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
cc.genes <- readxl::read_excel("data/regev_lab_cell_cycle_genes.xls")
s.genes <- na.omit(cc.genes$s.genes)
g2m.genes <- cc.genes$g2m.genes

df_samples <- readxl::read_excel("doc/190626_scRNAseq_info.xlsx") %>% as.data.frame()
colnames(df_samples) =  tolower(colnames(df_samples))
(keep = df_samples$tests %in% c("control",paste0("test",2:12)))
#===============================
# cell cycle B cells
#===============================
B_cells_MCL %<>% CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
B_meta.data = B_cells_MCL@meta.data[c("orig.ident","groups","Phase","S.Score","G2M.Score")]

jpeg(paste0(path,"cell_cycle_NK.jpeg"), units="in", width=3, height=3,res=600)
ggscatter(B_meta.data, x = "S.Score", y = "G2M.Score",
          palette = c("#E7B800","#D6604D","#4393C3"),
          color = "Phase")
dev.off()

df_sample_B <- df_samples[keep,c("sample","group")]
(df_sample_B$samples = c(df_sample_B$group[1:8],df_sample_B$sample[9:nrow(df_sample_B)]))

for(i in 1:nrow(df_sample_B)){
        Phase =B_meta.data[B_meta.data$orig.ident %in% df_sample_B$sample[i],"Phase"]  
        Phase_score = table(Phase) %>% as.data.frame() %>% .$Freq
        df_sample_B[i,c("G1", "G2M", "S")] = Phase_score/sum(Phase_score)*100
}
df_sample_B = df_sample_B[complete.cases(df_sample_B),]
(Groups = as.character(unique(df_sample_B$group)))
df_sample_B$group = factor(df_sample_B$group, levels = groups)
#===============================
# cell cycle CD4, CD8, NK
#===============================
Idents(T_NK_cells) = "manual"
(T.NK.cells = as.character(unique(Idents(T_NK_cells))))
T_NK_list <- lapply(T.NK.cells, function(x) subset(T_NK_cells, idents = x))
T_NK_list %<>% lapply(function(x) CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE))
T_NK_meta.data_list = lapply(T_NK_list, function(x) x@meta.data[c("orig.ident","groups","Phase","S.Score","G2M.Score")])

df_sample <- df_samples[keep,c("sample","group")]
(df_sample$samples = c(df_sample$group[1:8],df_sample$sample[9:nrow(df_sample)]))
Phase_score_list <- list()

for(i in 1:length(T.NK.cells)){
        
        Phase_score_list[[i]] = data.table::copy(df_sample)
        for(k in 1:nrow(Phase_score_list[[i]])){
                Phase =T_NK_meta.data_list[[i]][T_NK_meta.data_list[[i]]$orig.ident %in% 
                                                     Phase_score_list[[i]]$sample[k],"Phase"]  
                Phase_score = table(Phase) %>% as.data.frame() %>% .$Freq
                Phase_score_list[[i]][k,c("G1", "G2M", "S")] = Phase_score/sum(Phase_score)*100

        }
        Phase_score_list[[i]] = Phase_score_list[[i]][complete.cases(Phase_score_list[[i]]),]
        Phase_score_list[[i]]$group = factor(Phase_score_list[[i]]$group, levels = groups)
}

(Groups = Groups[c(2:7,9:11)])

#' produce likert style stack barplot
#' http://rnotr.com/likert/ggplot/barometer/likert-plots/
ggbarplot.1 <- function(data,title){
        data1 <- reshape::melt(data, id="samples")
        data_mean <- aggregate(. ~ samples + variable, data = data1, FUN=mean)
        data_se <- aggregate(. ~ samples + variable, data = data1, FUN=sd)
        data2 <- cbind.data.frame(data_mean, data_se[,3])
        colnames(data2)[4] = "mean_se"
        data2$mean_se[is.na(data2$mean_se)]=0
        
        pal = c("#E7B800","#D6604D","#4393C3")
        mylevels = c("G1","S", "G2M")
        data2$col <- rep(pal,each=nrow(data2)/length(pal))
        data2$samples <- stringr::str_wrap(data2$samples, width = 40)
        data2$samples <- factor(data2$samples, levels = unique(data$samples))
        
        data2 = data2[order(data2$samples),]
        
        S_G2M <- na.omit(data2[data2$variable %in% c("S","G2M"),])
        S_G2M <- within(S_G2M,cumsum <- ave(value,samples,FUN=cumsum))
        
        G1 <- na.omit(data2[data2$variable %in% "G1",])
        G1 <- G1[rev(rownames(G1)),]

        ggplot() + geom_bar(data=S_G2M, aes(x = samples, y=value, fill=col), colour="black", 
                            position="stack", stat="identity") +
                geom_bar(data=G1, aes(x = samples, y=-value, fill=col), colour="black", 
                         position="stack", stat="identity") +
                geom_errorbar(data=S_G2M, aes(x = samples, ymin=cumsum-mean_se, 
                                              ymax=cumsum+mean_se), width=.2)+
                geom_errorbar(data=G1, aes(x = samples, ymin= -value-mean_se, 
                                              ymax=-value+mean_se), width=.2)+                
                geom_hline(yintercept = 0, color =c("white")) +
                scale_fill_identity("Cell cycle phase", labels = mylevels, breaks=pal, guide="legend") + 
                ggthemes::theme_fivethirtyeight() + 
                labs(title=title, y="",x="") +
                theme(plot.title = element_text(size=14, hjust=0.5)) +
                theme(axis.text.y = element_text(hjust=0)) +
                theme(legend.position = "bottom")
        }

#' produce stack ggplot barchart
ggbarplot.2 <- function(data,title){
        data1 <- reshape::melt(data, id="samples")
        data_mean <- aggregate(. ~ samples + variable, data = data1, FUN=mean)
        data_se <- aggregate(. ~ samples + variable, data = data1, FUN=sd)
        data2 <- cbind.data.frame(data_mean, data_se[,3])
        colnames(data2)[4] = "mean_se"
        data2$mean_se[is.na(data2$mean_se)]=0
        
        pal = c("#E7B800","#D6604D","#4393C3")
        mylevels = c("G1","S", "G2M")
        data2$col <- rep(pal,each=nrow(data2)/length(pal))
        data2$samples <- stringr::str_wrap(data2$samples, width = 40)
        data2$samples <- factor(data2$samples, levels = unique(data$samples))
        
        data2 = data2[order(data2$samples),]
        data3 <- within(data2,cumsum <- ave(value,samples,FUN=cumsum))
        
        ggplot(data=data3, aes(x = samples, y=value, fill=col)) + 
                geom_bar(colour="black", position="stack", stat="identity") +
                geom_errorbar(aes(ymin=cumsum-mean_se, ymax=cumsum+mean_se), width=.2) +
                geom_hline(yintercept = 0, color =c("white")) +
                scale_fill_identity("Cell cycle phase", labels = mylevels, breaks=pal, guide="legend") + 
                ggthemes::theme_fivethirtyeight() + 
                labs(title=title, y="",x="") +
                theme(plot.title = element_text(size=14, hjust=0.5)) +
                theme(axis.text.y = element_text(hjust=0)) +
                theme(legend.position = "bottom")
}

# https://felixfan.github.io/stacking-plots-same-x/
# B cells
for(i in 1:length(Groups)){
        index = df_sample_B$group %in% c("Normal", "Untreated", Groups[i])
        data = df_sample_B[index,c(ifelse(i==1,"sample","samples"),"G1","S","G2M")]
        colnames(data)[1] = "samples"
        
        title= ifelse(i==1,
                      "Percentage of G2M, S, and G1 phase cells in \nNormal and Untreated patients' B and MCL cells",
                      paste0("Percentage of G2M, S, and G1 cells in \nNormal, Untreated and patient ",Groups[i],"'s B and MCL cells"))
        
        figure2 <- ggbarplot.2(data, title)
        B_path <- paste0(path,"B_Barchart/")
        if(!dir.exists(B_path)) dir.create(B_path, recursive = T)
        jpeg(paste0(B_path,"B_cell_cycle_",Groups[i],".jpeg"), units="in", width=10, height=7,res=600)
        print(figure2)
        dev.off()
}
for(i in 1:length(phases)){
        res<- t.test(df_sample_B[df_sample_B$samples %in% "Normal",phases[i]],
               df_sample_B[df_sample_B$samples %in% "Untreated",phases[i]],)
        print(res)
}
phases = c("G1","G2M","S")
for(i in 1:length(phases)){
        res<- t.test(df_sample_B[df_sample_B$samples %in% 
                                         c("Pt-17-LN-C1",
                                           "Pt-17-C2","Pt-17-C7"),phases[i]])
        print(res)
}

# T cells =============================================
for(i in 1:length(T.NK.cells)){
        for(k in 1:length(Groups)){
                index = Phase_score_list[[i]]$group %in% c("Normal", "Untreated", Groups[k])
                data = Phase_score_list[[i]][index,c(ifelse(k==1,"sample","samples"),"G1","S","G2M")]
                colnames(data)[1] = "samples"
                
                title= ifelse(k==1,
                              paste0("Percentage of G2M, S, and G1 phase cells in \nNormal and Untreated patients' ",T.NK.cells[i]),
                              paste0("Percentage of G2M, S, and G1 cells in \nNormal, Untreated and patient ",Groups[k],"'s ",T.NK.cells[i]))
                
                figure2 <- ggbarplot.2(data, title)
                T_path <- paste0(path,"T_Barchart/",T.NK.cells[i],"/")
                if(!dir.exists(T_path)) dir.create(T_path, recursive = T)
                jpeg(paste0(T_path,"Cell_cycle_",Groups[k],"_",T.NK.cells[i],".jpeg"),
                     units="in", width=10, height=7,res=600)
                print(figure2)
                dev.off()
        }
}

df_sample_T = Phase_score_list[[2]]
phases = c("G1","G2M","S")
for(i in 1:length(phases)){
        res<- t.test(df_sample_T[df_sample_T$samples %in% "Normal",phases[i]],
                     df_sample_T[df_sample_T$samples %in% "Untreated",phases[i]],)
        print(res)
}
for(i in 1:length(phases)){
        res<- t.test(df_sample_T[df_sample_T$samples %in% c("Pt-17-C2",
                                                          "Pt-17-C7",
                                                          "Pt-17-C31"),phases[i]])
        print(res)
}

}