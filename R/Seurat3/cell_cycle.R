########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(SingleR)
library(dplyr)
library(plyr)
library(tidyr)
library(kableExtra)
library(magrittr)
library(harmony)
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
(load(file = "data/MCL_V3_Harmony_43_20190610.Rda"))
(load(file = "data/T_NK_cells_43_20190611.Rda"))
(load(file = "data/B_cells_MCL_43_20190619.Rda"))
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
object %<>% CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
T_NK_cells %<>% CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
B_cells_MCL %<>% CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

meta.data = object@meta.data[c("orig.ident","groups","S.Score","G2M.Score","Phase", "singler1sub")]
meta.data = cbind.data.frame(meta.data, object@reductions$tsne@cell.embeddings)

T_NK_meta.data = T_NK_cells@meta.data[c("orig.ident","groups","Phase")]
B_meta.data = B_cells_MCL@meta.data[c("orig.ident","groups","Phase")]

jpeg(paste0(path,"cell_cycle.jpeg"), units="in", width=3, height=3,res=600)
ggscatter(test, x = "S.Score", y = "G2M.Score",
          palette = c("#808080","#E7B800", "#FC4E07"),
          color = "Phase")
dev.off()

df_samples <- readxl::read_excel("doc/190429_scRNAseq_info.xlsx")
colnames(df_samples) =  tolower(colnames(df_samples))

(keep = df_samples$tests %in% c("control",paste0("test",2:12)))
df_sample_cc <- df_samples[keep,c("sample","group")]
df_sample_cc$samples = c(df_sample_cc$group[1:8],df_sample_cc$sample[9:nrow(df_sample_cc)])

for(i in 1:nrow(df_sample_cc)){
        Phase =meta.data[meta.data$orig.ident %in% df_sample_cc$sample[i],"Phase"]
        Phase_score = table(Phase) %>% as.data.frame() %>% .$Freq
        df_sample_cc[i,"G2_M_S_all"]=sum(Phase_score[2:3])/sum(Phase_score)*100

        Phase =B_meta.data[B_meta.data$orig.ident %in% df_sample_cc$sample[i],"Phase"]
        Phase_score = table(Phase) %>% as.data.frame() %>% .$Freq
        df_sample_cc[i,"G2_M_S_B"]=sum(Phase_score[2:3])/sum(Phase_score)*100
        
        Phase =T_NK_meta.data[T_NK_meta.data$orig.ident %in% df_sample_cc$sample[i],"Phase"]
        Phase_score = table(Phase) %>% as.data.frame() %>% .$Freq
        df_sample_cc[i,"G2_M_S_T"]=sum(Phase_score[2:3])/sum(Phase_score)*100
        }
df_sample_cc
(groups = unique(df_sample_cc$group))
df_sample_cc$group = factor(df_sample_cc$group, levels = groups)
Groups = c("Normal","Untreated","Pt-17","Pt-25")

# https://felixfan.github.io/stacking-plots-same-x/
for(i in 2:length(Groups)){
        index = df_sample_cc$group %in% c("Normal", "Untreated", Groups[i])
        suppressWarnings(g1 <- ggbarplot(df_sample_cc[index,], x = "samples", 
                        y = "G2_M_S_B", xlab = "",ylab= "",title="B and MCL cells",
                        fill = "group",xtickslab.rt= 45,
                        add = c("mean_se"), palette = c("#00AFBB", "#E7B800","#FC4E07"),
                        position = position_dodge())+
                                theme(axis.title.x = element_blank(), 
                                      axis.text.x = element_blank(),
                                      plot.title = element_text(hjust = 0.5)))
        suppressWarnings(g2 <- ggbarplot(df_sample_cc[index,], x = "samples", 
                        y = "G2_M_S_T", xlab = "",ylab= "",title="T and NK cells",
                        fill = "group",xtickslab.rt= 45,
                        add = c("mean_se"), palette = c("#00AFBB", "#E7B800","#FC4E07"),
                        position = position_dodge())+
                                theme(plot.title = element_text(hjust = 0.5)))
        suppressWarnings(figure <- ggarrange(g1,g2,ncol = 1, nrow = 2,
                                             common.legend = TRUE, legend="right"))
        #g <- grid.draw(rbind(ggplotGrob(g1), ggplotGrob(g2), size = "last"))
        figure1 <- annotate_figure(figure,
                        top = text_grob("Proliferating cell percentage in B & MCL cells and T & NK cells",
                                        color = "black", face = "plain", size = 15),
                        left = text_grob("Proliferating cell percentage (G2+M+S %)", 
                                         color = "black", rot = 90),
                        #fig.lab = "Figure 1", fig.lab.face = "bold"
        )
        jpeg(paste0(path,"ccbar_",Groups[i],".jpeg"), units="in", width=10, height=7,res=600)
        print(figure1)
        dev.off()
}


