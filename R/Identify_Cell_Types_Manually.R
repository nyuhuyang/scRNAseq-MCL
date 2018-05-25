library(Seurat)
library(dplyr)
source("./R/Seurat_functions.R")

#====== 2.1 identify phenotype for each cluster  ==========================================
lnames = load(file = "./data/MCL_alignment.Rda")
lnames

Featureplot <- function(x){
    p <- FeaturePlot(object = MCL, 
                     reduction.use = "tsne",
                     features.plot = x, min.cutoff = NA, 
                     cols.use = c("lightgrey","blue"), pt.size = 0.5)
    return(p)
}
Adipocytes <- HumanGenes(MCL,c("SLC36A2","P2RX5","MYF5","UCP1","TRIP4","ASCC1"))
Endothelium <- HumanGenes(MCL,c("Cdh5","Pecam1","Flt1","Plvap","Kdr","ptprb",
                                        "Vwf","EMCN","Car4"))
Epithelium <- HumanGenes(MCL,c("KRT19","Epcam","KRT5",
                                       "MUC1","SCGB3A2","SCGB1A1","SCGB3A1","SFTPB","FOXJ1","Rpe65",
                                       "Rlbp1","Msln","Upk3b","Lrrn4"))
RPE <- HumanGenes(MCL,c("Rpe65","Rlbp1"))
Fibroblast <- HumanGenes(MCL,c("FGF1","FGF9","SFRP1"))
Hematopoietic <- HumanGenes(MCL,c("PTPRC","LAPTM5","SRGN"))
Myeloid <-  HumanGenes(MCL,c("PPBP","GNG11","HBA2","HBB","Cma1","Mcpt4","Tpsb2",
                                     "Cpa3","LYZ","S100A9","CD14","CCL2","FCGR3A","MS4A7","VMO1"))
Lymphoid <- HumanGenes(MCL,c("CD3G","CD3D","CD2","Cd19","CD79A","MS4A1",
                                     "GNLY","Ncr1","CCL5","KLRD1","NKG7"))
Melanocytes <- HumanGenes(MCL,c("Pmel","Mlana"))
Mesenchymal <- HumanGenes(MCL,c("Pdgfrb","Has2","Dcn"))
Myelinating_Schwann_cells <- HumanGenes(MCL,c("MBP","MPZ"))
Pericytes <- HumanGenes(MCL,c("Pdgfrb","Cspg4","Anpep","Rgs5",
                                      "Myh11","Mylk","Des","Vtn","Ifitm1"))
Smooth_muscle_cells <- HumanGenes(MCL,c("Acta2","Myh11"))
Stem_cell <- HumanGenes(MCL,c("POU5F1","FUT4","CD34","PROM1","ABCG2","Runx1","ATXN1",
                                      "Nes","NCAM","NGFR"))
Stromal_fibroblasts <- HumanGenes(MCL,c("DCN","COL6A1","TIMP3","PDGFRA"))
Neurons <- HumanGenes(MCL,c("Ihh","Gli1", "Ptch1", "Hhip"))
# Featureplot
Featureplot(Adipocytes) # Adipocytes
Featureplot(Endothelium) # Endothelial Cells
Featureplot(Epithelium) # Epithelium
Featureplot(c(RPE,Melanocytes,Myelinating_Schwann_cells)) # RPE, Melanocytes, Myelinating Schwann cells
Featureplot(Fibroblast) # Fibroblasts
Featureplot(c(Hematopoietic,Myeloid[7:9],Lymphoid[1:3])) # Hematopoietic cells
Featureplot(Myeloid) # Myeloid cells
Featureplot(Lymphoid) # Lymphoid cells
Featureplot(Mesenchymal) # Mesenchymal cells
Featureplot(Pericytes) # Pericytes
Featureplot(Smooth_muscle_cells)
Featureplot(Stem_cell)
Featureplot(Stromal_fibroblasts)
Featureplot(Neurons)

markers.to.plot <- c(Melanocytes,Myelinating_Schwann_cells,Endothelium[c(1:3,5,7)],
                     Hematopoietic[1:2],Pericytes[c(4,6:7)],Mesenchymal[c(1,4)],RPE,
                     Smooth_muscle_cells)
markers.to.plot <- unique(markers.to.plot)
DotPlot(MCL, genes.plot = rev(markers.to.plot),
        cols.use = c("blue","red"), x.lab.rot = T, plot.legend = F,
        dot.scale = 8, do.return = T)

# Rename ident
table(MCL@ident)
idents <- as.data.frame(table(MCL@ident))
old.ident.ids <- idents$Var1
new.cluster.ids <- c("Pericytes 0",
                     "RPE 1",
                     "Endothelial cells 2",
                     "Pericytes 3",
                     "Endothelial cells 4",
                     "RPE 5",
                     "Smooth muscle cells 6",
                     "Pericytes 7",
                     "Monocytes 8",
                     "schwann cells 9",
                     "Endothelial cells 10",
                     "T cells 11",
                     "Monocytes 12")

MCL@ident <- plyr::mapvalues(x = MCL@ident,
                                    from = old.ident.ids,
                                    to = new.cluster.ids)
markers.to.plot <- c(Hematopoietic[1:2], Lymphoid[1:2],Myelinating_Schwann_cells,
                     Myeloid[c(7,9)],Pericytes[c(4,6:7)],Mesenchymal[c(1,3)],
                     Smooth_muscle_cells, Endothelium[c(1:3,5,7)],RPE,Melanocytes)
#markers.to.plot <- unique(markers.to.plot)
DotPlot(MCL, genes.plot = rev(markers.to.plot),
        cols.use = c("blue","red"), x.lab.rot = T, plot.legend = F,
        dot.scale = 8, do.return = T)
# MCL <- RenameIdentBack(MCL)

#====== 2.2 dot Plots ==========================================
lnames = load(file = "./data/MCL_alignment.Rda")
lnames
table(MCL@ident)
idents <- as.data.frame(table(MCL@ident))
old.ident.ids <- idents$Var1
new.cluster.ids <- c("Pericytes",
                     "Retinal pigmented epithelium",
                     "Endothelial cells",
                     "Pericytes",
                     "Endothelial cells",
                     "Retinal pigmented epithelium",
                     "Smooth muscle cells",
                     "Pericytes",
                     "Monocytes",
                     "Myelinating schwann cells",
                     "Endothelial cells",
                     "T cells",
                     "Monocytes")

MCL@ident <- plyr::mapvalues(x = MCL@ident,
                                    from = old.ident.ids,
                                    to = new.cluster.ids)
markers.to.plot <- c(Melanocytes,Myelinating_Schwann_cells,Hematopoietic[1:2],
                     RPE,Smooth_muscle_cells,Endothelium[c(1:3,5,7)],
                     Pericytes[c(4,6:7)],Mesenchymal[c(1,3)])
markers.to.plot <- unique(markers.to.plot)
DotPlot(MCL, genes.plot = rev(markers.to.plot),
        cols.use = c("blue","red"), x.lab.rot = T, plot.legend = F,
        dot.scale = 8, do.return = T)

freq_table <- prop.table(x = table(MCL@ident, MCL@meta.data[, "conditions"]),
margin = 2)
barplot(height = freq_table)
freq_table
table(MCL@meta.data[, "conditions"])

#=====2.3 tsne plot=============================
lnames = load(file = "./data/MCL_alignment.Rda")
lnames
table(MCL@ident)
idents <- as.data.frame(table(MCL@ident))
old.ident.ids <- idents$Var1
new.cluster.ids <- c("0) Pericytes",
                     "1) Retinal pigmented epithelium",
                     "2) Endothelial cells",
                     "3) Pericytes",
                     "4) Endothelial cells",
                     "5) Retinal pigmented epithelium",
                     "6) Smooth muscle cells",
                     "7) Pericytes",
                     "8) Monocytes",
                     "9) Myelinating schwann cells",
                     "10) Endothelial cells",
                     "11) T cells",
                     "12) Monocytes")

MCL@ident <- plyr::mapvalues(x = MCL@ident,
                                    from = old.ident.ids,
                                    to = new.cluster.ids)
TSNEPlot(object = MCL, no.legend = TRUE, do.label = TRUE,
         do.return = TRUE, label.size = 5)+
  ggtitle("TSNE plot of major cell types")+
  theme(text = element_text(size=20),     #larger text including legend title							
        plot.title = element_text(hjust = 0.5)) #title in middle


#====== 2.4 Compare cell type changes across conditions  ==========================================
# the two patients profiled have very different composition
# Compare clusters for each dataset
SplitTSNEPlot(MCL, "conditions",do.label = F, 
              do.return = TRUE, no.legend = T )
