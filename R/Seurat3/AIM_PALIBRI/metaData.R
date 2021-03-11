########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
#
# ######################################################################
invisible(lapply(c("Seurat","dplyr","magrittr"), function(x) {
    suppressPackageStartupMessages(library(x,character.only = T))
}))
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
#======1.1 Setup the Seurat objects =========================
(load(file = "data/MCL_AIM_58_20201009.Rda"))
Idents(object) = "Doublets"
object %<>% subset(idents = "Singlet")

Idents(object) = "cell.types"
object %<>% subset(idents = c("HSC/progenitors","Nonhematopoietic_cells","Plasma_cells"), invert = T)

df <- list()
df[["cell.types"]] <- table(object$orig.ident, object$cell.types) %>% as.data.frame.matrix
df[["percentage.cell.types"]] <- table(object$orig.ident, object$cell.types) %>% prop.table(margin = 1)
df[["percentage.cell.types"]] = df[["percentage.cell.types"]]*100
openxlsx::write.xlsx(df, file =  paste0(path,"20210119_cell.types_orig.ident.xlsx"),
                    colNames = TRUE,row.names = TRUE,borders = "surrounding",colWidths = c(NA, "auto", "auto"))

#====== compare with phase I data ===============
Idents(object) = "phase"
MCL_41  <- subset(object, idents = "I")
meta.data = MCL_41@meta.data
rownames(meta.data) %<>% gsub("-1$","",.)

rm(object)
(load(file = "data/MCL_41_harmony_20191231.Rda"))
meta.data1 = object@meta.data
rownames(meta.data1) %<>% gsub("Pt-28-PB-C25D1","Pt-28-PB-C28D1",.)
rownames(meta.data1) %<>% gsub("Pt-13-BMA","Pt-13-BM",.)

cells = intersect(rownames(meta.data1), rownames(meta.data))
meta.data1 = meta.data1[cells,]
meta.data = meta.data[cells,]

table(rownames(meta.data1) == rownames(meta.data))
table(meta.data$cell.types, meta.data1$cell.types) %>% as.data.frame.matrix() %>%
    write.csv(file =  paste0(path,"20210118_cell.types_comparision.csv"))
