library(ggplot2)
library(Seurat)
library(dplyr)  
library(magrittr)
library(Matrix)
library(qlcMatrix)
#' prepare exp and tsne file
#' prepare exp and tsne file
PrepareShiny <- function(object, samples, Rshiny_path, split.by = "orig.ident",reduction = "tsne",
                         verbose = F,scale =NULL){
        if(is.null(samples)) stop("samples is not provided")
        if(is.null(Rshiny_path)) stop("Rshiny_path is not provided")
        Idents(object) <-  split.by
        avaible_samples <- samples %in% c("All_samples",object@meta.data[,split.by])
        if (!all(avaible_samples))
                stop(paste(paste(samples[!avaible_samples],collapse = " "),
                           "are not exist in the data."))
        exp <- list()
        tsne <- list()
        max_exp <- list()
        for (i in seq_along(samples)){
                sample <- samples[i]
                if(sample == "All_samples") {
                        single_object <- object
                } else single_object <- subset(object, idents = sample)
                #============== exp csv===============
                data <- GetAssayData(single_object)
                data <- as(data, "sparseMatrix")
                bad <- rowMax(data) == 0
                data = data[!bad,]
                if(!is.null(scale)){
                        range <- rowMax(data) - rowMin(data) # range <- apply(data,1,max) - apply(data,1,min)
                        data = sweep(data, 1, range,"/")*scale
                }
                
                if(verbose) {
                        print(sample)
                        print(format(object.size(data),units="MB"))
                }
                exp[[i]] = data
                #============== tsne csv===============
                tsne[[i]] = Embeddings(single_object, reduction = reduction)
                
                svMisc::progress(i/length(samples)*100)
        }
        names(exp) = samples
        names(tsne) = samples
        shiny_data_path <- paste0(Rshiny_path, "data/")
        if(!dir.exists(shiny_data_path)) dir.create(shiny_data_path, recursive = T)
        save(exp,tsne, file = paste0(shiny_data_path,basename(Rshiny_path),".Rda"))
}
(load(file = "data/MCL_V3_Harmony_43_20190627.Rda"))
(load(file = "data/MCL_41_harmony_20191231.Rda"))
#============== expression Rda ===============
df_samples <- readxl::read_excel("doc/191030_scRNAseq_info.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
object$sample = object$orig.ident
object$orig.ident %<>% plyr::mapvalues(from = unique(df_samples$sample),
                                       to = unique(df_samples$`sample name`))

object$orig.ident %<>% gsub("N01|N02|N03","Normal",.)
Idents(object) = "Doublets"
object %<>% subset(idents = "Singlet")
PrepareShiny(object, samples = samples, Rshiny_path = Rshiny_path, 
             verbose = T)
