library(ggplot2)
library(Seurat)
library(dplyr)  
library(magrittr)
library(Matrix)
library(qlcMatrix)
source("../R/Seurat3_functions.R")
#' prepare exp and tsne file
#' prepare exp and tsne file
PrepareShiny <- function(object, samples, Rshiny_path, split.by = "orig.ident",reduction = "tsne",
                         verbose = F,scale =NULL, assay = NULL){
        if(missing(object) | class(object) != "Seurat") stop("samples is not provided")
        if(missing(samples)) stop("samples is not provided")
        if(missing(Rshiny_path)) stop("Rshiny_path is not provided")
        assay = DefaultAssay(object) %||% assay
        Idents(object) <-  split.by
        avaible_samples <- samples %in% c("All_samples",object@meta.data[,split.by])
        if (!all(avaible_samples))
                stop(paste(paste(samples[!avaible_samples],collapse = " "),
                           "are not exist in the data."))
        max_exp <- list()
        if("All_samples" %in% samples) {
                single_object <- object
        } else single_object <- subset(object, idents = samples)
        max_exp = rowMax(single_object[[assay]]@data) %>% as.vector()
        max_exp = max_exp/log(2)
        names(max_exp) = rownames(single_object)
        
        exp <- list()
        tsne <- list()
        for (i in seq_along(samples)){
                sample <- samples[i]
                if(sample == "All_samples") {
                        single_object <- object
                } else single_object <- subset(object, idents = sample)
                #============== exp csv===============
                data <- GetAssayData(single_object)
                data <- as(data, "sparseMatrix")
                data = data/log(2)
                #bad <- rowMax(data) == 0
                #data = data[!bad,]
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
        save(exp,tsne,max_exp, file = paste0(shiny_data_path,basename(Rshiny_path),".Rda"))
}
(load(file = "data/MCL_41_harmony_20200225.Rda"))
#============== expression Rda ===============
object$orig.ident %<>% gsub("N01|N02|N03","Normal",.)
Idents(object) = "Doublets"
object %<>% subset(idents = "Singlet")
PrepareShiny(object, samples = samples, Rshiny_path = Rshiny_path, 
             verbose = T)
