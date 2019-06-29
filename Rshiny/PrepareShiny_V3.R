library(ggplot2)
library(Seurat)
library(dplyr)  
library(magrittr)
library(Matrix)
library(qlcMatrix)
#' prepare exp and tsne file
PrepareShiny <- function(object, samples, path, verbose = F,scale =NULL){
        
        Idents(object) <-  "orig.ident"
        object@meta.data$orig.ident %>% unique %>% sort %>% print
        
        exp <- lapply(samples, function(sample) {
                single_object <- subset(object, idents = sample)
                data <- GetAssayData(single_object)
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
                return(data)
        })
        names(exp) = samples
        message("express data complete. Access reductions data")
        #============== tsne csv===============
        tsne <- lapply(samples, function(sample) {
                single_object <- subset(object, idents = sample)
                single_object@reductions$tsne@cell.embeddings[colnames(single_object),]
        })
        names(tsne) = samples
        shiny_path <- paste0(path, "data/")
        if(!dir.exists(shiny_path)) dir.create(shiny_path, recursive = T)
        save(exp,tsne, file = paste0(shiny_path,basename(path),".Rda"))
}

(load(file = "data/B_cells_MCL_43_20190615.Rda"))
#============== expression Rda ===============
# remove cluster with less than 3 cells in Differential_analysis_B.R ======
Idents(B_cells_MCL) <- "X6_clusters"
B_cells_MCL_cluster4 <- subset(B_cells_MCL,idents = c(4,6))
PrepareShiny(B_cells_MCL_cluster4, samples, Rshiny_path, verbose = T,scale=2.5)
