library(ggplot2)
library(Seurat)
library(dplyr)  
library(magrittr)
library(Matrix)

#' prepare exp and tsne file
PrepareShiny <- function(object, samples, path, verbose = F){
        
        object <- SetAllIdent(object, id = "orig.ident")
        subset.object <- SubsetData(object, ident.use = samples)
        subset.object@meta.data$orig.ident %>% unique %>% sort %>% print
        
        exp <- lapply(samples, function(sample) {
                single_object <- SubsetData(subset.object, ident.use = sample)
                data <- single_object@data
                bad <- apply(data,1,max) == 0
                data = data[!bad,]
                range <- apply(data,1,max) - apply(data,1,min)
                scale_data = sweep(data, 1, range,"/")*2.5
                
                if(verbose) {
                    print(sample)
                    print(format(object.size(scale_data),units="MB"))
                }
                return(scale_data)
        })
        names(exp) = samples
        message("express data complete")
        #============== tsne csv===============
        tsne <- lapply(samples, function(sample) {
                single_object <- SubsetData(subset.object, ident.use = sample)
                single_object@dr$tsne@cell.embeddings
        })
        names(tsne) = samples
        shiny_path <- paste0(path, "data/")
        if(!dir.exists(shiny_path)) dir.create(shiny_path, recursive = T)
        save(exp,tsne, file = paste0(shiny_path,basename(path),".Rda"))
}


#' select 1/4 of cell from control
ScaleDown <- function(object, control=c("BH","DJ","MD","NZ")){
        
        normal_cells = lapply(control, function(x){
                rownames(object@meta.data)[(object@meta.data$orig.ident %in% x)]
        })
        set.seed(101)
        remove_normal_cells = lapply(normal_cells, function(x) sample(x, size = length(x)*3/4)) %>% unlist
        table(object@cell.names %in% remove_normal_cells)
        cell.use <- object@cell.names[!(object@cell.names %in% remove_normal_cells)]
        object <- SubsetData(object, cells.use = cell.use)
        object@meta.data$orig.ident = gsub(paste(control,collapse = "|"),"Normal",object@meta.data$orig.ident)
        return(object)
}

(load(file="data/MCL_Harmony_43_20190430.Rda"))
# =========== select 1/4 of cell from control ===============
object <- ScaleDown(object = object)
#============== expression Rda ===============
PrepareShiny(object, samples, path, verbose = T)
