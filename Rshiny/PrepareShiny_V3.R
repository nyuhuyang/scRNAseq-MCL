library(ggplot2)
library(Seurat)
library(dplyr)
library(magrittr)
library(Matrix)
library(qlcMatrix)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
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

#Prepare Rshiny from multiple Seurat objects
#' @param ... path of multiple Seurat objects
#' @orig.idents original identities, sample names
#' @Rshiny_path Rshiny folder path
#' @param split.by split objecty by. Colname of meta.data. Set to FALSE if want to use all
#' @example
#' (data_path = list.files("output/20200612/B-Correlation-pvalues",pattern = ".rds"))
# PrepareShiny.1(paste0("output/20200612/B-Correlation-pvalues/",data_path),
#               samples = samples, Rshiny_path = Rshiny_path,
#               split.by = "conditions", args = c("harmony","pca"),
#               reduction = "umap",verbose = T)
PrepareCorShiny <- function(..., orig.idents = NULL, Rshiny_path = NULL, verbose = F){
        if(is.null(orig.idents)) stop("samples is not provided")
        if(is.null(Rshiny_path)) stop("Rshiny_path is not provided")
        data_path <- c(...)
        if(!all(sapply(data_path, file.exists))) stop("Not all Rshiny_path has file exist")
        cor_list <- list()
        pvalue_list <- list()
        s = length(data_path)
        for(m in seq_along(data_path)){
                if(tools::file_ext(data_path[[m]]) == "Rda") load(file = data_path[[m]])
                if(tools::file_ext(data_path[m]) == "rds") cor_res = readRDS(file = data_path[m])
                cor = cor_res$r
                pvalue = cor_res$P
                M = lower.tri(cor_res$r, diag = T)
                cor[M] = 0
                pvalue[M] = 0
                cor_list[[m]] = as(cor, "sparseMatrix")
                pvalue_list[[m]] = as(pvalue, "sparseMatrix")
                if(verbose) Progress(m, s)
        }

        names(cor_list) = gsub(".*-","",data_path) %>% gsub("\\.rds","",.)
        names(pvalue_list) = gsub(".*-","",data_path) %>% gsub("\\.rds","",.)
        cor_list = cor_list[orig.idents]
        pvalue_list = pvalue_list[orig.idents]

        message("writing files")
        shiny_data_path <- paste0(Rshiny_path, "data/")
        if(!dir.exists(shiny_data_path)) dir.create(shiny_data_path, recursive = T)
        save(cor_list,pvalue_list, file = paste0(shiny_data_path,basename(Rshiny_path),".Rda"))
}

data_path = list.files("output/20200612/B-Correlation-pvalues",pattern = ".rds")
data_path = grep(paste(orig.idents,collapse = "|"), data_path, value = T)
(data_path = paste0("output/20200612/B-Correlation-pvalues/",data_path))
Rshiny_path <- "Rshiny/B_MCL_41_correlation1/"
PrepareCorShiny(data_path, orig.idents = orig.idents, Rshiny_path = Rshiny_path, verbose = T)


(load(file = "data/MCL_41_harmony_20200225.Rda"))
object = readRDS("data/MCL_41_B_20200225.rds")
#============== expression Rda ===============
object$orig.ident %<>% gsub("N01|N02|N03","Normal",.)
Idents(object) = "Doublets"
object %<>% subset(idents = "Singlet")

PrepareShiny(object, samples = samples, Rshiny_path = Rshiny_path,
             verbose = T)
