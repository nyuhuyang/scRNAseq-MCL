#remove NA columns and NA rows, remove duplicate Gene_name
CleanUp <- function(df){
    
    rm_NA_col <- df[which(df[,1] == "Gene_name"),] %>% is.na %>% as.vector
    df = df[,!rm_NA_col]
    rm_NA_row <- apply(df,1, function(x) all(is.na(x)))
    df = df[!rm_NA_row,]
    colnames(df) = df[which(df[,1] == "Gene_name"),] %>% as.character
    df = df[-which(df[,1] == "Gene_name"),]
    
    rm_col <- colnames(df) %in% c("Gene_id","biotype")
    df = df[,!rm_col]
    df = RemoveDup(df)
    
    return(df)
}

#remove duplicate rownames with lower rowsumns
#' @param mat input as data.frame with gene name
#' @export mat matrix with gene as rownames, no duplicated genes
RemoveDup <- function(mat){
    gene_id <- as.matrix(mat[,1])
    mat <- mat[,-1]
    if(!is.matrix(mat)) mat <- sapply(mat,as.numeric)
    rownames(mat) <- 1:nrow(mat)
    mat[is.na(mat)] = 0
    mat <- cbind(mat, "rowSums" = rowSums(mat))
    mat <- mat[order(mat[,"rowSums"],decreasing = T),]
    gene_id <- gene_id[as.numeric(rownames(mat))]
    remove_index <- duplicated(gene_id)
    mat <- mat[!remove_index,]
    rownames(mat) <- gene_id[!remove_index]
    return(mat[,-ncol(mat)])
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
