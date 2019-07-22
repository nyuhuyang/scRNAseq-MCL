
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


#' modify FeaturePlot to return TSNEplot like ggplot results
FeaturePlot.1 <- function (object, features, dims = c(1, 2), cells = NULL, 
                           cols = ifelse(test = c(blend,blend), yes = c("#ff0000", "#00ff00"), no = c("lightgrey", "blue")),
                           pt.size = NULL, order = FALSE, min.cutoff = NA, 
                           max.cutoff = NA, reduction = NULL, split.by = "orig.ident", shape.by = NULL, 
                           slot = "data", blend = FALSE, blend.threshold = 0.5, label = FALSE, 
                           label.size = 4, repel = FALSE, ncol = NULL, combine = TRUE, 
                           coord.fixed = FALSE, by.col = TRUE) {
    no.right <- theme(axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(), 
                      axis.text.y.right = element_blank(), axis.title.y.right = element_text(face = "bold",size = 14, margin = margin(r = 7)))
    if (is.null(reduction)) {
        default_order <- c("umap", "tsne", "pca")
        reducs <- which(default_order %in% names(object@reductions))
        reduction <- default_order[reducs[1]]
    }
    if (length(x = dims) != 2 || !is.numeric(x = dims)) {
        stop("'dims' must be a two-length integer vector")
    }
    if (blend && length(x = features) != 2) {
        stop("Blending feature plots only works with two features")
    }
    if (blend) {
        if (length(x = cols) > 2) {
            warning("Blending feature plots only works with two colors; using first two colors", 
                    call. = FALSE, immediate. = TRUE)
        } else if (length(x = cols) < 2) {
            warning("Blended feature plots require two colors, using default colors", 
                    call. = FALSE, immediate. = TRUE)
            cols <- c("#ff0000", "#00ff00")
        }
    }
    if (blend && length(x = cols) != 2) {
        stop("Blending feature plots only works with two colors")
    }
    dims <- paste0(Key(object = object[[reduction]]), dims)
    cells <- cells %||% colnames(x = object)
    data <- FetchData(object = object, vars = c(dims, "ident", features), 
                      cells = cells, slot = slot)
    if (ncol(x = data) < 4) {
        stop("None of the requested features were found: ", 
             paste(features, collapse = ", "), " in slot ", slot, 
             call. = FALSE)
    } else if (!all(dims %in% colnames(x = data))) {
        stop("The dimensions requested were not found", call. = FALSE)
    }
    features <- colnames(x = data)[4:ncol(x = data)]
    min.cutoff <- mapply(FUN = function(cutoff, feature) {
        return(ifelse(test = is.na(x = cutoff), yes = min(data[,feature]), no = cutoff))
    }, cutoff = min.cutoff, feature = features)
    max.cutoff <- mapply(FUN = function(cutoff, feature) {
        return(ifelse(test = is.na(x = cutoff), yes = max(data[,feature]), no = cutoff))
    }, cutoff = max.cutoff, feature = features)
    check.lengths <- unique(x = vapply(X = list(features, min.cutoff, 
                                                max.cutoff), FUN = length, FUN.VALUE = numeric(length = 1)))
    if (length(x = check.lengths) != 1) {
        stop("There must be the same number of minimum and maximum cuttoffs as there are features")
    }
    brewer.gran <- ifelse(test = length(x = cols) == 1, yes = brewer.pal.info[cols, 
                                                                              ]$maxcolors, no = length(x = cols))
    data[, 4:ncol(x = data)] <- sapply(X = 4:ncol(x = data), 
                                       FUN = function(index) {
                                           data.feature <- as.vector(x = data[, index])
                                           min.use <- Seurat:::SetQuantile(cutoff = min.cutoff[index - 
                                                                                          3], data.feature)
                                           max.use <- Seurat:::SetQuantile(cutoff = max.cutoff[index - 
                                                                                          3], data.feature)
                                           data.feature[data.feature < min.use] <- min.use
                                           data.feature[data.feature > max.use] <- max.use
                                           if (brewer.gran == 2) {
                                               return(data.feature)
                                           }
                                           data.cut <- if (all(data.feature == 0)) {
                                               0
                                           } else {
                                               as.numeric(x = as.factor(x = cut(x = as.numeric(x = data.feature), 
                                                                                breaks = brewer.gran)))
                                           }
                                           return(data.cut)
                                       })
    colnames(x = data)[4:ncol(x = data)] <- features
    rownames(x = data) <- cells
    data$split.by <- if (is.null(x = split.by)) {
        RandomName()
    } else {
        switch(EXPR = split.by, ident = Idents(object = object)[cells], 
               object[[split.by, drop = TRUE]][cells])
    }
    if (!is.factor(x = data$split.by)) {
        data$split.by <- factor(x = data$split.by)
    }
    if (!is.null(x = shape.by)) {
        data[, shape.by] <- object[[shape.by, drop = TRUE]]
    }
    xlims <- c(floor(x = min(data[, dims[1]])), ceiling(x = max(data[, dims[1]])))
    ylims <- c(floor(min(data[, dims[2]])), ceiling(x = max(data[, dims[2]])))
    if (blend) {
        ncol <- 4
        color.matrix <- BlendMatrix(two.colors = cols, col.threshold = blend.threshold)
        colors <- list(color.matrix[, 1], color.matrix[1, ], 
                       as.vector(x = color.matrix))
    }
    plots <- vector(mode = "list", length = ifelse(test = blend, 
                                                   yes = 4, 
                                                   no = length(x = features) ))
    for (j in 1:length(x = features)) {
            feature <- features[j]
            if (blend) {
                cols.use <- as.numeric(x = as.character(x = data.plot[, 
                                                                      feature])) + 1
                cols.use <- colors[[j]][sort(x = unique(x = cols.use))]
            } else {
                cols.use <- NULL
            }
            
            plot <- SingleDimPlot(data = data[, c(dims,"ident", feature, "split.by", shape.by)], 
                                   dims = dims, col.by = feature, 
                                  order = order, pt.size = pt.size, cols = cols.use, 
                                  shape.by = shape.by, label = FALSE) + scale_x_continuous(limits = xlims) + 
                scale_y_continuous(limits = ylims) +  cowplot::theme_cowplot()
        
            if (label) {
                plot <- LabelClusters(plot = plot, id = "ident", 
                                       repel = repel, size = label.size)
            }
            if (!is.null(x = split.by)) {
                plot <- plot + Seurat:::FacetTheme() + facet_wrap(facets = vars(!!sym(x = "split.by")),
                                                                  ncol = if (length(x = group.by) > 1 || is.null(x = ncol)) {
                                                                      length(x = unique(x = data[, "split.by"]))
                                                                  } else {
                                                                      ncol
                                                                  })
            }
            if (!blend) {
                plot <- plot + guides(color = NULL)
                cols.grad <- cols
                if (length(x = cols) == 1) {
                    plot <- plot + scale_color_brewer(palette = cols)
                } else if (length(x = cols) > 1) {
                    if (all(data[, feature] == data[,feature][1])) {
                        warning("All cells have the same value (", 
                                data[1, feature], ") of ", feature, 
                                ".")
                        if (data[1, feature][1] == 0) {
                            cols.grad <- cols[1]
                        }
                        else {
                            cols.grad <- cols
                        }
                    }
                    plot <- suppressMessages(expr = plot + scale_color_gradientn(colors = cols.grad, 
                                                                                 guide = "colorbar"))
                }
            }
            if (coord.fixed) {
                plot <- plot + coord_fixed()
            }
            plot <- plot + ylab(feature)
            plots[[j]] <- plot
        }
    if (blend) {
        blend.legend <- BlendMap(color.matrix = color.matrix)
        for (ii in 1:length(x = levels(x = data$split))) {
            suppressMessages(
                expr = plots <- append(x = plots, 
                           values = list(blend.legend + 
                              scale_y_continuous(sec.axis = dup_axis(name = ifelse(test = length(x = levels(x = data$split)) > 1, 
                                                yes = levels(x = data$split)[ii], no = "")), 
                                                expand = c(0, 0)) + labs(x = features[1], 
                                                 y = features[2], title = if (ii == 1) {
                                                     paste("Color threshold:", blend.threshold)
                                                 } else {
                                                     NULL
                                                 }) + no.right), after = 4 * ii - 1))
        }
    }
    plots <- Filter(f = Negate(f = is.null), x = plots)
    if (combine) {
        legend <- if (blend) {
            "none"
        } else {
            split.by %iff% "none"
        }
            plots <- CombinePlots(plots = plots, ncol = 1, 
                                  legend = legend, nrow = length(features))
    }
    return(plots)
}

FeaturePlot.1(object,features = c("CCND1","CCND2"),split.by = "orig.ident")

# re-calculate the cdr3 frequency in data.frame
Frequency <- function(df, key = "cdr3",top=NULL,remove.na =T,
                      remove.dup =T, verbose = FALSE){
        
        #colnames(df) = tolower(colnames(df))
        if(any(colnames(df) %in% "frequency")){
                df = df[,-which(colnames(df) %in% "frequency")]
        }
        cdr3_table <- table(df[,key]) %>% as.data.frame
        if(verbose) table(cdr3_table$Freq) %>% print
        Sum <- sum(cdr3_table$Freq)
        if(colnames(df)[1] == "ngene"|colnames(df)[1] == "nGene"){
        cdr3_table$frequency = cdr3_table$Freq/Sum
        }
        cdr3_table = cdr3_table[order(cdr3_table$Freq,decreasing = T),]
        if(remove.na) cdr3_table = cdr3_table[cdr3_table$Var1!="na",]
        colnames(cdr3_table)[1] = "cdr3"
        df_new <- left_join(df, cdr3_table,by = "cdr3")
        if(remove.dup) df_new = df_new[!duplicated(df_new$cdr3),]
        df_new$Freq[is.na(df_new$Freq)]=0
        df_new$frequency[is.na(df_new$frequency)]=0
        df_new = df_new[order(df_new$freq,decreasing = T),]
        if(is.null(top)) top = nrow(df_new)
        top = min(top,nrow(df_new))
        return(df_new[1:top,])
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

#' select 10% of cell from control
#' Seurat 3
ScaleDown.1 <- function(object, control=c("BH","DJ","MD","NZ"), scale.down=0.1){
    
    normal_cells = lapply(control, function(x){
        rownames(object@meta.data)[(object@meta.data$orig.ident %in% x)]
    })
    set.seed(101)
    remove_normal_cells = lapply(normal_cells, function(x) {
        sample(x, size = length(x)*(1-scale.down))
        }) %>% unlist
    table(colnames(object) %in% remove_normal_cells)
    cell.use <- colnames(object)[!(colnames(object) %in% remove_normal_cells)]
    object <- object[, cell.use]
    object@meta.data$orig.ident = gsub(paste(control,collapse = "|"),"Normal",object@meta.data$orig.ident)
    
    return(object)
}


#' Generate Trends of top 20 TCR clones in longitudinal samples
#' @param TCR_data TCR data.frame report. must contain columns like c("orig.ident","cdr3","frequency")
#' @param order.by character vector indicates sample longitudinal order
#' @param group character, the summary of all samples. For output figure titile only
#' @param size ggplot point size
#' @param top integer, include top N TCR clones from each time points
#' @param x.margin integer, blank margin on x axis
#' @example order.by = c("Pt17_C03_TCRB","Pt17_C07_TCRB","Pt17_C31_TCRB")
#' TCR_Longitud_Plot(Pt_17_TCR,order.by, "Pt-17 bulk")
TCR_Longitud_Plot <- function(TCR_data, order.by, group,color=NULL,
size=2,colors = singler.colors,x.margin=0.125,
do.return = TRUE, do.print = FALSE,do.log =TRUE, top=20){
    
    colnames(TCR_data) = tolower(colnames(TCR_data))
    TCR_data$orig.ident = as.character(TCR_data$orig.ident)
    TCR_data <- split(TCR_data, TCR_data$orig.ident) %>%
    lapply(function(x) Frequency(x,top=top)) %>% bind_rows
    keep_col <- c("orig.ident","cdr3","frequency")
    TCR_data = TCR_data[order(TCR_data$frequency,decreasing = T),]
    TCR_data <- TCR_data[,keep_col]
    # fill up 0
    TCR_data %<>% spread(key="orig.ident", value="frequency",fill = 0)
    TCR_data %<>% gather(key="orig.ident",value="frequency",-cdr3)
    
    TCR_data$orig.ident = factor(TCR_data$orig.ident,levels = order.by)
    TCR_data$samples <- as.integer(TCR_data$orig.ident)
    g1 <- ggplot(data=TCR_data, aes(x=samples, y=frequency, group=cdr3,colour=cdr3)) +
    geom_line(size = size) +
    geom_point(size = size)+
    #scale_fill_manual(values = colors)+
    theme_minimal()+
    xlab("longitudinal samples")+
    scale_x_continuous(limits=c(min(TCR_data$samples)-x.margin,
    max(TCR_data$samples)+x.margin),
    breaks=(min(TCR_data$samples)):(max(TCR_data$samples)),
    labels=order.by)+
    scale_y_continuous(labels = scales::percent)+
    ggtitle(paste("Trends of top",top,"TCR clones in longitudinal",group,
    "samples"))+
    theme(text = element_text(size=15),
    legend.position="none",
    plot.title = element_text(hjust = 0.5),
    axis.text.x  = element_text(angle=70, vjust=0.5))
    if(do.log) g1 = g1 + scale_y_log10()
    if(do.print){
        path <- paste0("output/",gsub("-","",Sys.Date()),"/")
        if(!dir.exists(path)) dir.create(path, recursive = T)
        jpeg(paste0(path,group,"_TCR.jpeg"), units="in", width=10, height=7,res=600)
        print(g1)
        dev.off()
    }
    if(do.return) return(g1)
}


#' Generate Trends of top 20 TCR clones in longitudinal samples
#' @param TCR_data TCR data.frame report. must contain columns like c("orig.ident","cdr3","frequency")
#' @param order.by length =2character vector indicates pariwise sample at x-axis and y-axis
#' @param top integer, include top N TCR clones from each time points
#' @param size ggplot point size
#' @param margin integer, blank margin on both x and y axis
#' @example TCRPairPlot(Pt_17_TCR,c("Pt17_C03_TCRB","Pt17_C07_TCRB"))
TCRPairPlot <- function(TCR_data, order.by, size=2,margin=0.1,
do.return = TRUE, do.print = FALSE,top=20){
    colnames(TCR_data) = tolower(colnames(TCR_data))
    TCR_data = TCR_data[(TCR_data$orig.ident %in% order.by),]
    TCR_data <- split(TCR_data, TCR_data$orig.ident) %>%
    lapply(function(x) Frequency(x,top = top)) %>% bind_rows
    keep_col <- c("orig.ident","cdr3","frequency")
    TCR_data = TCR_data[order(TCR_data$frequency,decreasing = T),]
    TCR_data <- TCR_data[,keep_col]
    # fill up 0
    TCR_data %<>% spread(key="orig.ident", value="frequency",fill = 0)
    TCR_data = TCR_data[,1:3]
    colnames(TCR_data)[2:3] =c("sample1","sample2")
    TCR_data$counts = as.numeric(TCR_data$sample1>0)+
    as.numeric(TCR_data$sample2>0)
    
    g1 <- ggplot(data=TCR_data, aes(x=sample1, y=sample2,
    group=counts, color=counts)) +
    geom_point(size = 2)+
    scale_fill_manual(values = singler.colors)+
    theme_minimal()+
    xlab(paste("frequency in",order.by[1]))+
    ylab(paste("frequency in",order.by[2]))+
    scale_x_sqrt(labels = scales::percent,expand=c(0,0),
    limits=c(0,max(TCR_data$sample1)*(1+margin)))+
    scale_y_sqrt(labels = scales::percent,expand=c(0,0),
    limits=c(0,max(TCR_data$sample2)*(1+margin)))+
    ggtitle(paste("Top",top,"TCR clones between",
    paste(order.by,collapse = " and ")))+
    theme(text = element_text(size=15),
    legend.position="none",
    plot.title = element_text(hjust = 0.5),
    axis.text.x  = element_text(angle=70, vjust=0.5))
    
    if(do.print){
        path <- paste0("output/",gsub("-","",Sys.Date()),"/")
        if(!dir.exists(path)) dir.create(path, recursive = T)
        jpeg(paste0(path,paste0(order.by,collapse = "_"),"_pairTCR.jpeg"), units="in", width=10, height=7,res=600)
        print(g1)
        dev.off()
    }
    if(do.return) return(g1)
}


#' Combine TSNEPlot and Feature Plot into one figure.
#' @param features Vector of features to plot, parameter pass to FeaturePlot
#' @param no.legend remove legend
#' @param title add ggplot title
#' @param do.print save jpeg file
#' @param unique.name save jpeg file with unique name
#' @param do.return return plot
TSNEPlot_FeaturePlot <- function(object,dims = c(1, 2),cells = NULL,cols = NULL,pt.size = NULL,
                       reduction = NULL,group.by = NULL,split.by = NULL,shape.by = NULL,
                       features = NULL,
                       order = NULL,label = FALSE,label.size = 4,repel = FALSE,
                       cells.highlight = NULL,cols.highlight = 'red',sizes.highlight = 1,
                       na.value = 'grey50',ncol = NULL,title = NULL,
                       no.legend = F,do.print = F,do.return = T,unique.name=F,
                       width=10, height=7,...) {
    if(unique.name) {
        VarName <- unique(object$orig.ident)
        VarName <- paste(VarName[1:min(5,length(VarName))],collapse = "_")
    } else VarName <- deparse(substitute(object))
    VarName = paste0(VarName,"_",group.by %||% FindIdentLabel(object))
    plots <- TSNEPlot(object = object,dims = dims, cells = cells, cols = cols, 
                      pt.size = pt.size, reduction = 'tsne', group.by = group.by, split.by = split.by, 
                      shape.by = shape.by, order = order, label = label, label.size = label.size, 
                      repel = repel, cells.highlight = cells.highlight, cols.highlight = cols.highlight, 
                      sizes.highlight = sizes.highlight, na.value = na.value, combine = combine, 
                      ncol = ncol, ...)
    if(no.legend) {
        plots = plots + NoLegend()
        L = ""
    } else L = "_Legend"
    featureplot <- FeaturePlot(subset_object, features = features, split.by = split.by)
    if (combine) {
        plots <- CombinePlots(
        plots = plots,align= align,
        ncol = if (!is.null(x = split.by) && length(x = group.by) > 1) {
            1
        } else {
            ncol
        },
        ...
        )
    }
    if(!is.null(title)){
        plots = plots + ggtitle(title)+
        theme(plot.title = element_text(hjust = 0.5,size=15,face = "plain"))
    }
    if(no.legend) {
        plots = plots + NoLegend()
        L = ""
    } else L = "_Legend"
    if(do.print) {
        path <- paste0("output/",gsub("-","",Sys.Date()),"/")
        if(!dir.exists(path)) dir.create(path, recursive = T)
        jpeg(paste0(path,"TSNEPlot_",VarName,L,".jpeg"),
        units="in", width=width, height=height,res=600)
        print(plots)
        dev.off()
    }
    if(do.return) return(plots)
}
