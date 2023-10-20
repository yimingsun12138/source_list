###################################################
## Project: sc_multiomics                        ##
## Script Purpose: collect the function I wrote  ##
## Data: 2022.07.17                              ##
## Author: Yiming Sun                            ##
###################################################

# Seurat related ----------------------------------------------------------

#' SCT integrate function in Seurat
#' @param object.list Seurat object list to be integrated
#' @param nfeature variable genes number for each seurat object
My_SCT_Integrate <- function(object.list,nfeature=3000){
  require(Seurat)
  for (i in 1:length(object.list)) {
    object.list[[i]] <- SCTransform(object.list[[i]], verbose = T)
  }
  object.features <- SelectIntegrationFeatures(object.list = object.list, nfeatures = nfeature)
  object.list <- PrepSCTIntegration(object.list = object.list, anchor.features = object.features, verbose = T)
  object.anchors <- FindIntegrationAnchors(object.list = object.list, normalization.method = "SCT", anchor.features = object.features, verbose = T)
  object.integrated <- IntegrateData(anchorset = object.anchors, normalization.method = "SCT", verbose = T)
  DefaultAssay(object.integrated) <- 'integrated'
  object.integrated <- RunPCA(object.integrated, npcs = 50, verbose = TRUE)
  return(object.integrated)
}

#' Function to plot gene expression from a seurat object
#' @param genes a vector of genes
#' @param object a seurat object or an express matrix
#' @param object.class the type of object."seurat" or "matrix"
#' @param assay the assay to pull data from
#' @param projection if your object is just a matrix, this is the coordinate 
#' projection
#' @param embedding In Seurat object, what is the Embedding name
#' @param reduction_key prefix of the embedding colname
#' @param plot.title Plot title
#' @param guide.name Legend title (ggplot2)
#' @param plot.type Whether you are plotting 1 gene, the average of the genes 
#' listed, or individual panels for each of the genes listed (currently uses 
#' ggplot2 facets, cowplot grids coming)
#' @param w A vector of weights for gene expression averages. For example, 
#' cluster membership scores. 
#' Automatically computes a weighted mean if not NULL.
#' @param scaled Whether or not to scale gene expression across cells
#' @param color.palette A vector of colors, passed to scale_color_gradientn()
#' @param aspectratio The desired aspect ratio of each plot/panel
#' @param point.size custom point size
#' @param trim For plotting large datasets, trims the uppermost and lowermost
#' percentiles, eg c(0.01, 0.99) trims the highest and lowest percentile, so
#' the color scale is not driven by outlier points
#' @param lims X / Y limits of coordinates. subset your projection space
#' @param print Whether the plot is passed through "print"
#' @param num.panel.rows For facet plots, how many rows of facets?
geneSetAveragePlot <- function(genes, object, object.class = "seurat", assay = 'RNA', 
                               projection = NULL, embedding = 'umap', reduction_key = 'UMAP',
                               plot.title = "Gene Expression", guide.name = "Gene expression", 
                               plot.type = c("single", "averages", "panels"), w = NULL, scaled = FALSE, 
                               color.palette = c("#B0A1AF","#441F26","#7D56AA","#A176CF","#DF9418","#E5C252"), 
                               aspectratio = 1.25, point.size = 0.3, trim = NULL, lims = NULL, print = TRUE, 
                               num.panel.rows = NULL) {
  #load library
  require(Seurat)
  require(ggplot2)
  require(dplyr)
  require(reshape2)
  
  gene.ids <- genes
  # Object class handling
  
  object.class <- tolower(object.class)
  
  if (object.class == "seurat") {
    DefaultAssay(object) <- assay
    rna.mat <- expm1(object[[assay]]@data)
    umap <- Seurat::Embeddings(object,embedding)
  }
  if (object.class == 'matrix') {
    rna.mat <- object
    if (is.null(projection)) {
      stop("Must provide a projection with matrix class")
    }
    umap <- projection
  }
  
  # Check again
  if (!all(gene.ids %in% rownames(rna.mat))) {
    print(gene.ids[!gene.ids %in% rownames(rna.mat)])
  }
  
  # Make submatrix
  gene.exp <- rna.mat[gene.ids, ]
  dims <- grep(reduction_key,colnames(umap),value = T)
  dim1 <- dims[1]
  dim2 <- dims[2]
  
  umap <- umap[,c(dim1,dim2)]
  
  # Trim
  if (! is.null(trim)) {
    gene.exp <- trimQuantiles(gene.exp,cuts = trim)
  }
  
  # Handle plot types. Make plot data frames
  plot.type <- match.arg(plot.type)
  
  if(plot.type == "single" || (length(gene.ids) == 1)) {
    gene.exp <- log1p(gene.exp)
    if(scaled == TRUE){
      gene.exp <- rowScale(gene.exp)
    }
    pl.df <- data.frame(umap,geneexpr = gene.exp)
    plot.title <- genes
  } 
  
  if (plot.type == "averages" & (length(gene.ids) > 1)) {
    if (is.null(w)) {
      gene.exp <- colMeans(gene.exp)
    } else {
      gene.exp <- apply(
        X = gene.exp, 
        MARGIN = 2,
        FUN = weighted.mean,
        w = w
      )
    }
    
    gene.exp <- log1p(gene.exp)
    
    if(scaled == TRUE){
      gene.exp <- rowScale(gene.exp)
    }
    
    pl.df <- data.frame(umap,geneexpr = gene.exp)
    
  }
  
  if(plot.type == "panels" & length(gene.ids) >1 ) {
    
    gene.exp <- log1p(gene.exp)
    
    if(scaled == TRUE){
      gene.exp <- rowScale(gene.exp)
    }
    
    gene.exp <- t(as.matrix(gene.exp))
    colnames(gene.exp) <- genes
    
    cat(sprintf("%s genes; %s cells\n",dim(gene.exp)[2],dim(gene.exp)[1]))
    
    pl.df2 <- as.data.frame(cbind(umap,gene.exp))
    pl.df <- melt(pl.df2,id.vars = c(dim1,dim2),variable.name = "Gene",value.name = "geneexpr")
    
    if (is.null(num.panel.rows) ) {
      num.panel.rows <- pmax(1,floor(length(genes)/3))
    }
    
  }
  
  # Plot
  gg <- ggplot(pl.df, aes_string(x = dim1, y = dim2, color = 'geneexpr')) + 
    theme_classic() + 
    scale_color_gradientn(colours = color.palette, name = guide.name) +
    labs(title = plot.title) +
    theme(plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(colour = "black", fill=NA, size=1))
  
  gg <- gg + geom_point(size = point.size, shape = 16) + theme(aspect.ratio = aspectratio)
  
  if (plot.type == "panels") {
    gg <- gg + facet_wrap(~Gene, nrow = num.panel.rows) + theme(axis.line = element_blank())
  }
  
  if (!is.null(lims)) {
    x.lims <- lims[[1]]
    y.lims <- lims[[2]]
    gg <- gg + xlim(x.lims) + ylim(y.lims)
  }
  
  if (print == TRUE) { 
    print(gg)
  } else {
    return(gg)
  }
}

#' origin from seurat DotPlot function
#' @param features Input vector of features, or named list of feature vectors if feature-grouped panels are desired
#' @param col.min Minimum scaled average expression threshold
#' @param col.max Maximum scaled average expression threshold
#' @param dot.min The fraction of cells at which to draw the smallest dot (default is 0). All cell groups with less than this expressing the given gene will have no dot drawn
#' @param dot.scale the maximum size of dots
#' @param group.by same as group by
#' @param scale.by Scale the size of the points by 'size' or by 'radius'
#' @param scale.min Set lower limit for dot scaling, use NA for default
#' @param scale.max Set upper limit for dot scaling, use NA for default
#' @param return_dot_plot whether to return the data.frame to generate the plots or just generate the plots
#' @param data_plot give the data.frame used to generate the plots and draw the fig directly
#' @param scale whether to scale the expression matrix
my_dotplot <- function (object, assay = NULL, features, cols = c("lightgrey","blue"), 
                        col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, scale = FALSE, 
                        group.by = NULL, scale.by = "radius", scale.min = NA, scale.max = NA, 
                        return_data_plot = FALSE, data_plot = NULL) {
  require(Seurat)
  #default setting modified from Seurat::Dotplot
  idents <- NULL
  scale.func <- switch(EXPR = scale.by, size = scale_size, 
                       radius = scale_radius, 
                       stop("'scale.by' must be either 'size' or 'radius'"))
  #all settled
  
  if(is.null(data_plot)){
    PercentAbove <- function(x, threshold) {
      return(length(x = x[x > threshold]) / length(x = x))
    }
    assay <- assay %||% DefaultAssay(object = object)
    DefaultAssay(object = object) <- assay
    feature.groups <- NULL
    if (is.list(features) | any(!is.na(names(features)))) {
      feature.groups <- unlist(x = sapply(X = 1:length(features), 
                                          FUN = function(x) {
                                            return(rep(x = names(x = features)[x], each = length(features[[x]])))
                                          }))
      if (any(is.na(x = feature.groups))) {
        warning("Some feature groups are unnamed.", call. = FALSE, 
                immediate. = TRUE)
      }
      features <- unlist(x = features)
      names(x = feature.groups) <- features
    }
    cells <- unlist(x = CellsByIdentities(object = object, idents = idents))
    data.features <- FetchData(object = object, vars = features, cells = cells)
    data.features$id <- if (is.null(x = group.by)) {
      Idents(object = object)[cells, drop = TRUE]
    }
    else {
      object[[group.by, drop = TRUE]][cells, drop = TRUE]
    }
    if (!is.factor(x = data.features$id)) {
      data.features$id <- factor(x = data.features$id)
    }
    id.levels <- levels(x = data.features$id)
    data.features$id <- as.vector(x = data.features$id)
    
    data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
      data.use <- data.features[data.features$id == ident, 1:(ncol(x = data.features) - 1), drop = FALSE]
      avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
        return(mean(x = expm1(x = x)))
      })
      pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, threshold = 0)
      return(list(avg.exp = avg.exp, pct.exp = pct.exp))
    })
    
    names(x = data.plot) <- unique(x = data.features$id)
    
    data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
      data.use <- as.data.frame(x = data.plot[[x]])
      data.use$features.plot <- rownames(x = data.use)
      data.use$id <- x
      return(data.use)
    })
    data.plot <- do.call(what = "rbind", args = data.plot)
    if (!is.null(x = id.levels)) {
      data.plot$id <- factor(x = data.plot$id, levels = id.levels)
    }
    if (length(x = levels(x = data.plot$id)) == 1) {
      scale <- FALSE
      warning("Only one identity present, the expression values will be not scaled",
              call. = FALSE, immediate. = TRUE)
    }
    avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot),
                             FUN = function(x) {
                               data.use <- data.plot[data.plot$features.plot == x, "avg.exp"]
                               if (scale) {
                                 data.use <- scale(x = data.use,center = TRUE)
                                 data.use <- MinMax(data = data.use, min = col.min, max = col.max)
                               }
                               else {
                                 data.use <- log1p(x = data.use)
                               }
                               return(data.use)
                             })
    avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
    
    data.plot$avg.exp.scaled <- avg.exp.scaled
    data.plot$features.plot <- factor(x = data.plot$features.plot, levels = features)
    data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
    data.plot$pct.exp <- data.plot$pct.exp * 100
    
    color.by <- "avg.exp.scaled"
    if (!is.na(x = scale.min)) {
      data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
    }
    if (!is.na(x = scale.max)) {
      data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
    }
    if (!is.null(x = feature.groups)) {
      data.plot$feature.groups <- factor(x = feature.groups[data.plot$features.plot], 
                                         levels = unique(x = feature.groups))
    }
    if(return_data_plot){
      return(data.plot)
    } else{
      plot <- ggplot(data = data.plot, mapping = aes_string(x = "features.plot", y = "id")) + 
        geom_point(mapping = aes_string(size = "pct.exp", color = color.by)) + 
        scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) + 
        theme(axis.title.x = element_blank(), 
              axis.title.y = element_blank()) + 
        guides(size = guide_legend(title = "Percent Expressed")) + 
        labs(x = "Features", y = "Identity") + 
        theme_cowplot()
      if (!is.null(x = feature.groups)) {
        plot <- plot + facet_grid(facets = ~feature.groups, scales = "free_x", 
                                  space = "free_x", switch = "y") + 
          theme(panel.spacing = unit(x = 1, units = "lines"), strip.background = element_blank())
      }
      
      plot <- plot + scale_colour_gradientn(colours = cols)
      plot <- plot + guides(color = guide_colorbar(title = "Average Expression"))
      return(plot)
    }
  } else{
    data.plot <- data_plot
    color.by <- "avg.exp.scaled"
    plot <- ggplot(data = data.plot, mapping = aes_string(x = "features.plot", y = "id")) + 
      geom_point(mapping = aes_string(size = "pct.exp", color = color.by)) + 
      scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) + 
      theme(axis.title.x = element_blank(), 
            axis.title.y = element_blank()) + 
      guides(size = guide_legend(title = "Percent Expressed")) + 
      labs(x = "Features", y = "Identity") + 
      theme_cowplot()
    if (!is.null(x = data.plot$feature.groups)) {
      plot <- plot + facet_grid(facets = ~feature.groups, scales = "free_x", 
                                space = "free_x", switch = "y") + 
        theme(panel.spacing = unit(x = 1, units = "lines"), strip.background = element_blank())
    }
    
    plot <- plot + scale_colour_gradientn(colours = cols)
    plot <- plot + guides(color = guide_colorbar(title = "Average Expression"))
    return(plot)
  }
}

#' from greenleaf
#' get uwot model from seurat object
#' @param object the seurat object to pull uwot model from
#' @param assayUsed the assay used for umap projection in object
#' @param reductionUsed the dim-reduction object in seurat
getUwotModelFromSeurat <- function(object, 
                                   assayUsed = DefaultAssay(object), 
                                   reductionUsed = "pca", ...){
  
  require(uwot)
  
  paramL <- object@commands[[paste0("RunUMAP.",assayUsed,".",reductionUsed)]]@params
  
  if (is.null(paramL[["dims"]])) logger.error("Invalid dims")
  reqParamNames <- c(
    "n.neighbors", "n.components", "metric", "n.epochs", "learning.rate", "min.dist",
    "spread", "set.op.mix.ratio", "local.connectivity", "repulsion.strength",
    "negative.sample.rate", "a", "b", "uwot.sgd", "seed.use"
  )
  
  # assign NULL to missing parameters
  for (pn in reqParamNames){
    if (!is.element(pn, names(paramL))) paramL[pn] <- list(NULL)
  }
  
  X <- Embeddings(object[[reductionUsed]])[, paramL[["dims"]]]
  
  # make sure to use the same 'random' numbers
  
  if (!is.null(paramL$seed.use)) {
    set.seed(seed = paramL$seed.use)
  }
  
  umapRes <- umap(
    X = X,
    n_neighbors = as.integer(paramL$n.neighbors),
    n_components = as.integer(paramL$n.components),
    metric = paramL$metric,
    n_epochs = paramL$n.epochs,
    learning_rate = paramL$learning.rate,
    min_dist = paramL$min.dist,
    spread = paramL$spread,
    set_op_mix_ratio = paramL$set.op.mix.ratio,
    local_connectivity = paramL$local.connectivity,
    repulsion_strength = paramL$repulsion.strength,
    negative_sample_rate = paramL$negative.sample.rate,
    a = paramL$a,
    b = paramL$b,
    fast_sgd = paramL$uwot.sgd,
    ret_model=TRUE, # this is the important part
    ...
  )
  
  return(umapRes)
}

#' from greenleaf
#' project new data on the already existed UMAP
#' @param X_scaled the scaled express matrix to be projected on object embedding
#' @param object the seurat object that has already been visualized by seurat_UMAP
#' @param umap_model the uwot model used in object, set NULL to pull from object
#' @param assayUsed the default assay of object
#' @param missing_gene any missing genes in X? set TRUE and the missing gene loading value in PCA will be removed(better not)
projectMatrix_SeuratUMAP <- function(X_scaled, 
                                     object, 
                                     umap_model = NULL,
                                     assayUsed = DefaultAssay(object), 
                                     missing_gene = FALSE) {
  require(Seurat)
  require(dplyr)
  require(uwot)
  
  #check input
  paramL_pca <- object@commands[[paste0("RunPCA.",assayUsed)]]@params
  reductionUsed <- paramL_pca[['reduction.name']]
  paramL_umap <- object@commands[[paste0("RunUMAP.",assayUsed,".",reductionUsed)]]@params
  
  cat("PCA projection ...\n")
  X <- X_scaled
  features <- rownames(object@reductions[[reductionUsed]]@feature.loadings)
  
  if(missing_gene){
    char <- length(features)
    features <- dplyr::intersect(features,rownames(X))
    char <- char - length(features)
    char <- paste(as.character(char),'genes missing!',sep = ' ')
    print(char)
  } else{
    if (any(!(features %in% rownames(X)))){
      stop("Could not find all features in X_scaled for the PC projection")
    }
  }
  
  X <- t(X[features,])
  
  projM <- Loadings(object, reduction = reductionUsed)
  projM <- projM[features,]
  pcaCoord_proj <- X %*% projM
  
  pcaCoord_orig <- t(GetAssayData(object[[assayUsed]], slot="scale.data")[features,]) %*% projM
  
  cat("Retrieving UMAP model ...\n")
  if (is.null(umap_model)) {
    umapRes <- getUwotModelFromSeurat(object, assayUsed=assayUsed, reductionUsed=reductionUsed)
  } else {
    umapRes <- umap_model
  }
  
  cat("UMAP projection ...\n")
  umapCoord_orig <- uwot::umap_transform(pcaCoord_orig[ ,paramL_umap[["dims"]]], umapRes)
  colnames(umapCoord_orig) <- colnames(Embeddings(object[[paramL_umap[["reduction.name"]]]]))
  rownames(umapCoord_orig) <- rownames(pcaCoord_orig)
  
  umapCoord_proj <- uwot::umap_transform(pcaCoord_proj[ ,paramL_umap[["dims"]]], umapRes)
  
  rownames(umapCoord_proj) <- rownames(pcaCoord_proj)
  colnames(umapCoord_proj) <- colnames(umapCoord_orig)
  
  res <- list(
    pcaCoord_proj=pcaCoord_proj,
    umapCoord_proj=umapCoord_proj,
    pcaCoord_orig=pcaCoord_orig,
    umapCoord_orig=umapCoord_orig
  )
  
  return(res)
}

#' add gene module score
#' @param seu.obj the seurat object
#' @param assay the assay to pull data from
#' @param features the gene list of the module
#' @param meta.var the col you want to store module score in the seurat meta.data
#' @param scale whether to scale the express, if FALSE, return log1p(average counts), if TRUE, return scaled score
#' @param center whether to center the express
My_add_module_score <- function(seu.obj,assay = 'RNA',features,meta_var,scale = TRUE,center = TRUE) {
  #load library
  require(Seurat)
  require(dplyr)
  
  #extract express and gene list
  express_matrix <- expm1(seu.obj[[assay]]@data)
  valid_features <- dplyr::intersect(features,rownames(express_matrix))
  if(length(valid_features) == 0){
    stop('no feature intersect with express matrix!')
  } else{
    temp <- length(features) - length(valid_features)
    temp <- paste0(as.character(temp),' genes missing!')
    print(temp)
  }
  features <- valid_features
  
  #calculate mean express
  if(length(features) == 1){
    warning('only 1 gene in the module!')
    avg.express <- express_matrix[features,]
  } else{
    avg.express <- express_matrix[features,]
    avg.express <- colMeans(avg.express)
  }
  
  #normalize
  avg.express <- log1p(avg.express)
  avg.express <- scale(avg.express,center = center,scale = scale)
  
  #return
  seu.obj@meta.data[,as.character(meta_var)] <- avg.express
  return(seu.obj)
}

#' pre-process of seurat object
#' @param object the seurat object
#' @param assay the assay to use
#' @param reduction.name name of the reduction object to store the runpca output or to be used for findneighbours and runumap
#' @param nfeatures variable feature number by vst
#' @param npcs PCA dims
#' @param preprocess whether pre-process(normalize, scale and pca) or later-process(cluster and umap)
#' @param dim_to_use use 1:dim_to_use dims to find neighbors and run umap
#' @param resolution the resolution used to find cluster
#' @param group.by same as param group.by in dimplot
#' @param label whether show label on dimplot
#' @param variable.feature set the variable features of the seurat object, default is NULL
#' @param vars.to.regress var name in metadata, regress out the influence brought by this parameter
my_process_seurat <- function(object,assay='RNA',reduction.name='pca',variable.feature = NULL,
                              nfeatures=2000,vars.to.regress = NULL,npcs=50,
                              preprocess=TRUE,dim_to_use=30,resolution=0.8,
                              group.by='seurat_clusters',label = TRUE){
  require(Seurat)
  DefaultAssay(object) <- assay
  if(preprocess){
    object <- NormalizeData(object,normalization.method = "LogNormalize",scale.factor = 10000)
    if(is.null(variable.feature)){
      object <- FindVariableFeatures(object,selection.method = "vst",nfeatures=nfeatures)
    } else{
      VariableFeatures(object) <- variable.feature
    }
    object <- ScaleData(object,vars.to.regress = vars.to.regress)
    object <- RunPCA(object,npcs = npcs,reduction.name = reduction.name)
    print(ElbowPlot(object,ndims = npcs,reduction=reduction.name))
    return(object)
  } else{
    object <- FindNeighbors(object,dims = 1:dim_to_use,reduction=reduction.name)
    object <- FindClusters(object,resolution = resolution)
    object <- RunUMAP(object,dims = 1:dim_to_use,reduction=reduction.name)
    p <- DimPlot(object,group.by = group.by,label = label,repel = TRUE)
    print(p)
    return(object)
  }
}

#' integrate seurat object using harmony.
#' default method for integrating seurat objects in harmony scales express matrix all together, here we scale separately.
#' @param named_seurat_list named list of seurat objects for integration.
#' @param assay the seurat assay used for integration, note that the list of seurat objects need to contain this assay all, and the assay should contain all the variable_feature.
#' @param variable_feature variable features to do scale and pca.
#' @param var_to_regress_list named list of var names to regress out while scaling, names should be the same as named_seurat_list.
#' @param npcs total dims while doing pca, deafult is 50.
#' @param reference_loading dataset name used to provide the pca_loading, default is the dataset with the most cells.
#' @param integration_var the var in the metadata of all seurat objects, used to separate dataset.
#' @param harmony_input_dim the pca embedding dims used to do harmony, default is 50.
#' @param max.iter.harmony the maximum rounds of harmony.
#' @param lambda Ridge regression penalty parameter. Specify for each variable in vars_use. Default lambda=1. Lambda must be strictly positive. Smaller values result in more aggressive correction.
#' @param theta Diversity clustering penalty parameter. Specify for each variable in vars_use Default theta=2. theta=0 does not encourage any diversity. Larger values of theta result in more diverse clusters.
#' @param sigma Width of soft kmeans clusters. Default sigma=0.1. Sigma scales the distance from a cell to cluster centroids. Larger values of sigma result in cells assigned to more clusters. Smaller values of sigma make soft kmeans cluster approach hard clustering.
#' @param reference_dataset (Advanced Usage) Defines reference dataset(s). Cells that have batch variables values matching reference_dataset will not be moved.
#' @param UMAP_dim the pca embedding dims used to create UMAP, default is 30.
#' @param resolution seurat find cluster resolution.
my_harmony_integration <- function(named_seurat_list,
                                   assay='RNA',
                                   variable_feature,
                                   var_to_regress_list,
                                   npcs=50,reference_loading=NULL,
                                   integration_var,
                                   harmony_input_dim=50,
                                   max.iter.harmony=10,
                                   lambda=1,theta=2,sigma=0.1,
                                   reference_dataset=NULL,
                                   UMAP_dim=30,resolution=1,...){
  require(Seurat)
  require(dplyr)
  require(harmony)
  
  for (i in names(named_seurat_list)) {
    DefaultAssay(named_seurat_list[[i]]) <- assay
  }
  
  #check data
  for (i in names(named_seurat_list)) {
    temp <- Seurat::GetAssay(named_seurat_list[[i]],assay=assay)
    temp <- rownames(temp)
    if(sum(variable_feature %in% temp) < length(variable_feature)){
      temp <- paste0('variable feature missing in ',as.character(i),'!')
      stop(temp)
    }
  }
  print('variable feature check pass!')
  
  temp <- base::lapply(named_seurat_list,FUN = function(x){
    return(colnames(x))
  })
  temp <- unlist(temp)
  if(sum(duplicated(temp)) > 0){
    stop('duplicated cell name in datasets!')
  }
  print('data check pass!')
  #scale_data
  for (i in names(named_seurat_list)) {
    named_seurat_list[[i]] <- my_process_seurat(object = named_seurat_list[[i]],assay = assay,
                                                variable.feature = variable_feature,
                                                vars.to.regress = var_to_regress_list[[i]],
                                                npcs = npcs,preprocess = TRUE)
  }
  
  #get PCA loading
  if(is.null(reference_loading)){
    temp <- base::lapply(named_seurat_list,FUN = function(x){
      return(dim(x)[2])
    })
    temp <- unlist(temp)
    names(temp) <- names(named_seurat_list)
    reference_loading <- names(which.max(temp))
  }
  
  PCA_loading <- named_seurat_list[[reference_loading]]@reductions$pca@feature.loadings
  
  #create meta_data
  named_meta_data <- list()
  for (i in names(named_seurat_list)) {
    temp <- data.frame(cell_id=colnames(named_seurat_list[[i]]),dataset=named_seurat_list[[i]]@meta.data[,integration_var])
    named_meta_data[[i]] <- temp
  }
  named_meta_data <- do.call(rbind,named_meta_data)
  colnames(named_meta_data) <- c('cell_id','dataset')
  rownames(named_meta_data) <- named_meta_data$cell_id
  
  #reduce dim using PCA loading
  named_PCA_matrix <- list()
  for (i in names(named_seurat_list)) {
    temp <- Seurat::GetAssay(named_seurat_list[[i]],assay=assay)
    temp <- temp@scale.data
    named_PCA_matrix[[i]] <- t(temp[rownames(PCA_loading),]) %*% PCA_loading
  }
  named_PCA_matrix <- do.call(rbind,named_PCA_matrix)
  rownames(named_PCA_matrix) <- rownames(named_meta_data)
  colnames(named_PCA_matrix) <- colnames(PCA_loading)
  
  #run harmony
  if(harmony_input_dim > npcs){
    harmony_input_dim <- npcs
  }
  named_PCA_matrix <- named_PCA_matrix[,1:harmony_input_dim]
  
  harmonized_PCA_matrix <- harmony::HarmonyMatrix(
    data_mat = named_PCA_matrix,
    meta_data = named_meta_data,
    vars_use = 'dataset',
    do_pca = FALSE,
    return_object = FALSE,
    max.iter.harmony = max.iter.harmony,
    plot_convergence = FALSE,
    lambda = lambda,theta=theta,sigma=sigma,
    reference_values=reference_dataset,...
  )
  
  #create seurat object
  gene_list <- list()
  for (i in names(named_seurat_list)) {
    temp <- Seurat::GetAssay(named_seurat_list[[i]],assay=assay)
    gene_list[[i]] <- rownames(temp)
  }
  
  gene_list <- Reduce(intersect,gene_list)
  temp <- length(gene_list)
  temp <- paste0(as.character(temp),' genes intersect across datasets!')
  print(temp)
  
  seurat_obj <- base::lapply(named_seurat_list,FUN = function(x){
    return(Seurat::GetAssay(x,assay=assay)@counts[gene_list,])
  })
  
  cell_list <- base::lapply(named_seurat_list,FUN = function(x){
    return(colnames(x))
  })
  cell_list <- unlist(cell_list)
  
  seurat_obj <- do.call(cbind,seurat_obj)
  colnames(seurat_obj) <- cell_list
  rownames(seurat_obj) <- gene_list
  
  seurat_obj <- CreateSeuratObject(counts = seurat_obj,project = 'integration',assay = "integration",
                                   meta.data = named_meta_data,min.cells = 0,min.features = 0)
  
  seurat_obj[['pca']] <- CreateDimReducObject(embeddings = harmonized_PCA_matrix,global = TRUE,assay = 'integration')
  if(UMAP_dim > harmony_input_dim){
    UMAP_dim <- harmony_input_dim
  }
  seurat_obj <- FindNeighbors(object = seurat_obj,reduction = 'pca',dims = 1:UMAP_dim)
  seurat_obj <- FindClusters(object = seurat_obj,resolution=resolution)
  seurat_obj <- RunUMAP(object = seurat_obj,reduction = 'pca',dims = 1:UMAP_dim)
  return(seurat_obj)
}

#' Using low dimension graph( like PCA, NMF, CCA and so on ) to create mutual nearest neighbor graph between data and query.
#' Return the probability label matrix of query calculated by mnn graph and the corresponding label in data.
#' @param data The reference seurat object.
#' @param query The query seurat object.
#' @param reference_var The colname of the meta.data table storing the reference label.
#' @param reduction The low dimension graph to be used for calculation.
#' @param knn The number of nearest neighbors to be selected between data and query, small knn results a locally prediction of reference label.
my_low_dim_anchor_construction <- function(data,query,reference_var,reduction='pca',knn=30){
  
  require(Seurat)
  require(dplyr)
  
  #create nearest neighbour graph
  query_nearest <- My_FindKNN(data = data@reductions[[reduction]]@cell.embeddings,
                              query = query@reductions[[reduction]]@cell.embeddings,
                              k = knn)
  
  print('query knn graph created!')
  
  data_nearest <- My_FindKNN(data = query@reductions[[reduction]]@cell.embeddings,
                             query = data@reductions[[reduction]]@cell.embeddings,
                             k = knn)
  
  print('data knn graph created!')
  
  #create mutual nearest neighbour
  temp <- do.call(rbind,base::lapply(colnames(query),function(x){
    return(as.numeric(colnames(data) %in% query_nearest$nn.cell[x,]))
  }))
  
  rownames(temp) <- colnames(query)
  colnames(temp) <- colnames(data)
  temp <- as.data.frame(temp)
  nearest_matrix <- temp
  
  temp <- do.call(cbind,base::lapply(colnames(data),function(x){
    return(as.numeric(colnames(query) %in% data_nearest$nn.cell[x,]))
  }))
  rownames(temp) <- colnames(query)
  colnames(temp) <- colnames(data)
  temp <- as.data.frame(temp)
  
  nearest_matrix <- nearest_matrix+temp-1
  nearest_matrix[nearest_matrix <= 0] <- Inf
  
  #filter cells
  nearest_matrix <- do.call(rbind,base::lapply(rownames(query_nearest$nn.cell),FUN = function(x){
    return(as.numeric(nearest_matrix[x,query_nearest$nn.cell[x,]]))
  }))
  rownames(nearest_matrix) <- rownames(query_nearest$nn.cell)
  nearest_matrix <- as.data.frame(nearest_matrix)
  
  #correct query_nearest dist matrix
  query_nearest$nn.dists <- (query_nearest$nn.dists)*nearest_matrix
  rownames(query_nearest$nn.dists) <- rownames(query_nearest$nn.cell)
  print('MNN graph created!')
  
  #create anchor
  group_list <- c(unique(as.character(data@meta.data[,reference_var])),'unknown')
  
  anchor_matrix <- do.call(rbind,base::lapply(colnames(query),function(x){
    nn <- query_nearest$nn.cell[x,]
    nn <- as.character(data@meta.data[nn,reference_var])
    nn.dists <- as.numeric(query_nearest$nn.dists[x,])
    nn.dists <- exp(-nn.dists)
    if(sum(nn.dists) == 0){
      nn <- factor('unknown',levels = group_list)
      nn <- as.numeric(table(nn))
      names(nn) <- group_list
      return(nn)
    }else{
      nn_list <- base::lapply(group_list,FUN = function(label){
        if(sum(nn == label) == 0){
          return(0)
        }else{
          temp <- which(nn == label)
          return(sum(nn.dists[temp])/sum(nn.dists))
        }
      })
      nn_list <- unlist(nn_list)
      nn_list <- as.numeric(nn_list)
      names(nn_list) <- group_list
      return(nn_list)
    }
  }))
  
  colnames(anchor_matrix) <- group_list
  rownames(anchor_matrix) <- colnames(query)
  anchor_matrix <- as.data.frame(anchor_matrix)
  print('create anchor done!')
  gc()
  return(anchor_matrix)
}

#' Predict query cell label based on MNN graph between query and reference( data ).
#' @param data The reference seurat object.
#' @param query The query seurat object.
#' @param reference_var The colname of the meta.data table storing the reference label.
#' @param reduction The low dimension graph to be used for calculation.
#' @param mnn The number of nearest neighbors to be selected between data and query, small mnn results a locally prediction of reference label.
#' @param knn Create knn graph on query to smooth the mnn prediction, how many neighbors to be selected?
#' @param return_query Whether to return query seurat object?
#' @param Iteration Using query data knn graph to smooth the mnn prediction, how many times to iterate?
my_MNN_label_transfer <- function(data,query,reference_var,reduction='pca',mnn=30,knn=50,Iteration=5,return_query=TRUE){
  
  require(Seurat)
  require(dplyr)
  
  #create anchor matrix
  anchor_matrix <- my_low_dim_anchor_construction(data = data,query = query,reference_var = reference_var,reduction = reduction,knn = mnn)
  
  #create knn graph
  nearest <- My_FindKNN(data = query@reductions[[reduction]]@cell.embeddings,
                        query = query@reductions[[reduction]]@cell.embeddings,
                        k = knn+1)
  nearest$nn.idx <- nearest$nn.idx[,-1]
  nearest$nn.dists <- nearest$nn.dists[,-1]
  nearest$nn.cell <- nearest$nn.cell[,-1]
  
  #modify knn graph
  nearest_matrix <- do.call(rbind,base::lapply(colnames(query),FUN = function(x){
    temp <- rep(0,length(colnames(query)))
    names(temp) <- colnames(query)
    dist_list <- nearest$nn.dists[x,]
    if(var(dist_list) == 0){
      temp[nearest$nn.cell[x,]] <- rep(x=1/knn,times = knn)
    }else if(sum(exp(-(dist_list^2)/(2*var(dist_list)))) == 0){
      temp[nearest$nn.cell[x,]] <- rep(x=1/knn,times = knn)
    }else{
      temp[nearest$nn.cell[x,]] <- exp(-(dist_list^2)/(2*var(dist_list)))/sum(exp(-(dist_list^2)/(2*var(dist_list))))
    }
    return(temp)
  }))
  
  colnames(nearest_matrix) <- colnames(query)
  rownames(nearest_matrix) <- colnames(query)
  nearest_matrix <- as.matrix(nearest_matrix)
  
  #predict
  predict_matrix <- anchor_matrix
  if(Iteration){
    for (i in 1:Iteration) {
      predict_matrix <- nearest_matrix %*% as.matrix(predict_matrix[colnames(query),])
    }
  }
  rownames(predict_matrix) <- colnames(query)
  colnames(predict_matrix) <- colnames(anchor_matrix)
  gc()
  
  if(return_query){
    label <- base::lapply(colnames(query),FUN = function(x){
      return(which.max(predict_matrix[x,]))
    })
    label <- unlist(label)
    label <- colnames(predict_matrix)[label]
    query$predict_label <- label
    
    colnames(predict_matrix) <- paste(colnames(predict_matrix),'score',sep = '_')
    predict_matrix <- predict_matrix[colnames(query),]
    query@meta.data <- cbind(query@meta.data,predict_matrix)
    return(query)
  }else{
    return(predict_matrix)
  }
}

#' Predict query cell label based on KNN graph between query and reference( data ).
#' @param data The reference seurat object.
#' @param query The query seurat object.
#' @param reference_var The colname of the meta.data table storing the reference label.
#' @param reduction The low dimension graph to be used for calculation.
#' @param knn How many nearest neighbors to be selected in data seurat object?
#' @param return_query Whether to return query seurat object?
my_KNN_label_transfer <- function(data,query,reference_var,reduction='pca',knn=5,return_query=TRUE){
  
  #load package
  require(Seurat)
  require(dplyr)
  
  #create anchor matrix
  group_list <- unique(as.character(data@meta.data[,reference_var]))
  anchor_matrix <- do.call(rbind,base::lapply(colnames(data),FUN = function(x){
    return(as.numeric(as.character(data@meta.data[x,reference_var]) == group_list))
  }))
  rownames(anchor_matrix) <- colnames(data)
  colnames(anchor_matrix) <- group_list
  anchor_matrix <- as.matrix(anchor_matrix)
  print('create anchor matrix done!')
  
  #create knn graph
  nearest <- My_FindKNN(data = data@reductions[[reduction]]@cell.embeddings,
                        query = query@reductions[[reduction]]@cell.embeddings,
                        k = knn)
  
  #modify knn graph
  nearest_matrix <- do.call(rbind,base::lapply(colnames(query),FUN = function(x){
    temp <- rep(0,length(colnames(data)))
    names(temp) <- colnames(data)
    dist_list <- nearest$nn.dists[x,]
    temp[nearest$nn.cell[x,]] <- exp(-dist_list)/sum(exp(-dist_list))
    return(temp)
  }))
  rownames(nearest_matrix) <- colnames(query)
  colnames(nearest_matrix) <- colnames(data)
  nearest_matrix <- as.matrix(nearest_matrix)
  print('create knn graph done!')
  
  #predict
  predict_matrix <- nearest_matrix %*% anchor_matrix
  rownames(predict_matrix) <- colnames(query)
  colnames(predict_matrix) <- group_list
  gc()
  print('predict done!')
  
  #return
  if(return_query){
    label <- base::lapply(colnames(query),FUN = function(x){
      return(which.max(predict_matrix[x,]))
    })
    label <- unlist(label)
    label <- colnames(predict_matrix)[label]
    query$predict_label <- label
    
    colnames(predict_matrix) <- paste(colnames(predict_matrix),'score',sep = '_')
    predict_matrix <- predict_matrix[colnames(query),]
    query@meta.data <- cbind(query@meta.data,predict_matrix)
    return(query)
  }else{
    return(predict_matrix)
  }
}

#' Find label transfer anchor from low dim graph, modified from Seurat::FindTransferAnchors.
#' @param reference reference seurat object.
#' @param query query seurat object.
#' @param ref_reduction name of DimReducObject used for anchor creating in reference.
#' @param query_reduction name of DimReducObject used for anchor creating in query.
#' @param ref_assay name of assay used for anchor creating in reference.
#' @param query_assay name of assay used for anchor creating in query.
#' @param l2.norm whether perform L2 normalization on the cell embedding.
#' @param dims which dimensions to use from the DimReducObject to specify the neighbor search space.
#' @param features features used to create dimreduction object in reference and query, set NULL as default.
#' @param k.anchor how many neighbors (k) to use when finding anchors.
#' @param k.score how many neighbors (k) to use when scoring anchors.
#' @param nn.method method for nearest neighbor finding. Options include: rann, annoy.
#' @param n.trees more trees gives higher precision when using annoy approximate nearest neighbor search.
#' @param eps error bound on the neighbor finding algorithm (from RANN or RcppAnnoy).
#' @param verbose print progress bars and output.
my_FindTransferAnchors <- function(reference,query,ref_reduction = 'pca',query_reduction = 'pca',
                                   ref_assay = 'RNA',query_assay = 'RNA',l2.norm = TRUE,dims = 1:30,features = NULL,
                                   k.anchor = 5,k.score = 30,nn.method = "annoy",n.trees = 50,eps = 0,verbose = TRUE){
  require(Seurat)
  require(dplyr)
  
  #check input
  if((dim(reference@reductions[[ref_reduction]]@cell.embeddings)[2] < max(dims)) | (dim(query@reductions[[query_reduction]]@cell.embeddings)[2] < max(dims))){
    stop("dims exceed!")
  }
  
  #diet query and reference
  reference <- Seurat::DietSeurat(object = reference,counts = TRUE,data = TRUE,scale.data = TRUE,
                                  features = NULL,assays = ref_assay,dimreducs = ref_reduction)
  query <- Seurat::DietSeurat(object = query,counts = TRUE,data = TRUE,scale.data = TRUE,
                              features = NULL,assays = query_assay,dimreducs = query_reduction)
  
  #rename cells
  reference <- Seurat::RenameCells(object = reference,new.names = paste(colnames(reference),'reference',sep = '_'))
  query <- Seurat::RenameCells(object = query,new.names = paste(colnames(query),'query',sep = '_'))
  
  #combine reference and query dimreduction object
  reference_embedding <- reference@reductions[[ref_reduction]]@cell.embeddings[,dims]
  query_embedding <- query@reductions[[query_reduction]]@cell.embeddings[,dims]
  combined_embedding <- CreateDimReducObject(embeddings = as.matrix(rbind(reference_embedding,query_embedding)),assay = ref_assay,key = 'ProjectPC_')
  
  #create combined seurat object
  combined_object <- merge(x = Seurat::DietSeurat(object = reference,counts = TRUE,data = TRUE,scale.data = FALSE,features = NULL,assays = ref_assay,dimreducs = NULL,graphs = NULL),
                           y = Seurat::DietSeurat(object = query,counts = TRUE,data = TRUE,scale.data = FALSE,features = NULL,assays = query_assay,dimreducs = NULL,graphs = NULL))
  combined_object[['pcaproject']] <- combined_embedding
  reduction <- 'pcaproject'
  
  #l2.norm
  if(l2.norm){
    combined_object <- Seurat::L2Dim(object = combined_object,reduction = reduction)
    reduction <- paste0(reduction,".l2")
  }
  
  #prepare anchor finding
  precomputed.neighbors <- list(ref.neighbors = NULL,query.neighbors = NULL)
  reduction.2 <- character()
  
  #find anchor
  anchors <- Seurat:::FindAnchors(object.pair = combined_object,assay = c(ref_assay,query_assay),
                                  cells1 = colnames(x = reference),cells2 = colnames(x = query),reduction = reduction,
                                  reduction.2 = reduction.2,internal.neighbors = precomputed.neighbors,
                                  dims = 1:length(x = dims),k.anchor = k.anchor,k.filter = NA,
                                  k.score = k.score,nn.method = nn.method,n.trees = n.trees,
                                  nn.idx1 = NULL, nn.idx2 = NULL,eps = eps,verbose = verbose)
  
  anchor.set <- new(Class = "TransferAnchorSet",object.list = list(combined_object),
                    reference.cells = colnames(x = reference),query.cells = colnames(x = query),
                    anchors = anchors,anchor.features = features,command = NULL)
  return(anchor.set)
}

#' dimplot modified from Seurat dimplot, only need the embedding matrix and the meta data.
#' @param embedding the embedding matrix used to create the dimplot.
#' @param meta_data additional cell-level metadata to add to the Seurat object.
#' @param group.by name of one or more metadata columns to group (color) cells by.
#' @param split.by name of a metadata column to split plot by.
#' @param label whether to label the clusters.
#' @param repel repel labels.
#' @param ... extra parameters passed to DimPlot.
my_dimplot <- function(embedding,meta_data,group.by = NULL,split.by = NULL,label = TRUE,repel = TRUE,...){
  
  require(Seurat)
  
  #check input
  if(sum(!(rownames(embedding) == rownames(meta_data))) > 0){
    stop('embedding and meta_data not match!')
  }
  
  embedding <- as.matrix(embedding)
  meta_data <- as.data.frame(meta_data)
  
  #create counts
  counts <- matrix(ncol = dim(embedding)[1],nrow = 5,data = 1)
  rownames(counts) <- c('gene1','gene2','gene3','gene4','gene5')
  colnames(counts) <- rownames(embedding)
  
  #create seurat object
  counts <- Seurat::CreateSeuratObject(counts = counts,project = 'temp',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)
  temp <- Seurat::CreateDimReducObject(embeddings = embedding,assay = 'RNA',key = 'UMAP_')
  counts[['UMAP']] <- temp
  
  #use dimplot
  p <- DimPlot(object = counts,group.by = group.by,split.by = split.by,label = label,repel = repel,...)
  return(p)
}

#' do wilcoxon test on seurat object.
#' @param seu.obj the seurat object on which to do the wilcox test.
#' @param assay which assay to use in the seu.obj? Must contain the slot 'data'.
#' @param ident.1 the first cell group set, used as numerator in log2FC.
#' @param ident.2 the second cell group set, used as denominator in log2FC.
#' @param group.by by which meta data to re-group the cells?
#' @param min.express above which should we consider the gene is expressed in the cell?
#' @param log2fc_thresholf limit testing to genes which show, on average, at least X-fold difference (log2-scale) between the two groups of cells.
#' @param pct.1 only test genes that are detected in a minimum fraction of pct.1 cells in ident.1 group.
#' @param features genes to be tested, default is to use all genes.
#' @param only.pos only care about positive markers.
#' @param workers how many cores to be used for computation? Default use all cores.
#' @param future.globals.maxSize the max size of objects used in paralleled computation, default is 2GB.
my_seurat_marker_wilcox_test <- function(seu.obj,
                                         assay = 'RNA',
                                         ident.1,
                                         ident.2,
                                         group.by,
                                         min.express = 0,
                                         log2fc_thresholf = 0.25,
                                         pct.1 = 0.1,
                                         features = NULL,
                                         only.pos = TRUE,
                                         workers = NULL,
                                         future.globals.maxSize = 2*(1024^3)){
  
  #library
  require(Seurat)
  
  #check input
  if(class(seu.obj) != 'Seurat'){
    stop('seu.obj must be a seurat object!')
  }
  if(!(assay %in% names(seu.obj@assays))){
    stop('assay not found!')
  }
  if(!(group.by %in% colnames(seu.obj@meta.data))){
    stop('group not found in meta data!')
  }
  
  #modify input
  group_list <- base::unique(as.character(seu.obj@meta.data[,group.by]))
  if(sum(ident.1 %in% group_list) < length(ident.1)){
    stop('ident.1 not found in the group list!')
  }
  if(sum(ident.2 %in% group_list) < length(ident.2)){
    stop('ident.2 not found in the group list!')
  }
  
  express_matrix <- SeuratObject::GetAssayData(object = seu.obj,slot = 'data',assay = assay)
  if(is.null(features)){
    features <- rownames(express_matrix)
  }else{
    if(sum(!(features %in% rownames(express_matrix))) > 0){
      stop('features not found in the expression matrix!')
    }
  }
  
  if(only.pos == TRUE){
    alternative <- 'greater'
  }else if(only.pos == FALSE){
    alternative <- 'two.sided'
    #when finding negative markers, pct.1 is not suitable for filtering.
    pct.1 <- -1
  }else{
    stop('wrong input of the parameter only.pos!')
  }
  
  #get cell list
  cells.1 <- colnames(seu.obj)[seu.obj@meta.data[,group.by] %in% ident.1]
  cells.2 <- colnames(seu.obj)[seu.obj@meta.data[,group.by] %in% ident.2]
  mat.1 <- express_matrix[features,cells.1,drop = FALSE]
  mat.1 <- expm1(mat.1)
  mat.2 <- express_matrix[features,cells.2,drop = FALSE]
  mat.2 <- expm1(mat.2)
  
  #do wilcox test
  wilcox_out <- my_DF_wilcox_test(mat1 = mat.1,mat2 = mat.2,alternative = alternative,paired = FALSE,workers = workers,future.globals.maxSize = future.globals.maxSize)
  wilcox_out <- wilcox_out[,c('log2FC','pval','fdr','mean1','mean2'),drop = FALSE]
  
  temp <- rowSums(mat.1[rownames(wilcox_out),,drop = FALSE] > min.express)
  wilcox_out$pct.1 <- temp/ncol(mat.1)
  
  temp <- rowSums(mat.2[rownames(wilcox_out),,drop = FALSE] > min.express)
  wilcox_out$pct.2 <- temp/ncol(mat.2)
  
  #return
  wilcox_out <- wilcox_out[abs(wilcox_out$log2FC) > log2fc_thresholf,,drop = FALSE]
  wilcox_out <- wilcox_out[wilcox_out$pct.1 > pct.1,,drop = FALSE]
  if(only.pos == TRUE){
    wilcox_out <- wilcox_out[wilcox_out$log2FC > 0,,drop = FALSE]
  }
  wilcox_out$fdr <- p.adjust(p = wilcox_out$pval,method = 'fdr')
  gc()
  message('Recommend: using p.adjust(p = pval,method = \'fdr\') to modify the fdr if you want to filter more genes by yourself.')
  return(wilcox_out)
}

#' find cell type specific markers by doing wilcox test between the target cell type and each one of the remaining cell types.
#' @param seu.obj the seurat object on which to do the wilcox test.
#' @param assay which assay to use in the seu.obj? Must contain the slot 'data'.
#' @param ident.1 target cell type that you would like to know its specific markers.
#' @param group.by by which meta data to re-group the cells?
#' @param min.express above which should we consider the gene is expressed in the cell?
#' @param log2fc_thresholf limit testing to genes which show, on average, at least X-fold difference (log2-scale) between the two groups of cells.
#' @param pct.1 only test genes that are detected in a minimum fraction of pct.1 cells in ident.1 group.
#' @param features genes to be tested, default is to use all genes.
#' @param overlap_ratio above which ratio when a gene is calculated as marker gene, should we consider it as the specific marker gene? Default overlap_ratio = 1 means the specific marker must be called as marker in all rounds of wilcox tests.
#' @param workers how many cores to be used for computation? Default use all cores.
#' @param future.globals.maxSize the max size of objects used in paralleled computation, default is 2GB.
my_seurat_find_specific_marker <- function(seu.obj,
                                           assay = 'RNA',
                                           ident.1,
                                           group.by,
                                           min.express = 0,
                                           log2fc_thresholf = 0.25,
                                           pct.1 = 0.1,
                                           features = NULL,
                                           overlap_ratio = 1,
                                           workers = NULL,
                                           future.globals.maxSize = 2*(1024^3)){
  
  #library
  require(Seurat)
  
  #check input
  if(class(seu.obj) != 'Seurat'){
    stop('seu.obj must be a seurat object!')
  }
  if(!(assay %in% names(seu.obj@assays))){
    stop('assay not found!')
  }
  if(!(group.by %in% colnames(seu.obj@meta.data))){
    stop('group not found in meta data!')
  }
  if(!is.null(features)){
    if(sum(features %in% rownames(seu.obj@assays[[assay]]@data)) < length(features)){
      stop('features not found in the expression matrix!')
    }
  }
  
  #modify input
  group_list <- base::unique(as.character(seu.obj@meta.data[,group.by]))
  if(sum(ident.1 %in% group_list) < length(ident.1)){
    stop('ident.1 not found in the group list!')
  }
  if(length(group_list) == length(ident.1)){
    stop('no other groups to compare!')
  }
  group_list <- group_list[!(group_list %in% ident.1)]
  
  #calculate marker
  marker_list <- base::lapply(X = group_list,FUN = function(x){
    temp <- suppressMessages(my_seurat_marker_wilcox_test(seu.obj = seu.obj,assay = assay,ident.1 = ident.1,ident.2 = x,group.by = group.by,min.express = min.express,log2fc_thresholf = log2fc_thresholf,pct.1 = pct.1,features = features,only.pos = TRUE,workers = workers,future.globals.maxSize = future.globals.maxSize))
    char <- paste(x,'claculate done!',nrow(temp),'marker detected.',sep = ' ')
    message(char)
    return(rownames(temp))
  })
  gc()
  
  marker_list <- base::unlist(marker_list)
  if(length(marker_list) == 0){
    warning('no marker found!')
    return(NULL)
  }
  unique_marker_list <- base::unique(marker_list)
  
  #filter markers
  marker_list <- table(marker_list)[unique_marker_list]
  marker_list <- unique_marker_list[marker_list >= length(group_list)*overlap_ratio]
  char <- length(unique_marker_list) - length(marker_list)
  char <- paste(as.character(char),'gene filted due to overlap ratio.')
  message(char)
  gc()
  return(marker_list)
}

# liger related -----------------------------------------------------------

#' suggest k and lambda for NMF in liger
#' @param object the liger object
#' @param lambda the range of lambda to be tested
#' @param k the range of k to be tested
#' @param parl_num the number of cores to be used
#' @param library_path your library path for .libPaths()
#' @param pic_path the path to store the picture results
#' some default parameters, usually no need to change
MySuggestLiger <- function(object,lambda,k,parl_num,library_path,pic_path,
                           knn_k = 20,quantiles = 50,min_cells = 20,do.center = FALSE,
                           max_sample = 1000,refine.knn = TRUE,eps = 0.9,resolution = 0.4,min_dist = 0.1){
  require(rliger)
  require(parallel)
  
  test_lambda <- c()
  test_k <- c()
  for(i in lambda){
    for(j in k){
      test_lambda <- append(test_lambda,i)
      test_k <- append(test_k,j)
    }
  }
  
  cl <- makeCluster(parl_num)
  clusterExport(cl,c('library_path','pic_path','test_lambda','test_k','object','knn_k','quantiles','min_cells','do.center','max_sample','refine.knn','eps','resolution','min_dist'),envir = environment())
  clusterEvalQ(cl,.libPaths(library_path))
  clusterEvalQ(cl,library(rliger))
  alignment_score <- parLapply(cl = cl,1:length(test_lambda),function(x){
    lambda <- test_lambda[x]
    k <- test_k[x]
    object <- optimizeALS(object,lambda = lambda,k = k)
    #Quantile Normalization and Joint Clustering
    object <- quantile_norm(object,knn_k = knn_k,quantiles = quantiles,min_cells = min_cells,do.center = do.center,max_sample = max_sample,refine.knn = refine.knn,eps = eps)
    object <- louvainCluster(object,resolution = resolution)
    object <- runUMAP(object,distance = 'cosine',min_dist = min_dist)
    all.plots <- plotByDatasetAndCluster(object, axis.labels = c('UMAP 1', 'UMAP 2'), return.plots = T)
    char <- paste('lambda',as.character(lambda),'k',as.character(k),sep = '_')
    char <- paste(pic_path,char,sep = '/')
    char <- paste(char,'pdf',sep = '.')
    pdf(file = char,width = 12,height = 6)
    print(all.plots[[1]] + all.plots[[2]])
    dev.off()
    return(0)
  })
  stopCluster(cl)
  return('done')
}

# utility -----------------------------------------------------------------

#' homologous gene express convert
#' @param express_matrix express matrix need to be a table.
#' @param anno the annotation table.
#' @param filter_anno whether to filter the duplicated rows in col_1 and col_2 of anno,
#' useful when same gene in specie 1 (convert from) can be converted to several different gene names in species 2 (to converte to).
#' @param future.globals.maxSize The max size of objects used in paralleled computation, default is 2GB.
#' @param workers How many cores to be used for computation? Default use all cores.
#' note the annotation need 4 cols, which are gene name and gene id for the species to be converted and the species to convert to
My_Convert_Homology_Gene_ID <- function(express_matrix,anno,filter_anno = TRUE,future.globals.maxSize = 2*(1024^3),workers = NULL){
  require(dplyr)
  require(Matrix)
  require(Seurat)
  require(future.apply)
  
  options(future.globals.maxSize = future.globals.maxSize)
  if(is.null(workers)){
    workers <- future::availableCores()
  }
  anno <- as.data.frame(as.matrix(anno),stringsAsFactors = FALSE)
  gene_pair <- paste(anno[,1],anno[,2],sep = '_')
  duplicated_pair <- unique(gene_pair[duplicated(gene_pair)])
  anno <- anno[!(gene_pair %in% duplicated_pair),]
  
  if(filter_anno){
    #filter duplicate from col_1 and col_2 in anno
    accept_gene <- unique(anno[duplicated(anno[,1]),1])
    accept_gene <- anno[!(anno[,1] %in% accept_gene),1]
    anno <- anno[anno[,1] %in% accept_gene,]
    accept_gene <- unique(anno[duplicated(anno[,2]),2])
    accept_gene <- anno[!(anno[,2] %in% accept_gene),2]
    anno <- anno[anno[,2] %in% accept_gene,]
  }
  
  #create express matrix
  express_matrix <- SeuratObject::as.sparse(express_matrix)
  print('create sparse matrix done!')
  gc()
  
  #filter genes not in anno and express matrix
  gene <- rownames(express_matrix)
  gene <- dplyr::intersect(gene,c(as.character(anno[,1]),as.character(anno[,2])))
  express_matrix <- express_matrix[gene,]
  anno <- anno[(as.character(anno[,1]) %in% gene) | (as.character(anno[,2]) %in% gene),]
  
  #modify the annotation
  anno[is.na(anno[,3]) | (anno[,3] == ''),3] <- as.character(anno[is.na(anno[,3]) | (anno[,3] == ''),4])
  anno <- anno[(!(is.na(anno[,3]))) & (!(anno[,3] == '')),]
  gene_converted <- unique(anno[,3])
  print('modify annotation done!')
  gc()
  
  #convert
  col_name <- colnames(express_matrix)
  plan(multisession,workers = workers)
  express_matrix <- do.call(rbind,future_lapply(gene_converted,function(x){
    raw_gene <- dplyr::union(anno[anno[,3] == x,1],anno[anno[,3] == x,2])
    raw_gene <- raw_gene[(!is.na(raw_gene)) & (!(raw_gene == ''))]
    raw_gene <- dplyr::intersect(rownames(express_matrix),raw_gene)
    if(length(raw_gene) == 0){
      stop('something wrong when converting!')
    } else if(length(raw_gene) == 1){
      return(as.numeric(express_matrix[raw_gene,]))
    } else{
      return(as.numeric(colSums(express_matrix[raw_gene,])))
    }
  }))
  plan(sequential)
  gc()
  rownames(express_matrix) <- gene_converted
  colnames(express_matrix) <- col_name
  print('convert done!')
  
  #output
  express_matrix <- SeuratObject::as.sparse(express_matrix)
  gc()
  return(express_matrix)
}

#' Scale/zscore a vector of values. Compatible with dplyr::group_by() 
#' @param x a vector of numeric values
scale_this <- function(x) {
  stdev <- sd(x, na.rm=TRUE)
  if(stdev == 0) {
    return(0)
  } else {
    return((x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE))
  }
}

#' Scale rows of a matrix (row z scores)
#' @param mat a numeric matrix
rowScale <- function(mat) {
  
  if (is.null(nrow(mat))) {
    return(scale_this(mat))
  }
  
  rn = rownames(mat)
  cn = colnames(mat)
  
  out = do.call(rbind, lapply(1:nrow(mat), function(x) {
    s = scale_this(mat[x, ])
  }))
  
  rownames(out) = rn
  colnames(out) = cn
  
  return(out)
}

#' Trim the extreme quantiles of a numeric matrix or vector
#' Makes large heatmap plots more interpretable, since outlier values rescale color values
#' for the whole heatmap
#' @param mat a numeric matrix or vector to trim
#' @cuts a length-2 vector of quantiles, low and high, to trim from the object 
#' (they will be replaced by min and max respectively)
trimQuantiles <- function(x, cuts = c(0.1,0.99)) {
  
  if((
    any(is.null(x)) || 
    any(is.na(x)) ||
    any(is.infinite(x))
  )) {
    warning("Object contains non-numeric values.")
  }
  
  if(sum(cuts >= 1) == 2 & sum(cuts <= 100) ==2) {
    cuts <- cuts / 100
  }
  
  quantile.cuts <- quantile(x, cuts, na.rm = T)
  x[x < quantile.cuts[1]] <- quantile.cuts[1]
  x[x > quantile.cuts[2]] <- quantile.cuts[2]
  
  return(x)
}

#' Return the proportion of different cells.
#' @param meta_data the meta data table, each row is a different cell, each col is a different feature.
#' @param group.by group by which col in meta data table?
#' @param split.by split by which col in meta data table?
My_Cell_Proportion <- function(meta_data,group.by,split.by){
  
  require(dplyr)
  
  #check input
  if(!((group.by %in% colnames(meta_data)) & (split.by %in% colnames(meta_data)))){
    stop('missing colnames in meta_data!')
  }
  
  #group list and split list
  meta_data[,group.by] <- as.factor(as.character(meta_data[,group.by]))
  meta_data[,split.by] <- as.factor(as.character(meta_data[,split.by]))
  group_list <- levels(meta_data[,group.by])
  split_list <- levels(meta_data[,split.by])
  full_list <- list(group=rep(group_list,times = length(split_list)),split=rep(split_list,each = length(group_list)))
  
  #calculate
  proportion_matrix <- do.call(rbind,base::lapply(1:length(full_list$group),FUN = function(x){
    numerator <- sum((as.character(meta_data[,group.by]) == full_list$group[x]) & (as.character(meta_data[,split.by]) == full_list$split[x]))
    denominator <- sum((as.character(meta_data[,split.by]) == full_list$split[x]))
    if(denominator == 0){
      stop('denominator can not be zero, something wrong!')
    }
    return(c(numerator/denominator,full_list$group[x],full_list$split[x]))
  }))
  colnames(proportion_matrix) <- c('Proportion',group.by,split.by)
  proportion_matrix <- as.data.frame(proportion_matrix,stringsAsFactors = FALSE)
  proportion_matrix$Proportion <- as.numeric(proportion_matrix$Proportion)
  gc()
  return(proportion_matrix)
}

#' create a KNN object describe the distance between query and data
#' the same parameter as RANN::nn2
My_FindKNN <- function (data,query,k=min(10,nrow(data)),treetype = c("kd","bd"),
                        searchtype = c("standard", "priority", "radius"),radius = 0,eps = 0) {
  require(RANN)
  nearest <- RANN::nn2(data = data,query = query,k = k,treetype = treetype,searchtype = searchtype,radius = radius,eps = eps)
  nearest$nn.cell <- base::apply(nearest$nn.idx,2,function(idx){return(rownames(data)[idx])})
  rownames(nearest$nn.idx) <- rownames(query)
  rownames(nearest$nn.dists) <- rownames(query)
  rownames(nearest$nn.cell) <- rownames(query)
  return(nearest)
}

#' A function to sample nearest neighbors of a particular cell
#' 
#' Takes a cell ID or index, and a KNN object (output from My_FindKNN), 
#' and returns a specified number of cells that are nearest neighbors of that 
#' cell. Cell selection can be biased according to a function of the cell-cell
#' distances (FUN). If no cell sampled, return 'unknown'.
#' @param cell A cell ID or idx
#' @param knn.object An S4 object with names nn.idx, nn.cells, and nn.dists, providing NN info
#' @param N The number of cells to sample
#' @param FUN A quoted string, which evaluates to an expression for weighting the probability of randomly selecting cells. Allows you to weight by distance for example. Make sure FUN(Inf) = 0.
#' @param output String, either 'name' or not, which specifies whether cell IDs or numerical cell index to be used
#' @param replace same as replace parameter in sample function(allow duplicate?)
sampleNearestNeighbors <- function(cell, 
                                   knn.object, 
                                   N=50, 
                                   FUN="(1/exp(dists))*(1/sum(1/exp(dists)))",
                                   output='name',replace=FALSE) {
  
  
  if (output == 'name') {
    cells <- knn.object[['nn.cell']][cell, ]
  }
  else {
    cells <- knn.object[['nn.idx']][cell, ]
  }
  
  dists <- knn.object[['nn.dists']][cell, ]
  
  if(sum(dists == Inf) == length(dists)){
    res <- rep(x = 'unknown',N)
  } else{
    res <- sample(
      x = cells, 
      size = N, 
      replace = replace, 
      prob = eval(parse(text = FUN))
    )
  }
  
  return(res)
  
}

#' Grab genes from a particular gene ontology
#' Wrapped around 'clusterProfiler'
#' 
#' @param x GO path ID
#' @param OrgDb e.g. 'org.Hs.eg.db' database
#' @param ont ontology to use
#' @param keytype symbol? Entrez?
getGOgeneSet <- function(x, 
                         OrgDb = "org.Hs.eg.db", 
                         ont = "BP", 
                         keytype = "SYMBOL") {
  require(clusterProfiler)
  goObj <- clusterProfiler:::get_GO_data(OrgDb=OrgDb, ont=ont, keytype=keytype)
  mapObj <- goObj$PATHID2EXTID
  
  if (!is.element(x, names(mapObj))){
    tt <- names(goObj$PATHID2NAME)
    names(tt) <- goObj$PATHID2NAME
    x <- tt[x]
  }
  
  return(sort(unique(mapObj[[x]])))
}

#' Given a list of GO terms, filter the terms with ancestor terms in the list.
#' @param GO_id A list (vector) of GO terms id.
#' @param GO_ontology GO ontology, must be BP CC or MF.
filter_child_GO_term <- function(GO_id,GO_ontology = 'BP'){
  
  #load required packages
  require(GO.db)
  
  #check parameter
  if(!(GO_ontology %in% c('BP','CC','MF'))){
    stop('GO_ontology must be BP CC or MF!')
  }
  
  #get GO terms that have no ancestors in the GO_id list
  GO_list <- base::lapply(X = GO_id,FUN = function(x){
    GO_db <- paste0('GO',GO_ontology,'ANCESTOR')
    GO_db <- base::get(GO_db)
    GO_db <- GO_db[[x]]
    if(sum(GO_id %in% GO_db) > 0){
      return(NULL)
    }else{
      return(x)
    }
  })
  GO_list <- base::unlist(GO_list)
  
  #return
  return(GO_list)
}

#' return confusion matrix based on original label and predicted label.
#' @param ori the original label vector.
#' @param prd the predicted label vector.
my_confusion_matrix <- function(ori, prd){
  require(dplyr)
  cross.validation.filt <- tibble(ori = ori, prd = prd) %>% dplyr::count(ori, prd) %>% tidyr::spread(key = prd, value = n)
  cross.validation.filt[is.na(cross.validation.filt)] <- 0
  cross.validation.filt <- cross.validation.filt %>% tidyr::gather(key = "prd", value = "value", -ori) %>% filter(value > 0)
  return(cross.validation.filt)
}

#' do unpaired wilcoxon test on sparse matrix using presto, modified from ArchR:::.sparseMatWilcoxon.
#' @param mat1 the first feature-by-sample matrix, used as numerator in log2FC.
#' @param mat2 the second feature-by-sample matrix, used as denominator in log2FC.
my_sparseMatWilcoxon <- function(mat1,mat2){
  
  #require package
  require(presto)
  require(Matrix)
  
  #check input
  n1 <- ncol(mat1)
  n2 <- ncol(mat2)
  if(n1 != n2){
    warning('unequal col number in mat1 and mat2!')
  }
  if(n1 == 0 | n2 == 0){
    stop('col number can not be zero!')
  }
  if(nrow(mat1) != nrow(mat2)){
    stop('mat1 and mat2 must have identical features!')
  }
  mat1 <- as(mat1,'sparseMatrix')
  mat2 <- as(mat2,'sparseMatrix')
  
  #filter features
  m1 <- Matrix::rowSums(mat1,na.rm = TRUE)
  m2 <- Matrix::rowSums(mat2,na.rm = TRUE)
  temp <- m1+m2
  temp <- c(temp > 0)
  char <- as.character(length(m1) - sum(temp))
  char <- paste(char,'gene filted due to zero value',sep = ' ')
  message(char)
  
  mat1 <- mat1[temp,,drop = FALSE]
  mat2 <- mat2[temp,,drop = FALSE]
  
  #df analysis using presto
  df <- presto::wilcoxauc(X = cbind(mat1,mat2),y = c(rep('Top',times = n1),rep('Bot',times = n2)))
  df <- df[which(df$group == "Top"),,drop = FALSE]
  
  m1 <- Matrix::rowSums(mat1,na.rm = TRUE)
  m2 <- Matrix::rowSums(mat2,na.rm = TRUE)
  offset <- 1
  log2FC <- log2((m1 + offset)/(m2 + offset)) + log2(n2/n1)
  
  #return
  out <- data.frame(log2FC = log2FC,
                    fdr = df$padj,
                    pval = df$pval,
                    mean1 = Matrix::rowMeans(mat1,na.rm = TRUE),
                    mean2 = Matrix::rowMeans(mat2,na.rm = TRUE),
                    n1 = ncol(mat1),
                    n2 = ncol(mat2),
                    auc = df$auc)
  return(out)
}

#' do wilcoxon test on matrix using stats::wilcox.test.
#' @param mat1 the first feature-by-sample matrix, used as numerator in log2FC.
#' @param mat2 the second feature-by-sample matrix, used as denominator in log2FC.
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.
#' @param paired a logical indicating whether you want a paired test.
#' @param workers how many cores to be used for computation? Default use all cores.
#' @param future.globals.maxSize the max size of objects used in paralleled computation, default is 2GB.
my_DF_wilcox_test <- function(mat1,
                              mat2,
                              alternative = 'two.sided',
                              paired = FALSE,
                              workers = NULL,
                              future.globals.maxSize = 2*(1024^3),
                              ...){
  
  #require package
  require(Matrix)
  require(future)
  require(future.apply)
  
  #check input
  n1 <- ncol(mat1)
  n2 <- ncol(mat2)
  if(n1 == 0 | n2 == 0){
    stop('col number can not be zero!')
  }
  if(nrow(mat1) != nrow(mat2)){
    stop('mat1 and mat2 must have identical features!')
  }
  if((paired == TRUE) & (n1 != n2)){
    stop('paired wilcoxon test must have same sample numbers!')
  }
  mat1 <- as.matrix(mat1)
  mat2 <- as.matrix(mat2)
  
  #filter features
  m1 <- Matrix::rowSums(x = mat1,na.rm = TRUE)
  m2 <- Matrix::rowSums(x = mat2,na.rm = TRUE)
  temp <- m1+m2
  temp <- c(temp > 0)
  char <- as.character(length(m1) - sum(temp))
  char <- paste(char,'gene filted due to zero value',sep = ' ')
  message(char)
  
  mat1 <- mat1[temp,,drop = FALSE]
  mat2 <- mat2[temp,,drop = FALSE]
  
  #future package parallel
  options(future.globals.maxSize = future.globals.maxSize)
  if(is.null(workers)){
    workers <- future::availableCores()
  }
  
  #df analysis
  plan(multisession,workers = workers)
  p_val <- future.apply::future_lapply(X = 1:nrow(mat1),FUN = function(i){
    p <- stats::wilcox.test(x = as.numeric(mat1[i,]),
                            y = as.numeric(mat2[i,]),
                            alternative = alternative,
                            paired = paired,
                            ...)
    return(p$p.value)
  })
  plan(sequential)
  p_val <- unlist(p_val)
  
  m1 <- Matrix::rowSums(x = mat1,na.rm = TRUE)
  m2 <- Matrix::rowSums(x = mat2,na.rm = TRUE)
  offset <- 1
  log2FC <- log2((m1 + offset)/(m2 + offset)) + log2(n2/n1)
  
  #return
  out <- data.frame(log2FC = log2FC,
                    fdr = p.adjust(p = p_val,method = 'fdr'),
                    pval = p_val,
                    mean1 = Matrix::rowMeans(x = mat1,na.rm = TRUE),
                    mean2 = Matrix::rowMeans(x = mat2,na.rm = TRUE),
                    n1 = ncol(mat1),
                    n2 = ncol(mat2))
  rownames(out) <- rownames(mat1)
  gc()
  return(out)
}

#' do student t test on matrix using stats::t.test.
#' @param mat1 the first feature-by-sample matrix, used as numerator in log2FC.
#' @param mat2 the second feature-by-sample matrix, used as denominator in log2FC.
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.
#' @param paired a logical indicating whether you want a paired test.
#' @param workers how many cores to be used for computation? Default use all cores.
#' @param future.globals.maxSize the max size of objects used in paralleled computation, default is 2GB.
my_DF_t_test <- function(mat1,
                         mat2,
                         alternative = 'two.sided',
                         paired = FALSE,
                         workers = NULL,
                         future.globals.maxSize = 2*(1024^3),
                         ...){
  
  #require package
  require(Matrix)
  require(future)
  require(future.apply)
  
  #check input
  n1 <- ncol(mat1)
  n2 <- ncol(mat2)
  if(n1 == 0 | n2 == 0){
    stop('col number can not be zero!')
  }
  if(nrow(mat1) != nrow(mat2)){
    stop('mat1 and mat2 must have identical features!')
  }
  if((paired == TRUE) & (n1 != n2)){
    stop('paired t test must have same sample numbers!')
  }
  mat1 <- as.matrix(mat1)
  mat2 <- as.matrix(mat2)
  
  #filter features
  m1 <- Matrix::rowSums(x = mat1,na.rm = TRUE)
  m2 <- Matrix::rowSums(x = mat2,na.rm = TRUE)
  temp <- m1+m2
  temp <- c(temp > 0)
  char <- as.character(length(m1) - sum(temp))
  char <- paste(char,'gene filted due to zero value',sep = ' ')
  message(char)
  
  mat1 <- mat1[temp,,drop = FALSE]
  mat2 <- mat2[temp,,drop = FALSE]
  
  #future package parallel
  options(future.globals.maxSize = future.globals.maxSize)
  if(is.null(workers)){
    workers <- future::availableCores()
  }
  
  #df analysis
  plan(multisession,workers = workers)
  p_val <- future.apply::future_lapply(X = 1:nrow(mat1),FUN = function(i){
    p <- stats::t.test(x = as.numeric(mat1[i,]),
                       y = as.numeric(mat2[i,]),
                       alternative = alternative,
                       paired = paired,
                       ...)
    return(p$p.value)
  })
  plan(sequential)
  p_val <- unlist(p_val)
  
  m1 <- Matrix::rowSums(x = mat1,na.rm = TRUE)
  m2 <- Matrix::rowSums(x = mat2,na.rm = TRUE)
  offset <- 1
  log2FC <- log2((m1 + offset)/(m2 + offset)) + log2(n2/n1)
  
  #return
  out <- data.frame(log2FC = log2FC,
                    fdr = p.adjust(p = p_val,method = 'fdr'),
                    pval = p_val,
                    mean1 = Matrix::rowMeans(x = mat1,na.rm = TRUE),
                    mean2 = Matrix::rowMeans(x = mat2,na.rm = TRUE),
                    n1 = ncol(mat1),
                    n2 = ncol(mat2))
  rownames(out) <- rownames(mat1)
  gc()
  return(out)
}

#' create little axis for ggplot scatter plot.
#' @param x_range the range for x axis, must be a numeric vector like c(1,2).
#' @param y_range the range for y axis.
#' @param ratio ratio for scaling the axis.
#' @param margin_value margin value for new axis and title label.
#' @param x_label x axis title label.
#' @param y_label y axis title label.
little_axis <- function(x_range,y_range,ratio = 0.2,margin_value = 1,
                        x_label = 'UMAP_1',y_label = 'UMAP_2'){
  #load package
  require(ggplot2)
  
  #check input
  if(!(length(x_range == 2) & class(x_range) == 'numeric' & x_range[2] > x_range[1])){
    stop('x_range input error!')
  }
  if(!(length(y_range) == 2 & class(y_range) == 'numeric' & y_range[2] > y_range[1])){
    stop('y_range input error!')
  }
  if(!(class(ratio) == 'numeric')){
    stop('ratio must be numeric!')
  }
  if(!(class(margin_value) == 'numeric')){
    stop('margin_value must be numeric!')
  }
  
  #coordinate
  x_lower <- x_range[1] - margin_value
  y_lower <- y_range[1] - margin_value
  x_upper <- x_range[2]
  y_upper <- y_range[2]
  
  x_line_length <- ratio*(x_upper - x_lower) + x_lower
  y_line_length <- ratio*(y_upper - y_lower) + y_lower
  x_mid_point <- ratio*(x_upper - x_lower)/2 + x_lower
  y_mid_point <- ratio*(y_upper - y_lower)/2 + y_lower
  
  #axis data
  axis_data <- data.frame(x = c(x_lower,x_line_length,x_lower,x_lower),
                          y = c(y_lower,y_lower,y_lower,y_line_length),
                          group = c('x','x','y','y'))
  
  #label data
  label_data <- data.frame(x = c(x_mid_point,x_lower - margin_value),
                           y = c(y_lower - margin_value,y_mid_point),
                           angle = c(0,90),
                           label = c(x_label,y_label))
  
  #line data plot
  line_data <- ggplot2::geom_line(data = axis_data,aes(x = x,y = y,group = group),
                                  arrow = arrow(length = unit(0.1,'inches'),ends = 'last',type = 'closed'))
  
  #label data plot
  label_plot_data <- ggplot2::geom_text(data = label_data,aes(x = x,y = y,angle = angle,label = label),fontface = 'bold')
  
  #return
  return(list(line_data,label_plot_data))
}

#' a simple function to reorder the hclust tree.
#' @param cl.obj a hclust object for reordering.
#' @param order.list the order you want to display the hclust tree, notice that the tree can not intersect with itself, so carefully design your order!
my_order_hclust <- function(cl.obj,order.list){
  
  #load library
  #stats only!
  
  #check input
  if(class(cl.obj) != 'hclust'){
    stop('cl.obj must be a hclust object!')
  }
  if(sum(order.list %in% cl.obj$labels) < length(order.list) | length(order.list) != length(cl.obj$labels)){
    stop('order.list error!')
  }
  
  #get sequence
  seq_num <- base::lapply(X = cl.obj$labels,FUN = function(x){
    return(which(order.list == x) - 1)
  })
  seq_num <- 2^(unlist(seq_num))
  
  #reorder
  cl.dend <- as.dendrogram(cl.obj)
  cl.dend <- reorder(x = cl.dend,wts = seq_num,agglo.FUN = sum)
  
  #return
  return(cl.dend)
}

#' This function randomly splits a list into equal-sized sublists by shuffling elements and evenly distributing them among sublists.
#' @param x A given list or vector.
#' @param n How many sublists to split into (at least 2 sublists).
random_equal_split <- function(x,n){
  
  #check param
  if(n < 2){
    stop('At least split into 2 groups!')
  }
  if(length(x) < n){
    stop('Too much groups to split!')
  }
  
  #generate tags
  tag_list <- c()
  for (i in 1:(n-1)) {
    temp <- base::rep(x = paste0('group',i),times = round(length(x)/n))
    tag_list <- base::append(x = tag_list,values = temp)
  }
  temp <- base::rep(x = paste0('group',n),times = length(x) - (n-1)*round(length(x)/n))
  tag_list <- base::append(x = tag_list,values = temp)
  
  #random permutation
  tag_list <- base::sample(x = tag_list,replace = FALSE)
  
  #grouped list
  grouped_list <- base::lapply(X = paste0('group',1:n),FUN = function(label){
    idx <- which(tag_list == label)
    return(x[idx])
  })
  names(grouped_list) <- paste0('group',1:n)
  
  #return
  return(grouped_list)
}
