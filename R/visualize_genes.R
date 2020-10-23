#' This function visualizes the results of the clustering and analysis using spatial plots and heatmaps
#'
#' @param counts non-negative numeric matrix containing gene counts, 
#' rows correspond to spots, columns correspond to genes
#' @param ids data frame or matrix assigning spatial coordinates to the spots
#' @param img image of ST slide, read using EBImage::readImage()
#' @param clustering numeric vector, same order as rownames(counts), assigning each spot a cluster
#' @param genelist data frame containing genes of interest (e.g. output of filter_genes)
#' @param filepath character, directory in which to save plots; if not specified, plots will not be saved to disk
#' @param plot.params plotting parameters as returned by plot_adjustment()
#' @return list containing 4 plots (ggplot objects):\cr
#' 1) spatial.cluster - overlay of clustering on ST image if image was supplied; NULL otherwise\cr
#' 2) spatial.cluster.legend - spatial clustering, without ST image in background, but with legend\cr
#' 3) heatmap.plt - simplified heatmap, depicting mean expression of all interesting genes in the different clusters\cr
#' 4) heatmap.full.plt - heatmap visualizing expression of all interesting genes in all spots\cr
#' @export
visualize_genes <- function(counts, ids, img = NULL, 
                            clustering, genelist, filepath = NULL, 
                            plot.params = list(nx = 35, ny = 33, ox = 0, oy = 0)){
    # remove genes with 0 variance for the heatmaps
    if(any(apply(counts, 2, var) == 0)){
	    counts <- counts[,-which(apply(counts, 2, var) == 0)]
    }

  # create spatial plots of clustering
    spatial.cluster.nobg <- spatial_plot(
      rownames(counts), 
      ids, 
      clustering, 
      NULL, 
      "discrete", 
      plot.params
      )
    
    if(!is.null(img))
        spatial.cluster <- spatial_plot(
          rownames(counts), 
          ids, 
          clustering, 
          img, 
          "discrete", 
          plot.params
          )

    all.genes <- genelist[which(genelist %in% colnames(counts))]

    # create z-scores if counts are raw or lasso data
    # don't if they were normalized with sctransform
    if(!any(counts < 0)){
      heatmap.counts <- apply(counts, 2, function(x){(x-mean(x)) / sd(x)})
      counts <- heatmap.counts
    }else{
	    heatmap.counts <- counts
    }
    
    # create data frame containing average gene expression for each cluster
    heatmap.df <- matrix("", ncol = 3, nrow = length(unique(clustering))*length(all.genes))
    for(i in 1:length(unique(clustering))){
        cl <- unique(clustering)[i]
        heatmap.df[((i-1)*length(all.genes)+1):(i*length(all.genes)),] <- cbind(
          cl, 
          all.genes, 
          colMeans(counts[which(clustering == as.numeric(cl)), all.genes, drop = F])
          )
    }
    rm(counts)

    # order count matrix by clustering for better heatmap
    heatmap.counts <- heatmap.counts[order(clustering),]
    # data frame containing gene expression of all genes in all spots
    heatmap.full.df <- matrix("", 
                              nrow = length(all.genes)*nrow(heatmap.counts), 
                              ncol = 3
                              )
    for(i in 1:length(all.genes)){
        g <- all.genes[i]
         for(j in 1:nrow(heatmap.counts)){
             heatmap.full.df[(i-1)*nrow(heatmap.counts)+j,] <- c(
								 as.character(j), 
								 g, 
								 as.character(heatmap.counts[j, g])
	          )
         }
    }

    # hierarchical clustering and reordering of genes
    gene.cluster <- hclust(dist(t(heatmap.counts[,all.genes]), "euclidean"), "complete")
    gene.order <- gene.cluster$order
    genes <- all.genes[gene.order]

    # hierarchical clustering and reordering of spots
    spot.cluster <- hclust(dist(heatmap.counts[,all.genes], "euclidean"), "complete")
    spot.order <- spot.cluster$order
    spots <- rownames(heatmap.counts)[spot.order]
    
    # prepare data frames for plotting
    heatmap.df <- as.data.frame(heatmap.df)
    colnames(heatmap.df) <- c("cluster", "gene", "expression")
    heatmap.df$expression <- as.numeric(as.character(heatmap.df$expression))
    heatmap.df$cluster <- as.numeric(as.character(heatmap.df$cluster))

    heatmap.full.df <- as.data.frame(heatmap.full.df)
    colnames(heatmap.full.df) <- c("spot", "gene", "expression")
    heatmap.full.df$expression <- as.numeric(as.character(heatmap.full.df$expression))
    heatmap.full.df$spot <- as.numeric(as.character(heatmap.full.df$spot))

    # plot heatmaps
    heatmap.plt <- ggplot(heatmap.df, aes(x=cluster, y=gene, fill = expression, col = expression)) + 
	    		geom_tile(na.rm = TRUE) + 
			theme(axis.text.y = element_blank()) + 
			xlab("cluster") + 
			ylab("gene") +
    			scale_fill_gradient2(low = "blue", mid = "black", high = "yellow") +
		        scale_color_gradient2(low = "blue", mid = "black", high = "yellow") +	
			scale_y_discrete(limits = genes)
    
    heatmap.full.plt <- ggplot(heatmap.full.df, 
			       aes(x=spot, y=gene, fill = expression, col = expression)) + 
			geom_tile(na.rm = TRUE) + 
			theme(axis.text.y = element_blank()) + 
			xlab("spot") + 
			ylab("gene") +
    			scale_fill_gradient2(low = "blue", mid = "black", high = "yellow") + 
			scale_color_gradient2(low = "blue", mid = "black", high = "yellow") +
			scale_y_discrete(limits = genes) +
			scale_x_discrete(limits = spots)

    # create filepath ...
    if(!is.null(filepath) && !dir.exists(filepath)){
        dir.create(filepath, recursive = TRUE)
    }
    # ... and save plots
    if(!is.null(filepath)){
        pdf(paste(filepath, "clustering.pdf", sep = "/"), width = 12, height = 12)
        
        if(exists("spatial.cluster")) plot(spatial.cluster)
        plot(spatial.cluster.nobg)
        plot(heatmap.plt)
        plot(heatmap.full.plt)
        dev.off()
    }
    
    # this variable exists only if background image was passed to function
    if(!exists("spatial.cluster")) spatial.cluster <- NULL
    
    cluster.plots <- list(spatial = spatial.cluster, spatial.no_img = spatial.cluster.nobg, heatmap = heatmap.plt, heatmap.full = heatmap.full.plt)
    return(cluster.plots)
}
