#' create plots for visualizing the enrichment and clustering results
#' 
#' @param counts non-negative numeric matrix. rows correspond to spots, 
#' columns correspond to genes
#' @param ids data frame with rownames corresponding to barcodes / spot names, assigning
#' spatial coordinates to each spot; coordinates are contained in columns X and Y and conform to
#' the same standards as the output of process_input()
#' @param img EBImage image object the information will be plotted on; default NULL
#' @param clustering.list output of cluster_spots
#' @param filepath directory in which results should be stored
#' @param filename character, filename for the clustering and enrichment plots; only filename, the path is given in filepath
#' @param spatial.file.prefix character, filename_prefix for several PDFs containing spatial enrichment plots; will be store in filepath/
#' @param spot.cluster.names character vector assigning each cluster in clustering.list$spot.clustering a name; for better visualization 
#' (only really useful on second run, when clustering has been visualized and tissue can be named properly). Ordering of the names needs to
#' match a numerical ordering of the clusters.
#' @return NULL, plots to files
#' @export

plot_enrichment_clustering <- function(counts, ids, img, clustering.list, filepath, filename, spatial.file.prefix = NULL, spot.cluster.names = NULL, nx = 35, ny = 33, ox = -1000/70, oy = 1000/32){
    theme_transparent <- theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_blank(), # get rid of legend bg
    legend.box.background = element_blank(), # get rid of legend panel bg
    legend.key = element_blank(), # get rid of key legend fill, and of the surrounding
    legend.position = "none",
    axis.line = element_blank(), # adding a black line for x and y axis
    axis.text = element_blank(),
    axis.title = element_blank(),
    #plot.title = element_blank(), 
    axis.ticks = element_blank()
    )

    spot.enrichments <-  clustering.list$spot.enrichments
    enrichment.mat <- clustering.list$enrichment.mat 
    spot.clustering <- clustering.list$spot.clustering 
    term.clustering <- clustering.list$term.clustering 
    term.clusters <- clustering.list$term.clusters

    if(!dir.exists(filepath)){
        flag <- dir.create(filepath, recursive = TRUE)
        if(!flag){
            stop(paste("Could not create directory", filepath))
        }
    }
    filename <- paste(filepath,filename,sep="/")
    if(is.null(spatial.file.prefix)){
        spatial.file.prefix <- "terms_spatial"
    }
    spatial.file.prefix <- paste(filepath, spatial.file.prefix, sep = "/")

    term.cluster.pattern <- matrix(0, ncol = length(unique(spot.clustering)), nrow = length(term.clusters))
     
    rownames(term.cluster.pattern) <- names(term.clusters)
    colnames(term.cluster.pattern) <- 1:ncol(term.cluster.pattern)
    for(term in names(term.clusters)){
        for(clust in colnames(term.cluster.pattern)){
            term.cluster.pattern[term, clust] <- mean(apply(enrichment.mat[term.clusters[[term]], which(spot.clustering == as.numeric(clust)), drop = F], 2, mean))
        }
    }

    plot.cols <- FALSE
    if(length(unique(spot.clustering)) == length(spot.cluster.names) && length(spot.cluster.names) <= 12){
        plot.cols <- TRUE
        clust.cols <- brewer.pal(n = length(spot.cluster.names), name = "Paired")
        col.vec <- sapply(spot.clustering, function(x){clust.cols[x]})
        names(clust.cols) <- spot.cluster.names
        names(col.vec) <- sapply(spot.clustering, function(x){names(clust.cols)[x]})
    }

    term.col.vec <- rep("grey", length(term.clustering))
    last.col <- "grey"
    for(i in unique(term.clustering)){
        if(last.col == "grey"){
            term.col.vec[which(term.clustering == i)] <- "black"
            last.col <- "black"
        }else{
            last.col <- "grey"
        }
    }

    if(!exists("col.vec")){
        col.vec <- colnames(enrichment.mat)
        names(col.vec) <- spot.clustering
        clust.cols <- brewer.pal(n = length(unique(spot.clustering)), name = "Paired")
        names(clust.cols) <- unique(spot.clustering)
    }

    pdf(filename)
    plot(spatial_plot(colnames(enrichment.mat), ids, names(col.vec), "", NULL))
    plot(spatial_plot(colnames(enrichment.mat), ids, names(col.vec), "", img))
    if(!plot.cols){
        heatmap3(enrichment.mat, Colv=NA, Rowv = NA, col = c("blue", "yellow"), scale = "none", RowSideColors = term.col.vec)
    }else{
        heatmap3(enrichment.mat, Colv = NA, Rowv = NA, col = c("blue", "yellow"), scale = "none",
                ColSideColors = col.vec, RowSideColors = term.col.vec, legendfun = function()showLegend(legend = names(clust.cols), col = clust.cols))
    }
    image(t(term.cluster.pattern), axes = F, xlab = "spot cluster", ylab = "term cluster")
    axis(side = 4,
     at = (1:length(term.clusters)-1)/(length(term.clusters) -1),
     labels = names(term.clusters),
     las = 1)
    axis(side = 3,
     at = (1:ncol(term.cluster.pattern) - 1)/(ncol(term.cluster.pattern)-1),
     labels = names(clust.cols),
     las = 1)
    dev.off()


    # this plotting needs to be generalized using spatial_plot!

    img <- EBImage::resize(img, w =1000, h = 1000)
    enrichment.mat <- spot.enrichments$term.matrix
    ids <- ids[colnames(enrichment.mat), ]
    ids[, 2] <- ids[, 2] * 1000/35 - 1000/70
    ids[, 1] <- ids[, 1] * 1000/33 - 1000/66

    for(n in names(term.clusters)){
        pdf(paste(spatial.file.prefix,"_",n,".pdf",sep=""))
        for(t in term.clusters[[n]]){
            temp <- data.frame(x = -ids[, 2], y = ids[, 1], value = -log10(enrichment.mat[t,]))
            temp[temp$value == 0, "value"] <- NA
            gob <- rasterGrob(img)
            p <- ggplot(temp) + 
                ggtitle(Term(t), subtitle = t) + ylim(0, 1000) + xlim(-1000, 0) +
                annotation_custom(gob, -1000, 0, 0, 1000) +
                geom_point(aes(x=x, y=y, size =value, col = value), na.rm = TRUE) +
                theme_transparent + theme(plot.title = element_text(size =24))
            plot(p)
        }
        dev.off()
    }
    
}
