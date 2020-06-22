#' @export
visualize_genes <- function(counts, ids, img, clustering, gene.info, filepath = NULL, nx = 35, ny = 33, ox = -1000/70, oy = 1000/32){
    spatial.cluster.legend <- spatial_plot(rownames(counts), ids, clustering, NULL, "discrete", nx, ny, ox, oy)
    spatial.cluster <- spatial_plot(rownames(counts), ids, clustering, img, "discrete", nx, ny, ox, oy)

    #all.genes <- unique(unlist(sapply(gene.info[[2]], function(x){names(x)}), use.names = F))
    all.genes <- rownames(gene.info[[2]])
    heatmap.mat <- matrix(0, nrow = length(all.genes), ncol = length(unique(clustering)))
    rownames(heatmap.mat) <- all.genes
    colnames(heatmap.mat) <- unique(clustering)

    heatmap.counts <- apply(counts, 2, function(x){(x-mean(x)) / sd(x)})
    counts <- heatmap.counts
    
    heatmap.df <- matrix("", ncol = 3, nrow = ncol(heatmap.mat)*length(all.genes))
    for(i in 1:ncol(heatmap.mat)){
        cl <- colnames(heatmap.mat)[i]
        heatmap.df[((i-1)*length(all.genes)+1):(i*length(all.genes)),] <- cbind(cl, all.genes, colMeans(counts[which(clustering == as.numeric(cl)), all.genes, drop = F]))
    }

    heatmap.counts <- heatmap.counts[order(clustering),]
    heatmap.full.df <- matrix("", nrow = length(all.genes)*nrow(heatmap.counts), ncol = 3)
    for(i in 1:length(all.genes)){
        g <- all.genes[i]
         for(j in 1:nrow(heatmap.counts)){
             heatmap.full.df[(i-1)*nrow(heatmap.counts)+j,] <- c(as.character(j), g, as.character(heatmap.counts[j, g]))
         }
    }

    gene.cluster <- hclust(dist(t(heatmap.counts[,all.genes]), "euclidean"), "complete")
    gene.order <- gene.cluster$order
    genes <- all.genes[gene.order]
    
    heatmap.df <- as.data.frame(heatmap.df)
    colnames(heatmap.df) <- c("cluster", "gene", "expression")
    heatmap.df$expression <- as.numeric(as.character(heatmap.df$expression))
    heatmap.df$cluster <- as.numeric(as.character(heatmap.df$cluster))

    heatmap.full.df <- as.data.frame(heatmap.full.df)
    colnames(heatmap.full.df) <- c("spot", "gene", "expression")
    heatmap.full.df$expression <- as.numeric(as.character(heatmap.full.df$expression))
    heatmap.full.df$spot <- as.numeric(as.character(heatmap.full.df$spot))

    heatmap.plt <- ggplot(heatmap.df, aes(x=cluster, y=gene, fill = expression)) + geom_tile() + theme(axis.text.y = element_blank()) + xlab("cluster") + ylab("gene") +
    scale_fill_gradient(low = "blue", high = "yellow") + scale_y_discrete(limits = genes)
    heatmap.full.plt <- ggplot(heatmap.full.df, aes(x=spot, y=gene, fill = expression)) + geom_tile() + theme(axis.text.y = element_blank()) + xlab("spot") + ylab("gene") +
    scale_fill_gradient(low = "blue", high = "yellow") + scale_y_discrete(limits = genes)

    if(!is.null(filepath) && !dir.exists(filepath)){
        dir.create(filepath, recursive = TRUE)
    }

    if(!is.null(filepath)){
        pdf(paste(filepath, "clustering.pdf", sep = "/"), width = 10, height = 10)
        plot(spatial.cluster)
        plot(spatial.cluster.legend)
        plot(heatmap.plt)
        plot(heatmap.full.plt)
        dev.off()
    }

    cluster.plots <- list(spatial = spatial.cluster, spatial.no_img = spatial.cluster.legend, heatmap = heatmap.plt, heatmap.full = heatmap.full.plt)

    # if(!is.null(filepath)){
    #     pdf(paste(filepath, "spatial_genes.pdf", sep = "/"), width = 10, height = 10)
    #     for(g in all.genes){
    #         col.vec <- rep(0, nrow(counts))
    #         for(cl in unique(clustering)){
    #             col.vec[which(clustering == cl)] <- -log10(gene.info[[2]][[as.character(cl)]][g])
    #         }
    #         p <- spatial_plot(rownames(counts), ids, counts[,g,drop=F], img, mode = "continuous", nx, ny, ox, oy) + geom_point(aes(col = col.vec))
    #             plot(p)
    #     }
    #     dev.off()
    # }
    return(list(cluster.plots = cluster.plots))
}
