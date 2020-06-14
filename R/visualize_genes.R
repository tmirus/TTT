#' @export
visualize_genes <- function(counts, ids, img, clustering, gene.info, filepath = NULL, nx = 35, ny = 33, ox = -1000/70, oy = 1000/32){
    spatial.cluster.legend <- spatial_plot(rownames(counts), ids, clustering, NULL, "discrete", nx, ny, ox, oy)
    spatial.cluster <- spatial_plot(rownames(counts), ids, clustering, img, "discrete", nx, ny, ox, oy)

    print(str(gene.info[[2]]))
    all.genes <- unique(unlist(sapply(gene.info[[2]], function(x){names(x)}), use.names = F))
    print(all.genes)
    print(str(all.genes))
    heatmap.mat <- matrix(0, nrow = length(all.genes), ncol = length(unique(clustering)))
    rownames(heatmap.mat) <- all.genes
    colnames(heatmap.mat) <- unique(clustering)

    print("reduced heatmap df")
    heatmap.counts <- apply(counts, 2, function(x){(x-mean(x)) / sd(x)})
    heatmap.df <- c()
    heatmap.full.df <- c()

    for(cl in colnames(heatmap.mat)){
        heatmap.df <- rbind(heatmap.df,
                            cbind(cl, all.genes, log2(colMeans(counts[which(clustering == as.numeric(cl)), all.genes, drop = F])+1)))
    }

    print("full heatmap df")
    print(Sys.time())
    print(dim(heatmap.counts))
    heatmap.counts <- heatmap.counts[order(clustering),]
    for(g in all.genes){
         for(j in 1:nrow(heatmap.counts)){
             heatmap.full.df <- rbind(heatmap.full.df,
                                 c(j, g, heatmap.counts[j, g]))
         }
    }
    print(Sys.time())
    
    print("prepare reduced heatmap df")
    heatmap.df <- as.data.frame(heatmap.df)
    colnames(heatmap.df) <- c("cluster", "gene", "expression")
    heatmap.df$expression <- as.numeric(as.character(heatmap.df$expression))

    print("prepare full heatmap df")
    heatmap.full.df <- as.data.frame(heatmap.full.df)
    colnames(heatmap.full.df) <- c("spot", "gene", "expression")
    heatmap.full.df$expression <- as.numeric(as.character(heatmap.full.df$expression))

    print("plot heatmaps")
    heatmap.plt <- ggplot(heatmap.df, aes(x=cluster, y=gene, fill = expression)) + geom_tile() + theme(axis.text.y = element_blank()) + xlab("cluster") + ylab("gene") +
    scale_fill_gradient(low = "blue", high = "yellow")
    heatmap.full.plt <- ggplot(heatmap.full.df, aes(x=spot, y=gene, fill = expression)) + geom_tile() + theme(axis.text.y = element_blank()) + xlab("spot") + ylab("gene") +
    scale_fill_gradient(low = "blue", high = "yellow")

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

    cluster.plots <- list(spatial = spatial.cluster, spatial.no_img = spatial.cluster.legend, heatmap = heatmap.plt)

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