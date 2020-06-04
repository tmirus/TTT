#' @export
visualize_genes <- function(counts, ids, img, clustering, gene.info, filepath,nx = 35, ny = 33, ox = -1000/70, oy = 1000/32){
    spatial.cluster.legend <- spatial_plot(rownames(counts), ids, clustering, NULL, "discrete", nx, ny, ox, oy)
    spatial.cluster <- spatial_plot(rownames(counts), ids, clustering, img, "discrete", nx, ny, ox, oy)

    all.genes <- unique(unlist(sapply(gene.info[[2]], function(x){names(x)}), use.names = F))
    print(length(all.genes))

    heatmap.mat <- matrix(0, nrow = length(all.genes), ncol = length(unique(clustering)))
    rownames(heatmap.mat) <- all.genes
    colnames(heatmap.mat) <- unique(clustering)

    for(cl in colnames(heatmap.mat)){
        heatmap.mat[, cl] <- colMeans(counts[which(clustering == as.numeric(cl)), all.genes, drop = F])
    }

    if(!dir.exists(filepath)){
        dir.create(filepath, recursive = TRUE)
    }

    pdf(paste(filepath, "clustering.pdf", sep = "/"), width = 10, height = 10)
    plot(spatial.cluster)
    plot(spatial.cluster.legend)
    heatmap(heatmap.mat)
    dev.off()

    pdf(paste(filepath, "spatial_genes.pdf", sep = "/"), width = 10, height = 10)
    for(g in all.genes){
        col.vec <- rep(0, nrow(counts))
        for(cl in unique(clustering)){
            col.vec[which(clustering == cl)] <- -log10(gene.info[[2]][[as.character(cl)]][g])
        }
        p <- spatial_plot(rownames(counts), ids, counts[,g,drop=F], img, mode = "continuous", nx, ny, ox, oy) + geom_point(aes(col = col.vec))
        plot(p)
    }
    dev.off()
}