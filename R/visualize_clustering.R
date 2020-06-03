#' @export
visualize_clustering <- function(ids, clustering.info, img, nx = 35, ny = 33, ox = -1000/70, oy = 1000/32){
	p1 <- spatial_plot(rownames(clustering.info$scores), ids, clustering.info$clustering, img, mode = "discrete", nx, ny, ox, oy)
	df <- data.frame(clustering.info$tsne.embed$Y, clustering.info$clustering)
	colnames(df) <- c("X", "Y", "cluster")
	p2 <- ggplot(df) + geom_point(aes(x=X, y=Y, col = cluster)) + xlab("tSNE1") + ylab("tSNE2") + ggtitle("t-SNE embeddding coloured by clustering")
	return(list(spatial.plot = p1, tsne.plot = p2))
}