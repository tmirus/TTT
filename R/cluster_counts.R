#' cluster counts matrix based using louvain clustering based on UMAP and PCA
#'
#' @param counts non-negative numeric matrix containing gene counts,
#' rows correspond to spots, columns correspond to genes
#' @param pca.fraction numeric, 0 < pca.fraction < 1. Fraction of variance to be explained by 
#' principal components taken into account. default 0.8
#' @param umap.metric string, metric passed to uwot::umap for umap embedding. default 'cosine'
#' @param ..., parameters passed to uwot::umap function for fine-tuning
#' @return list containing two elements:\cr
#' 1) clustering - numeric clustering vector, contains clusters for all spots in the same order as rownames(counts)\cr
#' 2) umap.plot - ggplot object, UMAP plot of count matrix coloured by clustering
#' @export

cluster_counts <- function(counts, ids = NULL, pca.fraction = 0.8, umap.metric = "cosine", ...) {
	suppressMessages(library(igraph, quietly = TRUE))
	suppressMessages(library(FNN, quietly = TRUE))

	if(!is.numeric(pca.fraction) || pca.fraction <= 0 || pca.fraction >= 1){
		stop("Invalid parameter 'pca.fraction'. Must be numeric between 0 and 1.")
	}
	if(!umap.metric %in% c('cosine', 'euclidean', 'manhattan', 'hamming', 'categorical')){
		stop("Invalid umap metric. Must be one of 'cosine', 'euclidean', 'manhattan', 'hamming', 'categorical'.")
	}
  # calculate PCA and take as many components as are necessary to explain pca.fraction * 100 percent of variance
	counts.pca <- prcomp(counts)
    x <- counts.pca$x[, 1:min(which(cumsum(counts.pca$sdev**2) / sum(counts.pca$sdev**2) > pca.fraction))]
    cat("Using first ", ncol(x), " principal components\n", sep = "")

	# calculate umap embedding and cluster
  embed <- uwot::umap(x, metric = umap.metric, nn_method = "annoy", pca = NULL, ...)
  
  	if(!is.null(ids)){
  		embed.knn <- as.matrix(cbind(embed, ids[rownames(counts),]))
		embed.knn <- apply(embed.knn, 2, function(x){ (x - mean(x)) / sd(x)})
	}else{
		embed.knn <- embed
	}
  
	k <- 25
	knn.norm <- FNN::get.knn(as.matrix(embed), k = 50)
	knn.norm <- data.frame(
			       from = rep(1:nrow(knn.norm$nn.index), k),
			       to = as.vector(knn.norm$nn.index),
			       weight = 1/(1 + as.vector(knn.norm$nn.dist))
	)
	nw.norm <- igraph::graph_from_data_frame(knn.norm, directed = FALSE)
	nw.norm <- igraph::simplify(nw.norm)
	lc.norm <- igraph::cluster_louvain(nw.norm)
	clustering <- as.numeric(membership(lc.norm))

	# create umap plot
  df <- data.frame(embed, clustering)
  colnames(df) <- c("UMAP1", "UMAP2", "cluster")
  df$cluster <- as.factor(df$cluster)
  umap.plot <- ggplot(df, aes(x = UMAP1, y = UMAP2, col = cluster)) +
      geom_point()

	# return clustering and plot
  return(list(clustering = clustering, umap.plot = umap.plot))
}
