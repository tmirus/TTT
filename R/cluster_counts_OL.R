#' @export
cluster_counts_OL <- function(counts, nc = 6, n = 1000, alpha = 0.25, verbose = TRUE, ncores = 4){
	counts <- t(apply(counts, 2, function(x) {(x-mean(x)) / sd(x)}))
	#counts <- t(counts)

	base = exp(-alpha)
	if(verbose) cat("Create score matrix\n_____________\n\n")
	if(ncores == 1){
		score.mat <- matrix(0, ncol(counts), ncol(counts))
		for(i in 1:ncol(score.mat)){
			for(j in 1:ncol(score.mat)){
				#if(i == j) next
				score <- OrderedList::scoreRankings(rank(counts[, i]), rank(counts[, j]), n, base)
				score.mat[i,j] <- score
				score.mat[j,i] <- score
			}
		}
	}else{
		cl <- makeCluster(ncores)
		registerDoParallel(cl)
		score.mat <- foreach(i = 1:ncol(counts), .combine = 'cbind') %dopar% {
			sapply(1:ncol(counts), function(x){
					OrderedList::scoreRankings(rank(counts[, i]), rank(counts[, x]), n, base)
			})
		}
		stopCluster(cl)
	}

	colnames(score.mat) <- colnames(counts)
	rownames(score.mat) <- colnames(counts)
	clustering <- pam(score.mat, nc)$clustering

    tsne <- Rtsne(score.mat, check_duplicates = FALSE)
	df <- data.frame(tsne$Y, clustering)
	colnames(df) <- c("tSNE1", "tSNE2", "cluster")
	df$cluster <- as.factor(df$cluster)
	p <- ggplot(df, aes(x=tSNE1, y=tSNE2, col = cluster)) + geom_point()
	return(list(scores = score.mat, clustering = clustering, tsne = p))
}

