cluster_counts_OL <- function(counts, nc = 6, n = 1000, alpha = 0.25, verbose = FALSE){
	counts <- t(counts)
	score.mat <- matrix(0, ncol(counts), ncol(counts))
	colnames(score.mat) <- colnames(counts)
	rownames(score.mat) <- colnames(counts)

	base = exp(-alpha)
	if(verbose) cat("Create distance matrix\n__________\n\n")
	for(i in 1:ncol(score.mat)){
		for(j in 1:ncol(score.mat)){
			if(i == j) next
			score <- OrderedList::scoreRankings(rank(counts[, i]), rank(counts[, j]), n, base)
			score.mat[i,j] <- score
			score.mat[j,i] <- score
		}
	}

    d <- dist(score.mat, "manhattan")
    tree <- hclust(d, "single")

    score.mat <- score.mat[tree$order, tree$order]

    clustering <- kmeans(score.mat, nc)$cluster

    tsne <- Rtsne(score.mat, check_duplicates = FALSE)
	return(list(scores = score.mat, clustering = clustering, tsne.embed = tsne))
}

