#' This function uses the OrderedList package to calculate pairwise similarity scores for all spots and cluster on these scores using pam
#' 
#' @param counts non-negative numeric matrix containing gene counts, 
#' rows correspond to spots, columns correspond to genes
#' @param nc integer > 1, number of clusters the data should be divided in
#' @param n integer > 1, the number of genes OrderedList takes into account at both ends of the ranked lists
#' @param alpha numeric 0 < alpha; weight for genes in OrderedList. The larger alpha, less genes toward the middle of the lists are taken into account
#' @param ncores integer > 0, number of threads to use
#' @param .tsne logical, should tSNE plots be generated to visualize the clustering? default TRUE
#' @details For more details on n and alpha see OrderedList package
#' @return list containing 3 elements:\cr
#' 1) scores - score matrix containing pairwise OrderedList-scores for all spots\cr
#' 2) clustering - numeric clustering vector, contains clusters for all spots in the same order as rownames(counts)\cr
#' 3) tsne.scores - ggplot object, tSNE plot of score matrix coloured by clustering\cr
#' 4) tsne.counts - ggplot object, tSNE plot of count matrix coloured by clustering
#' @export
cluster_counts_OL <- function(counts, nc = 6, n = 1000, alpha = 0.25,  ncores = 4, score.mat = NULL, .tsne = TRUE){
	suppressMessages(library(parallel, quietly = TRUE))
	suppressMessages(library(doParallel, quietly = TRUE))
	suppressMessages(library(foreach, quietly = TRUE))
	# calculate similarity on z-scores
	counts <- t(apply(counts, 2, function(x) {(x-mean(x)) / sd(x)}))

	# needed for OrderedList
	base = exp(-alpha)

	if(is.null(score.mat) | !all(colnames(score.mat) == colnames(counts))){
		# serial version
		if(ncores == 1){
			# create and fill score matrix
			score.mat <- matrix(0, ncol(counts), ncol(counts))
			for(i in 1:ncol(score.mat)){
				for(j in 1:ncol(score.mat)){
					score <- OrderedList::scoreRankings(rank(counts[, i]), rank(counts[, j]), n, base)
					score.mat[i,j] <- score
					score.mat[j,i] <- score
				}
			}
		}else{
			# parallel version
			# initialize threading
			cl <- makeCluster(ncores)
			registerDoParallel(cl)

			# calculate scores between all spot and all other spots in parallel
			score.mat <- foreach(i = 1:ncol(counts), .combine = 'cbind') %dopar% {
				sapply(1:ncol(counts), function(x){
					OrderedList::scoreRankings(rank(counts[, i]), rank(counts[, x]), n, base)
				})
			}
			stopCluster(cl)
		}
		colnames(score.mat) <- colnames(counts)
		rownames(score.mat) <- colnames(counts)
	}
	
	# cluster using pam with nc starting points
	clustering <- cluster::pam(score.mat, nc)$clustering

	# create tSNE plots
	if(.tsne){
    tsne <- Rtsne::Rtsne(score.mat, check_duplicates = FALSE)
  	tsne.counts <- Rtsne::Rtsne(t(counts[order(apply(counts, 1, sd), decreasing = TRUE)[1:min(500, nrow(counts))],]), check_duplicates = F)
  	
  	df <- data.frame(tsne$Y, clustering)
  	colnames(df) <- c("tSNE1", "tSNE2", "cluster")
  	df$cluster <- as.factor(df$cluster)
  	
  	df.counts <- data.frame(tsne.counts$Y, clustering)
  	colnames(df.counts) <- c("tSNE1", "tSNE2", "cluster")
  	df.counts$cluster <- as.factor(df.counts$cluster)
  
  	p1 <- ggplot(df, aes(x=tSNE1, y=tSNE2, col = cluster)) + geom_point()
  	p2 <- ggplot(df.counts, aes(x = tSNE1, y = tSNE2, col = cluster)) + geom_point()
	}else{
	  p1 <- NULL
	  p2 <- NULL
	}

	return(list(scores = score.mat, clustering = clustering, tsne.scores = p1, tsne.counts = p2))
}

