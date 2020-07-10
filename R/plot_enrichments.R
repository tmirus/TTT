#' plot spatial heatmaps of enriched terms to pdf
#'
#' @param enrichments output of test_enrichment() function
#' @param ids barcode data frame assigning spatial coordinates to all spots; column names 'X' and 'Y'
#' @param counts non-negative numeric matrix containing gene expression for all spots (spots x genes)
#' @param clustering numeric vector assigning clusters to all spots. same order as spots in counts
#' @param img ST image (EBImage object) if available. Will be used as plot background. If no image is available, leave default NULL
#' @param plot.params list as created by plot_adjustment()
#' @param filepath character, directory in which the output file should be stored. If not specified, currecnt working directory is used.
#' @return plot to file "enrichment_spatial.pdf" in 'filepath'. Returns the data frame containing the enrichment p-values for each term for all clusters.
#' @export

plot_enrichments <- function(enrichments, ids, counts, clustering, img = NULL, plot.params = list(nx = 35, ny = 33, ox = 0, oy = 0), filepath = NULL){

    if(is.null(filepath)){
        filepath <- "./enrichment_spatial.pdf"
    }else{
	if(!dir.exists(filepath)) dir.create(filepath)
        filepath <- paste(filepath, "enrichment_spatial.pdf", sep = "/")
    }
    split_up_down <- TRUE

    enrichment.terms <- sapply(enrichments, function(x){
            terms <- c()
            if(!is.null(x$up)){
                terms <- c(terms, x$up$GO.ID)
            }
            if(!is.null(x$down)){
                terms <- c(terms, x$down$GO.ID)
            }
	    if(length(terms) == 0){
		    if(!is.null(x$GO.ID)){
			    terms <- c(terms, x$GO.ID)
			    split_up_down <- FALSE
		    }
	    }
            return(terms)
        })
    enrichment.terms <- unique(unlist(enrichment.terms))

    cat(length(enrichment.terms), "enriched terms total\n")

    # store p-value and up/down for each term in each cluster
    enrichment.mat <- matrix(1, nrow = length(enrichment.terms), ncol = length(unique(clustering)))
    rownames(enrichment.mat) <- enrichment.terms
    colnames(enrichment.mat) <- unique(clustering)

    enrichment.chars <- c()
    for(t in enrichment.terms){
        for(cl in unique(clustering)){
		if(split_up_down){
            if(t %in% enrichments[[cl]]$up$GO.ID){
                enrichment.mat[t, cl] <- -log10(as.numeric(enrichments[[cl]]$up[which(enrichments[[cl]]$up$GO.ID == t), "weightFisher"]))
                if(!enrichments[[cl]]$up[which(enrichments[[cl]]$up$GO.ID == t), "Term"] %in% enrichment.chars)
                    enrichment.chars <- c(enrichment.chars, enrichments[[cl]]$up[which(enrichments[[cl]]$up$GO.ID == t), "Term"])
            }else{
                if(t %in% enrichments[[cl]]$down$GO.ID){
                    enrichment.mat[t, cl] <- log10(as.numeric(enrichments[[cl]]$down[which(enrichments[[cl]]$down$GO.ID == t), "weightFisher"]))
                    if(!enrichments[[cl]]$down[which(enrichments[[cl]]$down$GO.ID == t), "Term"] %in% enrichment.chars)
                        enrichment.chars <- c(enrichment.chars, enrichments[[cl]]$down[which(enrichments[[cl]]$down$GO.ID == t), "Term"])
                }
            }
		}else{
			if(t %in% enrichments[[cl]]$GO.ID){
				enrichment.mat[t, cl] <- -log10(as.numeric(enrichments[[cl]][which(enrichments[[cl]]$GO.ID == t),"weightFisher"]))
				if(!enrichments[[cl]][which(enrichments[[cl]]$GO.ID == t), "Term"] %in% enrichment.chars)
					enrichments.chars <- c(enrichment.chars, enrichments[[cl]][which(enrichments[[cl]]$GO.ID == t), "Term"])
			}
		}
        }
    }
    names(enrichment.terms) <- enrichment.chars

    enrichment.mat[enrichment.mat == 1] <- NA

    pdf(filepath)
    for(i in 1:length(enrichment.terms)){
        t <- enrichment.terms[i]
        n <- names(enrichment.terms)[i]
        cluster.values <- rep(0, nrow(counts))
        for(cl in unique(clustering)){
            cluster.values[which(clustering == cl)] <- enrichment.mat[t, cl]
        }
        p <- spatial_plot(rownames(counts), ids, cluster.values, img, "discrete", plot.params) + 
            scale_colour_brewer(palette = "YlGnBu") + theme(plot.title = element_text(), legend.position = "right") + 
            ggtitle(n) + labs(col = "+/- log10(pVal)") 
        plot(p)
    }
    dev.off()
    return(enrichment.mat)
}
