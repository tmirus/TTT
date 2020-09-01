#' this function plot spatial patterns of enriched GO terms to pdf
#' 
#'@param enrichment.mat numeric matrix of enrichment p-values (terms x spots) as returned by test_spot_enrichment or filter_spot_enrichments
#'@param term.names character vector containing full term names for GO ids in enrichment.mat
#'@param ids barcode data frame assigning spots to x and y coordinates (see process_input)
#'@param plot.params list of parameters needed for good spatial visualization as returned by plot_adjustment
#'@param img EBImage image object the information will be plotted on; default NULL
#'@param filename character, name of the pdf; default "spot_enrichments.pdf"
#'@param output.path character specifying output directory; default "./"

plot_spot_enrichments <- function(enrichment.mat, term.names = NULL, ids, plot.params = list(nx = 35, ny = 33, ox = 0, oy = 0), img = NULL, filename = "spot_enrichments.pdf", output.path = "./"){
    if(is.null(term.names)){
        term.names <- rownames(enrichment.mat)
    }
    if(!dir.exists(output.path)){
        flag <- dir.create(output.path)
        if(!flag){
            stop("Could not find or create output directory")
        }
    }

    enrichment.sig <- apply(enrichment.mat, 1, function(x){mean(x[x>0])})
    enrichment.mat <- enrichment.mat[order(enrichment.sig),]
    
    pdf(paste0(output.path, "/", filename), width = 10, height = 10)
    enrichment.mat[enrichment.mat == 0] <- NA
    for(i in 1:nrow(enrichment.mat)){
        plot(
            spatial_plot(
                colnames(enrichment.mat),
                ids, 
                -log10(enrichment.mat[i,]),
                img,
                plot.params = plot.params,
                mode = "continuous"
            ) + 
            ggplot2::ggtitle(term.names[i], subtitle = rownames(enrichment.mat)[i]) +
            geom_point(aes(col = -log10(enrichment.mat[i,])), na.rm = TRUE) +
            theme(legend.position = "right") +
            labs(col = "-log10(pVal)", size = "-log10(pVal)")
        )
    }
    dev.off()
}