#' create spatial plots for the top n.genes genes in a given table and store them as pdf
#' 
#' @param gene.table data frame containing interesting genes, as returned by filter_genes (or entry differential.genes of analyze_clustering output)
#' @param counts non-negative numeric matrix containing gene expression for all spots (spots x genes)
#' @param ids barcode data frame assigning spatial coordinates to all spots; column names 'X' and 'Y'
#' @param img ST image (EBImage object) if available. Will be used as plotting background. If not available, leave default NULL.
#' @param plot.params list created by plot_adjustment; only needed if img is provided
#' @param n.genes integer, number of genes to be plotted (starting from the top of gene.table)
#' @param filepath character specifying the directory in which the output file ('spatial_genes.pdf') should be stored
#' @param colour character specifying the colour of the spots in the gene expression plots; default "black"
#' @return NULL, plots to file 'filepath'/spatial_genes.pdf
#' @export

spatial_gene_plots <- function(gene.table, counts, ids, img = NULL, plot.params = NULL, n.genes = 100, filepath = NULL, colour = "black"){
    if(is.null(filepath)){
        filepath <- "./spatial_genes.pdf"
    }else{
        if(!dir.exists(filepath)){
            dir.create(filepath)
        }
        filepath <- paste(filepath, "spatial_genes.pdf", sep = "/")
    }

    pdf(filepath, width = 10, height = 10)
    for(g in rownames(gene.table)[1:min(n.genes, nrow(gene.table))]){
        if(!is.null(img) && !is.null(plot.params)){
            p <- spatial_plot(
                barcodes = rownames(counts), 
                ids = ids, 
                cluster = counts[,g], 
                img = img, 
                mode = "continuous", 
                plot.params = plot.params, 
                spot.col = colour, 
                title = g
                )
        }else{
            p <- spatial_plot(
                barcodes = rownames(counts), 
                ids = ids, 
                cluster = counts[,g], 
                img = NULL, 
                mode = "continuous", 
                title = g
                )
        }
        plot(p)
    }
    dev.off()
}
