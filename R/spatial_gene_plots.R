#' create spatial plots for the top n.genes genes in a given table and store them as pdf
#' 
#' @param gene.table data frame containing interesting genes, as returned by filter_genes (or entry differential.genes of analyze_clustering output)
#' @param counts non-negative numeric matrix containing gene expression for all spots (spots x genes)
#' @param ids barcode data frame assigning spatial coordinates to all spots; column names 'X' and 'Y'
#' @param img ST image (EBImage object) if available. Will be used as plotting background if plot.params is provided. If not available, leave default NULL.
#' @param plot.params list created by plot_adjustment; only needed if img is provided
#' @param n.genes integer, number of genes to be plotted (starting from the top of gene.table); optional
#' @param filepath character specifying the directory in which the output file ('spatial_genes.pdf') should be stored
#' @param filename character specifying the name of the pdf file; default "spatial_genes.pdf"
#' @param colour character specifying the colour of the spots in the gene expression plots; default "black"
#' @return NULL, plots to file 'filepath'/spatial_genes.pdf
#' @export

spatial_gene_plots <- function(gene.table, counts, ids, img = NULL, plot.params = NULL, n.genes = NULL, filepath = NULL, filename = "spatial_genes.pdf", colour = "black"){
    # set output path
    if(is.null(filepath)){
        filepath <- paste0("./", filename)
    }else{
        if(!dir.exists(filepath)){
            dir.create(filepath)
        }
        filepath <- paste(filepath, filename, sep = "/")
    }
    
    # set n.genes
    if(is.null(n.genes)){
        n.genes <- nrow(gene.table)
    }else{
        n.genes <- min(n.genes, nrow(gene.table))
    }

    # output pdf
    pdf(filepath, width = 10, height = 10)
    # plot all genes in table (or up to specified number)
    for(g in gene.table$gene[1:n.genes]){
        # plot without background if no image 
        # or no plot.params available
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
            # plot with background if possible
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
