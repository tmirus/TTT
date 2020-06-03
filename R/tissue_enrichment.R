tissue_enrichment <- function(counts, ids, img, filepath, filename, n.cores = 4, spatial.file.prefix = NULL, db_dataset = 'mmusculus_gene_ensembl', 
    enrichment.list = NULL, term.cutoff = 10, spot.cutoff = 1, n.cluster.spots = 6, 
    n.cluster.terms = 6, verbose = TRUE, output_folder = NULL, method = 'normal',
    nx = 35, ny = 33, ox = -1000/70, oy = 1000/32, plot = TRUE){

    # create the enrichment data
    if(is.null(enrichment.list)){
        enrichment.list <- calculate_spot_enrichments(counts, n.cores, db_dataset)
    }
    
    clustering.list <- cluster_spots(counts, spot.enrichments, 
        term.cutoff, spot.cutoff, n.cluster.spots, 
        n.cluster.terms, n.cores, verbose, output_folder, 
        method)
    
    if(plot){
        plot_enrichment_clustering(counts, ids, img, clustering.list, filepath, 
        filename, spatial.file.prefix = NULL, spot.cluster.names = NULL,
        nx, ny, ox, oy)
    }
    

    return(list(enrichment.info = enrichment.list, clustering.info = clustering.list))
}