#' perform enrichment testing using topGO and fisher's exact test for a given gene list
#' 
#' @param all.genes character vector, the reference gene set. for example all genes in the experiment
#' @param genelist character vector, list of interesting genes
#' @param db biomaRt database for the data
#' @param go_ids table assigning GO ids to gene names, obtained using getBM
#' @param sig.level numeric, 0 < sig.level < 1. Significance cutoff, p-values above this threshold are discarded
#' @param max.terms numeric > 0, maximum number of terms to be considered
#' @param node.size ignore GO nodes with less than node.size annotated genes
#' @return table containing information about each significant GO term

enrichment_test <- function(all.genes, genelist, db = NULL, go_ids = NULL, sig.level = 0.05, max.terms = 5000, node.size = 10){
  suppressMessages(library(topGO, quietly = TRUE))
  
  interesting.genes <- genelist
  if(is.null(go_ids)){
    go_ids <- biomaRt::getBM(attributes=c('go_id', 'external_gene_name', 'namespace_1003'), filters = 'external_gene_name', values = all.genes, mart = db, verbose = F)
  }

  gene_2_go <- unstack(go_ids[, c(1,2)])

  interesting.genes <- interesting.genes[which(interesting.genes %in% go_ids[,2])]
  geneList <- factor(as.integer(all.genes %in% interesting.genes))

  names(geneList) <- all.genes
  if(!is.factor(geneList) || length(levels(geneList)) != 2){
	  return(NULL)
  }
  suppressMessages({
  GOdata <- new("topGOdata", ontology="BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = gene_2_go, nodeSize = node.size)
  })

  suppressMessages({
  fisher_result <- runTest(GOdata, algorithm = 'weight01', statistic = 'fisher', verbose = F)
  allGO <- usedGO(GOdata)
  all_res <- GenTable(GOdata, weightFisher = fisher_result, orderBy = 'weightFisher', topNodes = min(max.terms,length(score(fisher_result))))
  })

  all_res <- all_res[order(as.numeric(all_res$weightFisher)),]
  res.table.p <- all_res[which(as.numeric(all_res$weightFisher) <= sig.level),]
  return(res.table.p)
}
