#' process input data from ST pipeline to make sure they conform 
#' to the format the package expects
#' 
#' @param counts non-negative numeric matrix containing gene counts, 
#' rows correspond to spots, columns correspond to genes
#' @param ids data frame or matrix assigning spatial coordinates to the spots,
#' the spot 0,0 is in the right lower corner of the 
#' corresponding tissue slide and the number of spots along the X-axis is 
#' higher than the number of spots on the Y-axis.
#' optional if spatial information is contained in rownames of counts
#' @param img_path character, optional. file path to ST slide image
#' @param img_size, numeric value giving with and height the image will be scaled to
#' @param separate_by only necessary if spatial information is given as
#' in the output of ST Pipeline (Navarro et al.), e.g. 32x2; default x
#' @param force_counts logical, in case there are less genes than spots this needs to be TRUE; 
#' default FALSE
#' @param force_indices logical, in case the range of x and y indices is not the same as in the
#' ST barcode files
#' @param dup.sep character, optional; if duplicate gene names are distinguished by e.g. _1, _2
#' at the end, pass dup.sep="_" to try to combine these columns; default "_"
#' @param n.gene.cutoff integer > 0, number of genes that need to be expressed in a spots in order 
#' for the spot not to be removed from the data
#' @param verbose logical, default TRUE
#' @return list with three entries:\cr
#' 1) counts - count matrix
#' 2) ids - barcode data frame assigning spatial positions to spots
#' 3) img - image (EBImage object), if img_path is supplied
#' 4) counts.normalized - normalized counts matrix, if normalize is true; NULL else
#' @export

process_input <- function(counts, ids = NULL, img_path = NULL, img_size = 1000, separate_by = "x", force_counts = FALSE, force_indices = FALSE, dup.sep = "_", n.gene.cutoff = 100, verbose = TRUE){
  if(is.character(counts)){
    if(file.exists(counts)){
      counts <- as.matrix(read.table(counts, check.names = FALSE))
      if(!is.numeric(counts)){
        stop("An error occurred while reading the counts matrix. Try reading it manually and pass it as parameter.")
      }
    }
  }
  # check input parameters
  if(!is.matrix(counts)){
    stop("counts needs to be supplied as matrix")
  }
  if(!is.numeric(counts) || any(counts < 0)){
    stop("counts must be a non-negative numeric matrix")
  }
  if(nrow(counts) >= ncol(counts) && !force_counts){
    stop("There seem to be more spots than genes. Please check matrix
         orientation or set force_counts = TRUE to use the data as it is.")
  }
  if(n.gene.cutoff <= 0 || round(n.gene.cutoff) != n.gene.cutoff){
    stop("Invalid value for 'n.gene.cutoff'. Must be integer > 0.")
  }
  if(!is.null(ids)){
    if(! (is.matrix(ids) || is.data.frame(ids)) ){
      stop("ids must be a matrix or a data frame")
    }
    if(ncol(ids) != 2){
      stop("ids must contain two columns (X and Y)")
    }
    if(nrow(ids) < nrow(counts)){
      stop("ids does not have enough rows to contain spatial information for all spots in counts.")
    }
    if(!all(rownames(counts) %in% rownames(ids))){
      stop("One or more spots are missing spatial information.")
    }
  }else{
    if(! all(grep(pattern=separate_by, x = rownames(counts)) == 1:nrow(counts))){
      stop("Splitting character split_by does not seem to be present in all spot names.")
    }
  }
  if(!is.null(img_path)){
	  img <- read_image(img_path, hw = img_size)
  }else{
	  img <- NULL
  }
 
  # prepare either input ids or convert rownames to ids data frame
  if(!is.null(ids)){
    if(is.null(colnames(ids)) || !all(colnames(ids) %in% c("X", "Y"))){ 
    	stop("ids needs to have columns X and Y.") 
    }
  }else{
    # extract ids from rownames of counts and create a data frame
    ids.parts <- strsplit(rownames(counts), split = separate_by)
    coords1 <- as.numeric(sapply(ids.parts, function(x){x[[1]]}))
    coords2 <- as.numeric(sapply(ids.parts, function(x){x[[2]]}))
    
    if(any(is.na(c(coords1,coords2)))){
      stop("rownames of counts could not be converted to coordinates. Make sure the rownames
           conform to the pattern {coord1}{split_by}{coord2} or provide an ids matrix
           mapping barcodes to spatial positions")
    }
    if(any(c(coords1, coords2) < 0)){
      stop("Invalid ids extracted from rownames. Coordinates must be >= 0.")
    }
    
    r1 <- range(coords1)
    r2 <- range(coords2)
    if( (r1[2] - r1[1]) > (r2[2] - r2[1]) ){
      ids <- data.frame(Y = coords2, X = coords1)
    }else{
      ids <- data.frame(Y = coords1, X = coords2)
    }
    rownames(ids) <- rownames(counts)
  }
  
  # adjust coordinates for compatibility 
  # with the package if possible
  
  # deal with potential duplicates
  if(!is.null(dup.sep)){
    gene.names <- strsplit(colnames(counts), split = dup.sep)
    suppressWarnings(gene.names <- sapply(gene.names, function(x){
    if(!is.na(as.numeric(x[length(x)])) && !(nchar(x[length(x)]) > 1)){
        if(length(x) == 2){
        return(x[1])
      }else{
        return(paste(x[-length(x)], collapse = "_"))
      }
    }
    return(paste(x, collapse = "_"))
    }))
    colnames(counts) <- gene.names
  }
  
  
  # find all duplicate names, sum them up and remove superfluous columns
  if(sum(duplicated(colnames(counts))) > 0){
    dup.genes <- unique(colnames(counts)[duplicated(colnames(counts))])
    if(verbose) cat("Duplicates found for ", length(dup.genes), " genes.\nRemoving...\n")
    removal.indices <- c()
    for(g in dup.genes){
      gene.idx <- which(colnames(counts) == g)
      counts[,gene.idx[1]] <- rowSums(counts[, gene.idx])
      removal.indices <- c(removal.indices, gene.idx[-1])
    }
    counts <- counts[,-removal.indices]
  }
  
  # remove spots with exceptionally low number of expressed genes
  n.genes <- apply(counts, 1, function(x) sum(x>0))
  if(any(n.genes <= n.gene.cutoff)){
    counts <- counts[-which(n.genes < n.gene.cutoff),]
  }
  
  # normalization with sctransform
  
  counts.normalized <- normalize_counts(counts)
  
  return(list(counts = counts, ids = ids, img = img, counts.normalized = counts.normalized))
}
