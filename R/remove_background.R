#' remove background spots (not covered by tissue) from ST data based on tissue image.
#' Manual checking of the result and adjustment of the threshold is needed.
#' 
#' @param img image of the ST slide, cropped to the spot area. image read by EBimage::readImage()
#' @param nx number of spots shown in the image on the horizontal axis
#' @param ny number of spots shown in the image on the vertical axis
#' @param ids data frame assigning barcodes / spot names to spatial coordinates
#' @param counts non-negative numeric matrix containing gene counts, 
#' rows correspond to spots, columns correspond to genes
#' @param threshold relative brightness. spots above this threshold are discarded
#' @param plot.params list as returned by plot_adjustment(). optional but recommended for optimal results.
#' @param force.manual logical, forces the function to use manual threshold instead of attempting to find best threshold with Expectation Minimization
#' @details nx and ny depend on the cropping of the image. the area of measurements in an ST experiment is
#' @return list containing four entries\cr
#' 1) spots.to.keep - vector of barcodes / spot names of the spots that should not be removed
#' 2) spots.keep.clustering - vector of barcodes / spot names of the spots that should not be removed based on clustering in addition to image analysis
#' 3) image - image showing the areas that are retained.
#' typically framed by spots for which there exist no measurements. Depending on whether these are contained 
#' in the image or not, different values need to be chosen for nx, ny. If they are included, the 
#' default values of nx=35 and ny=33 are sensible. If not, it should be nx=33 and ny=31.\cr
#' The image is assumed to be oriented such that the rectangular 4x4 array of spots in one of the corners
#' is in the lower right corner and on the x-axis (left to right) there are more spots than on the y-axis.\cr
#' 4) clustering.tsne tSNE embedding of spots coloured by classification (background/tissue)\cr
#' 5) brightness.plot density plot of brightness values, fits for two groups and decision boundary for dividing tissue and background
#' @export

remove_background <- function(img = NULL, nx = 35, ny=33, ids, counts, threshold = 0.7, plot.params = NULL, force.manual = FALSE){
  if(!is.null(img)){
    # reduce to one pixel per spot
    img <- EBImage::channel(img, "gray")
    img <- EBImage::resize(img, 1000, 1000)
    
  
    if(!is.null(plot.params)){
      nx <- plot.params$nx
      ny <- plot.params$ny
      ox <- plot.params$ox
      oy <- plot.params$oy
      dx <- (1000-2*ox)/nx
      dy <- (1000-2*oy)/ny
  
      ox1 <- as.integer(ox)
      ox2 <- as.integer(1000-ox)
      oy1 <- as.integer(oy)
      oy2 <- as.integer(1000-oy)
  
      img <- img[ox1:ox2, oy1:oy2]
    }
    
    img.reduced <- EBImage::resize(img, w=nx,h=ny)
    names <- rownames(counts)
  
    # blur to reduce the risk of removing spots within tissue
    for(x in 1:nx){
      for(y in 1:ny){
        x.neighbours <- seq(x-1,x+1)
        y.neighbours <- seq(y-1,y+1)
        if(any(x.neighbours < 1 | x.neighbours > nx)) 
          x.neighbours <- x.neighbours[-which(x.neighbours < 1 | x.neighbours > nx)]
        if(any(y.neighbours < 1 | y.neighbours > ny)) 
          y.neighbours <- y.neighbours[-which(y.neighbours < 1 | y.neighbours > ny)]
        EBImage::imageData(img.reduced)[x,y] <- mean(EBImage::imageData(img.reduced)[x.neighbours, y.neighbours])
      }
    }
    
    EBImage::imageData(img.reduced) <- EBImage::imageData(img.reduced) / max(EBImage::imageData(img.reduced))
  
    # select the spots to remove
    brightness.values <- as.vector(EBImage::imageData(img.reduced))
    if(!force.manual){
    em.result <- mclust::Mclust(data = brightness.values, G = 2, modelNames = "V")
    
    # get solution for best split
    means <- em.result$parameters$mean
    sds <- em.result$parameters$variance$sigmasq
    pros <- em.result$parameters$pro
    
    f <- function(x){
      return(pros[1] * dnorm(x, means[1], sqrt(sds[1])) - pros[2] * dnorm(x, means[2], sqrt(sds[2])))
    }

    is.between <- function(x1, x2, y){
      if((y > x1 && y < x2) || (y > x2 && y < x1)){
        return(TRUE)
      }else{
        return(FALSE)
      }
    }
    
    x0 <- cmna::bisection(f, means[1], means[2],tol = 1e-5, m = 10000)
    if(!is.between(means[1], means[2], x0) || abs(f(x0) > 0.01)){
      x0 <- threshold
    }else{
      cat("Brightness Threshold determined by EM: ", x0, "\n")
    }
    
    dist1 <- pros[1] * dnorm(x = seq(0,1,by = 0.01), mean = means[1], sd = sqrt(sds[1]))
    dist2 <- pros[2] * dnorm(x = seq(0,1,by = 0.01), mean = means[2], sd = sqrt(sds[2]))
    brightness.df <- data.frame(
      brightness = rep(seq(0,1,by = 0.01), 3), 
      density = c(dist1, dist2, density(brightness.values, n = 101, from = 0, to = 1)$y),
      curve = c(rep("group1", 101), rep("group2", 101), rep("real", 101))
    )
    
    brightness.plot <- ggplot(brightness.df, 
                              aes(
                                x = brightness, 
                                y = density, 
                                group = curve, 
                                col = curve)
                              ) + 
      geom_line() +
      geom_vline(xintercept = x0, col = "black") +
      ggtitle("Brightness distribution, fits and threshold")
    }else{
	x0 <- threshold
	brightness.density <- as.data.frame(density(brightness.values, n = 101, from = 0, to = 1)[1:2])
    	brightness.plot <- ggplot(brightness.density, aes(x = x, y = y)) +
		geom_line() +
		ggtitle("Brightness distribution and threshold") +
		xlab("Brightness") +
		ylab("Density") +
		geom_vline(xintercept = x0, col = "black")
    }
    
    EBImage::imageData(img.reduced)[EBImage::imageData(img.reduced) > x0] <- 1
    to.remove <- t(matrix(EBImage::imageData(img.reduced) == 1, nrow=nx))
    
    # store in more compatible format
    spots.to.keep <- c()
    for(x in 1:nx){
      for(y in 1:ny){
        if(!to.remove[y,x]) spots.to.keep <- rbind(spots.to.keep, c(max(ids[,"X"])-x,max(ids[,"Y"])-y))
      }
    }
    spots.to.keep <- as.data.frame(spots.to.keep)
    colnames(spots.to.keep) <- c("X", "Y")
    
    # return the original image with blanks where spots were removed
    dx <- dim(img)[1] / nx
    dy <- dim(img)[2] / ny
    for(x in 1:nx){
      for(y in 1:ny){
        if(to.remove[y,x]){
          EBImage::imageData(img)[(as.integer((x-1)*dx+1)):as.integer(x*dx), (as.integer((y-1)*dy+1)):as.integer(y*dy)] <- 1
        }
      }
    }
    if(!is.null(ids) && !is.null(names)){
      spots <- sapply(names,
                      function(x){
                        if(length(intersect(which(spots.to.keep[,"X"] == ids[x,2]), which(spots.to.keep[,"Y"] == ids[x,1])))>0){
                          return(TRUE)
                        }else{
                          return(FALSE)
                        }
                      }
      )
      spots <- names(spots)[spots]
      spots.to.keep <- spots
    }
    ids.reduced <- ids[spots,]
    ids.reduced <- clean_ids(ids.reduced)
    spots.to.keep <- rownames(ids.reduced)
  
    # tSNE embedding
    tsne <- Rtsne::Rtsne(counts, check_duplicates = F)$Y
    # create comparatively large number of clusters
    clustering <- cluster::pam(tsne, 15)$clustering
    names(clustering) <- rownames(counts)
  
    # keep all clusters with a certain amount of spots in spots.to.keep
    clusters.to.keep <- c()
    for(cl in unique(clustering)){
      if(sum(names(clustering[clustering == cl]) %in% spots.to.keep) > 0.75 * length(which(clustering == cl))){
        clusters.to.keep <- c(clusters.to.keep, names(clustering[clustering == cl]))
      }
    }
  
    # vector containing background information by clustering
    background.vec <- rep("background", nrow(counts))
    names(background.vec) <- rownames(counts)
    background.vec[clusters.to.keep] <- "tissue"
    background.vec <- as.factor(background.vec)
    
    # vector containing background information by brightness alone
    background.temp <- rep("background", nrow(counts))
    names(background.temp) <- rownames(counts)
    background.temp[spots.to.keep] <- "tissue"
    background.temp <- as.factor(background.temp)
    
    # create plots
    df <- data.frame(tsne, background.vec, clustering, background.temp)
    colnames(df) <- c("tSNE1", "tSNE2", "background", "clustering", "background.temp")
  
    p <- ggplot(df, aes(x=tSNE1, y = tSNE2, col = background)) + geom_point() + labs(col = "class")
    p.bg <- ggplot(df, aes(x = tSNE1, y = tSNE2, col = as.factor(clustering))) + geom_point() + labs(col = "clustering")
    p.bg.temp <- ggplot(df, aes(x = tSNE1, y = tSNE2, col = background.temp)) + geom_point() + labs(col = "class\nby\nbrightness")
  
  
    # reduce ids to the tissue spots
    if(length(clusters.to.keep) > 0){
      ids.reduced <- ids[clusters.to.keep,]
      ids.reduced <- clean_ids(ids.reduced)
    }
    clusters.to.keep <- rownames(ids.reduced)
  }else{
	  # if no image available
	  # cluster all spots using pam on tSNE embedding
	  tsne <- Rtsne::Rtsne(counts, check_duplicates = F)$Y
	  clustering <- cluster::pam(tsne, 20)$clustering
	  clustering.coarse <- cluster::pam(counts, 2)$clustering

	  # try to find a cutoff for total library size and/or entropy
	  # to identify background spot
	  
	  # identify clusters with relatively low total expression per spot
	  clust.mean.exp <- sapply(unique(clustering), function(x){mean(rowSums(counts[which(clustering == x),]))})
	  clusters <- unique(clustering)[which(clust.mean.exp < summary(clust.mean.exp)[4])]

	  # if the majority of spots within one of the coarse clusters is in the clusters with low expression, remove it
	  if(length(which(clustering[which(clustering.coarse == 1)] %in% clusters)) / length(which(clustering.coarse == 1)) > 0.8 && length(which(clustering.coarse[which(clustering %in% clusters)] == 1)) / length(which(clustering %in% clusters)) > 0.9){
		  spots.to.keep <- rownames(counts)[-which(clustering.coarse == 1)]
	  }else{
		 if(length(which(clustering[which(clustering.coarse == 2)] %in% clusters)) / length(which(clustering.coarse == 2)) > 0.8 && length(which(clustering.coarse[which(clustering %in% clusters)] == 2)) / length(which(clustering %in% clusters)) > 0.9){
			 spots.to.keep <- rownames(counts)[-which(clustering.coarse == 2)]
		 } 
	  }
	  background.vec <- rep(0, nrow(counts))
	  names(background.vec) <- rownames(counts)
	  background.vec[spots.to.keep] <- 1
	  background.vec <- as.factor(!as.logical(background.vec))
	  df <- data.frame(tsne, background.vec)
	  colnames(df) <- c("tSNE1", "tSNE2", "background")
  	p <- ggplot(df, aes(x = tSNE1, y = tSNE2, col = background)) + geom_point()
  	p.bg <- NULL
  	p.bg.temp <- NULL

	  clusters.to.keep <- spots.to.keep
  }

  return(list(
    spots.to.keep = spots.to.keep, 
    spots.keep.clustering = clusters.to.keep, 
    image = img,
    clustering.tsne = p,
    fine.clustering.tsne = p.bg,
    temp.clustering.tsne = p.bg.temp,
    brightness.plot = brightness.plot
))
}
