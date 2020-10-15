#' create 2d spatial plots of ST data with optional image data
#' 
#' @param barcodes character vector, barcodes / names of the spots to be plotted
#' @param ids data frame with rownames corresponding to barcodes / spot names, assigning
#' spatial coordinates to each spot; coordinates are contained in columns X and Y and conform to
#' the same standards as the output of process_input()
#' @param cluster numeric vector or factor containing information to be plotted in the same order as barcodes; usually colour (e.g. for clustering) or size (e.g. total RNA per spot);
#' how this information is used is determined by the parameter mode
#' @param img EBImage image object the information will be plotted on; default NULL
#' @param mode character specifying how the information in 'cluster' should be visualized. Must be one of 'discrete' or 'continuous'.
#' "discrete" will be encoded in colors (e.g. clustering information), 'continuous' will be displayed by size (e.g. amount of RNA)
#' @param plot.params list of parameters needed for good spatial visualization as returned by plot_adjustment
#' @param spot.col character specifying the color of the spots if mode is "continuous"
#' @param title character, title of the plot"
#' @param indicator string, either "col" or "size", depending on wether information should be conveyed by spot size or colour.
#' If not specified the default for mode "discrete" will be "col" and the default for mode "continuous" will be "size". 
#' @return ggplot2 object (plot)
#' @export 

spatial_plot <- function(barcodes, ids, cluster, img=NULL, mode="discrete", plot.params = list(nx = 35, ny = 33, ox = 0, oy = 0), spot.col = "black", title = "", indicator = NULL, spot.size = NULL){
    if(!mode %in% c("discrete", "continuous")){
      stop("Invalid value for parameter 'mode'. Must be one of 'discrete' or 'continuous'.")
    }
    if(!is.null(indicator)){
      if(!indicator %in% c("col", "size")){
        stop("Invalid value for parameter 'indicator'. Must be one of 'col' or 'size'.")
      }
    }else{
      if(mode == "discrete"){
        indicator <- "col"
      }else{
        indicator <- "size"
      }
    }
  
    theme_transparent <- theme(
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        axis.line = element_blank(), # adding a black line for x and y axis
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 30),
        axis.ticks = element_blank()
    )

    nx <- plot.params$nx
    ny <- plot.params$ny
    ox <- plot.params$ox
    oy <- plot.params$oy
 
    # standardize image size
    if(!is.null(img)){
      img <- EBImage::resize(img,w = 1000,h=1000)
    }

    df <- c()
    # convert ids to image coordinates
    for(i in 1:length(barcodes)){
        temp <- c(ox+(ids[barcodes[i],"X"]-2)*(1000-ox)/nx,oy+(ids[barcodes[i],"Y"]-2)*(1000-oy)/ny, cluster[i])
        df <- rbind(df, temp)
    }
    df <- as.data.frame(df)
    colnames(df) <- c("X", "Y", "cluster")

    if(mode == "discrete"){
	    if(is.factor(cluster)){
		    df$cluster <- cluster
	    }else{
        	df$cluster <- factor(df$cluster, levels = c(unique(as.character(df$cluster))))
	    }
        if(indicator == "col"){
          p <- ggplot(
            df, 
            aes(
              x=-as.numeric(as.character(X)), 
              y=as.numeric(as.character(Y)),
              col=cluster
              )
            )
        }else{
          p <- ggplot(df, aes(x=-as.numeric(as.character(X)), y=as.numeric(as.character(Y)),size=cluster)) +
            theme(legend.position = "none")
        }
    }else{
        df$cluster <- as.numeric(as.character(df$cluster))
        df[which(df$cluster == 0),"cluster"] <- NA
        if(indicator == "size"){
          p <- ggplot(df,aes(x=-as.numeric(as.character(X)),y=as.numeric(as.character(Y)),size=cluster)) +
            theme(legend.position = "none")
        }else{
          p <- ggplot(df,aes(x=-as.numeric(as.character(X)),y=as.numeric(as.character(Y)),col=cluster))
        }
    }

    if(!is.null(img)){
        gob <- grid::rasterGrob(img)
        p <- p + annotation_custom(gob,-1000,0,0,1000)
    }
    #suppressWarnings({
    p <- p + 
      geom_point(na.rm = TRUE) + 
      ggtitle(title) + 
      xlab("") + 
      ylab("") + 
      theme_transparent + 
      coord_fixed(ratio = 1, xlim =c(-1000,0), ylim=c(0,1000))
	i#    })
    if(indicator == "size") p <- p + geom_point(col = spot.col)
    else if(!is.null(spot.size)){
	    p <- p + geom_point(na.rm = TRUE, size = 5)
    }

    return(p)
}
