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

spatial_plot <- function(barcodes, ids, cluster, img = NULL, img_size = NULL, mode="discrete", plot.params = NULL, spot.col = "black", title = "", indicator = NULL, spot.size = NULL){
    # check input
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
  
    # set theme for spatial tissue plots 
    theme_transparent <- theme(
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        axis.line = element_blank(), # adding a black line for x and y axis
        axis.text = element_blank(),
        axis.title = element_blank(),
        #plot.title = element_text(size = 24),
        axis.ticks = element_blank()
    )

    # plot adjustment parameters
    if(!is.null(plot.params)){	  
    nx <- plot.params$nx
    ny <- plot.params$ny
    ox <- plot.params$ox
    oy <- plot.params$oy
    }else{
	nx <- 33
    	ny <- 31
	ox <- 0
	oy <- 0
    }
 

    # standardize image size
    if(!is.null(img)){
      
	
      if(is.null(img_size)){
      	dim_x <- dim(img)[1]
      dim_y <- dim(img)[2]
	      img_size <- min(dim_x, dim_y) 
      }else{
	      dim_x <- img_size
	      dim_y <- img_size
      }
      
      img <- EBImage::resize(img, w = img_size, h = img_size)  
    }else{
	    if(is.null(img_size) || is.null(plot.params)){
	    dim_x <- 1000
	    dim_y <- 1000
	    
	    # ignore plot params because they might be wrong for image width
	    nx <- max(ids$X)
	    ny <- max(ids$Y)
	    ox <- 0
	    oy <- 0
	    }
    }
    plot.params <- list(nx = nx, ny = ny, ox = ox, oy = oy)

    # convert grid coordinates to image coordinates
    df <- c()
    for(i in 1:length(barcodes)){
        temp <- c(
          ox+(ids[barcodes[i],"X"]-2)*(dim_x)/nx,
          oy+(ids[barcodes[i],"Y"]-2)*(dim_y)/ny, 
          cluster[i]
          )
    	
        df <- rbind(df, temp)
    }
    df <- as.data.frame(df)
    
    colnames(df) <- c("X", "Y", "cluster")

    # plotting mode
    # discrete uses input variable as factor, e.g. clustering
    # or converts it to fact
    # information coded via color
    if(mode == "discrete"){
      # factor
	    if(is.factor(cluster)){
		    df$cluster <- cluster
	    }else{
        	df$cluster <- factor(df$cluster, 
        	                     levels = c(unique(as.character(df$cluster))))
	    }
      # create basic plot variable with color
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
          # if specified otherwise, create basic plot variable with size
          p <- ggplot(df, aes(x=-as.numeric(as.character(X)), y=as.numeric(as.character(Y)),size=cluster)) +
            theme(legend.position = "none")
        }
      # continuous mode for numeric values, e.g. gene expression, lib size
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

    # set background
    if(!is.null(img)){
        gob <- grid::rasterGrob(img)
        p <- p + annotation_custom(gob,-dim_x,0,0,dim_y)
    }
    
    p <- p + 
      geom_point(na.rm = TRUE) + 
      ggtitle(title) + 
      xlab("") + 
      ylab("") + 
      theme_transparent + 
      coord_fixed(ratio = 1, xlim =c(-dim_x,0), ylim=c(0,dim_y))
    # if color does not code information allow manual setting
    if(indicator == "size") p <- p + geom_point(col = spot.col)
    else if(!is.null(spot.size)){
      # maual setting of spot size via parameter
	    p <- p + geom_point(na.rm = TRUE, size = spot.size)
    }

    return(p)
}
