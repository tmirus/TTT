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
#' @return ggplot2 object (plot)
#' @export 

spatial_plot <- function(barcodes, ids, cluster, img=NULL, mode="discrete", plot.params = list(nx = 35, ny = 33, ox = -1000/70, oy = 1000/32)){
    theme_transparent <- theme(
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        axis.line = element_blank(), # adding a black line for x and y axis
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_blank(),
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
        temp <- c(ox+(ids[barcodes[i],2])*(1000-2*ox)/nx,oy+(ids[barcodes[i],1])*(1000-2*ny)/ny,cluster[i])
        df <- rbind(df,temp)
    }
    df <- as.data.frame(df)
    colnames(df) <- c("X","Y","cluster")

    if(mode == "discrete"){
        df$cluster <- factor(df$cluster)
        p <- ggplot(df,aes(x=-as.numeric(as.character(X)),y=as.numeric(as.character(Y)),col=cluster)) + xlim(-1000,0) + ylim(0,1000)
    }else if(mode == "continuous"){
        df$cluster <- as.numeric(as.character(df$cluster))
        df[which(df$cluster == 0),"cluster"] <- NA
        p <- ggplot(df,aes(x=-as.numeric(as.character(X)),y=as.numeric(as.character(Y)),size=cluster)) + xlim(-1000,0) + ylim(0,1000)
    }else{
        stop("Invalid mode.")
    }

    if(!is.null(img)){
        gob <- rasterGrob(img)
        p <- p + annotation_custom(gob,-1000,0,0,1000)+geom_point(na.rm=TRUE)+theme(legend.position = "none") + ggtitle(title)
    }else{
        p <- p + geom_point(na.rm = TRUE)+xlim(-1000,0)+ylim(0,1000)+ggtitle(title)
    }
    p <- p + xlab("X") + ylab("Y") + theme_transparent
    return(p)
}
