% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatial_plot.R
\name{spatial_plot}
\alias{spatial_plot}
\title{create 2d spatial plots of ST data with optional image data}
\usage{
spatial_plot(
  barcodes,
  ids,
  cluster,
  img = NULL,
  mode = "discrete",
  plot.params = list(nx = 35, ny = 33, ox = 0, oy = 0),
  spot.col = "black",
  title = "",
  indicator = NULL,
  spot.size = NULL
)
}
\arguments{
\item{barcodes}{character vector, barcodes / names of the spots to be plotted}

\item{ids}{data frame with rownames corresponding to barcodes / spot names, assigning
spatial coordinates to each spot; coordinates are contained in columns X and Y and conform to
the same standards as the output of process_input()}

\item{cluster}{numeric vector or factor containing information to be plotted in the same order as barcodes; usually colour (e.g. for clustering) or size (e.g. total RNA per spot);
how this information is used is determined by the parameter mode}

\item{img}{EBImage image object the information will be plotted on; default NULL}

\item{mode}{character specifying how the information in 'cluster' should be visualized. Must be one of 'discrete' or 'continuous'.
"discrete" will be encoded in colors (e.g. clustering information), 'continuous' will be displayed by size (e.g. amount of RNA)}

\item{plot.params}{list of parameters needed for good spatial visualization as returned by plot_adjustment}

\item{spot.col}{character specifying the color of the spots if mode is "continuous"}

\item{title}{character, title of the plot"}

\item{indicator}{string, either "col" or "size", depending on wether information should be conveyed by spot size or colour.
If not specified the default for mode "discrete" will be "col" and the default for mode "continuous" will be "size".}
}
\value{
ggplot2 object (plot)
}
\description{
create 2d spatial plots of ST data with optional image data
}
