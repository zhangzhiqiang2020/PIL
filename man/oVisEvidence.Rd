% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oVisEvidence.r
\name{oVisEvidence}
\alias{oVisEvidence}
\title{Function to visualise evidence for prioritised genes in a gene network}
\usage{
oVisEvidence(
  xTarget,
  g = NA,
  nodes = NULL,
  node.info = c("smart", "simple", "none"),
  neighbor.order = 1,
  neighbor.seed = TRUE,
  neighbor.top = NULL,
  largest.comp = TRUE,
  show = TRUE,
  colormap = "ggplot2",
  legend.position = "topleft",
  legend.horiz = FALSE,
  mtext.side = 3,
  verbose = TRUE,
  edge.width = NULL,
  vertex.size = NULL,
  vertex.size.nonseed = NULL,
  vertex.label.color = "steelblue",
  vertex.label.color.nonseed = NULL,
  ...
)
}
\arguments{
\item{xTarget}{an object of class "dTarget", "sTarget" or "eTarget"}

\item{g}{an object of class "igraph". If NA, the 'metag' will be used, which is part of the input object "xTarget". If provided, it must have a node attribute called 'priority'}

\item{nodes}{which node genes are in query. If NULL, the top gene will be queried}

\item{node.info}{tells the additional information used to label nodes. It can be one of "none" (only gene labeling), "smart" for (by default) using three pieces of information (if any): genes, 5-star ratings, and associated ranks (marked by an @ icon)}

\item{neighbor.order}{an integer giving the order of the neighborhood. By default, it is 1-order neighborhood}

\item{neighbor.seed}{logical to indicate whether neighbors are seeds only. By default, it sets to true}

\item{neighbor.top}{the top number of the neighbors with the highest priority. By default, it sets to NULL to disable this parameter}

\item{largest.comp}{logical to indicate whether the largest component is only retained. By default, it sets to true for the largest component being left}

\item{show}{logical to indicate whether to show the graph}

\item{colormap}{short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta), and "ggplot2" (emulating ggplot2 default color palette). Alternatively, any hyphen-separated HTML color names, e.g. "lightyellow-orange" (by default), "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}}

\item{legend.position}{the legend position. If NA, the legend is not shown}

\item{legend.horiz}{logical specifying the legend horizon. If TRUE, set the legend horizontally rather than vertically}

\item{mtext.side}{the side of marginal text. If NA, it is not shown}

\item{verbose}{logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display}

\item{edge.width}{the width of the edge. If NULL, the width edge is proportional to the 'weight' edge attribute (if existed)}

\item{vertex.size}{the size of each vertex. If null, each vertex has the size proportional to the degree of nodes}

\item{vertex.size.nonseed}{the size of each nonseed vertex. If null, each vertex has the size proportional to the degree of nodes}

\item{vertex.label.color}{the color of vertex labels}

\item{vertex.label.color.nonseed}{the color of nonseed vertex labels}

\item{...}{additional graphic parameters. See \url{http://igraph.org/r/doc/plot.common.html} for the complete list.}
}
\value{
a subgraph, an object of class "igraph".
}
\description{
\code{oVisEvidence} is supposed to visualise evidence for prioritised genes in a gene network. It returns an object of class "igraph".
}
\examples{
\dontrun{
## TNFRSF1A
oVisEvidence(xTarget, nodes="TNFRSF1A", neighbor.order=1, neighbor.seed=TRUE, neighbor.top=NULL, vertex.label.color="black", vertex.label.cex=0.7, vertex.label.dist=0.6, vertex.label.font=1, vertex.label.family="Arial", legend.position="bottomleft", legend.horiz=TRUE, newpage=FALSE)
## UBA52
oVisEvidence(xTarget, nodes="UBA52", neighbor.order=1, neighbor.seed=TRUE, neighbor.top=20, vertex.label.color="black", vertex.label.cex=0.7, vertex.label.dist=0.6, vertex.label.font=1, legend.position="bottomleft", legend.horiz=TRUE, newpage=FALSE)
}
}
\seealso{
\code{\link{oColormap}}
}
