% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oCombineNet.r
\name{oCombineNet}
\alias{oCombineNet}
\title{Function to combine networks from a list of igraph objects}
\usage{
oCombineNet(
  list_ig,
  combineBy = c("union", "intersect"),
  attrBy = c("intersect", "union"),
  keep.all.vertices = FALSE,
  verbose = TRUE
)
}
\arguments{
\item{list_ig}{a list of "igraph" objects or a "igraph" object}

\item{combineBy}{how to resolve edges from a list of "igraph" objects. It can be "intersect" for intersecting edges and "union" for unionising edges (by default)}

\item{attrBy}{the method used to extract node attributes. It can be "intersect" for intersecting node attributes (by default) and "union" for unionising node attributes}

\item{keep.all.vertices}{logical to indicate whether all nodes are kept when intersecting edges. By default, it sets to false}

\item{verbose}{logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display}
}
\value{
an object of class "igraph"
}
\description{
\code{oCombineNet} is supposed to combine networks from a list of igraph objects.
}
\note{
none
}
\examples{
\dontrun{
g1 <- oDefineNet(network="KEGG_environmental")
g2 <- oDefineNet(network="KEGG_organismal")
ls_ig <- list(g1, g2)
ig <- oCombineNet(ls_ig, combineBy='union', attrBy="intersect", verbose=TRUE)
}
}
\seealso{
\code{\link{oDefineNet}}
}
