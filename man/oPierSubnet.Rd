% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oPierSubnet.r
\name{oPierSubnet}
\alias{oPierSubnet}
\title{Function to identify a gene network from top prioritised genes}
\usage{
oPierSubnet(
  pNode,
  priority.quantile = 0.1,
  network = NA,
  STRING.only = NA,
  network.customised = NULL,
  subnet.significance = 0.01,
  subnet.size = NULL,
  test.permutation = FALSE,
  num.permutation = 100,
  respect = c("none", "degree"),
  aggregateBy = c("Ztransform", "fishers", "logistic", "orderStatistic"),
  verbose = TRUE,
  placeholder = NULL,
  guid = NULL
)
}
\arguments{
\item{pNode}{an object of class "pNode" (or "sTarget" or "dTarget")}

\item{priority.quantile}{the quantile of the top priority genes. By default, 10% of top prioritised genes will be used for network analysis. If NULL or NA, all prioritised genes will be used}

\item{network}{the built-in network. Currently two sources of network information are supported: the STRING database (version 10) and the Pathway Commons database (version 7). STRING is a meta-integration of undirect interactions from the functional aspect, while Pathways Commons mainly contains both undirect and direct interactions from the physical/pathway aspect. Both have scores to control the confidence of interactions. Therefore, the user can choose the different quality of the interactions. In STRING, "STRING_highest" indicates interactions with highest confidence (confidence scores>=900), "STRING_high" for interactions with high confidence (confidence scores>=700), "STRING_medium" for interactions with medium confidence (confidence scores>=400), and "STRING_low" for interactions with low confidence (confidence scores>=150). For undirect/physical interactions from Pathways Commons, "PCommonsUN_high" indicates undirect interactions with high confidence (supported with the PubMed references plus at least 2 different sources), "PCommonsUN_medium" for undirect interactions with medium confidence (supported with the PubMed references). For direct (pathway-merged) interactions from Pathways Commons, "PCommonsDN_high" indicates direct interactions with high confidence (supported with the PubMed references plus at least 2 different sources), and "PCommonsUN_medium" for direct interactions with medium confidence (supported with the PubMed references). In addition to pooled version of pathways from all data sources, the user can also choose the pathway-merged network from individual sources, that is, "PCommonsDN_Reactome" for those from Reactome, "PCommonsDN_KEGG" for those from KEGG, "PCommonsDN_HumanCyc" for those from HumanCyc, "PCommonsDN_PID" for those froom PID, "PCommonsDN_PANTHER" for those from PANTHER, "PCommonsDN_ReconX" for those from ReconX, "PCommonsDN_TRANSFAC" for those from TRANSFAC, "PCommonsDN_PhosphoSite" for those from PhosphoSite, and "PCommonsDN_CTD" for those from CTD. For direct (pathway-merged) interactions sourced from KEGG, it can be 'KEGG' for all, 'KEGG_metabolism' for pathways grouped into 'Metabolism', 'KEGG_genetic' for 'Genetic Information Processing' pathways, 'KEGG_environmental' for 'Environmental Information Processing' pathways, 'KEGG_cellular' for 'Cellular Processes' pathways, 'KEGG_organismal' for 'Organismal Systems' pathways, and 'KEGG_disease' for 'Human Diseases' pathways. 'REACTOME' for protein-protein interactions derived from Reactome pathways}

\item{STRING.only}{the further restriction of STRING by interaction type. If NA, no such restriction. Otherwide, it can be one or more of "neighborhood_score","fusion_score","cooccurence_score","coexpression_score","experimental_score","database_score","textmining_score". Useful options are c("experimental_score","database_score"): only experimental data (extracted from BIND, DIP, GRID, HPRD, IntAct, MINT, and PID) and curated data (extracted from Biocarta, BioCyc, GO, KEGG, and Reactome) are used}

\item{network.customised}{an object of class "igraph". By default, it is NULL. It is designed to allow the user analysing their customised network data that are not listed in the above argument 'network'. This customisation (if provided) has the high priority over built-in network}

\item{subnet.significance}{the given significance threshold. By default, it is set to NULL, meaning there is no constraint on nodes/genes. If given, those nodes/genes with p-values below this are considered significant and thus scored positively. Instead, those p-values above this given significance threshold are considered insigificant and thus scored negatively}

\item{subnet.size}{the desired number of nodes constrained to the resulting subnet. It is not nulll, a wide range of significance thresholds will be scanned to find the optimal significance threshold leading to the desired number of nodes in the resulting subnet. Notably, the given significance threshold will be overwritten by this option}

\item{test.permutation}{logical to indicate whether the permutation test is perform to estimate the significance of identified network with the same number of nodes. By default, it sets to false}

\item{num.permutation}{the number of permutations generating the null distribution of the identified network}

\item{respect}{how to respect nodes to be sampled. It can be one of 'none' (randomly sampling) and 'degree' (degree-preserving sampling)}

\item{aggregateBy}{the aggregate method used to aggregate edge confidence p-values. It can be either "orderStatistic" for the method based on the order statistics of p-values, or "fishers" for Fisher's method, "Ztransform" for Z-transform method, "logistic" for the logistic method. Without loss of generality, the Z-transform method does well in problems where evidence against the combined null is spread widely (equal footings) or when the total evidence is weak; Fisher's method does best in problems where the evidence is concentrated in a relatively small fraction of the individual tests or when the evidence is at least moderately strong; the logistic method provides a compromise between these two. Notably, the aggregate methods 'Ztransform' and 'logistic' are preferred here}

\item{verbose}{logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display}

\item{placeholder}{the characters to tell the location of built-in RDS files. See \code{\link{oRDS}} for details}

\item{guid}{a valid (5-character) Global Unique IDentifier for an OSF project. See \code{\link{oRDS}} for details}
}
\value{
a subgraph with a maximum score, an object of class "igraph". It has ndoe attributes: signficance, score, type, priority (part of the "pNode" object). If permutation test is enabled, it also has a graph attribute (combinedP) and an edge attribute (edgeConfidence)
}
\description{
\code{oPierSubnet} is supposed to identify maximum-scoring gene subnetwork from a graph with the node information on priority scores, both are part of an object of class "pNode". It returns an object of class "igraph".
}
\note{
The priority score will be first scaled to the range x=[0 100] and then is converted to pvalue-like significant level: 10^(-x). Next, \code{\link{oSubneterGenes}} is used to identify a maximum-scoring gene subnetwork that contains as many highly prioritised genes as possible but a few lowly prioritised genes as linkers. An iterative procedure of scanning different priority thresholds is also used to identify the network with a desired number of nodes/genes. Notably, the preferential use of the same network as used in gene-level prioritisation is due to the fact that gene-level affinity/priority scores are smoothly distributed over the network after being walked. In other words, the chance of identifying such a gene network enriched with top prioritised genes is much higher.
}
\examples{
\dontrun{
# find maximum-scoring subnet with the desired node number=50
subnet <- oPierSubnet(pNode, priority.quantile=0.1, subnet.size=50)
}
}
\seealso{
\code{\link{oPierSubnet}}
}
