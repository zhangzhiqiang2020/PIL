% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oDefineNet.r
\name{oDefineNet}
\alias{oDefineNet}
\title{Function to define a gene network}
\usage{
oDefineNet(
  network = c("STRING_highest", "STRING_high", "STRING_medium", "STRING_low",
    "PCommonsUN_high", "PCommonsUN_medium", "PCommonsDN_high", "PCommonsDN_medium",
    "PCommonsDN_Reactome", "PCommonsDN_KEGG", "PCommonsDN_HumanCyc", "PCommonsDN_PID",
    "PCommonsDN_PANTHER", "PCommonsDN_ReconX", "PCommonsDN_TRANSFAC",
    "PCommonsDN_PhosphoSite", "PCommonsDN_CTD", "KEGG", "KEGGc_CellularProcess",
    "KEGGc_EnvironmentalProcess", "KEGGc_GeneticProcess", "KEGGc_HumanDisease",
    "KEGGc_Metabolism", "KEGGc_OrganismalSystem", "KEGGs_Aging", 
    
    "KEGGs_AminoAcidMetabolism", "KEGGs_BiosynthesisOfOtherSecondaryMetabolites",
    "KEGGs_CancerOverview", "KEGGs_CancerSpecificTypes", "KEGGs_CarbohydrateMetabolism",
    "KEGGs_CardiovascularDisease", "KEGGs_CellGrowthAndDeath", "KEGGs_CellMotility",
    "KEGGs_CellularCommunityEukaryotes", "KEGGs_CirculatorySystem",
    "KEGGs_DevelopmentAndRegeneration", "KEGGs_DigestiveSystem",
    "KEGGs_DrugResistanceAntineoplastic", "KEGGs_EndocrineAndMetabolicDisease",
    "KEGGs_EndocrineSystem", "KEGGs_EnergyMetabolism", "KEGGs_EnvironmentalAdaptation", 

        "KEGGs_ExcretorySystem", "KEGGs_FoldingSortingAndDegradation",
    "KEGGs_GlycanBiosynthesisAndMetabolism", "KEGGs_ImmuneDisease", "KEGGs_ImmuneSystem",
    "KEGGs_InfectiousDiseaseBacterial", "KEGGs_InfectiousDiseaseParasitic",
    "KEGGs_InfectiousDiseaseViral", "KEGGs_LipidMetabolism",
    "KEGGs_MetabolismOfCofactorsAndVitamins", "KEGGs_MetabolismOfOtherAminoAcids",
    "KEGGs_MetabolismOfTerpenoidsAndPolyketides", "KEGGs_NervousSystem",
    "KEGGs_NeurodegenerativeDisease", "KEGGs_NucleotideMetabolism",
    "KEGGs_ReplicationAndRepair", 
     "KEGGs_SensorySystem",
    "KEGGs_SignalTransduction", "KEGGs_SignalingMoleculesAndInteraction",
    "KEGGs_SubstanceDependence", "KEGGs_Translation", "KEGGs_TransportAndCatabolism",
    "KEGGs_XenobioticsBiodegradationAndMetabolism", "REACTOME", "HuRI", "HuRI_all",
    "HuRI_union", "HuRI_litbm"),
  STRING.only = c(NA, "neighborhood_score", "fusion_score", "cooccurence_score",
    "coexpression_score", "experimental_score", "database_score", "textmining_score")[1],
  weighted = FALSE,
  verbose = TRUE,
  placeholder = NULL,
  guid = NULL
)
}
\arguments{
\item{network}{the built-in network. Currently two sources of network information are supported: the STRING database (version 10) and the Pathway Commons database (version 7). STRING is a meta-integration of undirect interactions from the functional aspect, while Pathways Commons mainly contains both undirect and direct interactions from the physical/pathway aspect. Both have scores to control the confidence of interactions. Therefore, the user can choose the different quality of the interactions. In STRING, "STRING_highest" indicates interactions with highest confidence (confidence scores>=900), "STRING_high" for interactions with high confidence (confidence scores>=700), "STRING_medium" for interactions with medium confidence (confidence scores>=400), and "STRING_low" for interactions with low confidence (confidence scores>=150). For undirect/physical interactions from Pathways Commons, "PCommonsUN_high" indicates undirect interactions with high confidence (supported with the PubMed references plus at least 2 different sources), "PCommonsUN_medium" for undirect interactions with medium confidence (supported with the PubMed references). For direct (pathway-merged) interactions from Pathways Commons, "PCommonsDN_high" indicates direct interactions with high confidence (supported with the PubMed references plus at least 2 different sources), and "PCommonsUN_medium" for direct interactions with medium confidence (supported with the PubMed references). In addition to pooled version of pathways from all data sources, the user can also choose the pathway-merged network from individual sources, that is, "PCommonsDN_Reactome" for those from Reactome, "PCommonsDN_KEGG" for those from KEGG, "PCommonsDN_HumanCyc" for those from HumanCyc, "PCommonsDN_PID" for those froom PID, "PCommonsDN_PANTHER" for those from PANTHER, "PCommonsDN_ReconX" for those from ReconX, "PCommonsDN_TRANSFAC" for those from TRANSFAC, "PCommonsDN_PhosphoSite" for those from PhosphoSite, and "PCommonsDN_CTD" for those from CTD. For direct (pathway-merged) interactions sourced from KEGG, it can be 'KEGG' for all, 'KEGGc_Metabolism' for pathways grouped into 'Metabolism', 'KEGGc_GeneticProcess' for 'Genetic Process' pathways, 'KEGGc_EnvironmentalProcess' for 'Environmental Process' pathways, 'KEGGc_CellularProcess' for 'Cellular Process' pathways, 'KEGGc_OrganismalSystem' for 'Organismal System' pathways, and 'KEGGc_HumanDisease' for 'Human Disease' pathways. For KEGG 'Environmental Process' pathways, it can be further into 'KEGGs_ImmuneSystem','KEGGs_EndocrineSystem','KEGGs_CirculatorySystem','KEGGs_DigestiveSystem','KEGGs_ExcretorySystem','KEGGs_NervousSystem','KEGGs_SensorySystem','KEGGs_DevelopmentAndRegeneration','KEGGs_Aging','KEGGs_EnvironmentalAdaptation'. 'REACTOME' for protein-protein interactions derived from Reactome pathways. 'HuRI' for human reference direct interactome (biophysical), while 'HuRI_union' for all HuRI identified so far, 'HuRI_litbm' for literature-derived direct interactions, and 'HuRI_all' for the union of 'HuRI_union' and 'HuRI_litbm'}

\item{STRING.only}{the further restriction of STRING by interaction type. If NA, no such restriction. Otherwide, it can be one or more of "neighborhood_score","fusion_score","cooccurence_score","coexpression_score","experimental_score","database_score","textmining_score". Useful options are c("experimental_score","database_score"): only experimental data (extracted from BIND, DIP, GRID, HPRD, IntAct, MINT, and PID) and curated data (extracted from Biocarta, BioCyc, GO, KEGG, and Reactome) are used}

\item{weighted}{logical to indicate whether edge weights should be considered. By default, it sets to false. If true, it only works for the network from the STRING database}

\item{verbose}{logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display}

\item{placeholder}{the characters to tell the placeholder of RDS files. See \code{\link{oRDS}} for details}

\item{guid}{a valid (5-character) Global Unique IDentifier for an OSF project. See \code{\link{oRDS}} for details}
}
\value{
an object of class "igraph"
}
\description{
\code{oDefineNet} is supposed to define a gene network sourced from the STRING database or the Pathway Commons database. It returns an object of class "igraph".
}
\note{
The input graph will treat as an unweighted graph if there is no 'weight' edge attribute associated with
}
\examples{
\dontrun{
# STRING (high quality)
g <- oDefineNet("STRING_high", placeholder=placeholder)
# STRING (high quality), with edges weighted 
g <- oDefineNet("STRING_high", weighted=TRUE, placeholder=placeholder)
# STRING (high quality), only edges sourced from experimental or curated data
g <- oDefineNet("STRING_high", STRING.only=c("experimental_score","database_score"), placeholder=placeholder)

# Pathway Commons 
g <- oDefineNet("PCommonsDN_medium", placeholder=placeholder)

# KEGG (all)
g <- oDefineNet("KEGG", placeholder=placeholder)
# KEGG category ('Organismal Systems')
g <- oDefineNet("KEGGc_OrganismalSystem", placeholder=placeholder)
# KEGG subcategory ('Immune Systems')
g <- oDefineNet("KEGGs_ImmuneSystem", placeholder=placeholder)

# REACTOME 
g <- oDefineNet("REACTOME", placeholder=placeholder)
}
}
\seealso{
\code{\link{oRDS}}, \code{\link{oIG2TB}}, \code{\link{oTB2IG}}
}
