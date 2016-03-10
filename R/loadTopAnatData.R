########################
#' @title Retrieving Bgee data for gene set enrichment based on gene expression in anatomical structures.
#' @description The mapping from genes to anatomical structures is based on the detected presence of expression of the genes in anatomical structures. The expression comes from different data types (RNA-seq, Affymetrix microarray, ESTs, or in-situ hybridizations)
#'
#' @field species A character of species name as listed from Bgee.
#' Options are:
#'          "Anolis_carolinensis",
#'          "Bos_taurus",
#'          "Caenorhabditis_elegans",
#'          "Danio_rerio",
#'          "Drosophila_melanogaster",
#'          "Gallus_gallus",
#'          "Gorilla_gorilla",
#'          "Homo_sapiens",
#'          "Macaca_mulatta",
#'          "Monodelphis_domestica",
#'          "Mus_musculus",
#'          "Ornithorhynchus_anatinus",
#'          "Pan_paniscus",
#'          "Pan_troglodytes",
#'          "Rattus_norvegicus",
#'          "Sus_scrofa",
#'          "Xenopus_tropicalis"
#'
#' @field datatype A character of data platform.
#' Options are:
#'          "rna_seq",
#'          "affymetrix"
#'          "est",
#'          "in_situ"
#'
#' @field expressiontype A character of expression type to be used for enrichment. Only presence of expression implemented for now. Over-expression will be implemented in the future 
#' On default is "present"
#'
#' @field stage A character of developmental stage. Stages can be chosen from the metastages used in Bgee, and available at ftp://lausanne.isb-sib.ch/pub/databases/Bgee/current/stage_association.txt
#' On default is NULL: takes data from all developmental stages (BilaDO:0000001)
#'
#' @field confidence A character indicating if all data should be retrieved, or only high quality data 
#' Options are:
#'          "all",
#'          "high_quality"
#' On default is "all"
#'
#' @examples
#' \dontrun{
#'   topAnatData <- loadTopAnatData(species = "Mus_musculus", datatype = "rna_seq")
#' }
#'
#'
#' @import 
#' @export 

## TO DO: expressiontype
## Present
## over-expressed

## TO DO: implement
loadTopAnatData <- function(species, datatype, expressiontype="present", confidence="all", stage=NULL){
	url <-  "http://"

	return(list(gene2anatomy = ...))

}

