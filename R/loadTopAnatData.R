########################
#' @title Retrieve data from Bgee to perform GO-like enrichment of anatomical terms, mapped to genes by expression patterns.
#'
#' @description This function loads a mapping from genes to anatomical structures based on calls of expression in anatomical structures. It also loads the structure of the anatomical ontology. 
#' 
#' @details The expression calls come from Bgee (\url{http://bgee.org}), that integrates different expression data types (RNA-seq, Affymetrix microarray, ESTs, or in-situ hybridizations) in multiple animal species. Expression patterns are based exclusively on curated "normal", healthy, expression data (e.g., no gene knock-out, no treatment, no disease), to provide a reference of normal gene expression. Anatomical structures are identified using IDs from the Uberon ontology (browsable at \url{http://www.ontobee.org/ontology/UBERON}).
#'
#' @param species A character indicating species to be used. Options are:
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
#' @param datatype A vector of characters indicating data type(s) to be used. To be chosen among:
#'          "rna_seq",
#'          "affymetrix"
#'          "est",
#'          "in_situ".
#' Default is c("rna_seq","affymetrix","est","in_situ")
#'
#' @param calltype A character of indicating the type of expression calls to be used for enrichment. Only calls for significant presence of expression are implemented. Over-expression calls, based on differential expression analysis, will be implemented in the future. Default is "present"
#'
#' @param stage A character indicating the targeted developmental stages for the analysis. Developmental stages can be chosen from the metastage ontology used in Bgee (available at \url{ftp://lausanne.isb-sib.ch/pub/databases/Bgee/current/stages.obo}). The ID, not the name of the metastage need to be used (prefix "BilaDO:"). Default is BilaDO:0000001, the root of the metastage ontology, meaning that expression data from all developmental stages will be used.
#'
#' @param confidence A character indicating if all expression calls should be retrieved, or only high quality expression calls. Options are:
#'          "all",
#'          "high_quality".
#' Default is "all"
#' 
#' @param species.specific A character indicating if species-specific anatomical structures from the UBERON ontology should be used.
#' Default is 'TRUE': species-specific structures are used.
#'
#' @return A list of 3 elements. First, the \code{gene2anatomy} list, mapping genes to anatomical structures based on expression calls. Second, the \code{organ.names} data frame, with the name corresonding to UBERON IDs, Third, the \code{organ.relationships} list, giving the relationships between anatomical structures in the UBERON ontology (based on parent-child "is_a" and "part_of" relationships)
#'
#' @author Julien Roux \email{julien.roux at unil.ch}.
#'
#' @examples
#'   topAnatData <- loadTopAnatData(species = "Mus_musculus", datatype = "rna_seq")
#'
#' @import topGO
#' @export

loadTopAnatData <- function(species, datatype=c("rna_seq","affymetrix","est","in_situ"), calltype="present", confidence="all", stage=NULL, species.specific=TRUE){

  ## TO DO: implement

  ## Check parameters. Return error if data type not present for species? Hard code this?

  ## To know all species, use: listBgeeSpecies? Is it onl species wiht affy and RNA-seq?

  ## Build URL from options
  url <-  "http://"

  ## Launch query to topAnat server

  

  ## Format results
  ## tapply, etc

  
  return(list(gene2anatomy = ..., organ.relationships = ..., organ.names = ...))
}
