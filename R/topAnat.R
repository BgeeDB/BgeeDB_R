########################
#' @title Load function allowing to perform GO-like enrichment of anatomical terms using the topGO package
#'
#' @description This function loads the functions necessary to use topGO on the Uberon ontology instead of the GO ontology
#'
#' @details To perform the enrichment test for expression in anatomical structures for each term of Uberon ontology (browsable at \url{http://www.ontobee.org/ontology/UBERON}), it is interesting to use the topGO package since it allows to propagate the mapping of gene to terms to parent terms, and it possesses a pannel of enrichment test and decorrelation methods
#'
#' @param ...
#'
#' @return ...
#'
#' @author Julien Roux \email{julien.roux@unil.ch}.
#'
#' @examples
#'   topAnat()
#'
#' @export

topAnat <- function(topAnatData){

  ## TO DO: implement
  ## Do we even need to load topAnatData? Maybe just check that not empty?

  ## Load topOBO functions

  return()
}

## TO DO: why topGO not in DESCRIPTION file?
##        warning when loading topGO: problem?
## TO DO: topAnat_functions.R: add functions and make them start with a . to stay hidden
