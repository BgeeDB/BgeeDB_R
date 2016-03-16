########################
#' @title Retrieve data from Bgee to perform GO-like enrichment of anatomical terms, mapped to genes by expression patterns.
#'
#' @description This function loads a mapping from genes to anatomical structures based on calls of expression in anatomical structures. It also loads the structure of the anatomical ontology.
#'
#' @details The expression calls come from Bgee (\url{http://bgee.org}), that integrates different expression data types (RNA-seq, Affymetrix microarray, ESTs, or in-situ hybridizations) in multiple animal species. Expression patterns are based exclusively on curated "normal", healthy, expression data (e.g., no gene knock-out, no treatment, no disease), to provide a reference of normal gene expression. Anatomical structures are identified using IDs from the Uberon ontology (browsable at \url{http://www.ontobee.org/ontology/UBERON}). The mapping from genes to anatomical structures includes only the evidence of expression in these specific structures, and not the expression in their substructures (i.e., expression data are not propagated). The retrieval of prpagated expression data will likely be implemented in the future, but meanwhile, it can be obtained using specialized packages like topGO, see the topAnat.R function.
#'
#' @param host URL to Bgee webservice. Chnage host to access development or archive versions of Bgee. Default is "http://bgee.org" to access current Bgee release.
#'
#' @param species A numeric indicating the NCBI taxonomic ID of the species to be used among species in Bgee. Species present in bgee v13 include: 
#'  6239 (Caenorhabditis elegans)
#'  7227 (Drosophila melanogaster)
#'  7955 (Danio rerio)
#'  8364 (Xenopus tropicalis)
#'  9031 (Gallus gallus)
#'  9258 (Ornithorhynchus anatinus)
#'  9544 (Macaca mulatta)
#'  9593 (Gorilla gorilla)
#'  9597 (Pan paniscus)
#'  9598 (Pan troglodytes)
#'  9606 (Homo sapiens)
#'  9823 (Sus scrofa)
#'  9913 (Bos taurus)
#' 10090 (Mus musculus)
#' 10116 (Rattus norvegicus)
#' 13616 (Monodelphis domestica)
#' 28377 (Anolis carolinensis)
## TO DO: how to display a bullet point list?
#'
#' @param datatype A vector of characters indicating data type(s) to be used. To be chosen among:
#'          "rna_seq",
#'          "affymetrix"
#'          "est",
#'          "in_situ".
#' Default includes all data type for a given species: c("rna_seq","affymetrix","est","in_situ")
#'
#' @param calltype A character of indicating the type of expression calls to be used for enrichment. Only calls for significant presence of expression are implemented ("expressed"). Over-expression calls, based on differential expression analysis, will be implemented in the future. Default is "expressed"
#'
#' @param stage A character indicating the targeted developmental stages for the analysis. Developmental stages can be chosen from the developmental stage ontology used in Bgee (available at \url{TO DO}). If a stage ID is given, the expression pattern mapped to this stage and all children developmental stages (substages) will be retrieved. Default is NULL, meaning that expression patterns of genes are retrieved regardless of the stage of expression. Potential examples of useful stage specifications include:
#' UBERON:00000.. life cycle ## useless, since will give the same results as the default stage=NULL
#'   UBERON:0000068 embryonic stage
#'     UBERON...
#'     UBERON...
#'     UBERON...
#'   UBERON:0000092 post-embryonic stage
#'     UBERON...
#'     UBERON...
#'     UBERON...
## TO DO: how to display a bullet point list?
#'
#' @param confidence A character indicating if only high quality expression calls should be retrieved. Options are "all" or "high_quality". Default is "all".
#'
#' @return A list of 3 elements. First, a \code{gene2anatomy} list, mapping genes to anatomical structures based on expression calls. Second, a \code{organ.names} data frame, with the name corresonding to UBERON IDs, Third, a \code{organ.relationships} list, giving the relationships between anatomical structures in the UBERON ontology (based on parent-child "is_a" and "part_of" relationships)
#'
#' @author Julien Roux \email{julien.roux at unil.ch}.
#'
#' @examples
#' \dontrun{
#'   myTopAnatData <- loadTopAnatData(species = "Mus_musculus", datatype = "rna_seq")
#' }
#' @export

loadTopAnatData <- function(species, datatype=c("rna_seq","affymetrix","est","in_situ"), calltype="expressed", confidence="all", stage=NULL){
  allSpecies <- c(6239, 7227, 7955, 8364, 9031, 9258, 9544, 9593, 9597, 9598, 9606, 9823, 9913, 10090, 10116, 13616, 28377)
  ## Species is the only compulsory parameter
  if( length(species) == 0 ) {
    stop("Problem: you need to specify a species.")
  } else if ( length(species) > 1 ){
    stop("Problem: only one species is allowed.")
  } else if ( sum(species %in% allSpecies) == 0 ){
    stop("Problem: the specified speciesId is not amongthe list of species in Bgee.")
  }

  ## Test if parameters are in the range of allowed parameters
  if ( !sum(datatype %in% c("rna_seq","affymetrix","est","in_situ")) %in% 1:4 ){
    stop("Problem: you need to specify at least one valid data type to be used among rna_seq, affymetrix, est and in_situ.")
  }
  if ( length(datatype) != sum(datatype %in% c("rna_seq","affymetrix","est","in_situ")) ){
    warning("Warning: you apparently specified a data type that is not among rna_seq, affymetrix, est and in_situ. Please check for typos.")
  }
  if (calltype != "expressed"){
    stop("Problem: no other call types than present / absent can be retrieved for now.")
  }
  if ((confidence != "all") and (confidence != "high_quality")){
    stop("Problem: the data confidence parameter specified is not among the allowed values (all or high_quality).")
  }

  
  ## First query: organ relationships
  organRelationshipsFileName <- paste0("topAnat_AnatEntitiesRelationships_", species, ".tsv")
  ## Check if file is already in cache
  if (file.exists(organRelationshipsFileName)){
    cat("\nThe organ relationships file is already in the working directory and will be used as is. Please delete and rerun the function if you want the data to be updated.\n")
  } else {
    cat("\nBuilding URLs to retrieve organ relationships from Bgee...\n")

    ## TO DO: build URL + query webservice
    ## save organRelationshipsFileName in current working directory 
  }

  ## Second query: organ names
  organNamesFileName <- paste0("topAnat_AnatEntitiesNames_", species, ".tsv");
  ## Check if file is already in cache
  if (file.exists(organNamesFileName)){
    cat("\nThe organ names file is already in the working directory and will be used as is. Please delete and rerun the function if you want the data to be updated.\n")
  } else {
    cat("\nBuilding URLs to retrieve organ names from Bgee...\n")

    ## TO DO: build URL + query webservice
    ## save organNamesFileName in current working directory 
  }

  ## Third query: gene to organs mapping
  gene2anatomyFileName <- paste0("topAnat_GeneToAnatEntities_", species, "_", toupper(calltype))
  ## If a stage is specified, add it to file name
  if ( !is.null(stage) ){
    gene2anatomyFileName <- paste0(gene2anatomyFileName, "_", gsub(":", "_", stage))
  }
  ## If all datatypes specified, not need to add anything to file name. Otherwise, specify datatypes in file name
  if ( sum(datatype %in% c("rna_seq","affymetrix","est","in_situ")) < 4 ){
    gene2anatomyFileName <- paste0(gene2anatomyFileName, "_", toupper(paste(sort(datatype), collapse="_")))
  }
  ## If high quality data needed, specify in file name. Otherwise not specified
  if ( confidence == "high_quality" ){
    gene2anatomyFileName <- paste0(gene2anatomyFileName, "_HIGH")
  }
  gene2anatomyFileName <- paste0(gene2anatomyFileName, ".tsv")

  ## Check if file is already in cache
  if (file.exists(gene2anatomyFileName)){
    cat("\nThe gene to organs mapping file is already in the working directory and will be used as is. Please delete and rerun the function if you want the data to be updated.\n")
  } else {
    cat("\nBuilding URLs to retrieve mapping of gene to organs from Bgee...\n")

    ## TO DO: build URL + query webservice
    ## save gene2anatomyFileName in current working directory 
  }

  ## TO DO: maybe an argument could be path to the data / path to store the data?

  ## TO DO: URL building:
  ## Build invariable part of the URL: host and species
  myurl <-  paste0(host, "/?page=top_anat&action=...&species=", species)
  ## Add data type
  for (type in toupper(datatype)){
    myurl <- paste0(myurl, "&data_type=", type)
  }
  ## Add data quality
  if ( confidence == "high_quality" ){
    myurl <- paste0(myurl, "&data_qual=", confidence)
  }

  ## Add call type
  if (calltype == "expressed"){
    myurl <- paste0(myurl, "&expr_type=EXPRESSED")
  }
  ## Add developmental stage
  myurl <- paste0(myurl, "&stage_id=", stage)
  cat(paste0("   URL successfully built (", myurl,")\n"))

  ## TO DO? Check parameters. Return error if data type not present for species? Hard code this?

  cat("   Submitting URL to Bgee webservice (can be long)...\n")
  ## TO DO: Launch query to topAnat server
  ## How to do this? Is it like downloading a file?
  download.file(file.path(myurl), destfile = getwd())
  ## TO DO: add fileName here?

  cat(paste0("   Got answer from Bgee webservice. Result files are written in \"", getwd(), "\"\n"))


  ## Process the data and build the final list to return
  cat("\nParsing the results... ")

  ## Relationships between organs
  tab <- read.table(organRelationshipsFileName, header=FALSE, sep="\t")
  organRelationships <- tapply(as.character(tab[,2]), as.character(tab[,1]), unique)

  ## Organ names
  organNames <- read.table(organNamesFileName, header = FALSE, sep="\t", row.names=1)
  names(organNames) <- organNames

  ## Mapping of genes to tissues
  if (file.info(gene2anatomyFileName)$size != 0) {
    tab <- read.table(gene2anatomyFileName, header=FALSE, sep="\t")
    gene2anatomy <- tapply(as.character(tab[,2]), as.character(tab[,1]), unique)
  } else {
    ## TO DO: issue warning. Actually, do it for all files!
  }

  cat("Done.\n")
  return(list(gene2anatomy = gene2anatomy, organ.relationships = organRelationships, organ.names = organNames))
}
