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
#' @param calltype A character of indicating the type of expression calls to be used for enrichment. Only calls for significant presence of expression are implemented ("expressed"). Over-expression calls, based on differential expression analysis, will be implemented in the future. Default is "expressed"
#'
#' @param stage A character indicating the targeted developmental stages for the analysis. Developmental stages can be chosen from the metastage ontology used in Bgee (available at \url{ftp://lausanne.isb-sib.ch/pub/databases/Bgee/current/stages.obo}). The ID, not the name of the metastage need to be used (prefix "BilaDO:"). Default is BilaDO:0000001, the root of the metastage ontology, meaning that expression data from all developmental stages will be used.
#' TO DO: should we use UBERON stages??
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
#'   myTopAnatData <- loadTopAnatData(species = "Mus_musculus", datatype = "rna_seq")
#'
#' @export

loadTopAnatData <- function(species, datatype=c("rna_seq","affymetrix","est","in_situ"), calltype="expressed", confidence="all", stage=NULL, species.specific=TRUE){
  if(length(species)==0) {
    stop("Problem: you did not specify a species")
  }

  cat("\nBuilding URL to retrieve data from Bgee...\n")
  ## TO DO: Build URL from options
  myurl <-  "http://http://bgee.unil.ch/?page=top_anat&action=...&...=..."
  ## Add species
  myurl <- paste0(myurl, "&species=", species)
  ## Add data type
  for (type in toupper(datatype)){
    myurl <- paste0(myurl, "&data_type=", type)
  }
  ## Add data quality
  myurl <- paste0(myurl, "&data_qual=", confidence)
  ## TO DO: works with "all" but not with "high_quality"

  ## Add call type
  if (calltype == "expressed"){
    myurl <- paste0(myurl, "&expr_type=EXPRESSED")
  }
  ## Add developmental stage
  myurl <- paste0(myurl, "&stage_id=", stage)
  cat(paste0("   URL successfully built (", myurl,")\n"))
 
  ## TO DO? Check parameters. Return error if data type not present for species? Hard code this?

  cat("\nSubmitting URL to Bgee webservice...\n")
  ## TO DO: Launch query to topAnat server
  ## How to do this? Is it like downloading a file?
  ## download.file(file.path(myurl), destfile = getwd())

  cat(paste0("   Got answer from Bgee webservice. Result files are written in \"", getwd(), "\"\n"))
  cat("\nParsing the results... ")

  
  ## TO DO: read and format results

  ##   ## If data are zipepd
  ##   temp <- list.files(pattern="*.tsv.zip")
  ##   mydata <- lapply(temp2, unzip)
  ##   mydata_all <- lapply(mydata, fread)
  ##   ## Use data.table=F to get data frame?

  ## Relationships between tissues
  organRelationshipsFileName <- paste0("topAnat_AnatEntitiesRelationships_", species, ".tsv")
  tab <- read.table(organRelationshipsFileName, header=FALSE, sep="\t")
  organRelationships <- tapply(as.character(tab[,2]), as.character(tab[,1]), unique)

  ## Tissue names
  organNamesFileName <- paste0("topAnat_AnatEntitiesNames_", species, ".tsv");
  organNames <- read.table(organNamesFileName, header = FALSE, sep="\t", row.names=1)
  names(organNames) <- organNames

  ## Mapping of genes to tissues  
  gene2anatomyFileName <- paste0("topAnat_GeneToAnatEntities_", species, "_", toupper(calltype), "_", gsub(":", "_", stage), "_", toupper(paste(sort(datatype), collapse="_")), "_", toupper(confidence), ".tsv")
  if (file.info(gene2anatomyFileName)$size != 0) {
    tab <- read.table(gene2anatomyFileName, header=FALSE, sep="\t")
    gene2anatomy <- tapply(as.character(tab[,2]), as.character(tab[,1]), unique)
  }

  ## TO DO: filter out species-specific terms (non UBERON only?)
  if (species.specific == FALSE){
    
  }

  cat("Done.\n")
  return(list(gene2anatomy = gene2anatomy, organ.relationships = organRelationships, organ.names = organNames))
}
