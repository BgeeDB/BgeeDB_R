########################
#' @title Retrieve data from Bgee to perform GO-like enrichment of anatomical terms, mapped to genes by expression patterns.
#'
#' @description This function loads a mapping from genes to anatomical structures based on calls of expression in anatomical structures. It also loads the structure of the anatomical ontology.
#'
#' @details The expression calls come from Bgee (\url{http://bgee.org}), that integrates different expression data types (RNA-seq, Affymetrix microarray, ESTs, or in-situ hybridizations) in multiple animal species. Expression patterns are based exclusively on curated "normal", healthy, expression data (e.g., no gene knock-out, no treatment, no disease), to provide a reference of normal gene expression.
#'
#' Anatomical structures are identified using IDs from the Uberon ontology (browsable at \url{http://www.ontobee.org/ontology/UBERON}). The mapping from genes to anatomical structures includes only the evidence of expression in these specific structures, and not the expression in their substructures (i.e., expression data are not propagated). The retrieval of propagated expression data will likely be implemented in the future, but meanwhile, it can be obtained using specialized packages such as topGO, see the \code{topAnat.R} function.
#'
#' @param species A numeric indicating the NCBI taxonomic ID of the species to be used. The species has to be among species in Bgee v13, which include:
#' \itemize{
#'   \item{6239 (Caenorhabditis elegans)}
#'   \item{7227 (Drosophila melanogaster)}
#'   \item{7955 (Danio rerio)}
#'   \item{8364 (Xenopus tropicalis)}
#'   \item{9031 (Gallus gallus)}
#'   \item{9258 (Ornithorhynchus anatinus)}
#'   \item{9544 (Macaca mulatta)}
#'   \item{9593 (Gorilla gorilla)}
#'   \item{9597 (Pan paniscus)}
#'   \item{9598 (Pan troglodytes)}
#'   \item{9606 (Homo sapiens)}
#'   \item{9823 (Sus scrofa)}
#'   \item{9913 (Bos taurus)}
#'   \item{10090 (Mus musculus)}
#'   \item{10116 (Rattus norvegicus)}
#'   \item{13616 (Monodelphis domestica)}
#'   \item{28377 (Anolis carolinensis)}
#' }
#' See the listBgeeSpecies() function to get an up-to-date list of species.
#'
#' @param datatype A vector of characters indicating data type(s) to be used. To be chosen among:
#' \itemize{
#'   \item{"rna_seq"}
#'   \item{"affymetrix"}
#'   \item{"est"}
#'   \item{"in_situ"}
#' }
#' By default all data type are included: \code{c("rna_seq","affymetrix","est","in_situ")}.
#' Including a data type that is not present in Bgee for a given species has no effect.
#'
#' @param calltype A character of indicating the type of expression calls to be used for enrichment.
#' Only calls for significant presence of expression are implemented ("expressed").
#' Over-expression calls, based on differential expression analysis, will be implemented in the future.
#'
#' @param stage A character indicating the targeted developmental stages for the analysis.
#' Developmental stages can be chosen from the developmental stage ontology used in Bgee
#' (available at \url{https://github.com/obophenotype/developmental-stage-ontologies}).
#' If a stage ID is given, the expression pattern mapped to this stage and all children
#' developmental stages (substages) will be retrieved.
#' Default is NULL, meaning that expression patterns of genes are retrieved regardless of the stage of expression.
#' This is equivalent to specifying stage="UBERON:0000104" (life cycle, the root of the stage ontology).
#' The most useful stages (going no deeper than level 3 of the ontology) include:
#' \itemize{
#'   \item{UBERON:0000068 (embryo stage)}
#'   \itemize{
#'     \item{UBERON:0000106 (zygote stage)}
#'     \item{UBERON:0000107 (cleavage stage)}
#'     \item{UBERON:0000108 (blastula stage)}
#'     \item{UBERON:0000109 (gastrula stage)}
#'     \item{UBERON:0000110 (neurula stage)}
#'     \item{UBERON:0000111 (organogenesis stage)}
#'     \item{UBERON:0007220 (late embryonic stage)}
#'     \item{UBERON:0004707 (pharyngula stage)}
#'   }
#'   \item{UBERON:0000092 (post-embryonic stage)}
#'   \itemize{
#'     \item{UBERON:0000069 (larval stage)}
#'     \item{UBERON:0000070 (pupal stage)}
#'     \item{UBERON:0000066 (fully formed stage)}
#'   }
#' }
#'
#' @param confidence A character indicating if only high quality expression calls should be retrieved.
#'  Options are "all" or "high_quality". Default is "all".
#'
#' @param host URL to Bgee webservice.
#' Change host to access development or archive versions of Bgee. Default is "\url{http://bgee.org}" to access current Bgee release.
#'
#' @param pathToData Path to the directory where the data files are stored / will be stored. Default is the working directory.
#'
#' @return A list of 3 elements:
#' \itemize{
#'   \item{A \code{gene2anatomy} list, mapping genes to anatomical structures based on expression calls.}
#'   \item{A \code{organ.names} data frame, with the name corresonding to UBERON IDs.}
#'   \item{A \code{organ.relationships} list, giving the relationships between anatomical structures
#'   in the UBERON ontology (based on parent-child "is_a" and "part_of" relationships).}
#' }
#'
#' @author Julien Roux \email{julien.roux at unil.ch}.
#'
#' @examples{
#'   myTopAnatData <- loadTopAnatData(species = "10090", datatype = "rna_seq")
#' }
#' @export

loadTopAnatData <- function(species, datatype=c("rna_seq","affymetrix","est","in_situ"), calltype="expressed",
                            confidence="all", stage=NULL, host="http://bgee.org", pathToData=getwd()){
  allSpecies <- c(6239, 7227, 7955, 8364, 9031, 9258, 9544, 9593, 9597, 9598, 9606, 9823, 9913, 10090, 10116, 13616, 28377)
  ## Species is the only compulsory parameter
  if( length(species) == 0 ) {
    stop("Problem: you need to specify a species.")
  } else if ( length(species) > 1 ){
    stop("Problem: only one species is allowed.")
  } else if ( sum(species %in% allSpecies) == 0 ){
    stop("Problem: the specified speciesId is not among the list of species in Bgee.")
  }

  ## Test if parameters are in the range of allowed parameters
  if ( !sum(datatype %in% c("rna_seq","affymetrix","est","in_situ")) %in% 1:4 ){
    stop("Problem: you need to specify at least one valid data type to be used among \"rna_seq\", \"affymetrix\", \"est\" and \"in_situ\".")
  }
  if ( length(datatype) != sum(datatype %in% c("rna_seq","affymetrix","est","in_situ")) ){
    cat("Warning: you apparently specified a data type that is not among \"rna_seq\", \"affymetrix\", \"est\" and \"in_situ\". Please check for typos.\n")
  }
  if ( calltype != "expressed" ){
    stop("Problem: no other call types than present / absent can be retrieved for now.")
  }
  if ( (confidence != "all") && (confidence != "high_quality") ){
    stop("Problem: the data confidence parameter specified is not among the allowed values (\"all\" or \"high_quality\").")
  }
  if ( !grepl("^http://", host) && !grepl("^https://", host) ){
    host <- paste0("http://", host)
  }
  if ( !grepl("/$", host) ){
    host <- paste0(host, "/")
  }
  if ( !file.exists(pathToData) ){
    stop("Problem: please specify a valid path to data files.")
  }
  if ( !grepl("/$", pathToData) ){
    pathToData <- paste0(pathToData, "/")
  }

  ## Set the internet.info to 2 to have less verbose output (only reports critical warnings)
  options(internet.info=2)
  ## Set the timeout option to 600 seconds to let some time to the server to send data (default is 60s)
  options(timeout = 600)
  
  ## First query: organ relationships
  organRelationshipsFileName <- paste0("topAnat_AnatEntitiesRelationships_", species, ".tsv")
  ## Check if file is already in cache
  if (file.exists(paste0(pathToData, organRelationshipsFileName))){
    cat("\nWarning: an organ relationships file was found in the working directory and will be used as is. Please delete and rerun the function if you want the data to be updated.\n")
  } else {
    cat("\nBuilding URLs to retrieve organ relationships from Bgee.........\n")
    myurl <-  paste0(host, "?page=dao&action=org.bgee.model.dao.api.ontologycommon.RelationDAO.getAnatEntityRelations&display_type=tsv&species_list=", species,"&attr_list=SOURCE_ID&attr_list=TARGET_ID")

    ## Query webservice
    cat(paste0("   URL successfully built (", myurl,")\n   Submitting URL to Bgee webservice (can be long)\n  "))
    download.file(myurl, destfile = paste0(pathToData, organRelationshipsFileName, ".tmp"))

    ## Read 5 last lines of file: should be empty indicating success of data transmission
    ## We cannot use a system call to UNIX command since some user might be on Windows
    tmp <- tail(read.table(paste0(pathToData, organRelationshipsFileName, ".tmp"), header=TRUE, sep="\t", comment.char="", blank.lines.skip=FALSE, as.is=TRUE), n=5)
    if ( length(tmp[,1]) == 5 && (sum(tmp[,1] == "") == 5 || sum(is.na(tmp[,1])) == 5) ){
      ## The file transfer was successful, we rename the temporary file
      file.rename(paste0(pathToData, organRelationshipsFileName, ".tmp"), paste0(pathToData, organRelationshipsFileName))
    } else {
      ## delete the temporary file
      file.remove(paste0(pathToData, organRelationshipsFileName, ".tmp"))
      stop(paste0("File ", organRelationshipsFileName, " is truncated, there may be a temporary problem with the Bgee webservice, or there was an error in the parameters."))
    }
    cat(paste0("   Got results from Bgee webservice. Files are written in \"", pathToData, "\"\n"))
  }

  ## Second query: organ names
  organNamesFileName <- paste0("topAnat_AnatEntitiesNames_", species, ".tsv");
  ## Check if file is already in cache
  if (file.exists(paste0(pathToData, organNamesFileName))){
    cat("\nWarning: an organ names file was found in the working directory and will be used as is. Please delete and rerun the function if you want the data to be updated.\n")
  } else {
    cat("\nBuilding URLs to retrieve organ names from Bgee.................\n")
    myurl <-  paste0(host, "?page=dao&action=org.bgee.model.dao.api.anatdev.AnatEntityDAO.getAnatEntities&display_type=tsv&species_list=", species,"&attr_list=ID&attr_list=NAME")

    ## Query webservice
    cat(paste0("   URL successfully built (", myurl,")\n   Submitting URL to Bgee webservice (can be long)\n  "))
    download.file(myurl, destfile = paste0(pathToData, organNamesFileName, ".tmp"))

    ## Read 5 last lines of file: should be empty indicating success of data transmission
    ## We cannot use a system call to UNIX command since some user might be on Windows
    tmp <- tail(read.table(paste0(pathToData, organNamesFileName, ".tmp"), header=TRUE, sep="\t", comment.char="", blank.lines.skip=FALSE, as.is=TRUE, quote = ""), n=5)
    if ( length(tmp[,1]) == 5 && (sum(tmp[,1] == "") == 5 || sum(is.na(tmp[,1])) == 5) ){
      ## The file transfer was successful, we rename the temporary file
      file.rename(paste0(pathToData, organNamesFileName, ".tmp"), paste0(pathToData, organNamesFileName))
    } else {
      ## delete the temporary file
      file.remove(paste0(pathToData, organNamesFileName, ".tmp"))
      stop(paste0("File ", organNamesFileName, " is truncated, there may be a temporary problem with the Bgee webservice, or there was an error in the parameters."))
    }
    cat(paste0("   Got results from Bgee webservice. Files are written in \"", pathToData, "\"\n"))
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
  if (file.exists(paste0(pathToData, gene2anatomyFileName))){
    cat("\nWarning: a gene to organs mapping file was found in the working directory and will be used as is. Please delete and rerun the function if you want the data to be updated.\n")
  } else {
    cat("\nBuilding URLs to retrieve mapping of gene to organs from Bgee...\n")
    myurl <-  paste0(host, "?page=dao&action=org.bgee.model.dao.api.expressiondata.ExpressionCallDAO.getExpressionCalls&display_type=tsv&species_list=", species, "&attr_list=GENE_ID&attr_list=ANAT_ENTITY_ID")

    ## Add data type: only if not all data types needed
    if ( sum(datatype %in% c("rna_seq","affymetrix","est","in_situ")) < 4 ){
      for (type in toupper(sort(datatype))){
        myurl <- paste0(myurl, "&data_type=", type)
      }
    }
    ## Add data quality
    if ( confidence == "high_quality" ){
      myurl <- paste0(myurl, "&data_qual=HIGH")
    }
    if ( !is.null(stage) ){
      myurl <- paste0(myurl, "&stage_id=", stage)
    }

    ## Query webservice
    cat(paste0("   URL successfully built (", myurl,")\n   Submitting URL to Bgee webservice (can be long)\n  "))
    download.file(myurl, destfile = paste0(pathToData, gene2anatomyFileName, ".tmp"))

    ## Read 5 last lines of file: should be empty indicating success of data transmission
    ## We cannot use a system call to UNIX command since some user might be on Windows
    tmp <- tail(read.table(paste0(pathToData, gene2anatomyFileName, ".tmp"), header=TRUE, sep="\t", comment.char="", blank.lines.skip=FALSE, as.is=TRUE), n=5)
    if ( length(tmp[,1]) == 5 && (sum(tmp[,1] == "") == 5 || sum(is.na(tmp[,1])) == 5) ){
      ## The file transfer was successful, we rename the temporary file
      file.rename(paste0(pathToData, gene2anatomyFileName, ".tmp"), paste0(pathToData, gene2anatomyFileName))
    } else {
      ## delete the temporary file
      file.remove(paste0(pathToData, gene2anatomyFileName, ".tmp"))
      stop(paste0("File ", gene2anatomyFileName, " is truncated, there may be a temporary problem with the Bgee webservice, or there was an error in the parameters."))
    }
    cat(paste0("   Got results from Bgee webservice. Files are written in \"", pathToData, "\"\n"))
  }

  ## Process the data and build the final list to return
  cat("\nParsing the results...............................................\n")

  ## Relationships between organs
  if (file.exists(paste0(pathToData, organRelationshipsFileName))){
    if (file.info(paste0(pathToData, organRelationshipsFileName))$size != 0) {
      tab <- read.table(paste0(pathToData, organRelationshipsFileName), header=TRUE, sep="\t", blank.lines.skip=TRUE, as.is=TRUE)
      organRelationships <- tapply(as.character(tab[,2]), as.character(tab[,1]), unique)
    } else {
      stop(paste0("File ", organRelationshipsFileName, " is empty, there may be a temporary problem with the Bgee webservice, or there was an error in the parameters."))
    }
  } else {
    stop(paste0("File ", organRelationshipsFileName, " not found. There may be a temporary problem with the Bgee webservice, or there was an error in the parameters."))
  }
  ## Organ names
  if (file.exists(paste0(pathToData, organNamesFileName))){
    if (file.info(paste0(pathToData, organNamesFileName))$size != 0) {
      organNames <- read.table(paste0(pathToData, organNamesFileName), header=TRUE, sep="\t", comment.char="", blank.lines.skip=TRUE, as.is=TRUE, quote = "")
      names(organNames) <- c("organId", "organName")
    } else {
      stop(paste0("File ", organNamesFileName, " is empty, there may be a temporary problem with the Bgee webservice, or there was an error in the parameters."))
    }
  } else {
    stop(paste0("File ", organNamesFileName, " not found. There may be a temporary problem with the Bgee webservice, or there was an error in the parameters."))
  }
  ## Mapping of genes to tissues
  if (file.exists(paste0(pathToData, gene2anatomyFileName))){
    if (file.info(paste0(pathToData, gene2anatomyFileName))$size != 0) {
      tab <- read.table(paste0(pathToData, gene2anatomyFileName), header=TRUE, sep="\t", blank.lines.skip=TRUE, as.is=TRUE)
      if(length(tab[,1]) != 0){
        gene2anatomy <- tapply(as.character(tab[,2]), as.character(tab[,1]), unique)
      } else {
        stop(paste0("There was no mapping of genes to anatomical structures found. Probably the parameters are too stringent, or this data type is absent in this species."))
      }
    } else {
      stop(paste0("File ", gene2anatomyFileName, " is empty, there may be a temporary problem with the Bgee webservice, or there was an error in the parameters. It is also possible that the parameters are too stringent and returned no data, please try to relax them."))
    }
  } else {
    stop(paste0("File ", gene2anatomyFileName, " not found. There may be a temporary problem with the Bgee webservice, or there was an error in the parameters. It is also possible that the parameters are too stringent and returned no data, please try to relax them."))
  }

  cat("\nAdding BGEE:0 as unique root of all terms of the ontology.........\n")
  ## There can be multiple roots among all the terms downloaded. We need to add one unique root for topGO to work: BGEE:0
  ## Add all organs from organNames that are not source (child / names of the list) in organsRelationship to the organsRelationship list (with value / target / parent = BGEE:0)
  missingParents <- organNames[!organNames[,1] %in% names(organRelationships), 1]
  ## Add new values
  organRelationships <- c(organRelationships, as.list(rep("BGEE:0", times=length(missingParents))))
  ## Add new keys
  names(organRelationships)[(length(organRelationships)-length(missingParents)+1):length(organRelationships)] = as.character(missingParents)
  ## Add BGEE:0	/ root to organNames
  organNames <- rbind(organNames, c("BGEE:0", "root"))

  cat("\nDone.\n")
  return(list(gene2anatomy = gene2anatomy, organ.relationships = organRelationships, organ.names = organNames))
}
