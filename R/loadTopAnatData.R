#' @title Retrieve data from Bgee to perform GO-like enrichment of anatomical terms, mapped to genes by expression patterns.
#'
#' @description This function loads a mapping from genes to anatomical structures based on calls of expression in anatomical structures. It also loads the structure of the anatomical ontology.
#'
#' @details The expression calls come from Bgee (\url{http://bgee.org}), that integrates different expression data types (RNA-seq, Affymetrix microarray, ESTs, or in-situ hybridizations) from multiple animal species. Expression patterns are based exclusively on curated "normal", healthy, expression data (e.g., no gene knock-out, no treatment, no disease), to provide a reference atlas of normal gene expression. Anatomical structures are identified using IDs from the Uberon ontology (browsable at \url{http://www.ontobee.org/ontology/UBERON}). The mapping from genes to anatomical structures includes only the evidence of expression in these specific structures, and not the expression in their substructures (i.e., expression data are not propagated). The retrieval of propagated expression data might be implemented in the future, but meanwhile, it can be obtained using specialized packages such as topGO, see the \code{topAnat.R} function.
#'
#' @param myBgeeObject An output object from Bgee$new().
#'
#' @param callType A character of indicating the type of expression calls to be used for enrichment. Only calls for significant detection of expression are implemented so far ("presence"). Differential expression calls, based on differential expression analysis, might be implemented in the future.
#'
#' @param stage A character indicating the targeted developmental stages for the analysis. Developmental stages can be chosen from the developmental stage ontology used in Bgee (available at \url{https://github.com/obophenotype/developmental-stage-ontologies}). If a stage is specified, the expression pattern mapped to this stage and all children developmental stages (substages) will be retrieved. Default is NULL, meaning that expression patterns of genes are retrieved regardless of the developmental stage displaying expression; this is equivalent to specifying stage="UBERON:0000104" (life cycle, the root of the stage ontology).
#' For information, the most useful stages (going no deeper than level 3 of the ontology) include:
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
#' @param confidence A character indicating if only high quality present calls should be retrieved. For Bgee releases prior to 14, options are "all" (default) or "high_quality". For Bgee release 14 and above, options are "silver" (default) and "gold".
#'
#' @return A list of 4 elements:
#' \itemize{
#'   \item{A \code{gene2anatomy} list, mapping genes to anatomical structures based on expression calls.}
#'   \item{A \code{organ.names} data frame, with the name corresponding to UBERON IDs.}
#'   \item{A \code{organ.relationships} list, giving the relationships between anatomical structures in the UBERON ontology (based on parent-child "is_a" and "part_of" relationships).}
#'   \item{The Bgee class object thta was used to retrieve the data.}
#' }
#'
#' @author Julien Roux, Julien Wollbrett
#'
#' @examples{
#'   bgee <- Bgee$new(species = "Bos_taurus", dataType = "rna_seq")
#'   myTopAnatData <- loadTopAnatData(bgee)
#' }
#'
#' @import utils digest
#' @export

loadTopAnatData <- function(myBgeeObject, callType="presence", confidence=NULL, stage=NULL){
  OLD_WEBSERVICE_VERSION = '13.2'

  ## check that fields of Bgee object are not empty
  if (length(myBgeeObject$speciesId) == 0 | length(myBgeeObject$topAnatUrl) == 0 | length(myBgeeObject$dataType) == 0 | length(myBgeeObject$pathToData) == 0 | length(myBgeeObject$sendStats) == 0){
    stop("ERROR: there seems to be a problem with the input Bgee class object, some fields are empty. Please check that the object is valid.")
  }
  if ( callType != "presence" ){
    stop("ERROR: no other call types than \"presence\" expression calls can be retrieved for now.")
  }
  if ( compareVersion(gsub("_", ".", myBgeeObject$release), OLD_WEBSERVICE_VERSION) > 0 ){
    if ( is.null(confidence) ){
      confidence = "silver"
    }
    if ( (confidence != "gold") && (confidence != "silver") ){
      stop(paste0("ERROR: the data confidence parameter specified is not among the allowed values. For Bgee ", myBgeeObject$release, " allowed values are \"silver\" or \"gold\".\nBy default \"silver\" quality is selected."))
    }
  } else {
    if(is.null(confidence)){
      confidence = "all"
    }
    if ( (confidence != "all") && (confidence != "high_quality") ){
      stop(paste0("ERROR: the data confidence parameter specified is not among the allowed values. For Bgee ", myBgeeObject$release, " allowed values are \"all\" or \"high_quality\".\nBy default \"all\" quality is selected."))
    }
  }

  ## Set the internet.info to 2 to have less verbose output (only reports critical warnings)
  options(internet.info = 2)
  ## Set the timeout option to 600 seconds to let some time to the server to send data (default is 60s)
  options(timeout = 600)


  ## First query: organ relationships
  organRelationshipsFileName <- paste0("topAnat_AnatEntitiesRelationships_", myBgeeObject$speciesId, ".tsv")
  ## Check if file is already in cache
  if (file.exists(file.path(myBgeeObject$pathToData, organRelationshipsFileName))){
    cat(paste0("\nNOTE: an organ relationships file was found in the download directory ", myBgeeObject$pathToData,
        ". Data will not be redownloaded.\n"))
  } else {
    cat("\nBuilding URLs to retrieve organ relationships from Bgee.........\n")
    myUrl <- myBgeeObject$topAnatUrl
    if(compareVersion(gsub("_", ".", myBgeeObject$release), OLD_WEBSERVICE_VERSION) > 0){
      myUrl <- paste0(myUrl, "?page=r_package&action=get_anat_entity_relations&display_type=tsv&species_list=", myBgeeObject$speciesId, "&attr_list=SOURCE_ID&attr_list=TARGET_ID&api_key=", myBgeeObject$apiKey, "&source=BgeeDB_R_package&source_version=", as.character(packageVersion("BgeeDB")))
    } else {
      myUrl <- paste0(myUrl, "?page=dao&action=org.bgee.model.dao.api.ontologycommon.RelationDAO.getAnatEntityRelations&display_type=tsv&species_list=", myBgeeObject$speciesId, "&attr_list=SOURCE_ID&attr_list=TARGET_ID&api_key=", myBgeeObject$apiKey, "&source=BgeeDB_R_package&source_version=", as.character(packageVersion("BgeeDB")))
    }
    ## Query webservice
    cat(paste0("   URL successfully built (", myUrl,")\n   Submitting URL to Bgee webservice (can be long)\n"))
    success <- download.file(myUrl, destfile = paste0(myBgeeObject$pathToData, "/", organRelationshipsFileName, ".tmp"))

    if (success == 0){
      ## Read 5 last lines of file: should be empty indicating success of data transmission
      ## We cannot use a system call to UNIX command since some user might be on Windows
      tmp <- tail(read.table(paste0(myBgeeObject$pathToData, "/", organRelationshipsFileName, ".tmp"), header=TRUE, sep="\t", comment.char="", blank.lines.skip=FALSE, as.is=TRUE), n=5)
      if ( length(tmp[,1]) == 5 && (sum(tmp[,1] == "") == 5 || sum(is.na(tmp[,1])) == 5) ){
        ## The file transfer was successful, we rename the temporary file
        file.rename(paste0(myBgeeObject$pathToData, "/", organRelationshipsFileName, ".tmp"), paste0(myBgeeObject$pathToData, "/", organRelationshipsFileName))
      } else {
        ## delete the temporary file
        file.remove(paste0(myBgeeObject$pathToData, "/", organRelationshipsFileName, ".tmp"))
        stop(paste0("File ", organRelationshipsFileName, " is truncated, there may be a temporary problem with the Bgee webservice, or there was an error in the parameters."))
      }
      cat(paste0("   Got results from Bgee webservice. Files are written in \"", myBgeeObject$pathToData, "\"\n"))
    } else {
      serverAnswer = try(getURL(myUrl))
      if (class(serverAnswer) == "try-error"){
        stop("ERROR: the query to the server was not successful. Is your internet connection working?\n")
      } else {
        stop(paste0("ERROR: the query to the server was not successful. The server returned the following answer:\n", serverAnswer))
      }
    }
  }

  ## Second query: organ names
  organNamesFileName <- paste0("topAnat_AnatEntitiesNames_", myBgeeObject$speciesId, ".tsv");
  ## Check if file is already in cache
  if (file.exists(file.path(myBgeeObject$pathToData, organNamesFileName))){
    cat(paste0("\nNOTE: an organ names file was found in the download directory ", myBgeeObject$pathToData,
               ". Data will not be redownloaded.\n"))

  } else {
    cat("\nBuilding URLs to retrieve organ names from Bgee.................\n")
    myUrl <- myBgeeObject$topAnatUrl
    if(compareVersion(gsub("_", ".", myBgeeObject$release), OLD_WEBSERVICE_VERSION) > 0){
      myUrl <- paste0(myUrl, "?page=r_package&action=get_anat_entities&display_type=tsv&species_list=", myBgeeObject$speciesId, "&attr_list=ID&attr_list=NAME&api_key=", myBgeeObject$apiKey, "&source=BgeeDB_R_package&source_version=", as.character(packageVersion("BgeeDB")))
    }else {
      myUrl <- paste0(myUrl, "?page=dao&action=org.bgee.model.dao.api.anatdev.AnatEntityDAO.getAnatEntities&display_type=tsv&species_list=", myBgeeObject$speciesId, "&attr_list=ID&attr_list=NAME&api_key=", myBgeeObject$apiKey, "&source=BgeeDB_R_package&source_version=", as.character(packageVersion("BgeeDB")))
    }

    ## Query webservice
    cat(paste0("   URL successfully built (", myUrl,")\n   Submitting URL to Bgee webservice (can be long)\n"))
    success <- download.file(myUrl, destfile = paste0(myBgeeObject$pathToData, "/", organNamesFileName, ".tmp"))

    if (success == 0){
      ## Read 5 last lines of file: should be empty indicating success of data transmission
      ## We cannot use a system call to UNIX command since some user might be on Windows
      tmp <- tail(read.table(paste0(myBgeeObject$pathToData, "/", organNamesFileName, ".tmp"), header=TRUE, sep="\t", comment.char="", blank.lines.skip=FALSE, as.is=TRUE, quote = ""), n=5)
      if ( length(tmp[,1]) == 5 && (sum(tmp[,1] == "") == 5 || sum(is.na(tmp[,1])) == 5) ){
        ## The file transfer was successful, we rename the temporary file
        file.rename(paste0(myBgeeObject$pathToData, "/", organNamesFileName, ".tmp"), paste0(myBgeeObject$pathToData, "/", organNamesFileName))
      } else {
        ## delete the temporary file
        file.remove(paste0(myBgeeObject$pathToData, "/", organNamesFileName, ".tmp"))
        stop(paste0("File ", organNamesFileName, " is truncated, there may be a temporary problem with the Bgee webservice, or there was an error in the parameters."))
      }
     cat(paste0("   Got results from Bgee webservice. Files are written in \"", myBgeeObject$pathToData, "\"\n"))
    } else {
      serverAnswer = try(getURL(myUrl))
      if (class(serverAnswer) == "try-error"){
        stop("ERROR: the query to the server was not successful. Is your internet connection working?\n")
      } else {
        stop(paste0("ERROR: the query to the server was not successful. The server returned the following answer:\n", serverAnswer))
      }
    }
  }

  ## Third query: gene to organs mapping
  gene2anatomyFileName <- paste0("topAnat_GeneToAnatEntities_", myBgeeObject$speciesId, "_", toupper(callType))
  ## If a stage is specified, add it to file name
  if ( !is.null(stage) ){
    gene2anatomyFileName <- paste0(gene2anatomyFileName, "_", gsub(":", "_", stage))
  }
  ## If all data types specified, no need to add anything to file name. Otherwise, specify data types in file name
  if ( sum(myBgeeObject$dataType %in% c("rna_seq","affymetrix","est","in_situ")) < 4 ){
    gene2anatomyFileName <- paste0(gene2anatomyFileName, "_", toupper(paste(sort(myBgeeObject$dataType), collapse="_")))
  }
  ## If high quality data needed, specify in file name. Otherwise not specified
  if(compareVersion(gsub("_", ".", myBgeeObject$release), OLD_WEBSERVICE_VERSION) > 0){
    gene2anatomyFileName <- paste0(gene2anatomyFileName, toupper(confidence))
  } else {
    if ( confidence == "high_quality" ){
      gene2anatomyFileName <- paste0(gene2anatomyFileName, "_HIGH")
    }
  }
  gene2anatomyFileName <- paste0(gene2anatomyFileName, ".tsv")

  ## Check if file is already in cache
  if (file.exists(file.path(myBgeeObject$pathToData, gene2anatomyFileName))){
    cat(paste0("\nNOTE: a gene to organs mapping file was found in the download directory ", myBgeeObject$pathToData,
               ". Data will not be redownloaded.\n"))

  } else {
    cat("\nBuilding URLs to retrieve mapping of gene to organs from Bgee...\n")
    myUrl <- myBgeeObject$topAnatUrl
    if(compareVersion(gsub("_", ".", myBgeeObject$release), OLD_WEBSERVICE_VERSION) > 0){
      myUrl <- paste0(myBgeeObject$topAnatUrl, "?page=r_package&action=get_expression_calls&display_type=tsv&species_list=", myBgeeObject$speciesId, "&attr_list=GENE_ID&attr_list=ANAT_ENTITY_ID&api_key=", myBgeeObject$apiKey, "&source=BgeeDB_R_package&source_version=", as.character(packageVersion("BgeeDB")))
    }else {
      myUrl <- paste0(myUrl, "?page=dao&action=org.bgee.model.dao.api.expressiondata.ExpressionCallDAO.getExpressionCalls&display_type=tsv&species_list=", myBgeeObject$speciesId, "&attr_list=GENE_ID&attr_list=ANAT_ENTITY_ID&api_key=", myBgeeObject$apiKey, "&source=BgeeDB_R_package&source_version=", as.character(packageVersion("BgeeDB")))
    }

    ## Add data type to file name: only if not all data types asked
    if ( sum(myBgeeObject$dataType %in% c("rna_seq","affymetrix","est","in_situ")) < 4 ){
      for (type in toupper(sort(myBgeeObject$dataType))){
        myUrl <- paste0(myUrl, "&data_type=", type)
      }
    }
    ## Add data quality
    if(compareVersion(gsub("_", ".", myBgeeObject$release), OLD_WEBSERVICE_VERSION) > 0){
      myUrl <- paste0(myUrl, "&data_qual=", toupper(confidence))
    } else {
      if(confidence == "high_quality"){
        myUrl <- paste0(myUrl, "&data_qual=HIGH")
      }
    }

    if ( !is.null(stage) ){
      myUrl <- paste0(myUrl, "&stage_id=", stage)
    }

    ## Query webservice
    cat(paste0("   URL successfully built (", myUrl,")\n   Submitting URL to Bgee webservice (can be long)\n"))
    success <- download.file(myUrl, destfile = paste0(myBgeeObject$pathToData, "/", gene2anatomyFileName, ".tmp"))

    if (success == 0){
      ## Read 5 last lines of file: should be empty indicating success of data transmission
      ## We cannot use a system call to UNIX command since some user might be on Windows
      tmp <- tail(read.table(paste0(myBgeeObject$pathToData, "/", gene2anatomyFileName, ".tmp"), header=TRUE, sep="\t", comment.char="", blank.lines.skip=FALSE, as.is=TRUE), n=5)
      if ( length(tmp[,1]) == 5 && (sum(tmp[,1] == "") == 5 || sum(is.na(tmp[,1])) == 5) ){
        ## The file transfer was successful, we rename the temporary file
        file.rename(paste0(myBgeeObject$pathToData, "/", gene2anatomyFileName, ".tmp"), paste0(myBgeeObject$pathToData, "/", gene2anatomyFileName))
      } else {
        ## delete the temporary file
        file.remove(paste0(myBgeeObject$pathToData, "/", gene2anatomyFileName, ".tmp"))
        stop(paste0("File ", gene2anatomyFileName, " is truncated, there may be a temporary problem with the Bgee webservice, or there was an error in the parameters."))
      }
      cat(paste0("   Got results from Bgee webservice. Files are written in \"", myBgeeObject$pathToData, "\"\n"))
    } else {
      serverAnswer = try(getURL(myUrl))
      if (class(serverAnswer) == "try-error"){
        stop("ERROR: the query to the server was not successful. Is your internet connection working?\n")
      } else {
        stop(paste0("ERROR: the query to the server was not successful. The server returned the following answer:\n", serverAnswer))
      }
    }
  }

  ## Process the data and build the final list to return
  cat("\nParsing the results.............................................\n")

  ## Relationships between organs
  if (file.exists(file.path(myBgeeObject$pathToData, organRelationshipsFileName))){
    if (file.info(file.path(myBgeeObject$pathToData, organRelationshipsFileName))$size != 0) {
      tab <- read.table(file.path(myBgeeObject$pathToData, organRelationshipsFileName), header=TRUE, sep="\t", blank.lines.skip=TRUE, as.is=TRUE)
      organRelationships <- tapply(as.character(tab$TARGET_ID), as.character(tab$SOURCE_ID), unique)
    } else {
      stop(paste0("File ", organRelationshipsFileName, " is empty, there may be a temporary problem with the Bgee webservice, or there was an error in the parameters."))
    }
  } else {
    stop(paste0("File ", organRelationshipsFileName, " not found. There may be a temporary problem with the Bgee webservice, or there was an error in the parameters."))
  }
  ## Organ names
  if (file.exists(file.path(myBgeeObject$pathToData, organNamesFileName))){
    if (file.info(file.path(myBgeeObject$pathToData, organNamesFileName))$size != 0) {
      organNames <- read.table(file.path(myBgeeObject$pathToData, organNamesFileName), header=TRUE, sep="\t", comment.char="", blank.lines.skip=TRUE, as.is=TRUE, quote = "")
    } else {
      stop(paste0("File ", organNamesFileName, " is empty, there may be a temporary problem with the Bgee webservice, or there was an error in the parameters."))
    }
  } else {
    stop(paste0("File ", organNamesFileName, " not found. There may be a temporary problem with the Bgee webservice, or there was an error in the parameters."))
  }
  ## Mapping of genes to tissues
  if (file.exists(file.path(myBgeeObject$pathToData, gene2anatomyFileName))){
    if (file.info(file.path(myBgeeObject$pathToData, gene2anatomyFileName))$size != 0) {
      tab <- read.table(file.path(myBgeeObject$pathToData, gene2anatomyFileName), header=TRUE, sep="\t", blank.lines.skip=TRUE, as.is=TRUE)
      if(length(tab$GENE_ID) != 0){
        gene2anatomy <- tapply(as.character(tab$ANAT_ENTITY_ID), as.character(tab$GENE_ID), unique)
      } else {
        stop("There was no mapping of genes to anatomical structures found. Probably the parameters are too stringent, or this data type is absent in this species. See listBgeeSpecies() for data types availability.")
      }
    } else {
      stop(paste0("File ", gene2anatomyFileName, " is empty, there may be a temporary problem with the Bgee webservice, or there was an error in the parameters. It is also possible that the parameters are too stringent and returned no data, please try to relax them."))
    }
  } else {
    stop(paste0("File ", gene2anatomyFileName, " not found. There may be a temporary problem with the Bgee webservice, or there was an error in the parameters. It is also possible that the parameters are too stringent and returned no data, please try to relax them."))
  }

  cat("\nAdding BGEE:0 as unique root of all terms of the ontology.......\n")
  ## There can be multiple roots among all the terms downloaded. We need to add one unique root for topGO to work: BGEE:0
  ## Add all organs from organNames that are not source (child / names of the list) in organsRelationship to the organsRelationship list (with value / target / parent = BGEE:0)
  ## Some organs are not present in the organRelationships file because they have no relations to other organs (not linked to a target).
  ## That's why we use all organs present in organNames to find missingParents
  missingParents <- setdiff(organNames[, 1], names(organRelationships))

  ## Add new values
  organRelationships <- c(organRelationships, as.list(rep("BGEE:0", times=length(missingParents))))
  ## Add new keys
  names(organRelationships)[(length(organRelationships) - length(missingParents) + 1):length(organRelationships)] = as.character(missingParents)
  ## Add BGEE:0	/ root to organNames
  organNames <- rbind(organNames, c("BGEE:0", "root"))

  ## Check if some organ names are missing, and add them if necessary
  missingNames <- setdiff(unique(unique(unlist(organRelationships, use.names = FALSE)), names(organRelationships)), organNames$ID)
  if (length(missingNames) > 0){
    cat(paste0("\nWARNING: some organs names appear to be missing. There might be some problem with the ontology data.\n"))
    organNames <- rbind(organNames, setNames(data.frame(missingNames, "?"), names(organNames)))
  }

  cat("\nDone.\n")
  return(list(gene2anatomy = gene2anatomy, organ.relationships = organRelationships, organ.names = organNames, bgee.object = myBgeeObject))
}
