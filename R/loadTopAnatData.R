########################
#' @title Retrieve data from Bgee to perform GO-like enrichment of anatomical terms, mapped to genes by expression patterns.
#'
#' @description This function loads a mapping from genes to anatomical structures based on calls of expression in anatomical structures. It also loads the structure of the anatomical ontology.
#'
#' @details The expression calls come from Bgee (\url{http://bgee.org}), that integrates different expression data types (RNA-seq, Affymetrix microarray, ESTs, or in-situ hybridizations) in multiple animal species. Expression patterns are based exclusively on curated "normal", healthy, expression data (e.g., no gene knock-out, no treatment, no disease), to provide a reference of normal gene expression.
#'
#' Anatomical structures are identified using IDs from the Uberon ontology (browsable at \url{http://www.ontobee.org/ontology/UBERON}). The mapping from genes to anatomical structures includes only the evidence of expression in these specific structures, and not the expression in their substructures (i.e., expression data are not propagated). The retrieval of propagated expression data will likely be implemented in the future, but meanwhile, it can be obtained using specialized packages such as topGO, see the \code{topAnat.R} function.
#'
#' @param species A character indicating the species to be used, in the form "Genus_species", or a numeric indicating the species NCBI taxonomic id. Only species with data in Bgee will work. See the listBgeeSpecies() function to get the list of species available in the Bgee release used.
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
#' Only calls for significant detection of expression are implemented ("presence").
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
#' @param confidence A character indicating if only high quality present calls should be retrieved.
#'  Options are "all" or "high_quality". Default is "all".
#'
#' @param release Bgee release number to download data from, in the form "Release.subrelease" or "Release_subrelease", e.g., "13.2" or 13_2". Will work for release >=13.2. By default, the latest relase of Bgee is used.
#'
#' @param pathToData Path to the directory where the data files are stored. By default the working directory is used.
#'
#' @return A list of 3 elements:
#' \itemize{
#'   \item{A \code{gene2anatomy} list, mapping genes to anatomical structures based on expression calls.}
#'   \item{A \code{organ.names} data frame, with the name corresponding to UBERON IDs.}
#'   \item{A \code{organ.relationships} list, giving the relationships between anatomical structures
#'   in the UBERON ontology (based on parent-child "is_a" and "part_of" relationships).}
#' }
#'
#' @author Julien Roux
#'
#' @examples{
#'   myTopAnatData <- loadTopAnatData(species = "Mus_musculus", datatype = "rna_seq")
#' }
#'
#' @import utils
#' @export

loadTopAnatData <- function(species, datatype=c("rna_seq","affymetrix","est","in_situ"), calltype="presence",
                            confidence="all", stage=NULL, release=NULL, pathToData=getwd()){

  cat("Querying Bgee to get release information...\n")
  allReleases <- .getRelease()
  if (length(release)==0) {
    release <- gsub("\\.", "_", allReleases$release[1])
  } else if (length(release)==1){
    # In case the release number is written with a dot
    release <- gsub("\\.", "_", release)
    # test if required release exists
    if (sum(allReleases$release == gsub("_", ".", release))!=1){
      stop("ERROR: The specified release number is invalid, or is not available for BgeeDB.")
    }
  } else {
    stop("ERROR: The specified release number is invalid.")
  }

  ## Specify host to be used
  host <- allReleases$TopAnat.URL[allReleases$release == gsub("_", ".", release)]
  if ( !grepl("/$", host) ){
    host <- paste0(host, "/")
  }

  ## Retrieve list of all species for queried release
  allSpecies <- listBgeeSpecies(release=release, allReleases=allReleases)

  ## Species is the only compulsory parameter
  if( length(species) == 0 ) {
    stop("ERROR: you need to specify a species.")
  } else if ( length(species) > 1 ){
    stop("ERROR: only one species is allowed.")
  } else if (grepl("^\\d+$", species)){
    ## of species was specified as a taxonomic ID
    if (sum(allSpecies$ID == species) != 1){
      stop(paste0("ERROR: The specified species Id is invalid, or not available in Bgee release ", release))
    } else {
      speciesId <- as.numeric(species)
      speciesName <- paste(allSpecies[allSpecies$ID == species, 2:3], collapse="_")
    }
  } else if (is.character(species)){
    speciesSplitted <- unlist(strsplit(species, split="_"))
    if (sum(allSpecies$GENUS == speciesSplitted[1] & allSpecies$SPECIES_NAME == speciesSplitted[2]) != 1){
      stop(paste0("ERROR: The specified species name is invalid, or not available in Bgee release ", release, "."))
    } else {
      speciesName <- species
      speciesId <- as.numeric(allSpecies$ID[allSpecies$GENUS == speciesSplitted[1] & allSpecies$SPECIES_NAME == speciesSplitted[2]])
    }
  }

  ## Test if parameters are in the range of allowed parameters
  if ( !sum(datatype %in% c("rna_seq","affymetrix","est","in_situ")) %in% 1:4 ){
    stop("ERROR: you need to specify at least one valid data type to be used among \"rna_seq\", \"affymetrix\", \"est\" and \"in_situ\".")
  }
  if ( length(datatype) != sum(datatype %in% c("rna_seq","affymetrix","est","in_situ")) ){
    cat("WARNING: you apparently specified a data type that is not among \"rna_seq\", \"affymetrix\", \"est\" and \"in_situ\". Please check for mistakes or typos.\n")
  }
  if ( calltype != "presence" ){
    stop("ERROR: no other call types than present expression calls can be retrieved for now.")
  }
  if ( (confidence != "all") && (confidence != "high_quality") ){
    stop("ERROR: the data confidence parameter specified is not among the allowed values (\"all\" or \"high_quality\").")
  }

  ## check path of folder to store cached files
  if(length(pathToData)==0) {
    pathToData <- paste0(getwd(), "/", speciesName, "_Bgee_", release)
  } else if (length(pathToData)==1){
    if ( !file.exists(pathToData) ){
      stop("ERROR: please specify a valid and existing path to store data files.")
    } else {
      pathToData <- paste0(pathToData, "/", speciesName, "_Bgee_", release)
    }
  } else {
    stop("ERROR: Invalid path for data files.")
  }
  ## create sub-folder with species name to store downloaded files
  if (!file.exists(pathToData)){
    dir.create(pathToData)
  }

  ## Set the internet.info to 2 to have less verbose output (only reports critical warnings)
  options(internet.info=2)
  ## Set the timeout option to 600 seconds to let some time to the server to send data (default is 60s)
  options(timeout = 600)

  ## First query: organ relationships
  organRelationshipsFileName <- paste0("topAnat_AnatEntitiesRelationships_", speciesId, ".tsv")
  ## Check if file is already in cache
  if (file.exists(file.path(pathToData, organRelationshipsFileName))){
    cat("\nWARNING: an organ relationships file was found in the working directory and will be used as is. Please delete and rerun the function if you want the data to be updated.\n")
  } else {
    cat("\nBuilding URLs to retrieve organ relationships from Bgee.........\n")
    myurl <- paste0(host, "?page=dao&action=org.bgee.model.dao.api.ontologycommon.RelationDAO.getAnatEntityRelations&display_type=tsv&species_list=", speciesId, "&attr_list=SOURCE_ID&attr_list=TARGET_ID")

    ## Query webservice
    cat(paste0("   URL successfully built (", myurl,")\n   Submitting URL to Bgee webservice (can be long)\n"))
    download.file(myurl, destfile = paste0(pathToData, "/", organRelationshipsFileName, ".tmp"))

    ## Read 5 last lines of file: should be empty indicating success of data transmission
    ## We cannot use a system call to UNIX command since some user might be on Windows
    tmp <- tail(read.table(paste0(pathToData, "/", organRelationshipsFileName, ".tmp"), header=TRUE, sep="\t", comment.char="", blank.lines.skip=FALSE, as.is=TRUE), n=5)
    if ( length(tmp[,1]) == 5 && (sum(tmp[,1] == "") == 5 || sum(is.na(tmp[,1])) == 5) ){
      ## The file transfer was successful, we rename the temporary file
      file.rename(paste0(pathToData, "/", organRelationshipsFileName, ".tmp"), paste0(pathToData, "/", organRelationshipsFileName))
    } else {
      ## delete the temporary file
      file.remove(paste0(pathToData, "/", organRelationshipsFileName, ".tmp"))
      stop(paste0("File ", organRelationshipsFileName, " is truncated, there may be a temporary problem with the Bgee webservice, or there was an error in the parameters."))
    }
    cat(paste0("   Got results from Bgee webservice. Files are written in \"", pathToData, "\"\n"))
  }

  ## Second query: organ names
  organNamesFileName <- paste0("topAnat_AnatEntitiesNames_", speciesId, ".tsv");
  ## Check if file is already in cache
  if (file.exists(file.path(pathToData, organNamesFileName))){
    cat("\nWARNING: an organ names file was found in the working directory and will be used as is. Please delete and rerun the function if you want the data to be updated.\n")
  } else {
    cat("\nBuilding URLs to retrieve organ names from Bgee.................\n")
    myurl <-  paste0(host, "?page=dao&action=org.bgee.model.dao.api.anatdev.AnatEntityDAO.getAnatEntities&display_type=tsv&species_list=", speciesId, "&attr_list=ID&attr_list=NAME")

    ## Query webservice
    cat(paste0("   URL successfully built (", myurl,")\n   Submitting URL to Bgee webservice (can be long)\n"))
    download.file(myurl, destfile = paste0(pathToData, "/", organNamesFileName, ".tmp"))

    ## Read 5 last lines of file: should be empty indicating success of data transmission
    ## We cannot use a system call to UNIX command since some user might be on Windows
    tmp <- tail(read.table(paste0(pathToData, "/", organNamesFileName, ".tmp"), header=TRUE, sep="\t", comment.char="", blank.lines.skip=FALSE, as.is=TRUE, quote = ""), n=5)
    if ( length(tmp[,1]) == 5 && (sum(tmp[,1] == "") == 5 || sum(is.na(tmp[,1])) == 5) ){
      ## The file transfer was successful, we rename the temporary file
      file.rename(paste0(pathToData, "/", organNamesFileName, ".tmp"), paste0(pathToData, "/", organNamesFileName))
    } else {
      ## delete the temporary file
      file.remove(paste0(pathToData, "/", organNamesFileName, ".tmp"))
      stop(paste0("File ", organNamesFileName, " is truncated, there may be a temporary problem with the Bgee webservice, or there was an error in the parameters."))
    }
    cat(paste0("   Got results from Bgee webservice. Files are written in \"", pathToData, "\"\n"))
  }

  ## Third query: gene to organs mapping
  gene2anatomyFileName <- paste0("topAnat_GeneToAnatEntities_", speciesId, "_", toupper(calltype))
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
  if (file.exists(file.path(pathToData, gene2anatomyFileName))){
    cat("\nWARNING: a gene to organs mapping file was found in the working directory and will be used as is. Please delete and rerun the function if you want the data to be updated.\n")
  } else {
    cat("\nBuilding URLs to retrieve mapping of gene to organs from Bgee...\n")
    myurl <-  paste0(host, "?page=dao&action=org.bgee.model.dao.api.expressiondata.ExpressionCallDAO.getExpressionCalls&display_type=tsv&species_list=", speciesId, "&attr_list=GENE_ID&attr_list=ANAT_ENTITY_ID")

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
    cat(paste0("   URL successfully built (", myurl,")\n   Submitting URL to Bgee webservice (can be long)\n"))
    download.file(myurl, destfile = paste0(pathToData, "/", gene2anatomyFileName, ".tmp"))

    ## Read 5 last lines of file: should be empty indicating success of data transmission
    ## We cannot use a system call to UNIX command since some user might be on Windows
    tmp <- tail(read.table(paste0(pathToData, "/", gene2anatomyFileName, ".tmp"), header=TRUE, sep="\t", comment.char="", blank.lines.skip=FALSE, as.is=TRUE), n=5)
    if ( length(tmp[,1]) == 5 && (sum(tmp[,1] == "") == 5 || sum(is.na(tmp[,1])) == 5) ){
      ## The file transfer was successful, we rename the temporary file
      file.rename(paste0(pathToData, "/", gene2anatomyFileName, ".tmp"), paste0(pathToData, "/", gene2anatomyFileName))
    } else {
      ## delete the temporary file
      file.remove(paste0(pathToData, "/", gene2anatomyFileName, ".tmp"))
      stop(paste0("File ", gene2anatomyFileName, " is truncated, there may be a temporary problem with the Bgee webservice, or there was an error in the parameters."))
    }
    cat(paste0("   Got results from Bgee webservice. Files are written in \"", pathToData, "\"\n"))
  }

  ## Process the data and build the final list to return
  cat("\nParsing the results.............................................\n")

  ## Relationships between organs
  if (file.exists(file.path(pathToData, organRelationshipsFileName))){
    if (file.info(file.path(pathToData, organRelationshipsFileName))$size != 0) {
      tab <- read.table(file.path(pathToData, organRelationshipsFileName), header=TRUE, sep="\t", blank.lines.skip=TRUE, as.is=TRUE)
      organRelationships <- tapply(as.character(tab$TARGET_ID), as.character(tab$SOURCE_ID), unique)
    } else {
      stop(paste0("File ", organRelationshipsFileName, " is empty, there may be a temporary problem with the Bgee webservice, or there was an error in the parameters."))
    }
  } else {
    stop(paste0("File ", organRelationshipsFileName, " not found. There may be a temporary problem with the Bgee webservice, or there was an error in the parameters."))
  }
  ## Organ names
  if (file.exists(file.path(pathToData, organNamesFileName))){
    if (file.info(file.path(pathToData, organNamesFileName))$size != 0) {
      organNames <- read.table(file.path(pathToData, organNamesFileName), header=TRUE, sep="\t", comment.char="", blank.lines.skip=TRUE, as.is=TRUE, quote = "")
    } else {
      stop(paste0("File ", organNamesFileName, " is empty, there may be a temporary problem with the Bgee webservice, or there was an error in the parameters."))
    }
  } else {
    stop(paste0("File ", organNamesFileName, " not found. There may be a temporary problem with the Bgee webservice, or there was an error in the parameters."))
  }
  ## Mapping of genes to tissues
  if (file.exists(file.path(pathToData, gene2anatomyFileName))){
    if (file.info(file.path(pathToData, gene2anatomyFileName))$size != 0) {
      tab <- read.table(file.path(pathToData, gene2anatomyFileName), header=TRUE, sep="\t", blank.lines.skip=TRUE, as.is=TRUE)
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
  missingParents <- organNames$ID[!organNames$ID %in% names(organRelationships)]
  ## Add new values
  organRelationships <- c(organRelationships, as.list(rep("BGEE:0", times=length(missingParents))))
  ## Add new keys
  names(organRelationships)[(length(organRelationships) - length(missingParents) + 1):length(organRelationships)] = as.character(missingParents)
  ## Add BGEE:0	/ root to organNames
  organNames <- rbind(organNames, c("BGEE:0", "root"))

  cat("\nDone.\n")
  return(list(gene2anatomy = gene2anatomy, organ.relationships = organRelationships, organ.names = organNames))
}
