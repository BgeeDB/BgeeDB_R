#' @title Retrieve Bgee calls from one species.
#'
#' @description This function loads the integrated expression calls from one species. These calls are 
#' equivalent to the one provided in the gene page of the Bgee website. The calls have been
#' generated using all datatypes present in Bgee (EST, In Situ, Affymetrix, bulk RNA-Seq and 
#' single-cell RNA-Seq).
#'
#' @param myBgeeObject A Reference Class Bgee object, notably specifying the targeted species 
#' and data type.
#'
#' @param conditonParameters Specify which condition parameter you are interested
#' in. It can be `anatEntity` or `allCondParams`. The `anatEntity` option retrieve calls generated 
#' using only anatomical entity and cell-type as a condition parameter. The `allCondParams` option retrieve calls generated taking
#' into consideration all condition parameters present in Bgee (anat, entity, cell-type, developmental stage,
#' sex and strain). By default the `anatEntity` option is selected.
#' 
#' @param advancedColumns Boolean allowing to specify if advancend columns are required. The default value
#' is `FALSE` meaning that only gene , condition parameters, call quality, FDR, expression rank and
#' expression score are retrieved. If `TRUE` then a lot more information used to generate the calls are
#' retrieved like the number of self and descendant observation for each datatype, or the rank and score
#' per datatype. 
#' 
#' @param geneIds List of genes for which expression calls have to be retrieved. It has to be ensembl
#' IDs if the Bgee genome source is Ensembl (e.g ENSG00000244734) or RefSeq IDs if the Bgee genome source
#' is RefSeq (e.g 734881) but not gene names (e.g HBB, Apoc1, etc.).
#' 
#' @param anatEntityIds List of anatomical entity IDs for which expression calls have to be retrieved.
#'
#' @return Return a data.table containing all calls from the specified species
#'
#' @author Julien Wollbrett, Andrea Komljenovic and Julien Roux.
#'
#' @examples{
#'   bgee <- Bgee$new(species = "Danio_rerio")
#'   callsHer1 <- getIntegratedCalls(bgee, geneIds="ENSDARG00000014722")
#' }
#'
#' @import data.table bread
#' @importFrom R.utils gunzip
#' @export
#' 
getIntegratedCalls <- function(myBgeeObject, conditonParameters="anatEntity", advancedColumns=FALSE,
    geneIds=NULL) {
  if (! conditionParameters %in% c("anatEntity", "allCondParams")) {
    stop("the value of the conditionParameters argument should be either anatEntity or callCondParams.")
  }
  if(compareVersion(a = gsub("_", ".", myBgeeObject$release), b = "14.0") < 0) {
    stop("Integrated expression calls are available starting from Bgee 14.0.")
  }
  if(compareVersion(a = gsub("_", ".", myBgeeObject$release), b = "15.0") < 0 & conditionParameters != "anatEntity") {
    stop("calls files containing all condition parameters are only available starting from Bgee 15.0.\nPlease choose ",
         "conditionParameters=\"anatEntity\" if you really want to retrieve integrated expression calls from Bgee ",
         gsub("_", ".", myBgeeObject$release))
  }
  message("Start to download calls generated using all Bgee datatypes for ",
    gsub(pattern = "_", replacement = " ", x = myBgeeObject$speciesName))
  # Start to generate the proper download URL
  callsUrl <- myBgeeObject$callsUrl
  # Replace SPECIES_NAME pattern by actual species name
  callsUrl <- gsub(pattern = "SPECIES_NAME", replacement = myBgeeObject$speciesName, x = callsUrl)
  # Now replace COLUMN_TYPE pattern by the type of columns to retrieve
  callsUrl <- gsub(pattern = "COLUMN_TYPE", replacement = ifelse(advancedColumns, "expr_advanced", "expr_simple"),
    x = callsUrl)
  # And finally replace _COND_PARAMS pattern by the condition parameters requested
  callsUrl <- gsub(pattern = "_COND_PARAMS", replacement = ifelse(conditionParameters == "anatEntity",
    "", "_all_conditions"), x = callsUrl)
  
  # download the file and uncompress it
  destFile <- file.path(myBgeeObject$pathToData, basename(callsUrl))
  uncompressedDestFile <- gsub(pattern = ".gz$", replacement = "", x = destFile)
  if (!file.exists(uncompressedDestFile)) {
    bgee_download_file(url = callsUrl, destfile = destFile, mode = 'wb')
    output <- gunzip(destFile)
  }
  
  # combine genes to search in order to retrieve all rows with one of the gene Id
  outputCalls <- NULL
  timeNeeded <- NULL
  if (is.null(geneIds) & is.null(anatEntityIds)) {
    print("load without filter")
    timeNeeded <- system.time(outputCalls <- fread(file = uncompressedDestFile, header = TRUE,
      quote = FALSE))
  } else {
    if (is.null(anatEntityIds) | !is.null(geneIds) & length(geneIds) < length(anatEntityIds)) {
      print("filter by genes with bread")
      geneSearchPattern <- paste(geneIds, collapse =  "|")
      timeNeeded <- system.time(outputCalls <- bread(file = uncompressedDestFile,
        patterns = geneSearchPattern, filtered_columns = "Gene ID"))
      if (!is.null(anatEntityIds)) {
        outputCalls <- outputCalls[outputCalls$`Anatomical entity ID` %in% anatEntityIds,]
      }
    } else if (is.null(geneIds) | !is.null(anatEntityIds) & length(anatEntityIds) < length(geneIds)) {
      print("filter by anatEntities with bread")
      anatEntitySearchPattern <- paste(anatEntityIds, collapse = "|")
      timeNeeded <- system.time(outputCalls <- bfilter(file = uncompressedDestFile,
        patterns = anatEntitySearchPattern, filtered_columns = "Anatomical entity ID"))
      if (!is.null(geneIds)) {
        outputCalls <- outputCalls[outputCalls$`Gene ID` %in% geneIds,]
      }
    }
  }
  message("Expression calls have been retrieved in ", count_time(timeNeeded["user.self"]
    + timeNeeded["user.child"]))
  return(outputCalls)
  
}

breadWithErrorManagement <- function (file, patterns, filtered_columns) {
  outputCalls <- NULL
  tryCatch( {
    outputCalls <- bread(file = uncompressedDestFile,
      patterns = geneSearchPattern, filtered_columns = "Gene ID")
    },
    error = function(e) {stop("can not run the bread function on your system. If you are ",
      "on Windows please first install 'RTools', 'Git' or 'Cygwin' to allow Unix commands ",
      "like \'grep\', \'sed\', \'wc\', \'awk\' or \'cut\'. If you do not want to install ",
      "such tools you can also run the getIntegratedCalls() function without geneIds and ",
      "anatEntityIds and filter afterwards.")})
}