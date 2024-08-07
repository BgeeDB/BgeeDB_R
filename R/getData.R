#' @title Retrieve Bgee RNA-seq or Affymetrix data.
#'
#' @description This function loads the quantitative expression data and presence calls 
#' for samples available from Bgee (rna_seq, affymetrix, sc_full_length).
#'
#' @param myBgeeObject A Reference Class Bgee object, notably specifying the targeted species 
#' and data type.
#'
#' @param experimentId Filter allowing to specify one or more ArrayExpress or GEO accession, e.g., 
#' GSE43721. Default is NULL: takes all available experiments for targeted species and data type.
#' 
#' @param sampleId Filter allowing to specify one or more sample ID. Depending on the selected 
#' datatype this sample IDs can correspond to Chip IDs (affymetrix) or RNA-Seq library IDs (rna_seq). 
#' Default is NULL: takes all available samples for targeted species and data type.
#'
#' @param anatEntityId Filter allowing to specify one or more anatomical entity IDs from the UBERON 
#' ontology (http://uberon.github.io/). Default is NULL: takes all available anatomical entities for 
#' targeted species and data type.
#'
#' @param stageId Filter allowing to specify one or more developmental stage IDs from Developmental 
#' Stage Ontology (https://github.com/obophenotype/developmental-stage-ontologies). Default is 
#' NULL: takes all available developmental stages for targeted species and data type.
#' 
#' @param cellTypeId Filter specific to single cell datatype (sc_full_length) allowing to specify 
#' one or more cell type IDs from the UBERON ontology (http://uberon.github.io/). Default is 
#' NULL: takes all available cell types for targeted species and data type. Available for Bgee 15.0 and after
#' 
#' @param sex Filter allowing to specify one or more sexes. Default is 
#' NULL: takes all available sexes for targeted species and data type. Available for Bgee 15.0 and after
#' 
#' @param strain Filter allowing to specify one or more strains. Default is 
#' NULL: takes all available strains for targeted species and data type. Available for Bgee 15.0 and after
#' 
#' @param withDescendantAnatEntities Allows to filter on the selected anatEntityId and all its descendants.
#' This functionality is available for Bgee 15.0 release and after
#'
#' @param withDescendantStages Allows to filter on the selected stageId and all its descendants.
#' This functionality is available for Bgee 15.0 release and after
#'
#' @param withDescendantCellTypes Allows to filter on the selected cellTypeId and all its descendants.
#' This functionality is available for Bgee 15.0 release and after
#'
#' @return Return a dataframe containing all Bgee processed expression data from the selected species 
#' and datatype using specified filters with operator AND.
#'
#' @author Julien Wollbrett, Andrea Komljenovic and Julien Roux.
#'
#' @examples{
#'   bgee <- Bgee$new(species = "Mus_musculus", dataType = "rna_seq")
#'   dataMouseGSE43721 <- getData(bgee, experimentId = "GSE43721")
#'   dataMouseVariousFilters <- getData(bgee, experimentId = c("GSE43721", "GSE36026"), 
#'                              anatEntityId = c("UBERON:0002107", "UBERON:0000956", "UBERON:0002048"))
#' }
#'
#' @import RSQLite
#' @importFrom R.utils gunzip
#' @export
#' 
getData <- function(myBgeeObject, experimentId = NULL, sampleId = NULL, 
    anatEntityId = NULL, stageId = NULL, cellTypeId = NULL, sex = NULL, strain = NULL, 
    withDescendantAnatEntities = FALSE, withDescendantStages = FALSE, 
    withDescendantCellTypes = FALSE) {
  .Deprecated(new = "getSampleRawData")
  check_object(myBgeeObject)
  # write a warning message if user tried to retrieve ontology terms with descendants for
  # a bgee release where this functionality was not yet implemented
  compare_version_number <- gsub("_", ".", myBgeeObject$release)
  if ((withDescendantAnatEntities | withDescendantStages | withDescendantCellTypes) &
      compareVersion(a = compare_version_number , b = "15.0") < 0) {
    message("withDescendant functionality is available only for Bgee 15.0",
            " release and after. Will not retrieve descendant of selected parameters.")
  }
  if (withDescendantAnatEntities & compareVersion(a = compare_version_number , b = "15.0") >= 0) {
    if(is.null(anatEntityId)) {
      warning("No anatomical entity was provided. Not possible to filter on descendant anatomical entities.")
    }
    anatEntityId <- c(anatEntityId, getDescendantAnatEntities(bgee = myBgeeObject, ids = anatEntityId))
  }
  if (withDescendantStages & compareVersion(a = compare_version_number , b = "15.0") >= 0) {
    if(is.null(stageId)) {
      warning("No developmental stage was provided. Not possible to filter on descendant developmental stages.")
    }
    stageId <- c(stageId, getDescendantStages(bgee = myBgeeObject, ids = stageId))
  }
  if (withDescendantCellTypes & compareVersion(a = compare_version_number , b = "15.0") >= 0) {
    if(is.null(cellTypeId)) {
      warning("No cell type was provided. Not possible to filter on descendant cell types.")
    }
    cellTypeId <- c(cellTypeId, getDescendantCellTypes(bgee = myBgeeObject, ids = cellTypeId))
  }
  check_condition_parameters(myBgeeObject = myBgeeObject, anatEntityId = anatEntityId, 
    stageId = stageId, cellTypeId = cellTypeId, sex = sex, strain = strain)
  import_data(myBgeeObject = myBgeeObject, experimentId = experimentId, sampleId = sampleId, 
    anatEntityId = anatEntityId, stageId = stageId, cellTypeId = cellTypeId, sex = sex, 
    strain = strain)
  return(query_data(myBgeeObject = myBgeeObject, experimentId = experimentId, sampleId = sampleId, 
    anatEntityId = anatEntityId, stageId = stageId, cellTypeId = cellTypeId, sex = sex, 
    strain = strain))
}
