#' @title Galaxy adapter to the `getAnnotation` function
#' 
#' @description
#' Adapter of the `getAnnotation` function allowing to run it on Galaxy.
#' 
#' @param species name of the species for which annotation have to be be retrieved (e.g Homo_sapiens)
#' @param dataType The datatype for which annotation have to be retrieved. To be chosen among:
#' \itemize{
#'   \item{"rna_seq: target only bulk RNA-seq data"}
#'   \item{"sc_full_length: target only full length single-cell RNA-seq data"}
#'   \item{"sc_droplet_based: target only droplet based single-cell RNA-seq data"}
#'   \item{"affymetrix: target only Affymtrix microarray data"}
#' }
#' @param experimentAnnotationPath path to the file containing experiment annotation
#' @param sampleAnnotationPath path to the file containing sample annotation
#' @export

getAnnotation_galaxy <- function(species = NULL, dataType = character(0),
                                 experimentAnnotationPath = NULL, sampleAnnotationPath = NULL) {
  if (is.null(experimentAnnotationPath) || is.null(sampleAnnotationPath)) {
    stop("path to output files should always be provided")
  }
  myBgeeObject <- createBgeeObject_galaxy(species = species, dataTypes = dataType)
  annotation <- getAnnotation(myBgeeObject = myBgeeObject)
  write.table(x = annotation$experiment.annotation, file = experimentAnnotationPath,
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(x = annotation$sample.annotation, file = sampleAnnotationPath,
              sep = "\t", quote = FALSE, row.names = FALSE)
}

#' @title Galaxy adapter to the `getSampleProcessedData` function
#' 
#' @description
#' This function loads the quantitative expression data and presence calls 
#' for samples available from Bgee. These data are available for all droplet based single-cell
#' RNA-seq, full length single-cell RNA-seq, bulk RNA-seq and affymetrix samples present in Bgee.
#' It is possible to filter the processed expression data using experiment IDs, chip IDs (affymetrix)
#' or library IDs (RNA-seq), anatomical entity IDs from the UBERON ontology, developmental stage IDs,
#' cell-type IDs, sex and/or strain. Celltype, anatomical entity and stage IDs are normalized and come from
#' ontologies, allowing to retrieve data annotated with the selected term or with the term and all its descendants.
#' This option allows for instance to retrieve all processed expression values coming from the brain or
#' any subpart of the brain.
#' 
#' @param species name of the species for which annotation have to be be retrieved (e.g Homo_sapiens)
#' @param dataType The datatype for which annotation have to be retrieved. To be chosen among:
#' \itemize{
#'   \item{"rna_seq: target only bulk RNA-seq data"}
#'   \item{"sc_full_length: target only full length single-cell RNA-seq data"}
#'   \item{"sc_droplet_based: target only droplet based single-cell RNA-seq data"}
#'   \item{"affymetrix: target only Affymtrix microarray data"}
#' }
#' @param experimentId Filter allowing to specify one or more ArrayExpress or GEO accession, e.g., 
#' GSE43721. Default is NULL: takes all available experiments for targeted species and data type.
#' @param sampleId Filter allowing to specify one or more sample ID. Depending on the selected 
#' datatype this sample IDs can correspond to Chip IDs (affymetrix) or RNA-Seq library IDs (rna_seq). 
#' Default is NULL: takes all available samples for targeted species and data type.
#' @param anatEntityId Filter allowing to specify one or more anatomical entity IDs from the UBERON 
#' ontology (http://uberon.github.io/). Default is NULL: takes all available anatomical entities for 
#' targeted species and data type.
#' @param stageId Filter allowing to specify one or more developmental stage IDs from Developmental 
#' Stage Ontology (https://github.com/obophenotype/developmental-stage-ontologies). Default is 
#' NULL: takes all available developmental stages for targeted species and data type.
#' @param cellTypeId Filter specific to single cell datatype (sc_full_length) allowing to specify 
#' one or more cell type IDs from the UBERON ontology (http://uberon.github.io/). Default is 
#' NULL: takes all available cell types for targeted species and data type. Available for Bgee 15.0 and after
#' @param sex Filter allowing to specify one or more sexes. Default is 
#' NULL: takes all available sexes for targeted species and data type. Available for Bgee 15.0 and after
#' @param strain Filter allowing to specify one or more strains. Default is 
#' NULL: takes all available strains for targeted species and data type. Available for Bgee 15.0 and after
#' @param withDescendantAnatEntities Allows to filter on the selected anatEntityId and all its descendants.
#' This functionality is available for Bgee 15.0 release and after
#' @param withDescendantStages Allows to filter on the selected stageId and all its descendants.
#' This functionality is available for Bgee 15.0 release and after
#' @param withDescendantCellTypes Allows to filter on the selected cellTypeId and all its descendants.
#' This functionality is available for Bgee 15.0 release and after
#' @param sampleProcessedOutputFile path to the file containing sample processed data.
#' @export

getSampleProcessedData_galaxy <- 
  function(species = NULL, dataType = character(0),
           experimentId = NULL, sampleId = NULL, anatEntityId = NULL, stageId = NULL,
           cellTypeId = NULL, sex = NULL, strain = NULL, withDescendantAnatEntities = FALSE,
           withDescendantStages = FALSE, withDescendantCellTypes = FALSE,
           sampleProcessedOutputFile = NULL) {
    if (is.null(sampleProcessedOutputFile)) {
      stop("path to the output file should always be provided")
    }
    myBgeeObject <- createBgeeObject_galaxy(species = species, dataTypes = dataType)
    sample_processed_data <-
      getSampleProcessedData(myBgeeObject = myBgeeObject,experimentId = experimentId,
                             sampleId = sampleId, anatEntityId = anatEntityId, stageId = stageId,
                             cellTypeId = cellTypeId, sex = sex, strain = strain,
                             withDescendantAnatEntities = withDescendantAnatEntities,
                             withDescendantStages = withDescendantStages,
                             withDescendantCellTypes = withDescendantCellTypes)
    write.table(x = sample_processed_data, file = sampleProcessedOutputFile, sep = "\t",
                quote = FALSE, row.names = FALSE)
  }

#' @title Run a GO-like enrichment of anatomical terms, mapped to genes by expression patterns
#'
#' @description This function GO-like enrichment of anatomical terms using the Uberon ontology.
#' 
#' @param species name of the species for which annotation have to be be retrieved (e.g Homo_sapiens)
#' @param dataTypes The datatypes for which calls will to be retrieved. To be chosen among:
#' \itemize{
#'   \item{"rna_seq: target only bulk RNA-seq data"}
#'   \item{"sc_full_length: target only full length single-cell RNA-seq data"}
#'   \item{"sc_droplet_based: target only droplet based single-cell RNA-seq data"}
#'   \item{"affymetrix: target only Affymtrix microarray data"}
#' }
#' @param foregroundGenes the list of genes for which TopAnat will find anatomical entities that
#' have over or under-represented expression using annotations for that gene set, compared to the
#' background genes
#'
#' @param backgroundGenes the list of genes you want to consider as the universe in your analysis
#'
#' @param nodeSize Minimum number of genes mapped to a node for it to be tested. Default is 10.
#' 
#' @param algorithm Decorrelation algorithm used to take into account the topology of the anatomical
#' ontology (default = weight). By default `classic` is used. The full list of algorithm can be
#' retrieved with topGO::whichAlgorithms()
#' 
#' @param statistics Character string specifing which test to use. By default `fisher`is used. The full list
#' of tests can be retrieved with topGO::whichTests()
#' 
#' @param confidence A character indicating if only high quality present calls should be
#' retrieved. For Bgee releases prior to 14, options are "all" (default) or "high_quality".
#' For Bgee release 14 and above, options are "silver" (default) and "gold".
#' 
#' @param resultFile The file where GO-like enrichment of anatomical terms is stored
#'
#' @author Julien Wollbrett
#'
#' @examples{
#' foregroundGenes <- c("ENSBTAG00000000011","ENSBTAG00000000014","ENSBTAG00000000016",
#' "ENSBTAG00000000026","ENSBTAG00000000039","ENSBTAG00000000040",
#' "ENSBTAG00000000042","ENSBTAG00000000050","ENSBTAG00000000056",
#' "ENSBTAG00000000064","ENSBTAG00000000067","ENSBTAG00000000071",
#' "ENSBTAG00000000072","ENSBTAG00000000080","ENSBTAG00000000081")
#' backgroundGenes <- c("ENSBTAG00000000011","ENSBTAG00000000014","ENSBTAG00000000016",
#' "ENSBTAG00000000026","ENSBTAG00000000039","ENSBTAG00000000040",
#' "ENSBTAG00000000042","ENSBTAG00000000050","ENSBTAG00000000056",
#' "ENSBTAG00000000064","ENSBTAG00000000067","ENSBTAG00000000071",
#' "ENSBTAG00000000072","ENSBTAG00000000080","ENSBTAG00000000081",
#' "ENSBTAG00000000084","ENSBTAG00000000091","ENSBTAG00000000099",
#' "ENSBTAG00000000111","ENSBTAG00000000123","ENSBTAG00000000132",
#' "ENSBTAG00000000153","ENSBTAG00000000162","ENSBTAG00000000163",
#' "ENSBTAG00000000169","ENSBTAG00000000179","ENSBTAG00000000197",
#' "ENSBTAG00000000199","ENSBTAG00000000202","ENSBTAG00000000203",
#' "ENSBTAG00000000204","ENSBTAG00000000213","ENSBTAG00000000215",
#' "ENSBTAG00000000223","ENSBTAG00000000224","ENSBTAG00000000225",
#' "ENSBTAG00000000236","ENSBTAG00000000250","ENSBTAG00000000251",
#' "ENSBTAG00000000252","ENSBTAG00000000253","ENSBTAG00000000261",
#' "ENSBTAG00000000274","ENSBTAG00000000277","ENSBTAG00000000279",
#' "ENSBTAG00000000285","ENSBTAG00000000286","ENSBTAG00000000287",
#' "ENSBTAG00000000289","ENSBTAG00000000297","ENSBTAG00000000305",
#' "ENSBTAG00000000312","ENSBTAG00000000328","ENSBTAG00000000335",
#' "ENSBTAG00000000341","ENSBTAG00000000343","ENSBTAG00000000354",
#' "ENSBTAG00000000355","ENSBTAG00000000356","ENSBTAG00000000365",
#' "ENSBTAG00000000372","ENSBTAG00000000379","ENSBTAG00000000380",
#' "ENSBTAG00000000382","ENSBTAG00000000396","ENSBTAG00000000404",
#' "ENSBTAG00000000405","ENSBTAG00000000406","ENSBTAG00000000411",
#' "ENSBTAG00000000425","ENSBTAG00000000434","ENSBTAG00000000435",
#' "ENSBTAG00000000438","ENSBTAG00000000448","ENSBTAG00000000451",
#' "ENSBTAG00000000454","ENSBTAG00000000456","ENSBTAG00000000457",
#' "ENSBTAG00000000459","ENSBTAG00000000462","ENSBTAG00000000469",
#' "ENSBTAG00000000470","ENSBTAG00000000484","ENSBTAG00000000497",
#' "ENSBTAG00000000501","ENSBTAG00000009707","ENSBTAG00000026266",
#' "ENSBTAG00000021992","ENSBTAG00000005353","ENSBTAG00000005333",
#' "ENSBTAG00000006424","ENSBTAG00000026972","ENSBTAG00000010799",
#' "ENSBTAG00000014614","ENSBTAG00000045757","ENSBTAG00000046332",
#' "ENSBTAG00000008394")
#' topAnat_galaxy(species = "Bos_taurus", stageId = "UBERON:0000092",
#'               foregroundGenes = foregroundGenes, backgroundGenes = backgroundGenes,
#'               resultFile = "result.tsv")
#' }
#'
#' @import topGO graph
#' @export
#' 
topAnat_galaxy <- function(species = NULL, dataTypes = character(0), stageId = NULL,
                           foregroundGenes = NULL, backgroundGenes = NULL, algorithm = "classic",
                           statistics = "fisher", resultFile = NULL, nodeSize = 10,
                           confidence = "silver") {
  myBgeeObject <- createBgeeObject_galaxy(species = species, dataTypes = dataTypes)
  myTopAnatData <- loadTopAnatData(myBgeeObject = myBgeeObject, stage = stageId)
  if (is.null(backgroundGenes) || is.null(foregroundGenes)) {
    stop("foreground and background can not be null")
  }
  # Ensure background genes are unique
  #TODO : should allow to retrieve all geneIDs of one species to automatically generate a background
  backgroundGenes <- unique(backgroundGenes)
  if (length(backgroundGenes[! foregroundGenes %in% backgroundGenes]) > 0) {
    stop("All foreground genes should be part of the background")
  }
  # create the geneList vector used as input of topAnat
  geneList <- rep(0, length(backgroundGenes))
  names(geneList) <- backgroundGenes
  geneList[names(geneList) %in% foregroundGenes] <- 1
  geneList <- as.factor(geneList)
  myTopAnatObject <- NULL
  tryCatch(
    {
      myTopAnatObject <- topAnat(topAnatData = myTopAnatData, geneList = geneList, nodeSize = nodeSize)
    },
    error = function(e) {
      stop("Did not manage to run topAnat. You probably did not provide proper gene list. It",
           "has to be gene IDs and not gene names.", conditionMessage(e))
    }
  )
  results <- runTest(object = myTopAnatObject, algorithm = algorithm, statistic = statistic)
  tableOver <- makeTable(myTopAnatData, myTopAnatObject, results)
  write.table(x = tableOver, file = resultFile, quote = FALSE, sep = "\t", row.names = FALSE)
}

#' @title Retrieve Bgee calls from one species.
#'
#' @description Loads the integrated expression calls from one species. These calls are 
#' equivalent to the one provided in the gene page of the Bgee website. The calls have been
#' generated using all datatypes present in Bgee (EST, In Situ, Affymetrix, bulk RNA-Seq and 
#' single-cell RNA-Seq).
#'
#' @param species name of the species for which annotation have to be be retrieved (e.g Homo_sapiens)
#' @param conditionParameters Specify which condition parameter you are interested
#' in. It can be `anatEntity` or `allCondParams`. The `anatEntity` option retrieve calls generated 
#' using only anatomical entity and cell-type as a condition parameter. The `allCondParams` option retrieve calls generated taking
#' into consideration all condition parameters present in Bgee (anat, entity, cell-type, developmental stage,
#' sex and strain). By default the `anatEntity` option is selected.
#' @param advancedColumns Boolean allowing to specify if advancend columns are required. The default value
#' is `FALSE` meaning that only gene , condition parameters, call quality, FDR, expression rank and
#' expression score are retrieved. If `TRUE` then a lot more information used to generate the calls are
#' retrieved like the number of self and descendant observation for each datatype, or the rank and score
#' per datatype. 
#' @param geneIds List of genes for which expression calls have to be retrieved. It has to be ensembl
#' IDs if the Bgee genome source is Ensembl (e.g ENSG00000244734) or RefSeq IDs if the Bgee genome source
#' is RefSeq (e.g 734881) but not gene names (e.g HBB, Apoc1, etc.).
#' @param anatEntityIds List of anatomical entity IDs for which expression calls have to be retrieved.
#' It has to be IDs (e.g UBERON:0000955) and not names (e.g brain)
#' @param callsOutputFile name of the file where intergrated calls will be stored
#' @export
#' 
getIntergratedCalls_galaxy <- function(species = NULL, conditionParameters = "anatEntity",
                                       advancedColumns = FALSE, geneIds = NULL,
                                       anatEntityIds = NULL, callsOutputFile = NULL) {
  myBgeeObject <- createBgeeObject_galaxy(species = species)
  integrated_calls <- getIntegratedCalls(myBgeeObject = myBgeeObject,
                                         conditionParameters = conditionParameters,
                                         advancedColumns = advancedColumns, geneIds = geneIds,
                                         anatEntityIds = =anatEntityIds)
  write.table(x = integrated_calls, file = callsOutputFile, quote = FALSE, sep = "\t",
              row.names = FALSE)
}

#' @title Retrieve Bgee processed expression values at cell level for single cell data.
#'
#' @description This function download the processed gene count sparse matrix of single cell experiments
#' in Bgee. These data are available for all droplet based single-cell RNA-seq and full length
#' single-cell RNA-seq. These files are stored in the H5AD format.
#'
#' @param species name of the species for which annotation have to be be retrieved (e.g Homo_sapiens)
#' @param experimentId Filter allowing to specify one ArrayExpress or GEO accession, e.g., 
#' GSE43721. Can not be null 
#' @param cellProcessedDataFile Name of the file where processed expression values at cell
#' level is stored
#' 
#' @author Julien Wollbrett
#' 
#' @export

downloadCellProcessedData_galaxy <- function(species = NULL, experimentId = NULL,
                                             cellProcessedDataFile = NULL) {
  myBgeeObject <- createBgeeObject_galaxy(species = species)
  destFile <- downloadCellProcessedFile(myBgeeObject = myBgeeObject,
                                        experimentId = experimentId)
  file.copy(from = destFile, to = cellProcessedDataFile)
  file.remove(destFile)
} 

#' @description
#' helper function to create a Bgee object for Galaxy wrappers. It allows to handle
#' errors in a cleaner way for Galaxy. Especially it lists the releases,
#' explains the pattern used for the species field and list available datatypes if
#' provided one does not exist
#' @noRd
#' @noMd
createBgeeObject_galaxy <- function(species = NULL, dataTypes = NULL) {
  # then check the species
  invisible(capture.output(bgeeSpecies <- listBgeeSpecies()))
  bgeeSpecies$formatedSpeciesName <- paste0(gsub(" ", "_", bgeeSpecies$GENUS),
                                            "_", gsub(" ", "_", bgeeSpecies$SPECIES_NAME))
  if (! species %in% bgeeSpecies$formatedSpeciesName) {
    stop("The species is not present in Bgee. the Species should follow the pattern",
    "`Genus_species` (e.g Homo_sapiens, Mus_musculus, Canis_lupus_familiaris)")
  }
  return(Bgee$new(species = species, dataType = dataTypes))
}

  
}
