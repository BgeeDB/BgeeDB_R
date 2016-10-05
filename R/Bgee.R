########################
#' @title Bgee Reference Class
#'
#' @description This is used to specify information at the beginning of a BgeeDB working session, for example, the targeted species and data type. An object of this class is then passed as argument to other functions of the package to provide these informations. See examples in vignette.
#'
#' @details Bgee (\url{http://bgee.org}) integrates different expression data types (RNA-seq, Affymetrix microarray, ESTs, and in-situ hybridizations) from multiple animal species. Expression patterns are based exclusively on curated "normal", healthy, expression data (e.g., no gene knock-out, no treatment, no disease), to provide a reference atlas of normal gene expression.
#'
#' @field species A character indicating the species to be used, in the form "Genus_species", or a numeric indicating the species NCBI taxonomic id. Only species with data in Bgee will work. See the listBgeeSpecies() function to get the list of species available in the Bgee release used.
#'
#' @field dataType A vector of characters indicating data type(s) to be used. To be chosen among:
#' \itemize{
#'   \item{"rna_seq"}
#'   \item{"affymetrix"}
#'   \item{"est"}
#'   \item{"in_situ"}
#' }
#' By default all data type are included: \code{c("rna_seq","affymetrix","est","in_situ")}. For download of quantitative expression data, a single data type should be chosen among "rna_seq" or 'affymetrix".
#'
#' @field pathToData Path to the directory where the data files are stored. By default the working directory is used.
#'
#' @field release Bgee release number to download data from, in the form "Release.subrelease" or "Release_subrelease", e.g., "13.2" or 13_2". Will work for release >=13.2. By default, the latest relase of Bgee is used.
#'
#' @field useApiKey A field specifying if users should be tracked with an API key for our internal usage statistics and for the management of queries to server. Default to TRUE. The API key created is a hash that does not allow formal identification of the user. If this is still not fitting your needs, useApiKey can be set to FALSE and the package will work anonymously. In this case, please be careful not to launch too many queries in parallel and try to reuse cached data files (see pathToData argument) as much as possible.
#'
#' @field quantitativeData A field specifying if a single type of quantitative expression data ("rna_seq" or "affymetrix") was specified and if it is available for targeted species, helping the package to know if it should proceed with the execution of getAnnotation() and getData() functions.
#'
#' @examples{
#'  bgee <- Bgee$new(species = "Mus_musculus", dataType = "rna_seq")
#'  bgee <- Bgee$new(species = "Mus_musculus")
#' }
#'
#' @import methods Biobase RCurl data.table
#' ## TO DO: what is methods? data.table needed here? RCurl still needed?
#' @export Bgee
#' @exportClass Bgee

Bgee <- setRefClass(
  "Bgee",

  fields = list(
    species = "character",
    speciesName = "character",
    speciesId = "numeric",
    dataType = "character",
    pathToData = "character",
    release = "character",
    annotationUrl = "character",
    experimentUrl = "character",
    allexperimentsUrl = "character",
    useApiKey = "boolean",
    quantitativeData = "boolean"
    ## TO DO: is "boolean" OK?
  ),

  methods = list(
    initialize = function(...) {
      callSuper(...)
      ## TO DO: what is this doing? Is it needed?

      ## check data type
      if (length(dataType) == 0) {
        stop("ERROR: You didn't specify a data type. Choose 'affymetrix' or 'rna_seq'.")
      } else if ((length(dataType) > 1) ||
                 (length(dataType) == 1 &&
                  dataType %in% c("rna_seq", "affymetrix") == "FALSE")) {
        stop("ERROR: Choose correct dataType argument: 'affymetrix' or 'rna_seq'.")
      }

      ## TO DO: Set to c("rna_seq","affymetrix","est","in_situ") by default
      ## TO DO: Initial code from loadTopAnatData.R. Check if fully redundant or integrate
      ## Test if parameters are in the range of allowed parameters
      # if ( !sum(dataType %in% c("rna_seq","affymetrix","est","in_situ")) %in% 1:4 ){
      #   stop("ERROR: you need to specify at least one valid data type to be used among \"rna_seq\", \"affymetrix\", \"est\" and \"in_situ\".")
      # }
      # if ( length(dataType) != sum(dataType %in% c("rna_seq","affymetrix","est","in_situ")) ){
      #   cat("WARNING: you apparently specified a data type that is not among \"rna_seq\", \"affymetrix\", \"est\" and \"in_situ\". Please check for mistakes or typos.\n")
      # }

      cat("Querying Bgee to get release information...\n")
      allReleases <- .getBgeeRelease()
      if (length(release) == 0) {
        release <<- gsub("\\.", "_", allReleases$release[1])
      } else if (length(release) == 1) {
        # In case the release number is written with a dot
        release <<- gsub("\\.", "_", release)
        # test if required release exists
        if (sum(allReleases$release == gsub("_", ".", release)) != 1) {
          stop("ERROR: The specified release number is invalid, or is not available for BgeeDB.")
        }
      } else {
        stop("ERROR: The specified release number is invalid.")
      }

      ## TO DO: Initial code from loadTopAnatData.R. Check if fully redundant or integrate
      # cat("Querying Bgee to get release information...\n")
      # allReleases <- .getBgeeRelease()
      # if (length(release)==0) {
      #   release <- gsub("\\.", "_", allReleases$release[1])
      # } else if (length(release)==1){
      #   # In case the release number is written with a dot
      #   release <- gsub("\\.", "_", release)
      #   # test if required release exists
      #   if (sum(allReleases$release == gsub("_", ".", release))!=1){
      #     stop("ERROR: The specified release number is invalid, or is not available for BgeeDB.")
      #   }
      # } else {
      #   stop("ERROR: The specified release number is invalid.")
      # }

      ## TO DO: store release.tsv at root of pathToData
      ##        if it is there, do not call .getBgeeRelease
      ##        (similar to listBgeeSpecies)

      ############## TO DO: move tis somewhere else? Next to other url sepcifications ############
      ## Specify host to be used
      host <- allReleases$TopAnat.URL[allReleases$release == gsub("_", ".", release)]
      if ( !grepl("/$", host) ){
        host <- paste0(host, "/")
      }
      ## TO DO: should that be specified in Bgee object ?
      ## replace by topAnat.url field or similar
      ########################

      ## First retrieve list of all species for queried release
      allSpecies <-
        listBgeeSpecies(release = release, allReleases = allReleases)

      ## check species argument
      if (length(species) == 0) {
        stop("ERROR: You didn't specify any species.")
      } else if (length(species) > 1) {
        stop("ERROR: only one species is allowed.")
      } else if (grepl("^\\d+$", species)) {
        ## of species was specified as a taxonomic ID
        if (sum(allSpecies$ID == species) != 1) {
          stop(
            paste0(
              "ERROR: The specified species Id is invalid, or not available in Bgee release ",
              release
            )
          )
        } else {
          speciesId <<- as.numeric(species)
          speciesName <<-
            paste(allSpecies[allSpecies$ID == species, 2:3], collapse = "_")
        }
      } else {
        speciesSplitted <- unlist(strsplit(species, split = "_"))
        if (sum(allSpecies$GENUS == speciesSplitted[1] &
                allSpecies$SPECIES_NAME == speciesSplitted[2]) != 1) {
          stop(
            paste0(
              "ERROR: The specified species name is invalid, or not available in Bgee release ",
              release,
              "."
            )
          )
        } else {
          speciesName <<- species
          speciesId <<-
            as.numeric(allSpecies$ID[allSpecies$GENUS == speciesSplitted[1] &
                                       allSpecies$SPECIES_NAME == speciesSplitted[2]])
        }
      }

      ## TO DO: Initial code from loadTopAnatData.R. Check if fully redundant or integrate
      # ## Retrieve list of all species for queried release
      # allSpecies <- listBgeeSpecies(release=release, allReleases=allReleases)
      #
      # ## Species is the only compulsory parameter
      # if( length(species) == 0 ) {
      #   stop("ERROR: you need to specify a species.")
      # } else if ( length(species) > 1 ){
      #   stop("ERROR: only one species is allowed.")
      # } else if (grepl("^\\d+$", species)){
      #   ## of species was specified as a taxonomic ID
      #   if (sum(allSpecies$ID == species) != 1){
      #     stop(paste0("ERROR: The specified species Id is invalid, or not available in Bgee release ", release))
      #   } else {
      #     speciesId <- as.numeric(species)
      #     speciesName <- paste(allSpecies[allSpecies$ID == species, 2:3], collapse="_")
      #   }
      # } else if (is.character(species)){
      #   speciesSplitted <- unlist(strsplit(species, split="_"))
      #   if (sum(allSpecies$GENUS == speciesSplitted[1] & allSpecies$SPECIES_NAME == speciesSplitted[2]) != 1){
      #     stop(paste0("ERROR: The specified species name is invalid, or not available in Bgee release ", release, "."))
      #   } else {
      #     speciesName <- species
      #     speciesId <- as.numeric(allSpecies$ID[allSpecies$GENUS == speciesSplitted[1] & allSpecies$SPECIES_NAME == speciesSplitted[2]])
      #   }
      # }
      #





      ## check data type availability for chosen species
      if (dataType == "rna_seq" &
          allSpecies$RNA_SEQ[allSpecies$ID == speciesId] == FALSE) {
        stop(
          "ERROR: The specified species name is not among the list of species with RNA-seq data in Bgee release ",
          release,
          ". See listBgeeSpecies() for details on data types availability for each species."
        )
      } else if (dataType == "affymetrix" &
                 allSpecies$AFFYMETRIX[allSpecies$ID == speciesId] == FALSE) {
        stop(
          "ERROR: The specified species name is not among the list of species with Affymetrix microarray data in Bgee release ",
          release,
          ". See listBgeeSpecies() for details on data types availability for each species."
        )
      }

      ## create URLs
      if (dataType == "rna_seq") {
        ## annotation file
        annotationUrl <<-
          as.character(allReleases$RNA.Seq.annotation.URL.pattern[allReleases$release == gsub("_", ".", release)])
        ## Data from specific experiment
        experimentUrl <<-
          as.character(allReleases$RNA.Seq.experiment.value.URL.pattern[allReleases$release == gsub("_", ".", release)])
        ## Data from all experiments
        allexperimentsUrl <<-
          as.character(allReleases$RNA.Seq.all.value.URL.pattern[allReleases$release == gsub("_", ".", release)])
      } else if (dataType == "affymetrix") {
        ## annotation file
        annotationUrl <<-
          as.character(allReleases$Affymetrix.annotation.URL.pattern[allReleases$release == gsub("_", ".", release)])
        ## Data from specific experiment
        experimentUrl <<-
          as.character(allReleases$Affymetrix.experiment.value.URL.pattern[allReleases$release == gsub("_", ".", release)])
        ## Data from all experiments
        allexperimentsUrl <<-
          as.character(allReleases$Affymetrix.all.value.URL.pattern[allReleases$release == gsub("_", ".", release)])
      }
      annotationUrl <<-
        gsub("SPECIESNAMEPATTERN", speciesName, annotationUrl)
      allexperimentsUrl <<-
        gsub("SPECIESNAMEPATTERN", speciesName, allexperimentsUrl)
      experimentUrl <<-
        gsub("SPECIESNAMEPATTERN", speciesName, experimentUrl)
      ## Note: one more substitution is needed here for the experiment id. This is done later in the get_data() function.

      ## check path of folder to store cached files
      if (length(pathToData) == 0) {
        pathToData <<- paste0(getwd(), "/", speciesName, "_Bgee_", release)
      } else if (length(pathToData) == 1) {
        if (!file.exists(pathToData)) {
          stop("ERROR: please specify a valid and existing path to store data files.")
        } else {
          pathToData <<-
            paste0(pathToData, "/", speciesName, "_Bgee_", release)
        }
      } else {
        stop("ERROR: Invalid path for data files.")
      }
      ## create sub-folder with species name to store downloaded files
      if (!file.exists(pathToData)) {
        dir.create(pathToData)
      }

      ## TO DO: Initial code from loadTopAnatData.R. Check if fully redundant or integrate
      # ## check path of folder to store cached files
      # if(length(pathToData)==0) {
      #   pathToData <- paste0(getwd(), "/", speciesName, "_Bgee_", release)
      # } else if (length(pathToData)==1){
      #   if ( !file.exists(pathToData) ){
      #     stop("ERROR: please specify a valid and existing path to store data files.")
      #   } else {
      #     pathToData <- paste0(pathToData, "/", speciesName, "_Bgee_", release)
      #   }
      # } else {
      #   stop("ERROR: Invalid path for data files.")
      # }
      # ## create sub-folder with species name to store downloaded files
      # if (!file.exists(pathToData)){
      #   dir.create(pathToData)
      # }

      ## TO DO: Set quantitative Data field to TRUE or FALSE to say if species/dataType is OK for getAnnotation and getData.

      ## TO DO: create API key here? add field to Bgee object?

    },
    get_annotation = function(...){
      stop("ERROR: this function is deprecated. Use getAnnotation() function instead.")
    },
    get_data = function(...){
      stop("ERROR: this function is deprecated. Use getData() function instead.")
    },
    format_data = function(...){
      stop("ERROR: this function is deprecated. Use formatData() function instead.")
    }
  )
)




