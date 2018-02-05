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
#' @field pathToData Path to the directory where the data files are stored. By default the working directory is used. If many analyses are launched in parallel, please consider re-using the cached data files instead of redownlaoding them for each analysis.
#'
#' @field release Bgee release number to download data from, in the form "Release.subrelease" or "Release_subrelease", e.g., "13.2" or 13_2". Work for release >=13.2. By default, the latest relase of Bgee is used.
#'
#' @field sendStats A field specifying whether monitoring of users is performed for our internal usage statistics. This is useful to improve the settings of our servers and to get reliable usage statistics (e.g., when asking for funding for Bgee). No identification of the users is attempted, nor possible. Default to TRUE. This option can be set to FALSE, notably if all data files are in cache and that users want to be able to work offline.
#'
#' @field quantitativeData A field specifying if a single type of quantitative expression data ("rna_seq" or "affymetrix") was specified and if it is available for targeted species, helping the package to know if it should proceed with the execution of getAnnotation() and getData() functions.
#'
#' @examples{
#'  bgee <- Bgee$new(species = "Mus_musculus", dataType = "rna_seq")
#'  bgee <- Bgee$new(species = "Mus_musculus")
#' }
#'
#' @import methods Biobase RCurl data.table
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
    allExperimentsUrl = "character",
    topAnatUrl = "character",
    sendStats = "logical",
    quantitativeData = "logical",
    apiKey = "character"
  ),

  methods = list(
    initialize = function(...) {
      callSuper(...)

      ## check data type
      if (length(dataType) == 0) {
        cat("\nNOTE: You did not specify any data type. The argument dataType will be set to c(\"rna_seq\",\"affymetrix\",\"est\",\"in_situ\") for the next steps.\n")
        dataType <<- c("rna_seq","affymetrix","est","in_situ")
      } else if ( !sum(dataType %in% c("rna_seq","affymetrix","est","in_situ")) %in% 1:4 ){
        stop("ERROR: you need to specify at least one valid data type to be used among \"rna_seq\", \"affymetrix\", \"est\" and \"in_situ\".")
      }
      if ( length(dataType) != sum(dataType %in% c("rna_seq","affymetrix","est","in_situ")) ){
        cat("\nWARNING: you apparently specified a data type that is not among \"rna_seq\", \"affymetrix\", \"est\" and \"in_situ\". This will be removed. Please check for typos.\n")
        dataType <<- dataType[dataType %in% c("rna_seq","affymetrix","est","in_situ")]
      }


      ## check path of folder to store cached files
      if (length(pathToData) == 0) {
        pathToData <<- getwd()
      } else if (length(pathToData) == 1) {
        if (!file.exists(pathToData)) {
          stop("ERROR: please specify a valid and existing path to store data files.")
        }
      } else {
        stop("ERROR: Invalid path for data files.")
      }


      ## Get release information
      cat("\nQuerying Bgee to get release information...\n")
      allReleases <- try(.getBgeeRelease(removeFile = FALSE), silent=TRUE)
      if (class(allReleases) == "data.frame"){
        file.rename(from=file.path(getwd(), 'release.tsv'),
                    to=file.path(pathToData, "release.tsv"))
      } else if (class(allReleases) == "try-error"){
        if (file.exists(file.path(pathToData, "release.tsv"))){
          cat(paste0("\nWARNING: BgeeDB could not access Bgee releases information from the internet, but a release information file was found in the download directory ",
                     pathToData, ". This release file will be used, but be warned that it may not be up to date!\n"))
          allReleases <- read.table(file.path(pathToData, "release.tsv"), header=TRUE, sep="\t")
        } else {
          stop("ERROR: BgeeDB could not access Bgee releases information. Is your internet connection working?")
        }
      }

      if (length(release) == 0) {
        release <<- gsub("\\.", "_", format(allReleases$release, drop0Trailing = F)[1])
      } else if (length(release) == 1) {
        # In case the release number is written with a dot
        release <<- as.character(gsub("\\.", "_", release))
        if(grepl("_",release)){
          # test if required release exists
          # TODO same code than in the listSpeciesRelease
          if (sum(as.numeric(allReleases$release) == as.numeric(gsub("_", ".", release))) != 1) {
            stop("ERROR: The specified release number is invalid, or is not available for BgeeDB.")
          }
        } else {
          stop("Bgee release number must follow the format \"Release.subrelease\" or \"Release_subrelease\", e.g., \"13.2\" or \"13_2\". See function listBgeeRelease() to see available releases.")
        }
      } else {
        stop("ERROR: The specified release number is invalid.")
      }


      ## Specify URL to be used for topAnat. Can be done for any species and data type
      topAnatUrl <<-  as.character(allReleases$TopAnat.URL[as.numeric(allReleases$release) == as.numeric(gsub("_", ".", release))])
      if ( !grepl("/$", topAnatUrl) ){
        topAnatUrl <<- paste0(topAnatUrl, "/")
      }


      ## Get species information
      if (file.exists(paste0(pathToData, "/species_Bgee_", release, ".tsv"))){
        cat(paste0("\nNOTE: the file describing Bgee species information for release ", release, " was found in the download directory ",
                   pathToData, ". Data will not be redownloaded.\n"))
        allSpecies <- read.table(paste0(pathToData, "/species_Bgee_", release, ".tsv"),
                                 header=TRUE,
                                 sep="\t",
                                 blank.lines.skip=TRUE,
                                 as.is=TRUE)
      } else {
        allSpecies <- listBgeeSpecies(release = release,
                                      allReleases = allReleases,
                                      removeFile = FALSE)
      }

      ## check species argument (compulsory, no default value)
      if (length(species) == 0) {
        stop("ERROR: You did not specify any species.")
      } else if (length(species) > 1) {
        stop("ERROR: only one species is allowed.")
      } else if (grepl("^\\d+$", species)) {
        ## if species was specified as a taxonomic ID
        if (sum(allSpecies$ID == species) != 1) {
          stop(paste0("ERROR: The specified species taxonomic Id is invalid, or not available in Bgee release ",
                      release, "."))
        } else {
          speciesId <<- as.numeric(species)
          speciesName <<- gsub(" ", "_", paste(allSpecies[allSpecies$ID == species, 2:3], collapse = "_"))
        }
      } else {
        speciesSplitted <- unlist(strsplit(species, split = "_"))
        if (sum(allSpecies$GENUS == speciesSplitted[1] &
                allSpecies$SPECIES_NAME == speciesSplitted[2]) != 1) {
          stop(paste0("ERROR: The specified species name is invalid, or not available in Bgee release ",
                      release, "."))
        } else {
          speciesName <<- gsub(" ", "_", species)
          speciesId <<- as.numeric(allSpecies$ID[allSpecies$GENUS == speciesSplitted[1] &
                                                   allSpecies$SPECIES_NAME == speciesSplitted[2]])
        }
      }


      ## Set quantitativeData field to FALSE by default
      quantitativeData <<- FALSE

      ## For quantitative data download, there should be only 1 data type specified: rna_seq or affymetrix
      if (length(dataType) == 1){
        ## If only RNA-seq data, check availability in species
        if (dataType == "rna_seq") {
          quantitativeData <<- allSpecies$RNA_SEQ[allSpecies$ID == speciesId]
        }
        ## If only Affymetrix data, check availability in species
        else if (dataType == "affymetrix") {
          quantitativeData <<- allSpecies$AFFYMETRIX[allSpecies$ID == speciesId]
        }

        ## if quantitative data download can be done, fill the URLs fields
        if (quantitativeData == TRUE){
          if (dataType == "rna_seq") {
            ## annotation file
            annotationUrl <<- as.character(
              allReleases$RNA.Seq.annotation.URL.pattern[as.numeric(allReleases$release) == as.numeric(gsub("_", ".", release))])
            ## Data from specific experiment
            experimentUrl <<- as.character(
              allReleases$RNA.Seq.experiment.value.URL.pattern[as.numeric(allReleases$release) == as.numeric(gsub("_", ".", release))])
            ## Data from all experiments
            allExperimentsUrl <<- as.character(
              allReleases$RNA.Seq.all.value.URL.pattern[as.numeric(allReleases$release) == as.numeric(gsub("_", ".", release))])
          } else if (dataType == "affymetrix") {
            ## annotation file
            annotationUrl <<- as.character(
              allReleases$Affymetrix.annotation.URL.pattern[as.numeric(allReleases$release) == as.numeric(gsub("_", ".", release))])
            ## Data from specific experiment
            experimentUrl <<- as.character(
              allReleases$Affymetrix.experiment.value.URL.pattern[as.numeric(allReleases$release) == as.numeric(gsub("_", ".", release))])
            ## Data from all experiments
            allExperimentsUrl <<- as.character(
              allReleases$Affymetrix.all.value.URL.pattern[as.numeric(allReleases$release) == as.numeric(gsub("_", ".", release))])
          }

          ## Regex substitution to get the correct URLs
          annotationUrl <<- gsub("SPECIESNAMEPATTERN", speciesName, annotationUrl)
          allExperimentsUrl <<- gsub("SPECIESNAMEPATTERN", speciesName, allExperimentsUrl)
          experimentUrl <<- gsub("SPECIESNAMEPATTERN", speciesName, experimentUrl)
          ## Note: one more substitution is needed here for the experiment id. This is done in the getData() function.
        }
      }


      ## create sub-folder with species name to store downloaded files
      pathToData <<- paste0(pathToData, "/", speciesName, "_Bgee_", release)
      if (!file.exists(pathToData)) {
        dir.create(pathToData)
      }


      ## sendStats field
      if (length(sendStats) == 0) {
        sendStats <<- TRUE
      } else if (length(sendStats) != 1 | (length(sendStats) == 1 & !is.logical(sendStats))) {
        cat("\nNOTE: You did not specify a valid value for the \"sendStats\" field (should be TRUE or FALSE). The field will be set to TRUE for the next steps.\n")
        sendStats <<- TRUE
      }


      ## Create a concatenated string made of several variables that should be unique to the user, using Sys.info(), Sys.getenv() and R.version variables.
      myUserInfo <- c(Sys.info()[c("sysname", "release", "version", "machine", "login")], Sys.getenv()[c("R_HOME", "R_LIBS_USER")], R.version)
      ## Use library digest to create a SHA512 hash, that will be used as API key
      apiKey <<- digest(myUserInfo, algo = "sha512")
      cat(paste0("\nAPI key built: ", apiKey, "\n"))

      ## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ## The API key allows to monitor usage of our package and to limit the number of simultaneous queries
      ## to the server from the same user. It is a secure hash that does not allow the identification of
      ## the user.
      ## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      ## Message to users. Can for example be used to encourage users to update their version of the package.
      messageToUsers <- as.character(allReleases$messageToUsers[as.numeric(allReleases$release) == as.numeric(gsub("_", ".", release))])
      if(!(is.na(messageToUsers) || is.null(messageToUsers))){
        cat(paste0("\nIMPORTANT INFORMATION: ", messageToUsers, "\n"))
      }

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




