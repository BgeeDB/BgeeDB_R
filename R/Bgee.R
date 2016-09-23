# Takes in unformatted data downloaded from Bgee and outputs a list of expression matrices, phenotype annotations, feature annotations, and calls.
.extract.data <- function(data, columns, calltype){
    # if multiple data frames (multiple experiments or chips)
    if(class(data) == "list"){
        calls <- lapply(data, function(x) .calling(x, calltype, columns[3]))
        expr <- lapply(calls,
                       function(x) {
                         # subset the data to keep relevant columns
                         xt <- x[, columns[1:3]]
                         # from sample and feature columns, create a matrix with features as rows and samples as columns
                         xtt <- xt %>% spread_(columns[1], columns[3])
                         rownames(xtt) <- xtt[,columns[2]]
                         # Remove feature column to keep only data
                         xtt[,-1, drop = FALSE]
                       }
                     )
        cat("Extracting features...\n")
        features <- mapply(.extract.data.feature, calls, expr, rep(list(columns), times=length(calls)))
        cat("Done...\n")
        cat("Extracting pheno...\n")
        phenos <- mapply(.extract.data.pheno, calls, rep(list(columns[1]), times=length(calls)))
        cat("Done...\n")
    } else {
        # if only a single dataframe
        calls <- .calling(data, calltype, columns[3])
        # subset the data to keep relevant columns
        xt <- calls[, columns[1:3]]
        xtt <- xt %>% spread_(columns[1], columns[3])
        rownames(xtt) <- xtt[,columns[2]]
        # Remove feature column to keep only data
        expr <- xtt[,-1, drop = FALSE]

        cat("Extracting features...\n")
        features <- .extract.data.feature( calls, expr, columns)
        cat("Done...\n")
        cat("Extracting pheno...\n")
        phenos <- .extract.data.pheno( calls, columns[1])
        cat("Done...\n")
    }
    return(list(assayData = expr, pheno = phenos, features = features, calls = calls))
}

# Extract feature data (probesets or genes)
.extract.data.feature <- function(calls, expr, columns){
    # RNA-seq
    if(length(columns) == 3){
        fdata <- calls[match(rownames(expr), calls[, columns[2]]), columns[2], drop = FALSE]
    }
    # Affymetrix, 4 columns
    else if(length(columns) == 4){
        fdata <- calls[match(rownames(expr), calls[, columns[2]]), columns[c(2,4)], drop = FALSE]
    }
    rownames(fdata) <- fdata[, columns[2]]
    fdata <- as(fdata, "AnnotatedDataFrame")
    return(fdata)
}

# Extract annotation of samples
.extract.data.pheno <- function(calls, column){
    phdata <- calls[, c(column, "Anatomical.entity.ID", "Anatomical.entity.name", "Stage.ID", "Stage.name")]
    phdata <- phdata[!duplicated(phdata[,column]), ]
    rownames(phdata) <- phdata[, column]
    phdata <- as.data.frame(phdata)
    metadata <- data.frame(labelDescription=colnames(phdata),
                           row.names=colnames(phdata))
    phenodata <- new("AnnotatedDataFrame",
                     data=phdata,
                     varMetadata=metadata
                     )
    return(phenodata)
}

.calling <- function(x, calltype, column){
    ## check datatype
    if(calltype == "present"){
        cat("keeping only present genes...\n")
        x[(x$Detection.flag == "absent"), column] <- NA
    } else if (calltype == "present high quality"){
        cat("keeping only present high quality genes...\n")
        x[which(x$Detection.flag == "absent" | x$Detection.quality == "poor quality"), column] <- NA
    }
    return(x)
}


########################
#' @title Retrieving the Bgee database data
#' @description A Reference Class to give annotation available on Bgee for particular species and the requested data (rna_seq, affymetrix)
#'
#' @details The expression data come from Bgee (http://.bgee.org), that integrates different expression data types (RNA-seq, Affymetrix microarray, ESTs, or in-situ hybridizations) in multiple animal species. Expression patterns are based exclusively on curated "normal", healthy, expression data (e.g., no gene knock-out, no treatment, no disease), to provide a reference of normal gene expression.
#' This Class retrieves annotation of experiments in Bgee database (get_annotation), downloads the data (get_data), and formats the data into expressionSet object (format_data). See examples and vignette.
#'
#' @field species A character indicating the species to be used, in the form "Genus_species", or a numeric indicating the species NCBI taxonomic id. Only species with quantitative expression data in Bgee will work (RNA-seq and Affymetrix microarray data). See the listBgeeSpecies() function to get the list of species available in the Bgee release used for these two data types.
#'
#' @field datatype A character of data platform. Quantitative expression levels can be downloaded for two data types:
#' \itemize{
#'      \item{"rna_seq"}
#'      \item{"affymetrix"}}
#'
#' @field experiment.id An ArrayExpress or GEO accession, e.g., GSE30617
#' On default is NULL: takes all available experiments for specified species and datatype.
#'
#' @field pathToData Path to the directory where the data files are stored. By default the working directory is used.
#'
#' @field release Bgee release number to download data from, in the form "Release.subrelease" or "Release_subrelease", e.g., "13.2" or 13_2". Will work for release >=13.2. By default, the latest relase of Bgee is used.
#'
#' @field calltype A character.
#'  \itemize{
#'    \item{"present"}
#'    \item{"present high quality"}
#'    \item{"all"}}
#' Retrieve intensities only for present (expressed) genes, present high quality genes, or all genes. The default is present.
#'
#' @field stats A character indicating what expression values should be used in the formatted data expressionSet object.
#'  \itemize{
#'    \item{"rpkm" for RNA-seq}
#'    \item{"counts" for RNA-seq}
#'    \item{"tpm" for RNA-seq (for Bgee release 14 and above only)}
#'    \item{"intensities" for Affymetrix microarrays}
#'    }
#'
#' @return
#' \itemize{
#'  \item{\code{get_annotation()} returns a list of the annotation of experiments for chosen species.}
#'  \item{\code{get_data()}, if experiment ID is empty, returns a list of experiments. If specified experiment ID, then returns the dataframe of the chosen experiment}
#'  \item{\code{format_data()}, if experiment ID is empty, returns a list of ExpressionSet objects. If specified experiment ID, then returns an ExpressionSet object}}
#'
#' @author Andrea Komljenovic and Julien Roux.
#'
#' @examples{
#'  bgee <- Bgee$new(species = "Mus_musculus", datatype = "rna_seq")
#'  annotation_bgee_mouse <- bgee$get_annotation()
#'  data_bgee_mouse <- bgee$get_data()
#'  data_bgee_mouse_gse30617 <- bgee$get_data(experiment.id = "GSE30617")
#'  gene.expression.mouse.rpkm <- bgee$format_data(data_bgee_mouse_gse30617,
#'  calltype = "present", stats = "rpkm")
#'  }
#'
#'
#' @import methods Biobase RCurl data.table
#' @importFrom dplyr %>%
#' @importFrom tidyr spread
#' @export Bgee
#' @exportClass Bgee

Bgee <- setRefClass("Bgee",

fields = list(
  species = "character",
  speciesName = "character",
  speciesId = "numeric",
  datatype = "character",
  experiment.id = "character",
  pathToData = "character",
  release = "character",
  calltype = "character",
  stats = "character",
  annotationUrl = "character",
  experimentUrl = "character",
  allexperimentsUrl = "character"
),

methods = list(

initialize=function(...) {
  callSuper(...)

  ## check datatype
  if(length(datatype)==0) {
    stop("ERROR: You didn't specify a data type. Choose 'affymetrix' or 'rna_seq'.")
  } else if ((length(datatype) > 1) || (length(datatype) == 1 && datatype %in% c("rna_seq", "affymetrix") == "FALSE")){
    stop("ERROR: Choose correct datatype argument: 'affymetrix' or 'rna_seq'.")
  }

  cat("Querying Bgee to get release information...\n")
  allReleases <- .getRelease()
  if (length(release)==0) {
    release <<- gsub("\\.", "_", allReleases$release[1])
  } else if (length(release)==1){
    # In case the release number is written with a dot
    release <<- gsub("\\.", "_", release)
    # test if required release exists
    if (sum(allReleases$release == gsub("_", ".", release))!=1){
      stop("ERROR: The specified release number is invalid, or is not available for BgeeDB.")
    }
  } else {
    stop("ERROR: The specified release number is invalid.")
  }

  ## First retrieve list of all species for queried release
  allSpecies <- listBgeeSpecies(release=release, allReleases=allReleases)

  ## check species argument
  if(length(species)==0) {
    stop("ERROR: You didn't specify any species.")
  } else if ( length(species) > 1 ){
    stop("ERROR: only one species is allowed.")
  } else if (grepl("^\\d+$", species)){
    ## of species was specified as a taxonomic ID
    if (sum(allSpecies$ID == species) != 1){
      stop(paste0("ERROR: The specified species Id is invalid, or not available in Bgee release ", release))
    } else {
      speciesId <<- as.numeric(species)
      speciesName <<- paste(allSpecies[allSpecies$ID == species, 2:3], collapse="_")
    }
  } else {
    speciesSplitted <- unlist(strsplit(species, split="_"))
    if (sum(allSpecies$GENUS == speciesSplitted[1] & allSpecies$SPECIES_NAME == speciesSplitted[2]) != 1){
      stop(paste0("ERROR: The specified species name is invalid, or not available in Bgee release ", release, "."))
    } else {
      speciesName <<- species
      speciesId <<- as.numeric(allSpecies$ID[allSpecies$GENUS == speciesSplitted[1] & allSpecies$SPECIES_NAME == speciesSplitted[2]])
    }
  }

  ## check data type availability for chosen species
  if(datatype == "rna_seq" & allSpecies$RNA_SEQ[allSpecies$ID == speciesId] == FALSE){
    stop("ERROR: The specified species name is not among the list of species with RNA-seq data in Bgee release ", release,". See listBgeeSpecies() for details on data types availability for each species.")
  } else if(datatype == "affymetrix" & allSpecies$AFFYMETRIX[allSpecies$ID == speciesId] == FALSE){
    stop("ERROR: The specified species name is not among the list of species with Affymetrix microarray data in Bgee release ", release,". See listBgeeSpecies() for details on data types availability for each species.")
  }

  ## create URLs
  if(datatype == "rna_seq"){
    ## annotation file
    annotationUrl <<- as.character(allReleases$RNA.Seq.annotation.URL.pattern[allReleases$release == gsub("_", ".", release)])
    ## Data from specific experiment
    experimentUrl <<- as.character(allReleases$RNA.Seq.experiment.value.URL.pattern[allReleases$release == gsub("_", ".", release)])
    ## Data from all experiments
    allexperimentsUrl <<- as.character(allReleases$RNA.Seq.all.value.URL.pattern[allReleases$release == gsub("_", ".", release)])
  } else if (datatype == "affymetrix"){
    ## annotation file
    annotationUrl <<- as.character(allReleases$Affymetrix.annotation.URL.pattern[allReleases$release == gsub("_", ".", release)])
    ## Data from specific experiment
    experimentUrl <<- as.character(allReleases$Affymetrix.experiment.value.URL.pattern[allReleases$release == gsub("_", ".", release)])
    ## Data from all experiments
    allexperimentsUrl <<- as.character(allReleases$Affymetrix.all.value.URL.pattern[allReleases$release == gsub("_", ".", release)])
  }
  annotationUrl <<- gsub("SPECIESNAMEPATTERN", speciesName, annotationUrl)
  allexperimentsUrl <<- gsub("SPECIESNAMEPATTERN", speciesName, allexperimentsUrl)
  experimentUrl <<- gsub("SPECIESNAMEPATTERN", speciesName, experimentUrl)
  ## Note: one more substitution is needed here for the experiment id. This is done later in the get_data() function.

  ## check path of folder to store cached files
  if(length(pathToData)==0) {
    pathToData <<- paste0(getwd(), "/", speciesName, "_Bgee_", release)
  } else if (length(pathToData)==1){
    if ( !file.exists(pathToData) ){
      stop("ERROR: please specify a valid and existing path to store data files.")
    } else {
      pathToData <<- paste0(pathToData, "/", speciesName, "_Bgee_", release)
    }
  } else {
    stop("ERROR: Invalid path for data files.")
  }
  ## create sub-folder with species name to store downloaded files
  if (!file.exists(pathToData)){
    dir.create(pathToData)
  }
},


get_annotation = function(...){

  annotation.file <- basename(annotationUrl)

  ## To get the names of experiment and sample files, we start frome the annotation.file
  if (datatype == "affymetrix"){
    annotation.experiments <- gsub("_chips", "", annotation.file)
  } else if (datatype == "rna_seq"){
    annotation.experiments <- gsub("_libraries", "", annotation.file)
  }
  annotation.experiments <- gsub("zip", "tsv", annotation.experiments)
  annotation.samples     <- gsub("_experiments", "", annotation.file)
  annotation.samples     <- gsub("zip", "tsv", annotation.samples)

  ## Check if file is already in cache. If so, skip download step
  if (file.exists(file.path(pathToData, annotation.experiments)) && file.exists(file.path(pathToData, annotation.samples))){
    cat("WARNING: annotation files for this species were found in the download directory and will be used as is. Please delete the files and rerun the function, if you want the data to be updated.\n")
  } else {
    cat("Downloading annotation files...\n")

    success <- download.file(annotationUrl,
                             destfile=file.path(pathToData, annotation.file),
                             mode='wb')
    if (success != 0){
      stop("ERROR: Download from FTP was not successful.")
    }
    unzip(file.path(pathToData, annotation.file), exdir=pathToData)
    cat("Saved annotation files in", pathToData, "folder.\n")
    ## Clean directory and remove .zip file
    file.remove(file.path(pathToData, annotation.file))
    ## Test if extracted files are OK
    if (!(file.exists(file.path(pathToData, annotation.experiments)) & file.exists(file.path(pathToData, annotation.samples)))){
      stop("ERROR: extraction of annotation files from downloaded zip file went wrong.")
    }
  }

  ## Read the 2 annotation files
  myanno <- list(sample_annotation=as.data.frame(fread(file.path(pathToData, annotation.samples))),
                 experiment_annotation=as.data.frame(fread(file.path(pathToData, annotation.experiments)))
                 )
  ## remove spaces in headers
  for (i in 1:length(myanno)){
    names(myanno[[i]]) <- make.names(names(myanno[[i]]))
  }
  return(myanno)
},


get_data = function(..., experiment.id = NULL){

    if (length(experiment.id) == 0){
        cat(paste0("The experiment is not defined. Hence taking all ", datatype, " experiments available for ", speciesName, ".\n"))

        all_expression_values <- basename(allexperimentsUrl)

        ## check if RDS file already in cache. If so, skip download step
        if (file.exists(paste0(pathToData, "/", datatype, "_all_experiments_expression_data.rds"))){
            cat("WARNING: expression data file (.rds file) was found in the download directory and will be used as is. Please delete and rerun the function if you want the data to be updated.\n")
            data_all <- readRDS(file = paste0(pathToData, "/", datatype, "_all_experiments_expression_data.rds"))
        } else {
            cat("Downloading expression data...\n")

            success <- download.file(allexperimentsUrl,
                                     destfile=file.path(pathToData, all_expression_values),
                                     mode='wb')
            if (success != 0){
              stop("ERROR: Download from FTP was not successful.")
            }
            cat("Saved expression data file in", pathToData, "folder.\n")
            cat("Unzipping file...\n")
            temp.files <- unzip(file.path(pathToData, all_expression_values), exdir=pathToData)
            #print(temp.files)
            mydata <- lapply(temp.files, unzip, exdir=pathToData)
            #print(mydata)
            data_all <- lapply(unlist(mydata, rec = TRUE), function(x) as.data.frame(suppressWarnings(fread(x))))

            ## remove spaces in headers
            for (i in 1:length(data_all)){
              names(data_all[[i]]) <- make.names(names(data_all[[i]]))
            }

            cat("Saving all data in .rds file...\n")
            saveRDS(data_all, file = paste0(pathToData, "/", datatype, "_all_experiments_expression_data.rds"))

            ## cleaning up downloaded files
            try(file.remove(temp.files))
            try(file.remove(unlist(mydata, rec = TRUE)))
            try(file.remove(file.path(pathToData, all_expression_values)))
        }
    } else if( length(experiment.id) == 1){
        experiment.id <<- experiment.id
        if (!grepl("^GSE\\d+$|^E-\\w+-\\d+.*$", experiment.id, perl = TRUE)){
            stop("The experiment.id field needs to be a GEO or ArrayExpress accession, e.g., 'GSE30617' or 'E-MEXP-2011'")
        } else {
            cat("Downloading expression data for the experiment", experiment.id, "\n")

            ## Since experiment Id is defined now, we can substitute it in the URL
            experimentUrl <<- gsub("EXPIDPATTERN", experiment.id, experimentUrl)
            temp.file <- file.path(pathToData, basename(experimentUrl))

            ## check if RDS file already in cache. If so, skip download step
            if (file.exists(paste0(pathToData, "/", datatype, "_", experiment.id, "_expression_data.rds"))){
                cat("WARNING: expression data file (.rds file) was found in the download directory for", experiment.id, "and will be used as is. Please delete and rerun the function if you want the data to be updated.\n")
                data_all <- readRDS(paste0(pathToData, "/", datatype, "_", experiment.id, "_expression_data.rds"))
            } else {
                cat("Downloading expression data...\n")

                success <- download.file(experimentUrl,
                                         destfile=temp.file,
                                         mode='wb')
                if (success != 0){
                  stop("ERROR: Download from FTP was not successful.")
                }
                cat("Saved expression data file in", pathToData, "folder.\n")
                cat(paste0("Unzipping ", temp.file," file...\n"))
                # Unzipping this file can give one expression data file or multiple ones (if multiple chip types used in experiment)
                mydata <- unzip(temp.file, exdir=pathToData)
                data_all <- lapply(mydata, function(x) as.data.frame(fread(x)))

                ## remove spaces in headers
                for (i in 1:length(data_all)){
                  names(data_all[[i]]) <- make.names(names(data_all[[i]]))
                }

                if (length(data_all) == 1){
                  data_all <- as.data.frame(data_all[[1]])
                }
                cat("Saving all data in .rds file...\n")
                saveRDS(data_all, file = paste0(pathToData, "/", datatype, "_", experiment.id, "_expression_data.rds"))

                ## cleaning up downloaded files
                try(file.remove(temp.file))
                try(file.remove(mydata))
            }
        }
    } else {
        stop("Please provide only one experiment ID. If you want to get all data for this species and datatype, leave experiment.id empty")
    }


    return(data_all)
    cat("Done.")
},


format_data = function(data, calltype = "all", stats = NULL){
  if (datatype == "affymetrix" & stats != "intensities"){
    cat("WARNING: For Affymetrix microarray data, stats parameter can only be set to \"intensities\", this will be used in the next steps.\n")
    stats <<- "intensities"
  } else if (datatype == "rna_seq" & length(stats) == 0){
    stop("Please specify the stats parameters. Should be set to \"rpkm\" or \"counts\"")
  } else if (datatype == "rna_seq" & !(stats %in% c('rpkm', 'counts', 'tpm'))){
    ## TO DO: add TPM to documentation
    stop("Choose whether data formatting should create a matrix of RPKMs or read counts, with stats option set as \"rpkm\" or \"counts\"")
  }

  if(!(calltype %in% c('present','present high quality','all'))){
    stop("Choose between displaying intensities for present genes, present high quality genes or all genes, e.g., 'present', 'present high quality', 'all' ")
  }

  if(length(data) == 1) data[[1]] else data

  if(stats  == "rpkm"){
    columns <- c("Library.ID", "Gene.ID", "RPKM")
    expr <- .extract.data(data, columns, calltype)
  } else if(stats  == "tpm"){
    ## TO DO: for Bgee 13, output error message for TPM
    columns <- c("Library.ID", "Gene.ID", "TPM")
    expr <- .extract.data(data, columns, calltype)
  } else if (stats == "counts"){
    columns <- c("Library.ID", "Gene.ID", "Read.count")
    expr <- .extract.data(data, columns, calltype)
  } else {
    cat("Extracting intensities...\n")
    columns <- c("Chip.ID", "Probeset.ID", "Log.of.normalized.signal.intensity", "Gene.ID")
    expr <- .extract.data(data, columns, calltype)
  }
  if(is.data.frame(expr$assayData)){
    # one data matrix

    ## sort objects to have samples in the same order
    expr$assayData <- expr$assayData[, order(names(expr$assayData))]
    expr$pheno <- expr$pheno[order(sampleNames(expr$pheno))]

    ## create ExpressionSet object
    eset <- new('ExpressionSet',
                exprs=as.matrix(expr$assayData),
                phenoData=expr$pheno,
                featureData=expr$features)
  } else if(is.list(expr$assayData)){
    # multiple data matrices
    eset <- mapply(function(x,y,z){
      ## sort objects to have samples in the same order
      x <- x[, order(names(x))]
      y <- y[order(sampleNames(y))]

      ## create ExpressionSet object
      new('ExpressionSet',
          exprs=as.matrix(x),
          phenoData=y,
          featureData=z)
    },
    expr$assayData, expr$pheno, expr$features)
  }
  return(eset)
}
))
