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
#' @field species A character indicating the species to be used, in the form "Genus_species", or a numeric indicatic spcies NCBI taxonomic idea. Only species with quantitative expression data in Bgee will work (RNA-seq and Affymetrix microarray data). See the listBgeeSpecies() function to get the list of species available in the Bgee release used for these two data types.
#'
#' @field datatype A character of data platform. Quantitative expression levels can be downloaded for two data types:
#' \itemize{
#'      \item{"rna_seq"}
#'      \item{"affymetrix"}}
#'
#' @field experiment.id An ArrayExpress or GEO accession, e.g., GSE30617
#' On default is NULL: takes all available experiments for specified species and datatype.
#'
#' @field pathToData Path to the directory where the data files are stored / will be stored. Default is a folder names from used species in working directory.
#'
#' @field release Bgee release number to download data from, in the form "Release.subrelease" or "Release_subrelease", e.g., "13.2" or 13_2". Will work for release >=13.2. By default, the latest relase of Bgee is used.
#'
#' @field data A dataframe of downloaded Bgee data.
#'
#' @field calltype A character.
#'  \itemize{
#'    \item{"present"}
#'    \item{"present high quality"}
#'    \item{"all"}}
#' Retrieve intensities only for present (expressed) genes, present high quality genes, or all genes. The default is present.
#'
#' @field stats A character. The expression values can be retrieved in RPKMs and raw counts:
#'  \itemize{
#'    \item{"rpkm"}
#'    \item{"counts"}
#'    \item{"intensities"}
#'    }
#'The default is RPKMs for RNA-seq and intensities for microarray.
#'
#' @return
#' \itemize{
#'  \item{\code{get_annotation()} returns a list of the annotation of experiments for chosen species.}
#'  \item{\code{get_data()}, if experiment ID is empty, returns a list of experiments. If specified experiment ID, then returns the dataframe of the chosen experiment}
#'  \item{\code{format_data()}, if experiment ID is empty, returns a list of ExpressionSet objects. If specified experiment ID, then returns an ExpressionSet object}}
#'
#' @author Andrea Komljenovic \email{andrea.komljenovic at unil.ch}.
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
  myurl = "character",
  fnames = "character"
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

  # TO DO: code below might also be used in loadTopAnat.R function
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
  } else if (is.character(species)){
    speciesSplitted <- unlist(strsplit(species, split="_"))
    if (sum(allSpecies$GENUS == speciesSplitted[1] & allSpecies$SPECIES_NAME == speciesSplitted[2]) != 1){
      stop(paste0("ERROR: The specified species name is invalid, or not available in Bgee release ", release, "."))
    } else {
      speciesName <<- species
      speciesId <<- as.numeric(allSpecies$ID[allSpecies$GENUS == speciesSplitted[1] & allSpecies$SPECIES_NAME == speciesSplitted[2]])
    }
  }

  ## check data type availability for chosen species
  #if( datatype == "rna_seq" & allSpecies$rna_seq[allSpecies$ID == speciesId] == FALSE){
  #  stop("ERROR: The specified species name is not among the list of species with RNA-seq data in Bgee release ", release,". See listBgeeSpecies() for details on data types availability for each species.")
  #} else if(datatype == "affymetrix" & allSpecies$affymetrix[allSpecies$ID == speciesId] == FALSE){
  #  stop("ERROR: The specified species name is not among the list of species with Affymetrix microarray data in Bgee release ", release,". See listBgeeSpecies() for details on data types availability for each species.")
  #}
  ## TO DO: uncomment code above when Fred has implemented webservice retrieval of data avaibility

  # Creating the FTP URL to get data and annotation
  myurl <<- paste0(allReleases$FTP.URL[allReleases$release == gsub("_", ".", release)], "download/processed_expr_values/", datatype, "/", speciesName, "/")
  ## list files on the FTP at this URL
  listFiles <- function(url) {
    tmpcon <- textConnection(getURL(url), "r")
    allFiles <- read.table(tmpcon)
    close(tmpcon)
    ## Keep only last column
    allFiles <- as.character(allFiles[, ncol(allFiles)])
    return(allFiles)
  }
  fnames <<- try(listFiles(myurl), silent=FALSE)
  if (class(fnames) == "try-error") {
    stop("ERROR: Connection to FTP was not successful.")
  } else if (length(fnames)==0){
    stop("ERROR: Connection to FTP was successful, but target directory seems empty.")
  }

  ## check path of folder to store cached files
  if(length(pathToData)==0) {
    pathToData <<- paste0(getwd(), "/", speciesName, "_Bgee_", release)
  } else if (length(pathToData)==1){
    if ( !file.exists(pathToData) ){
      stop("Problem: please specify a valid path to store data files.")
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

    ## Annotation file names
    if (datatype == "affymetrix"){
        annotation.experiments <- paste0(pathToData, "/", speciesName, "_Affymetrix_experiments.tsv")
        annotation.samples     <- paste0(pathToData, "/", speciesName, "_Affymetrix_chips.tsv")
    } else if (datatype == "rna_seq"){
        annotation.experiments <- paste0(pathToData, "/", speciesName, "_RNA-Seq_experiments.tsv")
        annotation.samples     <- paste0(pathToData, "/", speciesName, "_RNA-Seq_libraries.tsv")
    }

    ## Check if file is already in cache. If so, skip download step
    if (file.exists(annotation.experiments) && file.exists(annotation.samples)){
        cat("WARNING: annotation files for this species were found in the download directory and will be used as is. Please delete the files and rerun the function, if you want the data to be updated.\n")
    } else {
        cat("Downloading annotation files...\n")
        if (datatype == "affymetrix"){
            annotation.file <- paste0(speciesName, "_Affymetrix_experiments_chips.zip")
        } else if (datatype == "rna_seq"){
            annotation.file <- paste0(speciesName, "_RNA-Seq_experiments_libraries.zip")
        }
        if (sum(annotation.file %in% fnames) == 0){
            stop("WARNING. The annotation file was not found on the FTP repository.\n")
        }
        download.file(file.path(myurl, annotation.file),
                      destfile=file.path(pathToData, annotation.file),
                      mode='wb')
        unzip(paste0(pathToData, "/", annotation.file), exdir=pathToData)
        cat("Saved annotation files in", pathToData, "folder.\n")
        ## Clean directory
        file.remove(file.path(pathToData, annotation.file))
    }

    ## Read the 2 annotation files
    myanno <- list(sample_annotation=as.data.frame(fread(annotation.samples)), experiment_annotation=as.data.frame(fread(annotation.experiments)))
    ## remove spaces in headers
    for (i in 1:length(myanno)){
      names(myanno[[i]]) <- make.names(names(myanno[[i]]))
    }
    return(myanno)
},


get_data = function(..., experiment.id = NULL){

    if (length(experiment.id) == 0){
        cat("The experiment is not defined. Hence taking all", datatype, "experiments available for", speciesName, ".\n")

        ## expression data file name
        if (datatype == "affymetrix"){
          all_expression_values <- paste0(speciesName, "_Affymetrix_probesets.zip")
        } else if (datatype == "rna_seq"){
          all_expression_values <- paste0(speciesName, "_RNA-Seq_read_counts_RPKM.zip")
        }

        ## check if RDS file already in cache. If so, skip download step
        if (file.exists(paste0(pathToData, "/", datatype, "_all_experiments_expression_data.rds"))){
            cat("WARNING: expression data file (.rds file) was found in the download directory and will be used as is. Please delete and rerun the function if you want the data to be updated.\n")
            data_all <- readRDS(file = paste0(pathToData, "/", datatype, "_all_experiments_expression_data.rds"))
        } else {
            cat("Downloading expression data...\n")
            if (sum( all_expression_values %in% fnames) == 0){
                stop("WARNING. The expression data file was not found on the FTP repository.")
            }
            download.file(file.path(myurl, all_expression_values),
                          destfile=file.path(pathToData, all_expression_values),
                          mode='wb')
            cat("Saved expression data file in", pathToData, "folder.\n")
            cat("Unzipping file...\n")
            unzip(file.path(pathToData, all_expression_values), exdir=pathToData)
            if(datatype == "affymetrix"){
                temp.files <- list.files(path = pathToData, pattern=".*_probesets_.*.zip$")
            } else {
                temp.files <- list.files(path = pathToData, pattern=".*_RPKM_.*.zip$")
            }
            # print(temp.files)
            mydata <- lapply(file.path(pathToData, temp.files), unzip, exdir=pathToData)
            data_all <- lapply(unlist(mydata, rec = TRUE), function(x) as.data.frame(suppressWarnings(fread(x))))

            ## remove spaces in headers
            for (i in 1:length(data_all)){
              names(data_all[[i]]) <- make.names(names(data_all[[i]]))
            }

            cat("Saving all data in .rds file...\n")
            saveRDS(data_all, file = paste0(pathToData, "/", datatype, "_all_experiments_expression_data.rds"))
        }
    } else if( length(experiment.id) == 1){
      experiment.id <<- experiment.id
        if (!grepl("^GSE\\d+$|^E-\\w+-\\d+.*$", experiment.id, perl = TRUE)){
            stop("The experiment.id field needs to be a GEO or ArrayExpress accession, e.g., 'GSE30617' or 'E-MEXP-2011'")
        } else {
            cat("Downloading expression data for the experiment", experiment.id, "\n")

            ## expression data file name
            if (datatype == "affymetrix"){
                temp.file <- paste0(speciesName, "_Affymetrix_probesets_", experiment.id,".zip")
            } else if (datatype == "rna_seq"){
                temp.file <- paste0(speciesName, "_RNA-Seq_read_counts_RPKM_", experiment.id,".tsv.zip")
            }
            ## check if RDS file already in cache. If so, skip download step
            if (file.exists(paste0(pathToData, "/", datatype, "_", experiment.id, "_expression_data.rds"))){
                cat("WARNING: expression data file (.rds file) was found in the download directory for", experiment.id, "and will be used as is. Please delete and rerun the function if you want the data to be updated.\n")
                data_all <- readRDS(paste0(pathToData, "/", datatype, "_", experiment.id, "_expression_data.rds"))
            } else {
                cat("Downloading expression data...\n")
                if (sum( temp.file %in% fnames) == 0){
                    stop("WARNING. The expression data file for this experiment was not found on the FTP repository. Please use get_annotation() function to be sure it is avilable from Bgee.")
                }
                download.file(file.path(myurl, temp.file),
                              destfile=file.path(pathToData, temp.file),
                              mode='wb')
                cat("Saved expression data file in", pathToData, "folder.\n")
                cat("Unzipping file...\n")
                # Unzipping this file can give one expression data file or multiple ones (if multiple chip types used in experiment)
                mydata <- unzip(file.path(pathToData, temp.file), exdir=pathToData)
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
            }
        }
    } else {
        stop("Please provide only one experiment ID. If you want to get all data for this species and datatype, leave experiment.id empty")
    }
    ## cleaning up downloaded files
    if(datatype == "affymetrix"){
        try(file.remove(file.path(pathToData, list.files(path=pathToData,  pattern=".*_probesets.*.zip"))))
        try(file.remove(file.path(pathToData, list.files(path=pathToData,  pattern=".*_probesets.*.tsv"))))
    } else {
        try(file.remove(file.path(pathToData, list.files(path=pathToData,  pattern=".*_RPKM.*.zip"))))
        try(file.remove(file.path(pathToData, list.files(path=pathToData,  pattern=".*_RPKM.*.tsv"))))
    }
    return(data_all)
    cat("Done.")
},


format_data = function(data, calltype = "all", stats = NULL){
  if (length(stats) == 0){
    stop("Please specify the stats parameters. Should be set to \"intensities\" for Affymetrix data and to \"rpkm\" or \"counts\" for RNA-seq data.")
  } else {
    if (datatype == "affymetrix" & stats != "intensities"){
        stop("For Affymetrix microarray data, stats parameter should be set to \"intensities\"")
    } else if (datatype == "rna_seq" & !(stats %in% c('rpkm', 'counts'))){
        stop("Choose whether data formatting should create a matrix of RPKMs or read counts, with stats option set as \"rpkm\" or \"counts\"")
    }
  }
  if(!(calltype %in% c('present','present high quality','all'))){
    stop("Choose between displaying intensities for present genes, present high quality genes or all genes, e.g., 'present', 'present high quality', 'all' ")
  }

  if(length(data) == 1) data[[1]] else data

  if(stats  == "rpkm"){
    columns <- c("Library.ID", "Gene.ID", "RPKM")
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
    eset <- new('ExpressionSet',
                exprs=as.matrix(expr$assayData),
                phenoData=expr$pheno,
                featureData=expr$features)
  } else if(is.list(expr$assayData)){
    # multiple data matrices
    eset <- mapply(function(x,y,z){
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
