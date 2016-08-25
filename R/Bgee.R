.subsetting <- function(st, v1, v2 , calltype, column){
    if(class(st) == "list"){
        calls <- lapply(st, function(x) .calling(x, calltype, column))
        expr <- lapply(calls, function(x) {
            xt <- x[, v1]
            xtt <- xt %>% spread_(v1[1], v1[3])
            rownames(xtt) <- xtt[,v1[2]]
            xtt[,-1, drop = FALSE] })

        features <- mapply(.subsetting.feature, calls, expr, v1[2], v2 )
        phenos <- mapply(.subsetting.pheno, calls, v1[1] )
    } else {
        calls <- .calling(st, calltype, column)
        dt <- calls[, v1]
        dtt <- dt %>% spread_(v1[1], v1[3])
        rownames(dtt) <- dtt[,v1[2]]
        expr <- dtt[,-1,  drop = FALSE]
        cat("Subseting features...\n")
        features <- .subsetting.feature( calls, expr, v1[2], v2)
        cat("Done...\n")
        cat("Subseting pheno...\n")
        phenos <- .subsetting.pheno( calls, v1[1] )
        cat("Done...\n")
    }
    return(list(assayData = expr, pheno = phenos, features = features, calls = calls))
}
#this works
.subsetting.feature <- function(x, y, v1, v2){
    cat("building the feature data...\n")
    cat("building the feature data...\n")
    if(v1 == v2){
        fdata <- x[match(rownames(y), x[,v1]), v1, drop = FALSE]
        rownames(fdata) <- fdata[,v1]
        fdata <- as(fdata, "AnnotatedDataFrame")

    } else {
        fdata <- x[match(rownames(y), x[,v1]), c(v1,v2), drop = FALSE]
        rownames(fdata) <- fdata[,v1]
        fdata <- as(fdata, "AnnotatedDataFrame")
    }
    return(fdata)
}



.subsetting.pheno <- function(x, v4){
    cat("building the pheno data...\n")
    phdata <- x[, c(v4, "Anatomical entity ID", "Anatomical entity name", "Stage ID", "Stage name")]
    phdata <- phdata[!duplicated(phdata[,v4]),]
    rownames(phdata) <- phdata[,v4]
    phdata <- as(phdata,"AnnotatedDataFrame")
    return(phdata)
}


# this works
.calling <- function(x, calltype, column){
    ## check datatype
    if(calltype == "expressed"){
        cat("calling expressed...\n")
        x[(x$"Detection flag" == "absent"), column] <- NA
    } else if (calltype == "expressed high quality"){
        cat("calling expressed and high quality...\n")
        x[which(x$"Detection flag" == "absent" | x$"Detection quality" == "poor quality"), column] <- NA
    } else {
        cat("calling both present and absent calls with all type of quality...")
        x
    }
    return(x)
}


########################
#' @title Retrieving the Bgee database data
#' @description A Reference Class to give annotation available on Bgee for particular species and the requested data (rna_seq, affymetrix)
#'
#' @details The expression calls come from Bgee (http://r.bgee.org), that integrates different expression data types (RNA-seq, Affymetrix microarray, ESTs, or in-situ hybridizations) in multiple animal species. Expression patterns are based exclusively on curated "normal", healthy, expression data (e.g., no gene knock-out, no treatment, no disease), to provide a reference of normal gene expression.
#' This Class retrieves annotation of all experiments in Bgee database (get_annotation), downloading the data (get_data), and formating the data into expression matrix (format_data). See examples and vignette.
#'
#'
#'
#'
#' @field species A character of species name as listed from Bgee. The species are:
#' \itemize{
#'    \item{"Anolis_carolinensis"}
#'    \item{"Bos_taurus"}
#'    \item{"Caenorhabditis_elegans"}
#'    \item{"Danio_rerio"}
#'    \item{"Drosophila_melanogaster"}
#'    \item{"Gallus_gallus"}
#'    \item{"Gorilla_gorilla"}
#'    \item{"Homo_sapiens"}
#'    \item{"Macaca_mulatta"}
#'    \item{"Monodelphis_domestica"}
#'    \item{"Mus_musculus"}
#'    \item{"Ornithorhynchus_anatinus"}
#'    \item{"Pan_paniscus"}
#'    \item{"Pan_troglodytes"}
#'    \item{"Rattus_norvegicus"}
#'    \item{"Sus_scrofa"}
#'    \item{"Xenopus_tropicalis"}}
#'
#' Homo sapiens is default species.
#'
#'
#'
#' @field datatype A character of data platform. Two types of datasets can be downloaded:
#' \itemize{
#'      \item{"rna_seq"}
#'      \item{"affymetrix"}}
#' By default, RNA-seq data is retrieved from database.
#'
#'
#' @field experiment.id  A character.
#' On default is NULL: takes all available data for that species.
#' If GSE[0-9]+: takes specified experiment, eg. GSE30617.
#'
#' @field data A dataframe of downloaded Bgee data.
#'
#' @field calltype A character. There exist two types of expression calls in Bgee - present and absent.
#'  \itemize{
#'    \item{"expressed"}
#'    \item{"all"}}
#' User can retrieve only expressed (present) calls, or mixed (present and absent) calls. The default is expressed (present) calltype.
#'
#'
#' @field stats A character. The expression values can be retrieved in RPKMs and raw counts:
#'  \itemize{
#'    \item{"rpkm"}
#'    \item{"counts"}}
#'The default is RPKMs.
#'
#'
#'
#' @return
#' \itemize{
#'  \item{A \code{get_annotation()} list, lists the annotation of experiments for chosen species.}
#'  \item{A \code{get_data()}, if empty returns a list of experiments, if chosen experiment ID, then returns the dataframe of the chosen experiment; for chosen species}
#'  \item{A \code{format_data()}, transforms the data into matrix of expression values, e.g. RPKMs or raw counts}}
#'
#'
#' @author Andrea Komljenovic \email{andrea.komljenovic at unil.ch}.
#'
#' @examples{
#'  bgee <- Bgee$new(species = "Mus_musculus", datatype = "rna_seq")
#'  annotation_bgee_mouse <- bgee$get_annotation()
#'  data_bgee_mouse <- bgee$get_data()
#'  data_bgee_mouse_gse30617 <- bgee$get_data(experiment.id = "GSE30617")
#'  gene.expression.mouse.rpkm <- bgee$format_data(data_bgee_mouse_gse30617,
#'  calltype = "expressed", stats = "rpkm")
#'  }
#'
#'
#' @import methods
#' @importFrom dplyr %>%
#' @importFrom tidyr spread
#' @export Bgee
#' @exportClass Bgee

Bgee <- setRefClass("Bgee",

fields = list(
species="character",
datatype = "character",
experiment.id = "character",
data = "list",
calltype = "character",
stats = "character",
myurl = "character",
destdir = "character",
fnames = "character"),


methods = list(

initialize=function(...) {
    callSuper(...)

    ## Species is a compulsory parameter
    allSpecies <- c("Homo_sapiens", "Mus_musculus", "Danio_rerio", "Drosophila_melanogaster",
    "Caenorhabditis_elegans", "Pan_troglodytes", "Pan_paniscus",
    "Gorilla_gorilla", "Macaca_mulatta", "Rattus_norvegicus", "Bos_taurus",
    "Sus_scrofa", "Monodelphis_domestica", "Ornithorhynchus_anatinus",
    "Gallus_gallus", "Anolis_carolinensis", "Xenopus_tropicalis",
    "Pongo_pygmaeus", "Tetraodon_nigroviridis")
    # bothSpeciesAssays <- c("Homo_sapiens", "Caenorhabditis_elegans", "Mus_musculus")
    onlyAffymetrixSpecies <- c("Danio_rerio", "Drosophila_melanogaster")
    onlyRNAseqSpecies <- c("Pan_paniscus", "Pan_troglodytes", "Gorilla_gorilla", "Macaca_mulatta", "Rattus_norvegicus", "Bos_taurus",
    "Sus_scrofa", "Monodelphis_domestica", "Anolis_carolinensis", "Xenopus_tropicalis", "Tetraodon_nigroviridis",
    "Pongo_pygmaeus", "Gallus_gallus", "Ornithorhynchus_anatinus")

    if(length(species)==0) {
        stop("ERROR: You didn't specify species.")
    } else if ( length(species) > 1 ){
        stop("ERROR: only one species is allowed.")
    } else if ( sum(species %in% allSpecies) == 0 ){
        stop("ERROR: the specified speciesId is not among the list of species in Bgee. Maybe you did not specificy species name,
        but common name, or did not put an underscore between genus an species?
        Examples: 'Homo_sapiens', 'Mus_musculus', 'Drosophila_melanogaster', 'Caenorhabditis_elegans'.
        See listBgeeSpecies() for all species available.\n")
    }

    ## check datatype
    if(length(datatype)==0) {
        stop("ERROR: You didn't specify a data type. Choose 'affymetrix' or 'rna_seq'.")
    } else if ((length(datatype) > 1) || (length(datatype) == 1 && datatype %in% c("rna_seq", "affymetrix") == "FALSE")){
        stop("ERROR: Choose correct datatype argument: 'affymetrix' or 'rna_seq'.")
    }

    # check species type
    if( datatype == "rna_seq" && species %in% onlyAffymetrixSpecies){
        stop("ERROR: For this species there is no RNAseq data. Please change the datatype to 'affymetrix'.")
    } else if(datatype == "affymetrix" && species %in% onlyRNAseqSpecies){
        stop("ERROR: For this species there is no Affymetrix data. Please change the datatype to 'rna_seq'.")
    }

    # Creating the folder - common to get_data and get_annotation
    gdsurl <- 'ftp://ftp.bgee.org/current/download/processed_expr_values/%s/%s/'
    ## Built FTP URL for this datatype and species
    myurl <<- sprintf(gdsurl, datatype, species)
    ## list files in this folder
    fnames <<- try(.listDirectories(myurl), silent=TRUE)

    ## create a folder with species name to store downloaded files
    destdir <<- file.path(getwd(), species)
    if (!file.exists(destdir)){
        dir.create(destdir)
    }


},


get_annotation = function(...){

    ## Annotation file names
    if (datatype == "affymetrix"){
        annotation.experiments <- paste0(destdir, "/", species, "_Affymetrix_experiments.tsv")
        annotation.samples     <- paste0(destdir, "/", species, "_Affymetrix_chips.tsv")
    } else if (datatype == "rna_seq"){
        annotation.experiments <- paste0(destdir, "/", species, "_RNA-Seq_experiments.tsv")
        annotation.samples     <- paste0(destdir, "/", species, "_RNA-Seq_libraries.tsv")
    }

    ## Check if file is already in cache. If so, skip download step
    if (file.exists(annotation.experiments) && file.exists(annotation.samples)){
        cat("WARNING: annotation files for this species were found in the download directory and will be used as is.\n
        Please delete the files and rerun the function, if you want the data to be updated.\n")
    } else {
        cat("Downloading annotation files...\n")
        if (datatype == "affymetrix"){
            annotation.file <- paste0(species, "_Affymetrix_experiments_chips.zip")
        } else if (datatype == "rna_seq"){
            annotation.file <- paste0(species, "_RNA-Seq_experiments_libraries.zip")
        }
        if (sum(annotation.file %in% fnames) == 0){
            stop("WARNING. The annotation file was not found on the FTP repository.\n")
        }
        download.file(file.path(myurl, annotation.file),
        destfile=file.path(destdir, annotation.file),
        mode='wb')
        unzip(paste0(destdir, "/", annotation.file), exdir=destdir)
        cat("Saved annotation files in", species, "folder.\n")
        ## Clean directory
        file.remove(file.path(destdir, annotation.file))
    }

    ## Read the 2 annotation files
    myanno <- list(sample_annotation=as.data.frame(fread(annotation.samples)), experiment_annotation=as.data.frame(fread(annotation.experiments)))
    return(myanno)
},


get_data = function(..., experiment.id = NULL){

    if (length(experiment.id) == 0){
        cat("The experiment is not defined. Hence taking all", datatype, "available for", species, ".\n")

        ## expression data file name
        if (datatype == "affymetrix"){
            all_expression_values <- paste0(species, "_Affymetrix_probesets.zip")
        } else if (datatype == "rna_seq"){
            all_expression_values <- paste0(species, "_RNA-Seq_read_counts_RPKM.zip")
        }

        ## check if RDS file already in cache. If so, skip download step
        if (file.exists(paste0(destdir, "/", datatype, "_all_experiments_expression_data.rds"))){
            cat("WARNING: expression data file (.rds file) was found in the download directory and will be used as is.
            Please delete and rerun the function if you want the data to be updated.\n")
            data_all <- readRDS(file = paste0(destdir, "/", datatype, "_all_experiments_expression_data.rds"))
        } else {
            cat("Downloading expression data...\n")
            if (sum( all_expression_values %in% fnames) == 0){
                stop("WARNING. The expression data file was not found on the FTP repository.")
            }
            download.file(file.path(myurl, all_expression_values),
            destfile=file.path(destdir, all_expression_values),
            mode='wb')
            cat("Saved expression data file in", species, "folder.\n")
            cat("Unzipping file...\n")
            unzip(file.path(destdir, all_expression_values), exdir=destdir)
            if(datatype == "affymetrix"){
                temp.files <- list.files(path = destdir, pattern=".*_probesets_.*.zip$")
            } else {
                temp.files <- list.files(path = destdir, pattern=".*_RPKM_.*.zip$")
            }
            # print(temp.files)
            mydata <- lapply(file.path(destdir, temp.files), unzip, exdir=destdir)
            print(head(mydata))
            data_all <- lapply(unlist(mydata, rec = T), function(x) as.data.frame(suppressWarnings(fread(x))))

            cat("Saving all data in .rds file...\n")
            saveRDS(data_all, file = paste0(destdir, "/", datatype, "_all_experiments_expression_data.rds"))
        }
    } else if( length(experiment.id) == 1){
        if (!grepl("^GSE\\d+$|^E-\\w+-\\d+.*$", experiment.id, perl = TRUE)){
            stop("The experiment.id field need to be in GEO or ArrayExpress format, e.g., 'GSE30617' or 'E-MEXP-2011'")
        } else {
            cat("Downloading expression data for the experiment", experiment.id, "\n")

            ## expression data file name
            if (datatype == "affymetrix"){
                expression_values <- paste0(species, "_Affymetrix_probesets_", experiment.id,".zip")
            } else if (datatype == "rna_seq"){
                expression_values <- paste0(species, "_RNA-Seq_read_counts_RPKM_", experiment.id,".tsv.zip")
            }
            ## TO DO: check if RDS file already in cache. If so, skip download step
            if (file.exists(paste0(destdir, "/", datatype, "_", experiment.id, "_expression_data.rds"))){
                cat("WARNING: expression data file (.rds file) was found in the download directory for", experiment.id, ".
                These will be used as is. Please delete and rerun the function if you want the data to be updated.\n")
                data_all <- readRDS(paste0(destdir, "/", datatype, "_", experiment.id, "_expression_data.rds"))
            } else {
                cat("Downloading expression data...\n")
                if (sum( expression_values %in% fnames) == 0){
                    stop("WARNING. The expression data file for this experiment was not found on the FTP repository.")
                }
                download.file(file.path(myurl, expression_values),
                destfile=file.path(destdir, expression_values),
                mode='wb')
                cat("Saved expression data file in", species, "folder.\n")
                cat("Unzipping file...\n")
                mydata <- lapply(file.path(destdir, expression_values), unzip, exdir=destdir)
                data_all <- lapply(mydata, function(x) as.data.frame(fread(x)))

                data_all <- as.data.frame(data_all[[1]])
                cat("Saving all data in .rds file...\n")
                saveRDS(data_all, file = paste0(destdir, "/", datatype, "_", experiment.id, "_expression_data.rds"))
            }
        }
    } else {
        stop("Please provide only one experiment ID. If you want to get all data for this species and datatype, leave experiment.id empty")
    }
    ## cleaning up the files
    if(datatype == "affymetrix"){
        try(file.remove(file.path(destdir, list.files(path=destdir,  pattern=".*_probesets.*.zip"))))
        try(file.remove(file.path(destdir, list.files(path=destdir,  pattern=".*_probesets.*.tsv"))))
    } else {
        try(file.remove(file.path(destdir, list.files(path=destdir,  pattern=".*_RPKM.*.zip"))))
        try(file.remove(file.path(destdir, list.files(path=destdir,  pattern=".*_RPKM.*.tsv"))))
    }
    return(data_all)
    cat("Done.")
},


format_data = function(data, calltype, stats){

    if(!(stats %in% c('rpkm', 'counts', 'intensities'))) stop("Choose between RPKM, counts or affymetrix intensities,
    e.g. 'rpkm', 'counts', 'intensities' ")

    if(!(calltype %in% c('expressed','expressed high quality','none'))) stop("Choose between expressed, expressed high quality or none,
    e.g. 'expressed', 'expressed high quality', 'none' ")

    if(length(data) == 1) data[[1]] else data
    v2 <- "Gene ID"
    if(stats  == "rpkm"){
        v1 <- c("Library ID", "Gene ID", "RPKM")
        expr <- .subsetting(data, v1, v2, calltype, "RPKM")
    } else if (stats == "counts"){
        v1 <- c("Library ID", "Gene ID", "Read count")
        expr <- .subsetting(data, v1, v2, calltype, "Read count")
    } else {
        cat("Subseting intensities...\n")
        v1 <- c("Chip ID", "Probeset ID", "Log of normalized signal intensity")
        expr <- .subsetting(data, v1, v2, calltype, "Log of normalized signal intensity")
    }


    if(length(expr$assayData) == 1){
        eset2 <- new('ExpressionSet', exprs= as.matrix(expr$assayData),
        phenoData= expr$pheno,
        featureData = expr$features)
    } else {
        assaydata <- expr$assayData
        phenodata <- expr$pheno
        featdata <- expr$features

        eset2 <- mapply(function(x,y,z)
        new('ExpressionSet', exprs= as.matrix(x), phenoData= y, featureData = z),
        assaydata, phenodata, featdata)
    }
    cat("Done.\n")
    return(eset2)
}



))



