########################
#' @title Retrieving the Bgee database data
#' @description A Reference Class to give annotation available on Bgee for particular species and the requested data (rna_seq, affymetrix)
#'
#' @details The expression calls come from Bgee (http://bgee.org), that integrates different expression data types (RNA-seq, Affymetrix microarray, ESTs, or in-situ hybridizations) in multiple animal species. Expression patterns are based exclusively on curated "normal", healthy, expression data (e.g., no gene knock-out, no treatment, no disease), to provide a reference of normal gene expression.
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
#' @examples
#' \dontrun{
#' bgee <- Bgee$new(species = "Mus_musculus", datatype = "rna_seq")
#' annotation_bgee_mouse <- bgee$get_annotation()
#' data_bgee_mouse <- bgee$get_data()
#' data_bgee_mouse_gse30617 <- bgee$get_data(experiment.id = "GSE30617")
#' gene.expression.mouse.rpkm <- bgee$format_data(data_bgee_mouse_gse30617, calltype = "present", stats = "rpkm")
#' }
#'
#'
#' @import methods
#' @importFrom tidyr spread
#' @export Bgee
#' @exportClass Bgee



Bgee <- setRefClass("Bgee",

  fields = list(
        species="character",
        datatype = "character",
        experiment.id = "character",
        data = "data.frame",
        calltype = "character",
        stats = "character"),


  methods = list(


        initialize=function(...) {
          callSuper(...)

          if(length(species)==0) {
            cat("WARNING: You didn't specify species. Setting to Homo sapiens.\n")
            species <<- "Homo_sapiens"
          }
          if(length(datatype)==0) {
            cat("WARNING: You didn't specify a data type. Choose 'affymetrix' or 'rna_seq'. Setting to RNAseq for now.\n")
            datatype <<- "rna_seq"

          }},

          get_annotation = function(...){
            gdsurl <- 'ftp://ftp.bgee.org/current/download/processed_expr_values/%s/%s/'
            myurl <- sprintf(gdsurl,datatype,species)
            # first file is the annotation
            fnames <- try(.listDirectories(myurl),silent=TRUE)
            getwd()
            dir.create(file.path(getwd(), species))
            setwd(file.path(getwd(), species))
            distdir <- getwd()
            download.file(file.path(myurl,fnames[1]),
                        destfile=file.path(distdir,fnames[1]),
                        mode='wb')

            unzip( fnames[1])
            temp <- list.files(pattern="*.tsv$")
            cat("Saved files in ", species, " folder:\n")
            print(temp)
            myanno <- lapply(temp, function(x) as.data.frame(fread(x)))
            names(myanno) <- c("experiment_annotation", "sample_annotation")
            return(myanno)
          },

          get_data = function(..., experiment.id = NULL){
            gdsurl <- 'ftp://ftp.bgee.org/current/download/processed_expr_values/%s/%s/'
            myurl <- sprintf(gdsurl,datatype,species)
            fnames <- try(.listDirectories(myurl),silent=TRUE)
            distdir <- getwd()

            if(length(experiment.id) == 0){
              all_expression_values <- fnames[length(fnames)]


              cat("The experiment is not defined. Hence taking all ", datatype," available for ", species, "\n")

              download.file(file.path(myurl, all_expression_values),
                      destfile = file.path(distdir, all_expression_values),
                      mode = "wb")
              unzip(all_expression_values)
              cat("Unzipping files...\n")

              temp2 <- list.files(pattern="*.tsv.zip")
              mydata <- lapply(temp2, unzip)
              data_all <- lapply(mydata, function(x) as.data.frame(fread(x)))

              ## cleaning up the files
              file.remove(dir(path=distdir,  pattern="*.tsv.zip"))
              file.remove(dir(path=distdir,  pattern="*.tsv"))

              cat("Saving all data in .rds file...\n")

              saveRDS(data_all, file = "Bgee_all_experiments_expression_values.rds")
              return(data_all)
              cat("Done.")

            } else if(!grepl("^(GSE).*$|^(E-).*$", experiment.id, perl = TRUE)){

                  stop("The experiment needs to be empty (to download all data) or start with GSE/E- (for a specific experiment) e.g. 'GSE30617' or 'E-MEXP-2011'  ")

            } else {


              cat("Downloading the experiment id ", experiment.id, "\n")
              pk <- match(experiment.id, sapply(strsplit(sub(".*_", "", fnames), ".", fixed = TRUE), "[[", 1))
              gsedataset <- fnames[pk]

              download.file(file.path(myurl,gsedataset),
                      destfile = file.path(distdir, gsedataset),
                      mode = "wb")
              unzip(gsedataset)

              cat("Loading the data...")

              temp3 <- list.files(pattern="*.tsv$")
              kp <- match(experiment.id, sapply(strsplit(sub(".*_", "", temp3), ".", fixed = TRUE), "[[", 1))
              print(temp3[kp])
              gse <- as.data.frame(fread(temp3[kp]))
              return(gse)
              cat("Done.\n")
              }


            },

            format_data = function(data, calltype = "expressed", stats = "rpkm") {

                  ## warning messages
                  if(!(calls %in% c("expressed", "all"))) stop("Choose between only expressed(present) calls or all (present and absent).")
                  if(!(stats %in% c('rpkm', 'counts', 'tpm'))) stop("Choose between RPKM, counts or TPM, e.g. 'rpkm', 'counts', 'tpm' ")



                  l <- split(data, f = data$"Anatomical entity name")
                  cat("Selecting ", calls, " calls.\n")
                  if(calls == "expressed") lt <- lapply(l, function(x) x[which(x$"Detection flag" == "present"),]) else lt <- l

                  cat("Selecting ", stats, " values.\n")
                  cat("Transforming the data.\n")

                  if(stats == "rpkm"){
                    expr <- lapply(lt, function(x) x[, c("Library ID", "Gene ID", "RPKM")])
                    expr.final <- lapply(expr, function(x) x %>% spread("Library ID", "RPKM"))

                  } else{
                    expr <- lapply(lt, function(x) x[, c("Library ID", "Gene ID", "Read count")])
                    expr.final <- lapply(expr, function(x) x %>% spread("Library ID", "Read count"))

                  }
                  cat("Done.\n")
                  return(expr.final)
              }


          ))




