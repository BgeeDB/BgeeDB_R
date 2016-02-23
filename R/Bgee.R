########################
#' @title Bgee RC
#' @description A Reference Class to give annotation available on Bgee for particular species and the requested data (rna_seq, affymetrix)
#'
#'
#'
#'
#'
#' @field species A character of species name as listed from Bgee.
#' Options are:
#'          "Anolis_carolinensis",
#'          "Bos_taurus",
#'          "Caenorhabditis_elegans",
#'          "Danio_rerio",
#'          "Drosophila_melanogaster",
#'          "Gallus_gallus",
#'          "Gorilla_gorilla",
#'          "Homo_sapiens",
#'          "Macaca_mulatta",
#'          "Monodelphis_domestica",
#'          "Mus_musculus",
#'          "Ornithorhynchus_anatinus",
#'          "Pan_paniscus",
#'          "Pan_troglodytes",
#'          "Rattus_norvegicus",
#'          "Sus_scrofa",
#'          "Xenopus_tropicalis"
#'
#' @field datatype A character of data platform.
#' Options are:
#'          "rna_seq",
#'          "affymetrix"
#'
#' @field experiment.id  A character.
#' On default is NULL: takes all available data for that species
#' GSE[0-9]+: takes specified experiment, eg. GSE30617
#'
#' @import methods
#' @export Bgee
#' @exportClass Bgee



Bgee <- setRefClass("Bgee",

  fields = list(
        species="character",
        datatype = "character",
        experiment.id = "character"),


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
            fnames <- try(listDirectories(myurl),silent=TRUE)
            getwd()
            dir.create(paste(getwd(), species, sep = "/"))
            setwd(paste(getwd(), species, sep = "/"))
            distdir <- getwd()
            download.file(file.path(myurl,fnames[1]),
                        destfile=file.path(distdir,fnames[1]),
                        mode='wb')

            unzip( fnames[1])
            temp <- list.files(pattern="*.tsv$")
            print(temp)
            myanno<- lapply(temp, fread)
            return(myanno)
          },

          get_data = function(..., experiment.id = NULL){
            gdsurl <- 'ftp://ftp.bgee.org/current/download/processed_expr_values/%s/%s/'
            myurl <- sprintf(gdsurl,datatype,species)
            fnames <- try(listDirectories(myurl),silent=TRUE)

            if(length(experiment.id) == 0){
              all_expression_values <- fnames[length(fnames)]

            # if(is.null(experiment.id)){
              cat("The experiment is not defined. Hence taking all ", datatype," available for ", species, "\n")
              distdir <- getwd()
              download.file(file.path(myurl, all_expression_values),
                      destfile = file.path(distdir, all_expression_values),
                      mode = "wb")
              unzip(all_expression_values)
              cat("Unzipping files...\n")

              temp2 <- list.files(pattern="*.tsv.zip")
              mydata <- lapply(temp2, unzip)
              data_all <- lapply(mydata, fread)

              ## cleaning up the files
              file.remove(dir(path=distdir,  pattern="*.tsv.zip"))
              file.remove(dir(path=distdir,  pattern="*.tsv"))

              cat("Saving all data in .rds file...\n")

              saveRDS(data_all, file = "Bgee_all_experiments_expression_values.rds")
              return(data_all)
              cat("Done.")

            } else if(!grepl("^(GSE).*$", experiment.id, perl = TRUE)){

                  stop("The experiment needs to be empty (to download all data) or start with GSE (for a specific experiment) e.g. 'GSE30617'  ")

            } else {


              cat("Downloading the experiment id ", experiment.id, "\n")
              pk <- match(experiment.id, sapply(strsplit(sub(".*_", "", fnames), ".", fixed = TRUE), "[[", 1))
              gsedataset <- fnames[pk]

              # setwd(paste(getwd(), species, sep = "/"))
              distdir <- getwd()

              download.file(file.path(myurl,gsedataset),
                      destfile = file.path(distdir, gsedataset),
                      mode = "wb")
              unzip(gsedataset)

              cat("Loading the data...")

              require(data.table)
              temp3 <- list.files(pattern="*.tsv$")
              kp <- match(experiment.id, sapply(strsplit(sub(".*_", "", temp3), ".", fixed = TRUE), "[[", 1))
              print(temp3[kp])
              gse <- fread(temp3[kp])
              return(gse)
              cat("Done.\n")
              }
              #else {
              #  stop("The experiment needs to be empty (to download all data) or start with GSE (for a specific experiment) e.g. 'GSE30617'  ")
              #}


            }


          ))
