
.listDirectories <- function(url, version = FALSE) {
	tmpcon <- textConnection(getURL(url), "r")
	tx <- read.table(tmpcon)
	close(tmpcon)
	ifelse(version, tx1 <- paste0(levels(tx$V6), tx$V8)[1], tx1 <- as.character(tx[,ncol(tx)]))
	return(tx1)
}


#' @title List species with Affymetrix or RNA-seq data in the Bgee database
#' @description Returns a list of available genomes for different platforms in Bgee.
#' @param ... an empty parameter
#'
#'
#' @return A list of species with available Affymetrix or RNA-seq data in the Bgee
#'
#'@examples{
#'  listBgeeSpecies()
#' }
#'
#' @author Andrea Komljenovic \email{andrea.komljenovic at unil.ch}.
#' @import data.table RCurl
#' @export

listBgeeSpecies <- function(...){
	url_rnaseq <-  "ftp://ftp.bgee.org/current/download/processed_expr_values/rna_seq/"
	url_affymetrix <-  "ftp://ftp.bgee.org/current/download/processed_expr_values/affymetrix/"

  species.taxid <- c(9606, 10090, 7955, 7227, 6239, 9598, 9597, 9593, 9544, 10116, 9913, 9823, 13616, 9258, 9031, 28377, 8364, 9600, 99883)
  genus.name <- c("Homo", "Mus", "Danio", "Drosophila", "Caenorhabditis", "Pan", "Pan", "Gorilla", "Macaca", "Rattus", "Bos", "Sus", "Monodelphis", "Ornithorhynchus", "Gallus", "Anolis", "Xenopus", "Pongo", "Tetraodon")
  species.name <- c("sapiens", "musculus", "rerio", "melanogaster", "elegans", "troglodytes", "paniscus", "gorilla", "mulatta", "norvegicus", "taurus", "scrofa", "domestica", "anatinus", "gallus", "carolinensis", "tropicalis", "pygmaeus", "nigroviridis")
  common.name <- c("human", "mouse", "zebrafish", "fruitfly", "c.elegans", "chimpanzee", "bonobo", "gorilla", "macacque", "rat", "cow", "pig", "opossum", "platypus", "chicken", "anolis", "xenopus", "orangutan", "tetraodon")

  df.species.db <- data.frame(ID = species.taxid, GENUS = genus.name, SPECIES = species.name, COMMON = common.name)

  cat("Last update:" , .listDirectories(url_rnaseq, version = TRUE), "\n")
	return(list("Species in the Database"= df.species.db, RNAseq = .listDirectories(url_rnaseq), Affymetrix = .listDirectories(url_affymetrix)))

}



