########################
#' @title Formats results of the enrichment test on anatomical structures.
#'
#' @description This function loads the results from the topGO test and creates an output table with organ names,
#' fold enrichment and FDR. Data are sorted by p-value and only terms below the specified FDR cutoff are included.
#'
#' @param topAnatData A list produced by the function loadTopAnatData().
#'
#' @param topAnatObject An object produced by the function topAnat().
#'
#' @param results A result object, produced by the runtest() function of topGO.
#'
#' @param cutoff An FDR cutoff between 0 and 1. Only terms with FDR lower than this cutoff are included.
#' Default is 1, meaning that all terms are included.
#'
#' @return A data frame with significantly enriched anatomical structures, sorted by p-value.
#'
#' @author Julien Roux \email{julien.roux@unil.ch}.
#'
#' @examples{
#'  ## Launch topGO test on data loaded from Bgee
#'   myTopAnatData <- loadTopAnatData(species = "10090", datatype = "rna_seq")
#'   geneList <- c(0,1,0,1,0,1,0,1)
#'   names(geneList) <- c("gene1","gene3","gene3","gene4","gene5","gene6","gene7","gene8")
#'   myTopAnatObject <- topAnat(myTopAnatData, geneList)
#'   resFis <- runTest(myTopAnatObject, algorithm = 'elim', statistic = 'fisher')
#'  ## Format results
#'   tableOver <- makeTable(myTopAnatData, myTopAnatObject, resFis, 0.1)
#' }
#' @export

makeTable <- function(topAnatData, topAnatObject, results, cutoff=1){
  ## Perform some checks on the input data
  if(is.na(score(results)) || length(score(results)) == 0){
    stop("Problem: the results object is empty.")
  }
  if (!is.numeric(cutoff)){
    cutoff <- as.numeric(cutoff)
  }
  if(is.null(myTopAnatObject)){
    stop("Problem: the topAnatObject is empty.")
  }
  if( length(topAnatData$organ.names[,1]) == 0 ) {
    stop("Problem: the organ.names data frame of your topAnatData object is empty.")
  }

  ## retrieve p-values for the enrichment
  scores <- score(results)
  fdr <- p.adjust(p=scores, method = "fdr")
  topTerms <- sort(scores[fdr <= cutoff])
  topTerms <- as.data.frame(topTerms)

  if( nrow(topTerms) != 0 ){
    cat(paste0("\nBuilding the results table for the ", nrow(topTerms), " significant terms at FDR threshold of ", cutoff, "... "))
    odds <- termStat(topAnatObject, row.names(topTerms))
    foldEnrichment <- odds[2]/odds[3]

    # Rounding odds, P-values and FDR to 2 decimal places
    foldEnrichment <- format(foldEnrichment, digits=3)
    topTerms <- format(topTerms, digits=3)
    fdr[row.names(topTerms)] <- format(fdr[row.names(topTerms)], digits=3)

    topTerms <- cbind(odds, foldEnrichment, topTerms, fdr[row.names(topTerms)])
    names(topTerms) <- c("annotated", "significant", "expected", "foldEnrichment" , "pValue", "FDR")
    topTable <- merge(topAnatData$organ.names, topTerms, by.x=1, by.y=0)
    topTable <- topTable[order(as.numeric(topTable$pValue)), ]

    cat("Done\n\n")
    return(topTable)
  } else {
    cat("\nWarning: there was no significant term at a FDR threshold of", cutoff, "\n\n")
    return(NA)
  }
}
