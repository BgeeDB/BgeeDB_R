########################
#' @title Formats the results obtained from gene set enrichment test on anatomical structures
#'
#' @description This function loads the results from the test function (topGO package) and creates an output table with organ names, fold enrichment and FDR
#'
#' @param topAnatData A list produced by the function loadTopAnatData().
#'
#' @param topAnatObject An object produced by the function topAnat().
#'
#' @param results A result object, produced for example by the runtest() function of topGO
#'
#' @param cutoff An FDR cutoff between 0 and 1. Only terms with FDR lower than this cutoff will be in the output. Default is 1.
#'
#' @return A data frame with significantly enriched anatomical structures, sorted by FDR
#'
#' @author Julien Roux \email{julien.roux@unil.ch}.
#'
#' @examples
#'   \dontrun{
#'     ## Launch topGO test on ata loaded from Bgee
#'     resFis <- runTest(myTopAnatObject, algorithm = 'elim', statistic = 'fisher')
#'     tableOver <- makeTable(myTopAnatData, myTopAnatObject, resFis, 0.1)
#'   }
#' @export

makeTable <- function(topAnatData, topAnatObject, results, cutoff){
  ## TO DO: check input objects
  ## score(results)
  ## is.numeric(cutoff)
  ## termStat(topAnatObject, ...)
  ## topAnatData$organ.names not empty (see topAnat.R)

  ## retrieve p-values for the enrichment
  scores <- score(results)
  fdr <- p.adjust(p=scores, method = "fdr")
  topTerms <- sort(scores[fdr <= cutoff])
  topTerms <- as.data.frame(topTerms)

  if( nrow(topTerms) != 0 ){
    cat("\nBuilding the results table for the", nrow(topTerms), "significant terms at FDR threshold of", cutoff, "... ")
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
