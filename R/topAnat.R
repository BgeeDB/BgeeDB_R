########################
#' @title Produces an object allowing to perform GO-like enrichment of anatomical terms using the topGO package
#'
#' @description This function loads the functions necessary to use topGO on the Uberon ontology instead of the GO ontology. It then produces a topAnatObject, ready to use for gene set enrichment testing using functions from the topGO package.
#'
#' @details To perform the enrichment test for expression in anatomical structures for each term of Uberon ontology (browsable at \url{http://www.ontobee.org/ontology/UBERON}), it is interesting to use the topGO package since it allows to propagate the mapping of gene to terms to parent terms, and it possesses a pannel of enrichment tests and decorrelation methods.
#'
#' @param topAnatData a list produced by the function loadTopAnatData().
#'
#' @param geneList Vector indicating foreground and background genes. Names of the vector indicate teh background genes. Vector values are 0 (gene not in foreground) and 1 (gene in foreground).
#'
#' @param nodeSize Minimum number of genes mapped to a node for it to be tested. Default is 10.
#'
#' @param ... Additional parameters as passed to topGOdata object of topGO package.
#'
#' @return topAnatObject, a topGO-compatibe object ready for gene set enrichment testing. 
#'
#' @author Julien Roux \email{julien.roux@unil.ch}.
#'
#' @examples
#'   myTopAnatObject <- topAnat(myTopAnatData)
#' 
#' @import topGO
#' @export

topAnat <- function(topAnatData, geneList, nodeSize = 10, ... ){
  ## topAnatData is a list including a gene2anatomy list, an organ.relationships list and an organ.names data.frame

  ## TO DO: tests if topAnatData not empty
  ##        test if gene List is fine 
  ## if (length(geneList) > 0 & length(levels(geneList)) == 2 & length(geneList) >= 20) {

  ## TO DO: if geneList includes genes not present in topAnatData$gene2anatomy, restrict to these genes (and add warning)
  ## geneList <- factor(as.integer(names(gene2anatomy) %in% StringIDs))

  ## TO DO: report to the user how many genes in background / foreground
  
  topAnatObject <- makeTopAnatDataObject(
                                         parentMapping = topAnatData$organ.relationships,
                                         allGenes = geneList,
                                         nodeSize = nodeSize,
                                         gene2Nodes = topAnatData$gene2anatomy
                                         )
  ## TO DO: implement

 
  
  return(topAnatObject)
}

## TO DO: why topGO not in DESCRIPTION file?
##        warning when loading topGO: problem?
## TO DO: I added a . in front of topAnat functions: is it a good idea?

#################################################################
## Sets of functions given by Adrian Alexa (author of the topGO package, personnal communication) to make topGO work with another ontology than the Gene Ontology

.annFUN.gene2Nodes <- function(feasibleGenes = NULL, gene2Nodes) {
  ## Restrict the mappings to the feasibleGenes set  
  if(!is.null(feasibleGenes))
    gene2Nodes <- gene2Nodes[intersect(names(gene2Nodes), feasibleGenes)]

  ## Throw-up the genes which have no annotation
  if(any(is.na(gene2Nodes)))
    gene2Nodes <- gene2Nodes[!is.na(gene2Nodes)]

  gene2Nodes <- gene2Nodes[sapply(gene2Nodes, length) > 0]

  ## Get all the Terms and keep a one-to-one mapping with the genes
  allNodes <- unlist(gene2Nodes, use.names = FALSE)
  geneID <- rep(names(gene2Nodes), sapply(gene2Nodes, length))

  return(split(geneID, allNodes))
}
#################################################################

.buildGraph.topology <- function(knownNodes, parentMapping) {
  ## first build the lookUp table for the terms
  nodeLookUp <- new.env(hash = T, parent = emptyenv())

  ## warping functions for a easier acces to the lookUp table
  isNodeInDAG <- function(node) {
    return(exists(node, envir = nodeLookUp, mode = 'logical', inherits = FALSE))
  }
  setNodeInDAG <- function(node) {
    assign(node, TRUE, envir = nodeLookUp)
  }

  ## get the root of the ontology
  GENE.ONTO.ROOT <- setdiff(unique(unlist(parentMapping, use.names = FALSE)),
      names(parentMapping))

  ## we read all the database once and the access the list
  adjLookUP <- as.list(parentMapping)

  ## we use an environment of environments to store edges: (this way is faster)
  ## in the end we will coerce it to a list of list and build a graphNEL obj. 
  edgeEnv <- new.env(hash = T, parent = emptyenv())

  ## add the arc (u --> v) to edgeEnv of type :
  envAddEdge <- function(u, v) {
    assign(v, 0, envir = get(u, envir = edgeEnv))
  }

  ## recursivly build the induced graph starting from one node
  buildInducedGraph <- function(node) {
    ## if we have visited the node, there is nothing to do
    if(isNodeInDAG(node))
      return(1)

    ## we put the node in the graph and we get his parents
    setNodeInDAG(node)    # we visit the node
    assign(node, new.env(hash = T, parent = emptyenv()), envir = edgeEnv) # adj list

    if(node == GENE.ONTO.ROOT) 
      return(2)

    adjNodes <- adjLookUP[[node]]

    ## debuging point! should not happen!
    if(length(adjNodes) == 0)
      cat('\n There are no adj nodes for node: ', node, '\n')

    for(i in 1:length(adjNodes)) {
      x <- as.character(adjNodes[i])
      envAddEdge(node, x)
      buildInducedGraph(x)
    }

    return(0)
  }

  ## we start from the most specific nodes
  lapply(knownNodes, buildInducedGraph)

  ## now we must transform our env into a Graph structure
  ## for now we use lapply, later we can do it with eapply
  .graphNodes <- ls(edgeEnv)
  .edgeList <- eapply(edgeEnv,
                      function(adjEnv) {
                        aux <- as.list(adjEnv)
                        return(list(edges = match(names(aux), .graphNodes),
                                    weights = as.numeric(aux)))
                      })

  ## now we can build the graphNEL object
  return(new('graphNEL',
             nodes = .graphNodes,
             edgeL = .edgeList,
             edgemode = 'directed'))
}
#################################################################

makeTopAnatDataObject <- function(## the child-parrent relationship 
                                parentMapping,
                                ## a named numeric or factor, the names are the genes ID
                                allGenes, 
                                ## function to select the signif. genes
                                geneSelectionFun = NULL,
                                ## minimum node size
                                nodeSize = 0,
                                ## annotation data
                                gene2Nodes,
                                ## additional parameters
                                ...) {

  ## code from new()
  ClassDef <- getClass("topGOdata", where = topenv(parent.frame()))
  ## with R > 2.3.1, PACKAGE = "base" doesn't seem to work
  .Object <- .Call("R_do_new_object", ClassDef, PACKAGE = "base")
  ## .Object <- .Call("R_do_new_object", ClassDef)

  ## TO DO: - resolve namespace problem? Seem to work when called from worksheet
  ##        - how to access this function from base package?
  ##        - create a new prototype for topAnatData class instead of topGOdata?
  ##        - try without Rstudio, with ESS
  ##        - maybe just source the topGO functions?
  
  ## http://r.789695.n4.nabble.com/question-re-error-message-package-error-quot-functionName-quot-not-resolved-from-current-namespace-td4663892.html
  ## From ?.Call, 'PACKAGE' is I believe meant to name the DLL (libRantsImageRead in this case) rather than the R package.

  
  ## some checking
  if(is.null(names(allGenes)))
    stop("allGenes must be a named vector")

  if(!is.factor(allGenes) && !is.numeric(allGenes))
    stop("allGenes should be a factor or a numeric vector")

  .Object@allGenes <- names(allGenes)

  if(is.factor(allGenes)) {
    if(length(levels(allGenes)) != 2)
      stop("allGenes must be a factor with 2 levels")
    .Object@allScores <- factor(as.character(allGenes))
    .Object@geneSelectionFun <- function(x) {
      return(as.logical(as.integer(levels(x)))[x])
    }
  }
  else {
    .Object@allScores <- as.numeric(allGenes)

    ## function to select which genes are significant
    if(is.null(geneSelectionFun))
      warning("No function to select the significant genes provided!")
    .Object@geneSelectionFun <- geneSelectionFun
  }

  ## size of the nodes which will be pruned
  .Object@nodeSize = as.integer(nodeSize)

  ## this function is returning a list of terms from the specified ontology
  ## with each entry being a vector of genes
  cat("\nBuilding 'most specific' Terms .....")
  mostSpecificTerms <- .annFUN.gene2Nodes(feasibleGenes = .Object@allGenes, gene2Nodes = gene2Nodes)
  cat("\t(", length(mostSpecificTerms), "Terms found. )\n")

  ## build the graph starting from the most specific terms ...
  cat("\nBuild DAG topology ..........")
  g <- .buildGraph.topology(names(mostSpecificTerms), parentMapping)
  cat("\t(",  numNodes(g), "terms and", numEdges(g), "relations. )\n")

  ## probably is good to store the levels but for the moment we don't 
  .nodeLevel <- buildLevels(g, leafs2root = TRUE)

  ## annotate the nodes in the graph with genes
  cat("\nAnnotating nodes ...............")
  g <- mapGenes2GOgraph(g, mostSpecificTerms, nodeLevel = .nodeLevel) ## leafs2root

  ## select the feasible genes
  gRoot <- getGraphRoot(g)
  feasibleGenes <- ls(nodeData(g, n = gRoot, attr = "genes")[[gRoot]])
  cat("\t(", length(feasibleGenes), "genes annotated to the nodes. )\n")

  .Object@feasible <- .Object@allGenes %in% feasibleGenes

  cc <- .countsInNode(g, nodes(g))
  .Object@graph <- subGraph(names(cc)[cc >= .Object@nodeSize], g)

  .Object
}
