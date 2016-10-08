---
output: html_document
---
# BgeeDB, an R package for retrieval of curated expression datasets and for gene list enrichment tests
##### Andrea Komljenović, Julien Roux, Marc Robinson-Rechavi, Frédéric Bastian

```BgeeDB``` is a collection of functions to import data from the Bgee database (<http://bgee.org/>) directly into R, and to facilitate downstream analyses, such as gene set enrichment test based on expression of genes in anatomical structures. Bgee provides annotated and processed expression data and expression calls from curated wild-type healthy samples, from humans and many other animal species.
 
The package retrieves the annotation of RNA-seq or Affymetrix experiments integrated into the Bgee database, and downloads into R the quantitative data and expression calls produced by the Bgee pipeline. The package also allows to run GO-like enrichment analyses based on anatomical terms, where genes are mapped to anatomical terms by expression patterns, based on the ```topGO``` package. This is the same as the TopAnat web-service available at (<http://bgee.org/?page=top_anat#/>), but with more flexibility in the choice of parameters and developmental stages.

In summary, the BgeeDB package allows to: 
* 1. List annotation of RNA-seq and microarry data available the Bgee database
* 2. Download the processed gene expression data available in the Bgee database
* 3. Download the gene expression calls and use them to perform TopAnat analyses 

## Installation

In R:
``` {r}
source("https://bioconductor.org/biocLite.R")
biocLite("BgeeDB")
```

## How to use BgeeDB package

### Load the package
``` {r, message = FALSE, warning = FALSE}
library(BgeeDB)
```

### Running example: downloading and formatting RNA-seq data

#### List available species in Bgee
The ```listBgeeSpecies()``` function allows to retrieve available species in the Bgee database, and which data types are available for each species. 

``` {r}
listBgeeSpecies()
```

It is possible to list species from a specific release of Bgee with the ```release``` argument (see ```listBgeeRelease()``` function), and order the species according to a specific columns with the ```ordering``` argument. For example:

``` {r}
listBgeeSpecies(release = "13.2", order = 2)
```

#### Choose the species and data type to create a new Bgee object

In the following example we will choose mouse ("Mus\_musculus"), but any other species with RNA-seq or Affymetrix microarray data from the above list could be chosen. Species can also be specified using their NCBI taxonomic IDs. Specify that RNA-seq data are wanted with the ```dataType``` argument set to "rna\_seq". To download Affymetrix microarray data, set this argument to "affymetrix". 

``` {r}
bgee <- Bgee$new(species = "Mus_musculus", dataType = "rna_seq")
```

*Note*: It is possible to work with data from a specific release of Bgee by specifying the ```release``` argument, see ```listBgeeRelease()``` function. The functions of the package will store the downloaded files in a versioned folder in your working directory to allow faster later access to the data. The root folder where this versioned cache folder is created can be changed with the ```pathToData``` argument.

#### Retrieve the annotation of mouse RNA-seq datasets

The ```getAnnotation()``` function will output the list of RNA-seq experiments and libraries available in Bgee for mouse. 

``` {r}
annotation_bgee_mouse <- getAnnotation(bgee)
# list the first experiments and libraries
lapply(annotation_bgee_mouse, head)
```

#### Download the processed mouse RNA-seq data

The ```getData()``` will download RNA-seq data from all available mouse experiments in Bgee as a list. The downloaded files will be stored in the versioned folder created by the ```getAnnotation()``` function above.

``` {r}
# download all RNA-seq experiments from mouse
data_bgee_mouse <- getData(bgee)
# the number of experiments downloaded from Bgee
length(data_bgee_mouse)
# check the downloaded data
lapply(data_bgee_mouse, head)
# isolate the first experiment
data_bgee_experiment1 <- data_bgee_mouse[[1]]
```

*Note*: An additional column in the data frame including expression in the TPM unit is going to be available in Bgee release 14 (planned for the end of 2016). 

Alternatively, you can choose to download data from only one particular mouse experiment from Bgee:

``` {r}
# download data for GSE30617
data_bgee_mouse_gse30617 <- getData(bgee, experimentId = "GSE30617")
```

The result of the ```getData()``` function is, for each experiment, a data frame with the different samples listed in rows, one after the other. It is sometimes easier to work with data organized as a matrix, where different columns represent different samples. The ```formatData()``` function reformats the data into an ExpressionSet object including:
* An expression data matrix, with genes or probesets as rows, and samples as columns (```assayData``` slot). The ```stats``` argument allows to choose if the matrix should be filled with read counts, RPKMs (and soon TPMs) for RNA-seq data. For micoarray data the matrix is filled with log2 expression intensities.
* A data frame listing the samples and their anatomical structure and developmental stage annotation (```phenoData``` slot)
* For microarray data, the mapping from probesets to Ensembl genes (```featureData``` slot)

This function also allows to retain only actively expressed genes or probesets, with the ```callType``` argument set to "present" or "present high quality". Genes or probesets that are absent in a given sample are given ```NA``` values.

```{r}
# use only present calls and fill expression matric with RPKM values
gene.expression.mouse.rpkm <- formatData(bgee, data_bgee_mouse_gse30617, callType = "present", stats = "rpkm")
gene.expression.mouse.rpkm 
```

### Running example: TopAnat gene expression enrichment analysis

For some documentation on the TopAnat analysis, please refer to our publications, or to the web-tool page (<http://bgee.org/?page=top_anat#/>).

#### Download the data allowing to perform TopAnat analysis

The ```loadTopAnatData()``` function loads a mapping from genes to anatomical structures based on calls of expression in anatomical structures. It also loads the structure of the anatomical ontology and the names of anatomical structures. Below, we will choose zebrafish as targeted species:

```{r}
# Creating new Bgee class object
bgee <- Bgee$new(species = "Danio_rerio")
# Loading calls of expression
myTopAnatData <- loadTopAnatData(bgee)
```
TO DO: mention thta new Bgee object needs to be built

By default all data types available for the targeted species are used. This can be changed using the ```dataType``` argument. The data quality can be changed with the ```confidence``` argument. Finally, if you are interested in expression data coming from a particular developmental stage or a group of stages, please specify the ```stage``` argument. 

TO DO: rewrite this

```{r, eval=FALSE}
# Loading expression calls from affymetrix data made on embryonic samples only 
# Not to be run on this vignette
## bgee <- Bgee$new(species = "Danio_rerio", dataType="affymetrix")
## myTopAnatData <- loadTopAnatData(bgee, stage="UBERON:0000068")
# TO DO: add a comment to say not to run it

# Look at the data
str(myTopAnatData)
```

Similarly to the examples above, the downloaded data files are stored in a versioned folder that can be set with the ```pathToData``` argument. If you query again Bgee with the exact same parameters, these cached files will be read instead of querying the web-service. It is possible to work with data from a specific release of Bgee by specifying the ```release``` argument, see ```listBgeeRelease()``` function.

TO DO: update this paragraph

#### Prepare a topGO object allowing to perform TopAnat analysis

First we need to prepare a list of genes in the foreground and in the background. The input format is the same as the gene list required to build a ```topGOdata``` object in the ```topGO``` package: a vector with background genes as names, and 0 or 1 values depending if a gene is in the foreground or not. In this example we will look at genes, annotated with "spermatogenesis" in the Gene Ontology (using the ```biomaRt``` package). The background is set to all genes annotated to at least one Gene Ontology term.

```{r}
source("https://bioconductor.org/biocLite.R")
biocLite("biomaRt")
library(biomaRt)
ensembl <- useMart("ensembl")
ensembl <- useDataset("drerio_gene_ensembl", mart=ensembl)

# Foreground genes are those with GO annotation "spermatogenesis"
myGenes <- getBM(attributes= "ensembl_gene_id", filters=c("go_id"), values=list(c("GO:0007283")), mart=ensembl)

# Background are all genes with GO annotation
universe <- getBM(attributes= "ensembl_gene_id", filters=c("with_go_go"), values=list(c(TRUE)), mart=ensembl)

# Prepare the gene list vector 
geneList <- factor(as.integer(universe[,1] %in% myGenes[,1]))
names(geneList) <- universe[,1]
head(geneList)
summary(geneList == 1)

# Prepare the topGO object
myTopAnatObject <-  topAnat(myTopAnatData, geneList)
```

*Warning*: This can be long, especially if the gene list is large, since the anatomical ontology is large and expression calls will be propagated through the whole ontology (e.g., expression in the forebrain will also be counted as expression in parent structures such as the brain, nervous system, etc). Consider running a script in batch mode if you have multiple analyses to do.

#### Launch the enrichment test

For this step, see the vignette of the ```topGO``` package for more details, as you have to directly use the tests implemented in the ```topGO``` package, as shown in this example:

```{r}
results <- runTest(myTopAnatObject, algorithm = 'classic', statistic = 'fisher')
```

You can also choose one of the topGO decorrelation methods, for example the "weight" method, allowing to get less redundant results
```{r, eval=FALSE}
## Not to be run on this vignette
## results <- runTest(myTopAnatObject, algorithm = 'weight', statistic = 'fisher')

```

*Warning*: This can be long because of the size of the ontology. Consider running a script in batch mode if you have multiple analyses to do.

#### Format the table of results after an enrichment test for anatomical terms

The ```makeTable``` function allows to filter and format the test results, and calculate FDR values. 

```{r}
# Display results sigificant at a 5% FDR threshold
tableOver <- makeTable(myTopAnatData, myTopAnatObject, results, cutoff = 0.05)
# Display all results, sorted by p-value
tableOver <- makeTable(myTopAnatData, myTopAnatObject, results, cutoff = 1)
```

By default results are sorted by p-value, but this can be changed with the ```ordering``` parameter by specifying which column should be used to order the results (preceded by a "-" sign to indicate that ordering should be made in decreasing order). For example, it is often convenient to sort  significant structures by decreasing enrichment fold.

*Warning*: it is debated whether FDR correction is appropriate on enrichment test results, since tests on different terms of the ontologies are not independent. A nice discussion can be found in the vignette of the ```topGO``` package.
