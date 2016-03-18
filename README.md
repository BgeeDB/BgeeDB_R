# BgeeDB: an R package for datasets retrieval from Bgee database

```BgeeDB``` is a collection of functions to import data from the Bgee database (<http://bgee.org/>) directly into R, and to facilitate downstream analyses, such as gene set enrichment test based on expression of genes in anatomical structures.
 
The package retrieves the annotation of RNA-seq or Affymetrix experiments integrated into the Bgee database, and downloads into R the data reprocessed by the Bgee pipeline. Currently, Bgee database includes gene expression data from 17 species. The package also allows to run GO-like enrichment analyses based on anatomical terms, where genes are mapped to anatomical terms by expression patterns. This gives similar results as the TopAnat webservive available at (<http://bgee.org/?page=top_anat#/>), and is based on the ```topGO``` package.

This package allows: 
* 1. Listing annotation files gene of expression data available in the current version of Bgee database
* 2. Downloading the processed gene expression data available in the current version of Bgee database
* 3. Downloading the gene expression calls and annotations allowing to perform TopAnat analysis 

## Installation

### Install via install\_github()

``` {r}
# install the package
install.packages("devtools") # if you don't have devtools installed
library(devtools) 
install_github("wirawara/BgeeDB", build_vignettes=FALSE)
library("BgeeDB")
```

### Install via classic install

In the terminal:

    git clone https://github.com/wirawara/BgeeDB.git

In R:
``` {r}
# install the package
install.packages("./BgeeDB", repos = NULL, type="source")
library("BgeeDB")
```

## How to use BgeeDB package

### Running example for Mus musculus

#### List available species in current version of Bgee

``` {r}
library(BgeeDB)
listBgeeSpecies()
```

#### Choose the species and create a new Bgee species object

From the list above, please choose one species (e.g., "Homo\_sapiens", "Mus\_musculus",...) and platform ("rna\_seq" or "affymetrix").

``` {r}
bgee <- Bgee$new(species = "Mus_musculus", datatype = "rna_seq")
```

#### Retrieve annotation for Mus musculus RNA-seq data

The ```bgee$get_annotation()``` will output the list of experiments and libraries currently available in Bgee for RNA-seq of Mus musculus. The ```bgee$get_annotation()``` loads the annotation in R, but also creates the Mus musculus folder in your current path, where it saves the downloaded annotation locally, so you can use the annotation for later as well.

``` {r}
# the path where your folder with annotation will be saved. The folder is named after your chosen species.
getwd()
annotation_bgee_mouse <- bgee$get_annotation()
# head the experiments and libraries
lapply(annotation_bgee_mouse, head)
```

#### Download the processed RNA-seq data for Mus musculus

The ```bgee$get_data()``` will download read counts and RPKMs for Mus musculus from all available experiments in Bgee database as a list (see below). In case of downloaded data from all experiments for Mus musculus, ```bgee$get_data()``` will save the downloaded data in your current folder for later usage. 

``` {r}
# the path where your data will be saved. 
getwd()

# download all RPKMs and counts for Mus musculus
data_bgee_mouse <- bgee$get_data()

# the number of experiments downloaded from Bgee
length(data_bgee_mouse)
# check the data
sapply(data_bgee_mouse, head)
# see your first experiment
data_bgee_experiment1 <- data_bgee_mouse[[1]]
```

Alternatively, you can choose to download only one experiment from Bgee, as in the example below. The data is then saved in .tsv file in your current folder.

``` {r}
# download RPKMs and counts only for GSE30617 for Mus musculus
data_bgee_mouse_gse30617 <- bgee$get_data(experiment.id = "GSE30617")
```

The data from different samples will be listed in rows, one after the other. It is sometimes easier to work with data organized as a matrix, where different columns represent different samples. To transform the data into a matrix with genes in rows and samples in columns, you can use the ```bgee$format_data()``` function. This function also allows to filter out genes that are not called present in a given sample (gives NA values).

```{r}
# only present calls and rpkm values
gene.expression.mouse.rpkm <- bgee$format_data(data_bgee_mouse_gse30617, "present", "rpkm")
# only present calls and raw counts
gene.expression.mouse.counts <- bgee$format_data(data_bgee_mouse_gse30617, "present", "counts")
# all calls and raw counts
gene.expression.mouse.all.counts <- bgee$format_data(data_bgee_mouse_gse30617, "all", "counts")
```

#### Download the data allowing to perform GO-like enrichment test for anatomical terms

The ```loadTopAnatData()``` function loads a mapping from genes to anatomical structures based on calls of expression in anatomical structures. It also loads the structure of the anatomical ontology and the names of anatomical structures.

```{r}
myTopAnatData <- loadTopAnatData(species=10090, datatype=c("rna_seq","affymetrix","est","in_situ"))
# Loading calls of expression based on RNA-seq data only
myTopAnatData <- loadTopAnatData(species=10090, datatype="rna_seq")
# Loading high-quality calls only
myTopAnatData <- loadTopAnatData(species=10090, datatype=c("rna_seq","affymetrix","est","in_situ"), confidence="high")
# Loading calls observed in embryonic stages only
myTopAnatData <- loadTopAnatData(species=10090, datatype=c("rna_seq","affymetrix","est","in_situ"), stage="UBERON:0000068")

# Look at the data
lapply(myTopAnatData, head)
```

#### Prepare a topGO object allowing to perform GO-like enrichment test for anatomical terms, for Mus musculus

First we need to prepare a list of genes in the foreground and in the background. The input format is the same as the gene list required to build a ```topGOdata``` object in the ```topGO``` package: a vector with background genes as names, and 0 or 1 values depending if a gene is in the foreground or not. In this example we will look at genes, annotated with "spermatogenesis" in the Gene Ontology (using the ```biomaRt``` package). 

```{r}
source("https://bioconductor.org/biocLite.R")
biocLite("biomaRt")
library(biomaRt)
ensembl <- useMart("ensembl")
ensembl <- useDataset("mmusculus_gene_ensembl", mart=ensembl)

# Foreground genes are those with a GO annotation "spermatogenesis"
myGenes <- getBM(attributes= "ensembl_gene_id", filters=c("go_id"), values=list(c("GO:0007283")), mart=ensembl)

# Background are all genes with a GO annotation
universe <- getBM(attributes= "ensembl_gene_id", filters=c("with_go_go"), values=list(c(TRUE)), mart=ensembl)

# Prepares the gene list vector 
geneList <- factor(as.integer(universe[,1] %in% myGenes[,1]))
names(geneList) <- universe[,1]
head(geneList)
summary(geneList == 1)

# Prepares the topGO object
myTopAnatObject <-  topAnat(myTopAnatData, geneList)
```

*Warning*: This can be long, especially if the gene list is large, since the anatomical ontology is large and expression calls will be propagated through the whole ontology (e.g., expression in the forebrain will also be counted as expression in parent structures such as the brain, nervous system, etc). Consider running a script in batch mode if you have multiple analyses to do.

#### Launch an enrichment test for anatomical terms

This part is not dependent on the ```BgeeDB``` package, and you can readily use the ```topGO``` package and all its functionalities for this step. See the vignette of the ```topGO``` package for more details. For example:
```{r}
results <- runTest(myTopAnatObject, algorithm = 'classic', statistic = 'fisher')
# You can also use the topGO decorrelation methods, for example the "weight" method to get less redundant results
results <- runTest(myTopAnatObject, algorithm = 'weight', statistic = 'fisher')
```

*Warning*: This can be long because of the size of the ontology. Consider running a script in batch mode if you have multiple analyses to do.

#### Format the table of results after an enrichment test for anatomical terms

We built the ```makeTable``` function to filter and format the test results. Results are sorted by p-value, and FDR values are calculated. 

```{r}
# Display results sigificant at a 1% FDR threshold
tableOver <- makeTable(myTopAnatData, myTopAnatObject, results, 0.01)
# Display all results
tableOver <- makeTable(myTopAnatData, myTopAnatObject, results, 1)
```

*Warning*: it is debated if FDR correction is appropriate on enrichment test results, since tests on different terms of the ontologies are not independent. A nice discussion can be found in the vignette of the ```topGO``` package.

