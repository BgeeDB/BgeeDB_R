# install and load BgeeDB
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require("BgeeDB", quietly = TRUE))
  BiocManager::install("BgeeDB")
if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr")
if (!require("BgeeDB", quietly = TRUE))
  install.packages("ggfortify", repos = "http://cran.us.r-project.org")

# to update depending on wd
source("./R/retrieveExpression.R")
library(BgeeDB)
library(dplyr)
library(ggfortify)

############################# MAIN CODE RETRIEVE ONE-TO-ONE ORTHOLOGS ################################
####################################### AND GENERATE PCA #############################################
# retrieve one to one orthologs
orthologs_one_to_one <- listOrthologs(mandatorySpecies = c(7955, 7918, 8090, 7936, 8010, 7994, 8049), onlyOneToOne = TRUE)

# retrieve expression for each phylofish experiments present in Bgee
#TODO for now it is not possible to generate the expression data.frame by calling only once the function retrieveOrthologsExpression
# Could be useful to update the function to be able to retrieve expression with one single call like :
# expression_all <- retrieveOrthologsExpression(experimentId = c("SRP044782","SRP044781", "SRP044784", "SRP045099", "SRP045141", "SRP058863", "SRP058865"),
#   speciesIds = c(7918, 7955, 8090, 7936, 8010, 7994, 8049), orthologs = orthologs_one_to_one)

# For now we have to call per species :
expression_7918 <- retrieveOrthologsExpression(experimentId = "SRP044782", speciesIds = 7918, orthologs = orthologs_one_to_one)
expression_7955 <- retrieveOrthologsExpression(experimentId = "SRP044781", speciesIds = 7955, orthologs = orthologs_one_to_one)
expression_8090 <- retrieveOrthologsExpression(experimentId = "SRP044784", speciesIds = 8090, orthologs = orthologs_one_to_one)
expression_7936 <- retrieveOrthologsExpression(experimentId = "SRP045099", speciesIds = 7936, orthologs = orthologs_one_to_one)
expression_8010 <- retrieveOrthologsExpression(experimentId = "SRP045141", speciesIds = 8010, orthologs = orthologs_one_to_one)
expression_7994 <- retrieveOrthologsExpression(experimentId = "SRP058863", speciesIds = 7994, orthologs = orthologs_one_to_one)
expression_8049 <- retrieveOrthologsExpression(experimentId = "SRP058865", speciesIds = 8049, orthologs = orthologs_one_to_one)
# And then merge the results
expression_phylofish_bgee <- rbind(expression_7918, expression_7955, expression_8090,
    expression_7936, expression_7994, expression_8049)


# transform table result to be compatible with the input format used to generate the PCA plot
all_exression_transformed <- as.data.frame(matrix(nrow = 0, ncol = nrow(orthologs_one_to_one) + 3))
for (speciesId in unique(expression_phylofish_bgee$Species.ID)) {
  for (libraryId in unique(expression_phylofish_bgee$Library.ID[expression_phylofish_bgee$Species.ID == speciesId])) {

    transformed <- as.data.frame(t(expression_phylofish_bgee[expression_phylofish_bgee$Species.ID == speciesId & expression_phylofish_bgee$Library.ID== libraryId, c("TPM", "Gene.Family")]))
    colnames(transformed) <- as.integer(transformed["Gene.Family",])
    transformed <- transformed[!row.names(transformed) %in% "Gene.Family",]
    transformed$libraryId <- libraryId
    transformed$speciesId <- speciesId
    transformed$anatEntityName <- unique(expression_phylofish_bgee$Anatomical.entity.name[expression_phylofish_bgee$Species.ID == speciesId
        & expression_phylofish_bgee$Library.ID== libraryId])
    if(nrow(all_exression_transformed) == 0) {
      all_exression_transformed <- transformed
    } else {
      all_exression_transformed <- dplyr::bind_rows(all_exression_transformed, transformed)
    }
  }
}
#reorder columns
all_exression_transformed <- all_exression_transformed[, c("libraryId", "speciesId", "anatEntityName", seq_len(nrow(orthologs_one_to_one)))]

#remove gene family for which some species have NA TPM value
#TODO : check where these NA come from. Could be a bug
expression_dataset <- all_exression_transformed %>% select_if(~ !any(is.na(.)))

#filter on organ names
expression_filtered <- expression_dataset[expression_dataset$anatEntityName %in% c("\"brain\"", "\"testis\"", "\"heart\"", "\"intestine\""),]

# generate PCA
expression_filtered$speciesId <- as.character(expression_filtered$speciesId)

pca_res <- prcomp(expression_filtered[4:length(expression_filtered)], scale. = TRUE)

jpeg(file="~/Downloads/Species_PCA.jpeg")
autoplot(pca_res, data = expression_filtered, colour = 'speciesId')
dev.off()
jpeg(file="~/Downloads/Organ_PCA.jpeg")
autoplot(pca_res, data = expression_filtered, colour = 'anatEntityName')
dev.off()

############################# MAIN CODE 1-TO-TWO ORTHOLOGS ################################

## retrieve orthologs for all expected species
reference_species <- 7918
gene_duplicated_species <- c(7955)
one_to_many_orthologs <- listOrthologs(referenceSpecies = reference_species, mandatorySpecies = gene_duplicated_species)

# now only keep one to two orthologs
one_to_two_df <- formatOrthologsOneToN(orthologs = one_to_many_orthologs, n = 2, referenceSpecies = reference_species,
  geneDuplicatedSpecies = gene_duplicated_species)

# retrieve expression data for 1-to-2 orthologs only
expression_7955_one_to_two <- retrieveOrthologsExpression(experimentId="SRP044781", speciesIds = c(7955), orthologs = one_to_two_df)
expression_7918_one_to_two <- retrieveOrthologsExpression(experimentId = "SRP044782", speciesIds = c(7918), orthologs = one_to_two_df)
expression_gar_zebrafish_bgee <- rbind(expression_7918_one_to_two, expression_7955_one_to_two)


# transform table result to be compatible with the input format used to generate the PCA plot
expression_gar_zebrafish_bgee_transformed <- NULL
group_by_gene_family <- expression_gar_zebrafish_bgee %>% arrange(Gene.ID) %>% group_by(Gene.Family) %>% group_map(~ .x, .keep = TRUE)
for (i in seq_len(length(group_by_gene_family))) {
  gene_family_expression_transformed <- group_by_gene_family[[i]][,c("Library.ID", "Species.ID", "Anatomical.entity.name", "TPM")]
  colnames(gene_family_expression_transformed)[4] <- group_by_gene_family[[i]][1,"Gene.Family"]
  #TODO quick and dirty that only works for 1-to-2 orthologs.
  # It has to be improved to be compatible with all values of n
  gene_family_expression_transformed$geneDuplicate <- ifelse(duplicated(gene_family_expression_transformed[,1:3]), 2, 1)
  if(is.null(expression_gar_zebrafish_bgee_transformed)) {
    expression_gar_zebrafish_bgee_transformed <- gene_family_expression_transformed
  } else {
    expression_gar_zebrafish_bgee_transformed <- left_join(expression_gar_zebrafish_bgee_transformed, gene_family_expression_transformed, by = c("Library.ID", "Species.ID", "Anatomical.entity.name", "geneDuplicate"))
  }
}

#reorder columns
expression_gar_zebrafish_bgee_transformed <- as.data.frame(expression_gar_zebrafish_bgee_transformed[, c("Library.ID", "Species.ID", "Anatomical.entity.name", "geneDuplicate", seq_len(nrow(one_to_two_df)))])

#remove gene family for which some species have NA TPM value
expression_gar_zebrafish_dataset <- expression_gar_zebrafish_bgee_transformed %>% select_if(~ !any(is.na(.)))

#filter on organ names
expression_gar_zebrafish_filtered <- expression_gar_zebrafish_dataset[expression_gar_zebrafish_dataset$Anatomical.entity.name %in% c("\"brain\"", "\"testis\"", "\"heart\"", "\"intestine\""),]
# update type of Species.ID column from numeric to character
expression_gar_zebrafish_filtered$Species.ID <- as.character(expression_filtered$Species.ID)

# generate PCA
pca_gar_zebrafish_res <- prcomp(expression_gar_zebrafish_filtered[5:length(expression_gar_zebrafish_filtered)], scale. = TRUE)
jpeg(file="~/Downloads/Species_PCA.jpeg")
autoplot(pca_gar_zebrafish_res, data = expression_gar_zebrafish_filtered, colour = 'Species.ID')
dev.off()
jpeg(file="~/Downloads/Organ_PCA.jpeg")
autoplot(pca_gar_zebrafish_res, data = expression_gar_zebrafish_filtered, colour = 'Anatomical.entity.name')
dev.off()
