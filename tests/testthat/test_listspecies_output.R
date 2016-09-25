context("expect_output")

test_that("Distinct species", {
  species <- listBgeeSpecies()
	expect_that( species, is_a("data.frame") )
  expect_that( ncol(species), equals(8) )
  expect_message(message("Querying Bgee to get release information..."))
  expect_true( colnames(species) %in% c("ID", "GENUS",
                                        "SPECIES_NAME", "COMMON_NAME",
                                        "AFFYMETRIX", "EST",
                                        "IN_SITU", "RNA_SEQ"))
})
