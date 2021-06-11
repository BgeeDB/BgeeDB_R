context("expect_output")

test_that("Distinct species", {
  species <- listBgeeSpecies(release = "14_2")
	expect_that( species, is_a("data.frame") )
  expect_that( ncol(species), equals(8) )
  expect_message(message("Query to Bgee webservice successful!"))
  expect_true( all(colnames(species) %in% c("ID", "GENUS",
                                            "SPECIES_NAME", "COMMON_NAME",
                                            "AFFYMETRIX", "EST",
                                            "IN_SITU", "RNA_SEQ")))
})
