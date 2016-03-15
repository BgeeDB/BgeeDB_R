context("expect_output")

test_that("Distinct species", {

	species <- listBgeeSpecies()

  	expect_that( species, is_a("list") )
    expect_that( length(species), equals(2) )
    expect_that( unique(names(species) %in% c("rna_seq", "affymetrix")), is_true() )
})
