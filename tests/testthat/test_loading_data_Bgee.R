context("expect_output")


test_that("Gene expression files", {
  bgee <- Bgee$new(species = "Mus_musculus", dataType = "rna_seq")
  data_bgee_mouse <- getData(bgee)


  expect_that( data_bgee_mouse, is_a("list") )
  expect_message(message("Saving all data in .rds file..."))

})
