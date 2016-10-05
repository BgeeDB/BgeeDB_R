context("expect_output")


test_that("Single experiment gene expression files", {
  bgee <- Bgee$new(species = "Mus_musculus", dataType = "rna_seq")
  data_bgee_mouse_exp <- bgee$get_data(experiment.id  = "GSE30617")

  expect_that( data_bgee_mouse_exp, is_a("data.frame") )
})
