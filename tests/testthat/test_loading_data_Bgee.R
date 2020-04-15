context("expect_output")


test_that("Gene expression files", {
  bgee <- Bgee$new(species = "Drosophila_simulans", dataType = "rna_seq")
  data_bgee_mouse <- getData(bgee)


  expect_that( data_bgee_mouse, is_a("data.frame") )
  expect_message(message("Save data in local sqlite database"))
  expect_message(message("Load queried data. The query is : SELECT * from rna_seq"))

})
