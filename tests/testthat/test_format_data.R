context("expect_output")


test_that("Formatting gene expression files", {
  bgee <- Bgee$new(species = "Mus_musculus", datatype = "rna_seq")
  data_bgee_mouse_exp <- bgee$get_data(experiment.id  = "GSE30617")
  gene.expression <- bgee$format_data(data_bgee_mouse_exp, calltype = "expressed", stats = "rpkm")

  expect_that( gene.expression, is_a("data.frame") )
  expect_true( colnames(gene.expression)[1] == "Gene ID" )

})
