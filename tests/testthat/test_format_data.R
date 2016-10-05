context("expect_output")


test_that("Formatting gene expression files", {
  bgee <- Bgee$new(species = "Mus_musculus", dataType = "rna_seq")
  data_bgee_mouse_exp <- bgee$get_data(experiment.id  = "GSE30617")
  gene.expression <- bgee$format_data(data_bgee_mouse_exp, calltype = "present", stats = "rpkm")

  expect_true(class(gene.expression) == "ExpressionSet" )
})
