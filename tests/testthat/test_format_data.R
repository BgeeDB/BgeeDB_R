context("expect_output")


test_that("Formatting gene expression files", {
  bgee <- Bgee$new(species = "Mus_musculus", dataType = "rna_seq")
  data_bgee_mouse_exp <- getData(bgee, experimentId  = "GSE30617")
  gene.expression <- formatData(bgee, data_bgee_mouse_exp, callType = "present", stats = "rpkm")

  expect_true(class(gene.expression) == "ExpressionSet" )
})
