context("expect_output")


test_that("Formatting gene expression files", {
  bgee <- Bgee$new(species = "Mus_musculus", dataType = "rna_seq")
  data_bgee_mouse_exp <- getData(bgee, experimentId  = "GSE43721")
  gene.expression <- formatData(bgee, data_bgee_mouse_exp, callType = "present", stats = "tpm")
  expect_true(class(gene.expression) == "ExpressionSet" )
  gene.expression <- formatData(bgee, data_bgee_mouse_exp, callType = "present", stats = "counts")
  expect_true(class(gene.expression) == "ExpressionSet" )
  # not available since bgee_v15_1
  expect_error(formatData(bgee, data_bgee_mouse_exp, callType = "present", stats = "fpkm"),
               "Choose whether data formatting should create a matrix of TPMs or read counts, with stats option set as \"tpm\" or \"counts\"")
  # not available since bgee_v14
  expect_error(formatData(bgee, data_bgee_mouse_exp, callType = "present", stats = "rpkm"), 
               "Choose whether data formatting should create a matrix of FPKMs, TPMs or read counts, with stats option set as \"fpkm\", \"tpm\" or \"counts\"")
  #test formatting for files generated before bgee_v14
  bgee_v13 <- Bgee$new(species = "Mus_musculus", dataType = "rna_seq", release = "13_2")
  data_bgee_mouse_exp_v13 <- getData(bgee_v13, experimentId  = "GSE43721")
  expect_error(formatData(bgee_v13, data_bgee_mouse_exp_v13, callType = "present", stats = "tpm"), 
               "Choose whether data formatting should create a matrix of RPKMs or read counts, with stats option set as \"rpkm\" or \"counts\"")
  expect_error(formatData(bgee_v13, data_bgee_mouse_exp_v13, callType = "present", stats = "fpkm"), 
               "Choose whether data formatting should create a matrix of RPKMs or read counts, with stats option set as \"rpkm\" or \"counts\"")  
  gene.expression <- formatData(bgee_v13, data_bgee_mouse_exp_v13, callType = "present", stats = "counts")
  expect_true(class(gene.expression) == "ExpressionSet" )
  gene.expression <- formatData(bgee_v13, data_bgee_mouse_exp_v13, callType = "present", stats = "rpkm")
  expect_true(class(gene.expression) == "ExpressionSet" )

  
})



