context("expect_output")


test_that("Annotation files", {
  bgee <- Bgee$new(species = "Mus_musculus", dataType = "rna_seq")
  annotation_bgee_mouse <- getAnnotation()


  expect_that( annotation_bgee_mouse, is_a("list") )
  expect_that( length(annotation_bgee_mouse), equals(2) )
  expect_that( unique(names(annotation_bgee_mouse) %in% c("experiment.annotation", "sample.annotation")), is_true() )
})
