context("expect_output")


test_that("Annotation files", {
  bgee <- Bgee$new(species = "Mus_musculus", dataType = "rna_seq")
  annotation_bgee_mouse <- getAnnotation(bgee)
  expect_that( annotation_bgee_mouse, is_a("list") )
  expect_that( length(annotation_bgee_mouse), equals(2) )
  expect_true( unique(names(annotation_bgee_mouse) %in% c("experiment.annotation", "sample.annotation")))
  bgee <- Bgee$new(species = "Mus_musculus", dataType = "affymetrix")
  annotation_bgee_mouse <- getAnnotation(bgee)
  expect_that( annotation_bgee_mouse, is_a("list") )
  
  #test loading old version of annotation (bgee version < 14.0)
  bgee_v13 <- Bgee$new(species = "Mus_musculus", dataType = "rna_seq", release = "13.2")
  annotation_bgee_mouse_v13 <- getAnnotation(bgee_v13)
  expect_that( annotation_bgee_mouse_v13, is_a("list") )
  bgee_v13 <- Bgee$new(species = "Mus_musculus", dataType = "affymetrix", release = "13.2")
  annotation_bgee_mouse_v13 <- getAnnotation(bgee_v13)
  expect_that( annotation_bgee_mouse_v13, is_a("list") )
})
