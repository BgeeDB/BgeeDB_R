context("expect_output")


test_that("Annotation files", {
  bgee <- Bgee$new(species = "Mus_musculus", datatype = "rna_seq")
  annotation_bgee_mouse <- bgee$get_annotation()


  expect_that( annotation_bgee_mouse, is_a("list") )
  expect_that( length(annotation_bgee_mouse), equals(2) )
  expect_that( unique(names(annotation_bgee_mouse) %in% c("experiment_annotation", "sample_annotation")), is_true() )
  expect_message(message("Saved files "))
  # expect_true(dir.exists(file.path("..", "Mus_musculus")))
})
