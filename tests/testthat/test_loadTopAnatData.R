context("expect_output")


test_that("Loading of topAnatData files is working", {
  bgee <- Bgee$new(species="Mus_musculus", dataType="rna_seq")
  myTopAnatData <- loadTopAnatData(bgee, stage="UBERON:0000068")

  expect_that( myTopAnatData, is_a("list") )
  expect_that( length(myTopAnatData), equals(4) )
  expect_that( myTopAnatData$gene2anatomy, is_a("array") )
  expect_that( myTopAnatData$organ.relationships, is_a("list") )
  expect_that( myTopAnatData$organ.names, is_a("data.frame") )
  expect_that( myTopAnatData$bgee.object, is_a("Bgee") )

  ## TO DO: test that all organs in gene2anatomy are in organ.names
  ## expect_true( sum(myTopAnatData$gene2anatomy %in% myTopAnatData$organ.names[,1]) == length(myTopAnatData$gene2anatomy) )

  expect_message( message("Done.") )

})
