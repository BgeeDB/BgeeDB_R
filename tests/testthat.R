library(testthat)
library(BgeeDB)

Sys.setenv("R_TESTS" = "")
test_check("BgeeDB")
