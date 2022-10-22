
test_that("Locating Sample Table file in the package", {

  test_data_path <- system.file("extdata", "Template_sampleTable.csv",
                                package = "dPCP")

  expect_true(file.exists(test_data_path))

  })
