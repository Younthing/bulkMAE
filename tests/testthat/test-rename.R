test_that("bulkmae_rename maps recorded 0.3 names and rejects unknowns", {
  expect_identical(bulkmae_rename("run_gsva"), "score_gsva()")
  expect_identical(bulkmae_rename("run_gsva()"), "score_gsva()")
  expect_identical(bulkmae_rename("run_go_ora"), "enrich_go(method = \"ora\")")
  expect_error(bulkmae_rename("score_gsva"), "already a 0.4 name")
  expect_error(bulkmae_rename("not_a_bulkmae_function"), "not a recorded pre-0.4 name")
  expect_error(bulkmae_rename(c("run_gsva", "run_pca")), "single function name")

  map <- bulkmae_rename()
  expect_s3_class(map, "data.frame")
  expect_identical(names(map), c("old_name", "new_name"))
  expect_gt(nrow(map), 50L)
  expect_false(anyDuplicated(map$old_name) > 0L)
})
