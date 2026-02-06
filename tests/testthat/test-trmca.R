test_that("Single TMRCA call", {
  skip("Experimental function — under active development")

  set.seed(1234)
  tree = ape::rcoal(40, br = runif)
  states = ape::rTraitDisc(tree)

  out = .tmrca_onetree(tree, states)

  expect_equal(out$Node, c("56", "59", "62", "67", "76", "77"))

  expect_equal(out, structure(
    list(
      Node = c("56", "59", "62", "67", "76", "77"),
      Date = c(10.257, 9.755, 9.02, 6.414, 1.969, 1.345),
      n_tips = c(2L,
                 6L, 8L, 4L, 2L, 2L),
      tips = c(
        "t19; t22",
        "t12; t20; t27; t3; t33; t8",
        "t10; t11; t18; t25; t26; t28; t38; t5",
        "t31; t32; t35; t7",
        "t21; t39",
        "t34; t6"
      ),
      state = structure(
        c(2L, 1L, 1L, 2L,
          1L, 1L),
        levels = c("A", "B"),
        class = "factor"
      )
    ),
    row.names = c(NA,-6L),
    class = "data.frame"
    )
  )

})


test_that("Looped TMRCA call", {
  set.seed(1234)
  tree = ape::rcoal(40, br = runif)
  states = ape::rTraitDisc(tree)

  trees = rep(tree, 10)
  class(trees) = "multiPhylo"

  out = tmrca(trees, states)

  expect_equal(unique(out$Node), c("56", "59", "62", "67", "76", "77"))

})
