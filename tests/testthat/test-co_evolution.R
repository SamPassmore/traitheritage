## trait co-evolution

test_that("Simple test", {
  t = ape::read.tree(text = "((A,B),(C,D));")
  tree = ape::compute.brlen(t)

  trait = c("b", "a", "a", "a")
  names(trait) = t$tip.label


  distances = matrix(
    c(1, 2,
      2, 1,
      5, 6,
      6, 5),
    byrow = TRUE,
    ncol = 2,
    dimnames = list(t$tip.label, c("X", "Y"))
  )
  distance_matrix = as.matrix(dist(distances))

  cut_off = 2
  generation_time = 0.5

  result = trait_coevolution(tree, trait, distance_matrix, generation_time, cut_off)

  expect_equal(
    result$time_df$denominator_sum,
    c(0, 4))
  expect_equal(
    result$time_df$p_lang_dist,
    c(0, 0.25))
  expect_equal(
    result$pairs_df$dist_trait,
    c(T, T, F, F, F, F, T, T))
})
