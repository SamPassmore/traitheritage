# Test Distance

test_that("#1. Simple distance test", {

  t = ape::read.tree(text = "((A,B),(C,D));")
  t = ape::compute.brlen(t)

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

  result = distance_trait_heritage(
    t,
    distance_matrix = distance_matrix,
    generation_time = generation_time,
    cut_off = cut_off
  )

  expect_equal(result$clade_probability[2], 1)
  expect_equal(result$numerator_sum[2], 2)
  expect_equal(result$denominator_sum[2], 2)
})

test_that("#2. Simple distance test", {

  t = ape::read.tree(text = "((A),(B,C,D));")
  t = ape::compute.brlen(t)

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
  generation_time = 0.25

  result = distance_trait_heritage(
    t,
    distance_matrix = distance_matrix,
    generation_time = generation_time,
    cut_off = cut_off
  )

  expect_equal(result$clade_probability, c(NaN, NaN, NaN, 0.333333333333333))
  expect_equal(result$numerator_sum, c(0, 0, 0, 1))
  expect_equal(result$denominator_sum, c(0, 0, 0, 3))
})

test_that("#2. Simple Error test", {

  t = ape::read.tree(text = "((A),(B,C,D));")
  t = ape::compute.brlen(t)

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

  rownames(distance_matrix) = c("A", "X", "C", "D")

  cut_off = 2
  n_generations = 3

  expect_error(distance_trait_heritage(
    t,
    distance_matrix = distance_matrix,
    n_generations = n_generations,
    cut_off = cut_off
  ))
})
