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
  n_generations = 2

  result = distance_trait_heritage(
    t,
    distance_matrix = distance_matrix,
    n_generations = n_generations,
    cut_off = cut_off
  )

  expect_equal(result$clade_probability, 1)
  expect_equal(result$numerator_sum, 2)
  expect_equal(result$denominator_sum, 2)
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
  n_generations = 3

  result = distance_trait_heritage(
    t,
    distance_matrix = distance_matrix,
    n_generations = n_generations,
    cut_off = cut_off
  )

  expect_equal(result$clade_probability, c(NaN, 0.333333333333333))
  expect_equal(result$numerator_sum, c(0, 1))
  expect_equal(result$denominator_sum, c(0, 3))
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
