# Posterior Tests

test_that("#1. Simple posterior test", {

  t = ape::read.tree(text = "((A,B),(C,D));")
  t = ape::compute.brlen(t)

  trees = list(t, t, t, t)
  class(trees) = "multiPhylo"

  trait = c("a", "a", "b", "b")
  names(trait) = c("A", "B", "C", "D")

  result = posterior_trait_heritage(trees, trait, generation_time = 0.5)

  expect_equal(result$clade_probability, c(NaN, 1, NaN, 1, NaN, 1, NaN, 1))
})

test_that("#1. Second simple posterior test", {

  t = ape::read.tree(text = "((A,B),(C,D));")
  t = ape::compute.brlen(t)

  trees = list(t, t, t, t)
  class(trees) = "multiPhylo"

  trait = c("a", "b", "b", "b")
  names(trait) = c("A", "B", "C", "D")

  result = posterior_trait_heritage(trees, trait, generation_time =  0.5)

  expect_equal(result$clade_probability, c(NaN, 0.5, NaN, 0.5, NaN, 0.5, NaN, 0.5))
})
