# Posterior Tests

test_that("#1. Simple posterior test", {

  t = ape::read.tree(text = "((A,B),(C,D));")
  t = ape::compute.brlen(t)

  trees = list(t, t, t, t)
  class(trees) = "multiPhylo"

  trait = c("a", "a", "b", "b")
  names(trait) = c("A", "B", "C", "D")

  result = posterior_trait_heritage(trees, trait, n_generations = 2)

  expect_equal(result$clade_probability, c(1, 1, 1, 1))

})
