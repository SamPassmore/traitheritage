test_that("#5 Full pipeline simple test", {
  set.seed(123)

  tree = ape::read.tree(text = "(tA,(tB,(tC,tD)));")
  tree = ape::compute.brlen(tree)

  trait = c("b", "a", "a", "a")
  names(trait) = tree$tip.label

  ## Because A & B are singletons, they are not counted.
  ## Then CD are the same so the probability of a shared trait is 1
  result = permute_trait_heritage(tree, trait, generation_time = 0.2, n_permutations = 1)

  expect_equal(result$summary$generation, c(0.2, 0.4, 0.6, 0.8, 1.0))
  expect_equal(result$summary$numerator_sum, c(0, 0, 0, 1, 3))
  expect_equal(result$summary$denominator_sum, c(0, 1, 1, 3, 6))
  expect_equal(round(result$summary$clade_probability, 5), c(NaN, 0, 0, 0.33333, 0.5))
})


test_that("#6 Pipeline NA test", {
  t = ape::read.tree(text = "(tA,(tB,(tC,tD)));")
  t = ape::compute.brlen(t)

  trait = c("b", "a", NA, "a")
  names(trait) = t$tip.label

  expect_error(permute_trait_heritage(t, trait, generation_time = 0.2))

})


test_that("#6 Matching names test", {
  t = ape::read.tree(text = "(tA,(tB,(tC,tD)));")
  t = ape::compute.brlen(t)

  trait = c("b", "a", "a", "a")
  names(trait) = c("tA", "tB", "AnotherTaxa", "tD")

  expect_error(permute_trait_heritage(t, trait, generation_time = 0.2))

})

test_that("#6 No matches", {
  t = ape::read.tree(text = "(tA,(tB,(tC,tD)));")
  t = ape::compute.brlen(t)

  trait = c("b", "a", "c", "d")
  names(trait) = c("tA", "tB", "tC", "tD")

  result = permute_trait_heritage(t, trait, generation_time = 0.2, n_permutations = 1)

  expect_equal(result$summary$clade_probability, c(0, 0, 0, 0, 0))

})

test_that("#7 Correct Iterations", {
  t = ape::read.tree(text = "(tA,(tB,(tC,tD)));")
  t = ape::compute.brlen(t)

  trait = c("b", "a", "c", "d")
  names(trait) = c("tA", "tB", "tC", "tD")

  result = permute_trait_heritage(t, trait, generation_time = 0.2, n_permutations = 2)

  expect_equal(result$summary$iteration,
               c("p_1","p_1","p_1","p_1","p_1","p_2","p_2","p_2","p_2","p_2"))

})
