## Clade probabilities test

test_that("#1 Simple clade calculation works",{
  clade_states = c("A", "A")
  names(clade_states) = c("t1", "t2")
  expect_equal(.clade_probabilities(clade_states), list(numerator = 1, denominator = 1))
})

test_that("#2 Simple clade calculation works",{
  clade_states = c("A", "A", "B")
  names(clade_states) = c("t1", "t2", "t3")
  expect_equal(.clade_probabilities(clade_states), list(numerator = 1, denominator = 3))
})

test_that("#2 Many values clade calculation works",{
  clade_states = c("A", "A", "B", "C", "C")
  names(clade_states) = c("t1", "t2", "t3", "t4", "t5")
  expect_equal(.clade_probabilities(clade_states), list(numerator = 2, denominator = 10))
})

test_that("#3 A clade with one taxa should return all zeros",{
  clade_states = c("A")
  names(clade_states) = c("t1")
  expect_equal(.clade_probabilities(clade_states), list(numerator = 0, denominator = 0))
})


## Tree cut tests
test_that("#3 Simple tree cut test", {
  # make a fake tree
  t = ape::read.tree(text = "((A,B),(C,D));")
  t = ape::compute.brlen(t)

  expect_equal(.slice_tree(t, 0.33),
               structure(
                 list(
                   generation = c(
                     "g_0",
                     "g_0",
                     "g_0",
                     "g_0",
                     "g_0.33",
                     "g_0.33",
                     "g_0.33",
                     "g_0.33",
                     "g_0.66",
                     "g_0.66",
                     "g_0.66",
                     "g_0.66",
                     "g_0.99",
                     "g_0.99",
                     "g_0.99",
                     "g_0.99"
                   ),
                   taxa = c(
                     "A",
                     "B",
                     "C",
                     "D",
                     "A",
                     "B",
                     "C",
                     "D",
                     "A",
                     "B",
                     "C",
                     "D",
                     "A",
                     "B",
                     "C",
                     "D"
                   ),
                   clade = c(1, 2, 3, 4, 1, 2, 3, 4, 1, 1, 2, 2, 1, 1, 2, 2)
                 ),
                 row.names = c(NA,-16L),
                 class = c("data.table", "data.frame")
               ))
})

test_that("#4 Simple tree cut test", {
  # make a fake tree
  t = ape::read.tree(text = "((A),(B,(C,D)));")
  t = ape::compute.brlen(t)

  expect_equal(.slice_tree(t, 0.33),
               structure(
                 list(
                   generation = c(
                     "g_0",
                     "g_0",
                     "g_0",
                     "g_0",
                     "g_0.33",
                     "g_0.33",
                     "g_0.33",
                     "g_0.33",
                     "g_0.66",
                     "g_0.66",
                     "g_0.66",
                     "g_0.66",
                     "g_0.99",
                     "g_0.99",
                     "g_0.99",
                     "g_0.99"
                   ),
                   taxa = c(
                     "A",
                     "B",
                     "C",
                     "D",
                     "A",
                     "B",
                     "C",
                     "D",
                     "A",
                     "B",
                     "C",
                     "D",
                     "A",
                     "B",
                     "C",
                     "D"
                   ),
                   clade = c(1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 3, 1, 2, 2, 2)
                 ),
                 row.names = c(NA,-16L),
                 class =  c("data.table", "data.frame")
               ))
})

## Main function tests

test_that("#5 Full pipeline simple test", {
  tree = ape::read.tree(text = "((tA),(tB,(tC,tD)));")
  tree = ape::compute.brlen(tree)

  trait = c("b", "a", "a", "a")
  names(trait) = tree$tip.label

  ## Because A & B are singletons, they are not counted.
  ## Then CD are the same so the probability of a shared trait is 1
  expect_equal(trait_heritage(tree, trait, generation_time = 0.2),
               structure(
                 list(
                   generation = c("g_0", "g_0.2", "g_0.4", "g_0.6",
                                  "g_0.8"),
                   numerator_sum = c(0, 0, 1, 1, 3),
                   denominator_sum = c(0,
                                       0, 1, 1, 3),
                   clade_probability = c(NaN, NaN, 1, 1, 1)
                 ),
                 class = c("data.table", "data.frame"),
                 row.names = c(NA,-5L)
               ))
})


test_that("#6 Pipeline NA test", {
  t = ape::read.tree(text = "((tA),(tB,(tC,tD)));")
  t = ape::compute.brlen(t)

  trait = c("b", "a", NA, "a")
  names(trait) = t$tip.label

  expect_error(trait_heritage(t, trait, generation_time = 0.2))

})


test_that("#6 Matching names test", {
  t = ape::read.tree(text = "((tA),(tB,(tC,tD)));")
  t = ape::compute.brlen(t)

  trait = c("b", "a", "a", "a")
  names(trait) = c("tA", "tB", "AnotherTaxa", "tD")

  expect_error(trait_heritage(t, trait, generation_time = 0.2))

})

test_that("#6 No matches", {
  t = ape::read.tree(text = "((tA),(tB,(tC,tD)));")
  t = ape::compute.brlen(t)

  trait = c("b", "a", "c", "d")
  names(trait) = c("tA", "tB", "tC", "tD")

  result = trait_heritage(t, trait, generation_time = 0.2)

  expect_equal(result$clade_probability, c(NaN, NaN, 0, 0, 0))

})


