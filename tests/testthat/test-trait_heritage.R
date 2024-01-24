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

  expect_equal(.slice_tree(t, 2),
               structure(
                 list(
                   generation = c("g_0.5", "g_0.5", "g_0.5", "g_0.5"),
                   taxa = c("A", "B", "C", "D"),
                   clade = c(1, 1, 2, 2)
                 ),
                 row.names = c(NA,-4L),
                 class = "data.frame"
               ))
})

test_that("#4 Simple tree cut test", {
  # make a fake tree
  t = ape::read.tree(text = "((A),(B,(C,D)));")
  t = ape::compute.brlen(t)

  expect_equal(.slice_tree(t, 3),
               structure(
                 list(
                   generation = c(
                     "g_0.33",
                     "g_0.33",
                     "g_0.33",
                     "g_0.33",
                     "g_0.67",
                     "g_0.67",
                     "g_0.67",
                     "g_0.67"
                   ),
                   taxa = c("A", "B", "C",
                            "D", "A", "B", "C", "D"),
                   clade = c(1, 2, 3, 3, 1, 2, 2, 2)
                 ),
                 row.names = c(NA,-8L),
                 class = "data.frame"
               ))
})

## Main function tests

test_that("#5 Full pipeline simple test", {
  t = ape::read.tree(text = "((tA),(tB,(tC,tD)));")
  t = ape::compute.brlen(t)

  trait = c("b", "a", "a", "a")
  names(trait) = t$tip.label

  ## Because A & B are singletons, they are not counted.
  ## Then CD are the same so the probability of a shared trait is 1
  expect_equal(trait_heritage(t, trait, n_generations = 2)$clade_probability,
               1)

})

# plot(t)
# ape::axisPhylo()
# abline(v = 0.6)
# abline(v = 0.8)
#
# phyloregion::get_clades(t, cut = 1.1)
# phyloregion::get_clades(t, cut = 1.0)
# phyloregion::get_clades(t, cut = 0.8)
# phyloregion::get_clades(t, cut = 0.6)
# phyloregion::get_clades(t, cut = 0.4)
# phyloregion::get_clades(t, cut = 0.2)
# phyloregion::get_clades(t, cut = 0.0)




