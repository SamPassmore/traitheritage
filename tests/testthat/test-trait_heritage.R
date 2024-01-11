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


## Tree cut tests



test_that("#3 Simple tree cut test", {
  # make a fake tree
  t = ape::read.tree(text = "((A,B),(C,D));")
  t = ape::compute.brlen(t)

  expect_equal(slice_tree(t, 2),
               structure(
                 list(
                   generation = c("g_0", "g_0", "g_0", "g_0", "g_0.5",
                                  "g_0.5", "g_0.5", "g_0.5"),
                   taxa = c("A", "B", "C", "D", "A",
                            "B", "C", "D"),
                   clade = structure(
                     c(1L, 2L, 3L, 4L, 1L, 1L, 2L,
                       2L),
                     levels = c("1", "2", "3", "4"),
                     class = "factor"
                   )
                 ),
                 class = "data.frame",
                 row.names = c(NA,-8L)
               )
  )
})

test_that("#4 Simple tree cut test", {
  # make a fake tree
  t = ape::read.tree(text = "((A),(B,(C,D)));")
  t = ape::compute.brlen(t)

  expect_equal(slice_tree(t, 3),
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
                     "g_0.67",
                     "g_0.67",
                     "g_0.67",
                     "g_0.67"
                   ),
                   taxa = c("A", "B", "C", "D", "A", "B", "C", "D", "A", "B",
                            "C", "D"),
                   clade = structure(
                     c(1L, 2L, 3L, 4L, 1L, 2L, 3L, 3L,
                       1L, 2L, 2L, 2L),
                     levels = c("1", "2", "3", "4"),
                     class = "factor"
                   )
                 ),
                 class = "data.frame",
                 row.names = c(NA,-12L)
               ))
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




