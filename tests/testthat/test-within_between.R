
test_that("#1 Simple Within-Between test", {
  tree = ape::read.tree(text = "((tA),(tB,(tC,tD)));")
  tree = ape::compute.brlen(tree)

  trait = c("b", "a", "a", "a")
  names(trait) = tree$tip.label

  expect_equal(
    traitheritage::within_between(
      tree = tree,
      generation_time = 0.25,
      trait = trait
    ),
    list(
      v_within = structure(
        list(
          generation = c(
            "g_0",
            "g_0.25",
            "g_0.5",
            "g_0.75",
            "g_0",
            "g_0.25",
            "g_0.5",
            "g_0.75"
          ),
          trait = c("a", "a", "a", "a", "b", "b", "b", "b"),
          op = c(0, 0, 1, 3, 0, 0, 0, 0),
          tp = c(3, 3, 3, 3, 0, 0, 0, 0),
          v_within = c(0, 0, 0.333333333333333, 1, NaN, NaN, NaN, NaN)
        ),
        class = c("data.table", "data.frame"),
        row.names = c(NA, -8L)),
      v_between = structure(
        list(
          generation = c("g_0", "g_0.25", "g_0.5", "g_0.75"),
          trait_pairs = c("ab", "ab", "ab", "ab"),
          numerator = c(0, 0, 0, 0),
          denominator = c(3L, 3L, 3L, 3L),
          v_between = c(0, 0, 0, 0)
        ),
        class = c("data.table", "data.frame"),
        row.names = c(NA, -4L)
      )
    )
  )
})


test_that("#2 Simple Within-Between test", {
  tree = ape::read.tree(text = "((tA),(tB,(tC,tD)));")
  tree = ape::compute.brlen(tree)

  trait = c("b", "b", "a", "a")
  names(trait) = tree$tip.label

  # plot(tree)
  # tiplabels(pch = 19, col = factor(trait))
  # abline(v = c(0.25, 0.5, 0.75))

  expect_equal(
    traitheritage::within_between(
      tree = tree,
      generation_time = 0.25,
      trait = trait
    ),
    list(
      v_within = structure(
        list(
          generation = c(
            "g_0",
            "g_0.25",
            "g_0.5",
            "g_0.75",
            "g_0",
            "g_0.25",
            "g_0.5",
            "g_0.75"
          ),
          trait = c("a", "a", "a", "a", "b", "b", "b", "b"),
          op = c(0, 0, 1, 1, 0, 0, 0, 0),
          tp = c(1, 1, 1, 1, 1, 1, 1, 1),
          v_within = c(0, 0, 1, 1, 0, 0, 0, 0)
        ),
        class = c("data.table", "data.frame"),
        row.names = c(NA, -8L)
      ),
      v_between = structure(
        list(
          generation = c("g_0", "g_0.25", "g_0.5", "g_0.75"),
          trait_pairs = c("ab", "ab", "ab", "ab"),
          numerator = c(0, 0, 0, 2),
          denominator = c(4L, 4L, 4L, 4L),
          v_between = c(0, 0, 0, 0.5)
        ),
        class = c("data.table", "data.frame"),
        row.names = c(NA, -4L)
      )
    )
  )
})
