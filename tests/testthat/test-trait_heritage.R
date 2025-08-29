
test_that("#5 Full pipeline simple test", {
  tree = ape::read.tree(text = "(tA,(tB,(tC,tD)));")
  tree = ape::compute.brlen(tree)

  trait = c("b", "a", "a", "a")
  names(trait) = tree$tip.label

  ## Because A & B are singletons, they are not counted.
  ## Then CD are the same so the probability of a shared trait is 1
  out = trait_heritage(tree, trait, generation_time = 0.2)
  expect_equal(
    out$by_trait,
    structure(
      list(
        generation = c(0.2, 0.2, 0.4, 0.4, 0.6, 0.6, 0.8, 0.8, 1.0, 1.0),
        state = c("a", "b", "a", "b", "a", "b", "a", "b", "a", "b"),
        numerator_sum = c(0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 3.0, 0.0, 3.0, 0.0),
        denominator_sum = c(0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 3.0, 3.0, 6, 6),
        clade_probability = c(0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.5, 0.0)
      ),
      sorted = c("generation", "state"),
      row.names = c(NA, -10L),
      class = c("data.table", "data.frame")
    ))

  expect_equal(
    out$summary,
    structure(
      list(
        generation = c(0.2, 0.4, 0.6, 0.8, 1.0),
        numerator_sum = c(0, 1, 1, 3, 3),
        denominator_sum = c(0, 1, 1, 3, 6),
        clade_probability = c(NaN, 1, 1, 1, 0.5)
      ),
      row.names = c(NA, -5L),
      class = c("data.table", "data.frame"),
      sorted = "generation"
    )
  )
})

test_that("Simple test for trait_heritage_specific", {
  tree = ape::read.tree(text = "((tA:0.2,tX:0.2):0.8,(tB:0.65,(tC:0.3,tD:0.3):0.35):0.35);")

  trait = c("b", "b", "b", "a", "a")
  names(trait) = tree$tip.label

  ## Because A & B are singletons, they are not counted.
  ## Then CD are the same so the probability of a shared trait is 1
  out = trait_heritage_specific(tree, trait, generation_time = 0.2, condition = "a")
  expect_equal(
    out$by_trait,
    structure(
      list(
        generation = c(0.0, 0.2,0.4, 0.6, 0.8, 1.0),
        state = c("a", "a", "a","a","a", "a"),
        numerator_sum = c(0.0, 0.0, 1.0, 1.0, 1.0, 1.0),
        denominator_sum = c(0.0, 0.0, 1.0, 1.0, 3.0, 10),
        clade_probability = c(0.0, 0.0, 1.0, 1.0, 1/3, 0.1)
      ),
      sorted = c("generation", "state"),
      row.names = c(NA, -6L),
      class = c("data.table", "data.frame")
    ))
})

test_that("Reintroduce NAs in trait specific", {
  tree = ape::read.tree(text = "(((t1:0.21,t2:0.21):0.3,(t3:0.25, t4:0.25):0.25):0.5,(t5:0.41, t6:0.41):0.6);")

  trait = c("a", "a", "b", "b", "a", "a")
  names(trait) = tree$tip.label

  ## Because A & B are singletons, they are not counted.
  ## Then CD are the same so the probability of a shared trait is 1
  out = trait_heritage_specific(tree, trait, generation_time = 0.2, condition = "a")
  expect_equal(
    out$by_trait,
    structure(
      list(
        generation = c(0.0, 0.2,0.4, 0.6, 0.8, 1.0),
        state = c("a", "a", "a","a","a", "a"),
        numerator_sum = c(0.0, 0.0, 1.0, 1.0, 1.0, 1.0),
        denominator_sum = c(0.0, 0.0, 1.0, 1.0, 3.0, 10),
        clade_probability = c(0.0, 0.0, 1.0, 1.0, 1/3, 0.1)
      ),
      sorted = c("generation", "state"),
      row.names = c(NA, -6L),
      class = c("data.table", "data.frame")
    ))
})


test_that("#6 Pipeline NA test", {
  t = ape::read.tree(text = "(tA,(tB,(tC,tD)));")
  t = ape::compute.brlen(t)

  trait = c("b", "a", NA, "a")
  names(trait) = t$tip.label

  expect_error(trait_heritage(t, trait, generation_time = 0.2))

})


test_that("#6 Matching names test", {
  t = ape::read.tree(text = "(tA,(tB,(tC,tD)));")
  t = ape::compute.brlen(t)

  trait = c("b", "a", "a", "a")
  names(trait) = c("tA", "tB", "AnotherTaxa", "tD")

  expect_error(trait_heritage(t, trait, generation_time = 0.2))
})

test_that("#6 No matches", {
  t = ape::read.tree(text = "(tA,(tB,(tC,tD)));")
  t = ape::compute.brlen(t)

  trait = c("b", "a", "c", "d")
  names(trait) = c("tA", "tB", "tC", "tD")

  result = trait_heritage(t, trait, generation_time = 0.2)

  out_probability = result$summary

  expect_equal(out_probability$clade_probability, c(NaN, 0, 0, 0, 0))

})


test_that("#1. Complex test", {

  t_file = test_path("testdata", "20240617_CBES_bModelTest_adj-papua_n.mcct.trees")
  t_file = test_path("testdata", "STATE_140990000.nex")
  md_file = test_path("testdata", "language_metadata.csv")

  t = ape::read.nexus(t_file)
  md = read.csv(md_file)

  t = ape::keep.tip(t, md$PhyID)

  generation_time = 25

  trait = md$TNG_Branch
  names(trait) = md$PhyID

  result = trait_heritage(
    tree = t,
    trait = trait,
    generation_time = generation_time
  )

  expect_equal(round(result$summary$clade_probability[result$summary$generation == 62975], 8),  0.2292672)
  expect_equal(result$by_trait$numerator_sum[result$by_trait$generation == 62975], c(0, 78, 28, 190, 15, 741, 1, 990, 0))
  expect_equal(result$by_trait$denominator_sum[result$by_trait$generation == 62975], rep(8911, 9))
})


