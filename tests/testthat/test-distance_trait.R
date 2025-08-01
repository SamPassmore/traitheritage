# Test Distance

test_that("#1. Simple distance test", {

  t = ape::read.tree(text = "(tA,(tB,(tC,tD)));")
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
  generation_time = 0.2

  result = distance_trait_heritage(
    t,
    distance_matrix = distance_matrix,
    generation_time = generation_time,
    cut_off = cut_off
  )

  expect_equal(result$numerator_sum, c(0, 1, 1, 1, 2))
  expect_equal(result$denominator_sum, c(0, 1, 1, 3, 6))

})

test_that("#2. Simple distance test", {

  t = ape::read.tree(text = "(tA,(tB,(tC,tD)));")
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

  expect_equal(result$clade_probability, c(NaN, 1.000, 1/3, 1/3))
  expect_equal(result$numerator_sum, c(0, 1, 1, 2))
  expect_equal(result$denominator_sum, c(0, 1, 3, 6))
})

test_that("Same Results Test", {

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

  result1 = distance_trait_heritage(
    t,
    distance_matrix = distance_matrix,
    generation_time = generation_time,
    cut_off = cut_off
  )

  result2 = distance_trait_heritage(
    t,
    distance_matrix = distance_matrix,
    generation_time = generation_time,
    cut_off = cut_off
  )

  expect_equal(result1, result2)
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

test_that("Same Denominator", {

  t = ape::read.tree(text = "(tA,(tB,(tC,tD)));")
  t = ape::compute.brlen(t)

  trait = c("b", "a", "a", "a")
  names(trait) = t$tip.label

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
  generation_time = 0.1

  d_result = distance_trait_heritage(
    t,
    distance_matrix = distance_matrix,
    generation_time = generation_time,
    cut_off = cut_off
  )

  t_result = trait_heritage(t, trait, generation_time)

  expect_equal(d_result$denominator_sum, t_result$summary$denominator_sum)
})


test_that("#1. Complex distance test", {

  t_file = test_path("testdata", "20240617_CBES_bModelTest_adj-papua_n.mcct.trees")
  t_file = test_path("testdata", "STATE_140990000.nex")
  md_file = test_path("testdata", "language_metadata.csv")

  t = ape::read.nexus(t_file)
  md = read.csv(md_file)

  metadata_dist = md %>%
    dplyr::filter(!is.na(latitude_hybrid) & !is.na(longitude_hybrid))

  t = ape::keep.tip(t, metadata_dist$PhyID)

  ## Calculate distance on a sphere & convert units to Kms
  distance_matrix = geosphere::distm(metadata_dist[, c("longitude_hybrid", "latitude_hybrid")],
                                     fun = geosphere::distHaversine) / 1000
  dimnames(distance_matrix) = list(metadata_dist$PhyID, metadata_dist$PhyID)

  cut_off = 50
  generation_time = 25

  result = distance_trait_heritage(
    t,
    distance_matrix = distance_matrix,
    generation_time = generation_time,
    cut_off = cut_off
  )

  expect_equal(result$clade_probability[result$generation == 62975], 0.11911011)
  expect_equal(result$numerator_sum[result$generation == 62975], 953)
  expect_equal(result$denominator_sum[result$generation == 62975], 8001)
})
