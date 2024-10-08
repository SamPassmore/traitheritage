## trait co-evolution

test_that("Simple test", {
  t = ape::read.tree(text = "((A,B),(C,D));")
  tree = ape::compute.brlen(t)

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
  generation_time = 0.5

  result = trait_coevolution(tree, trait, distance_matrix, generation_time, cut_off)

  expect_equal(
    result$time_df$denominator_sum,
    c(0, 2))
  expect_equal(
    result$time_df$p_lang_dist,
    c(0, 0.5))
  expect_equal(
    result$pairs_df$dist_trait,
    c(T, T, F, F, F, F, T, T))
})

test_that("Complex test", {
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

  trait = metadata_dist$TNG
  names(trait) = metadata_dist$PhyID

  cut_off = 50
  generation_time = 25

  ce = trait_coevolution(
    t,
    trait,
    distance_matrix = distance_matrix,
    generation_time = generation_time,
    cut_off = cut_off
  )
})

