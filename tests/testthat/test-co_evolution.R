## trait co-evolution

test_that("Simple test: only langauge matches", {
  # tree = ape::read.tree(text = "(tA,(tB,(tC,tD)));")
  # tree = ape::compute.brlen(tree)

  tree = ape::read.tree(text =
                          "((tA:1.16,tB:1.16):0.44,((tC:0.20,tD:0.20):0.18,tE:0.39):1.22);")
  tree = ape::compute.brlen(tree) # makes the tree depth = 1

  trait = c("b", "b", "a", "a", "a")
  names(trait) = tree$tip.label


  distances = matrix(
    c(1, 2,
      1, 2,
      2, 1,
      5, 6,
      6, 5),
    byrow = TRUE,
    ncol = 2,
    dimnames = list(tree$tip.label, c("X", "Y"))
  )
  distance_matrix = as.matrix(dist(distances))

  cut_off = -1 # no distances are negative so this makes no distance pairs
  generation_time = 0.2

  result = trait_coevolution(tree, trait, distance_matrix, generation_time, cut_off)

  expect_equal(
    result$time_df$denominator_sum,
    c(0, 0, 2, 4, 4, 10, 10))
  expect_equal(
    result$time_df$nolang_dist,
    c(0, 0, 0, 0, 0, 0, 0)) # there should be no distance matches (cut_off = -1)
  expect_equal(
    rowSums(result$time_df[,c("lang_dist", "lang_nodist", "nolang_dist", "nolang_nodist")]),
    result$time_df$denominator_sum)
})


test_that("Simple test: specific language matches", {
  tree = ape::read.tree(text = "(((t1:0.15,t2:0.15):0.35,(t3:0.25, t4:0.25):0.25):0.5,(t5:0.35, t6:0.35):0.65);")

  trait = c("a", "a", "b", "b", "a", "a")
  names(trait) = tree$tip.label


  distances = matrix(
    c(1, 2,
      1, 2,
      2, 1,
      5, 6,
      6, 5,
      8, 8),
    byrow = TRUE,
    ncol = 2,
    dimnames = list(tree$tip.label, c("X", "Y"))
  )
  distance_matrix = as.matrix(dist(distances))

  cut_off = -1 # no distances are negative so this makes no distance pairs
  generation_time = 0.2

  result = trait_coevolution_specific(tree, trait, distance_matrix, generation_time, cut_off, condition = "a")

  expect_equal(
    result$time_df$denominator_sum,
    c(0, 1.0, 2.0, 7.0, 7.0, 15.0, 15.0)
    )
  expect_equal(
    result$time_df$nolang_dist,
    c(0, 0, 0, 0, 0, 0, 0)) # there should be no distance matches (cut_off = -1)

  expect_equal(
    result$time_df$lang_nodist,
    c(0, 1, 2, 2, 2, 6, 6)
  )

})

test_that("Same node heigts", {
  tree = ape::read.tree(text = "(((t1:0.15,t2:0.15):0.35,(t3:0.25, t4:0.25):0.25):0.5,(t5:0.15, t6:0.15):0.85);")

  trait = c("a", "a", "b", "b", "a", "a")
  names(trait) = tree$tip.label


  distances = matrix(
    c(1, 2,
      1, 2,
      2, 1,
      5, 6,
      6, 5,
      8, 8),
    byrow = TRUE,
    ncol = 2,
    dimnames = list(tree$tip.label, c("X", "Y"))
  )
  distance_matrix = as.matrix(dist(distances))

  cut_off = -1 # no distances are negative so this makes no distance pairs
  generation_time = 0.2

  # I know this will give a warning
  result = suppressWarnings(
    trait_coevolution_specific(tree, trait, distance_matrix, generation_time, cut_off, condition = "a")
  )


  expect_equal(
    result$time_df$denominator_sum,
    c(0, 2.0, 2.0, 7.0, 7.0, 15.0, 15.0)
  )
  expect_equal(
    result$time_df$nolang_dist,
    c(0, 0, 0, 0, 0, 0, 0)) # there should be no distance matches (cut_off = -1)

  expect_equal(
    result$time_df$lang_nodist,
    c(0, 2, 2, 2, 2, 6, 6)
  )

})



test_that("Complex test", {
  t_file = test_path("testdata", "20240617_CBES_bModelTest_adj-papua_n.mcct.trees")
  t_file = test_path("testdata", "STATE_140990000.nex")
  md_file = test_path("testdata", "language_metadata.csv")

  t = ape::read.nexus(t_file)
  md = read.csv(md_file)

  metadata_dist = md %>%
    dplyr::filter(!is.na(latitude_hybrid) & !is.na(longitude_hybrid) &
                    !is.na(TNG))

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

  expect_equal(
    tail(ce$time_df$denominator_sum, 1),
    choose(ape::Ntip(t), 2)
  )

  expect_equal(
    rowSums(ce$time_df[,c("lang_dist", "lang_nodist", "nolang_dist", "nolang_nodist")]),
    ce$time_df$denominator_sum)
})

