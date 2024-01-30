## Distance function

library(ape)

t = ape::read.tree(text = "((A,B),(C,D));")
t = ape::compute.brlen(t)

tree = t

trees = list(t, t, t, t)
class(trees) = "multiPhylo"

trait = c("a", "a", "b", "b")
names(trait) = c("A", "B", "C", "D")

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
n_generations = 2

distance_trait_heritage = function(tree, trait, n_generations, distance_matrix, cut_off){

  # 1. Calculate tree cuts
  clades = .slice_tree(tree, n_generations)

  # 2. Calculate trait distance
  ## This solution is taken from https://stackoverflow.com/a/76731093/1544746
  ## Which notes that it will not work for R < 4.0.0
  n = nrow(distance_matrix)
  s = seq_len(n) - 1L
  nms = dimnames(distance_matrix)
  dm = data.frame(value = distance_matrix[sequence(s, seq.int(1L, length(distance_matrix), n))],
                    row = gl(n, 1L, labels = nms[[1L]])[sequence(s, 1L)],
                    col = rep.int(gl(n, 1L, labels = nms[[2L]]), s))

  dm$trait = dm$value < cut_off

  # 3. For each generation, calculate the probability of a shared trait within each clade
  generations = unique(clades$generation)

  # Use lapply for generations and sapply for clade_sets
  generation_df = as.data.table(clades)
  # set key on data.table for speeding up subsetting
  setkey(generation_df, generation)

  output = lapply(generations, function(g) {
    generation_df = generation_df[generation_df$generation == g]
    clade_sets <- unique(generation_df$clade)

    clade_result <- sapply(clade_sets, function(cs) {
      clade_taxa = generation_df$taxa[generation_df$clade == cs]
      trait = dm$trait[dm$row %in% clade_taxa & dm$col %in% clade_taxa]

      if (length(trait) >= 1) {
        numerator = sum(trait)
        denominator = length(trait)
      } else {
        numerator = 0
        denominator = 0
      }
      list(numerator = numerator, denominator = denominator)
    })

    # Combine clade result with generation and clade information
    data.frame(
      generation = g,
      clade = clade_sets,
      numerator = unlist(clade_result[1, ]),
      denominator = unlist(clade_result[2, ]),
      stringsAsFactors = FALSE
    )
  })
  output = do.call(rbind, output)

  # 3. Calculate the probability of a shared trait for each generation
  output_dt = data.table::data.table(output)
  output_dt = output_dt[, list(numerator_sum = sum(numerator),
                               denominator_sum = sum(denominator)), by = generation]
  output_dt[, clade_probability := numerator_sum / denominator_sum]

  # Return
  data.frame(output_dt)
}
