library(profvis)
library(ape)
library(microbenchmark)

devtools::load_all()

tree = read.tree("../daru_phylogeny/processed_data/mtdna_ultra.nex")
metadata = read.csv('../daru_phylogeny/processed_data/daru_metadata.csv')

metadata = metadata[metadata$tips.x %in% tree$tip.label,]

# cuts for 25 year generations times
n_generations = max(ape::node.depth.edgelength(tree)) / 25

## Language spoken
trait = metadata$Language
names(trait) = metadata$tips.x

all(tree$tip.label %in% names(language))
length(tree$tip.label) == length(language)


profvis({
  ## Add argument tests to the function
  # 1. Calculate tree cuts
  clades = .slice_tree(tree, n_generations)

  # Add trait to splits
  clades = dplyr::full_join(clades, data.frame(taxa = names(trait), trait = trait), by = "taxa")

  # 2. For each generation, calculate the probability of a shared trait within each clade
  generations = unique(clades$generation)

  # Use lapply for generations and sapply for clade_sets
  generation_df = as.data.table(clades)
  output = lapply(generations, function(g) {
    generation_df = generation_df[generation == g]
    clade_sets <- unique(generation_df$clade)

    clade_result <- sapply(clade_sets, function(cs) {
      trait_states <- generation_df$trait[generation_df$clade == cs]
      names(trait_states) <- generation_df$taxa[generation_df$clade == cs]

      .clade_probabilities(clade_states = trait_states)
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
})


## Subsetting test
g = "g_7725"
generation_df = as.data.table(clades)

microbenchmark(
  clades[clades$generation == g, ],
  {
    generation_df[generation == g]
  },
  times = 10
)

## Lapply test

## Data-table is faster
f1_datatable = function(){
  generation_df = as.data.table(clades)
  lapply(generations, function(g) {
    generation_df = generation_df[generation == g]
    clade_sets <- unique(generation_df$clade)

    clade_result <- sapply(clade_sets, function(cs) {
      trait_states <- generation_df$trait[generation_df$clade == cs]
      names(trait_states) <- generation_df$taxa[generation_df$clade == cs]

      .clade_probabilities(clade_states = trait_states)
    })

    # Combine clade result with generation and clade information
    x = data.frame(
      generation = g,
      clade = clade_sets,
      numerator = unlist(clade_result[1, ]),
      denominator = unlist(clade_result[2, ]),
      stringsAsFactors = FALSE
    )
  })
}

f2_baseR = function(){
  lapply(generations, function(g) {
    generation_df <- clades[clades$generation == g, ]
    clade_sets <- unique(generation_df$clade)

    clade_result <- sapply(clade_sets, function(cs) {
      trait_states <- generation_df$trait[generation_df$clade == cs]
      names(trait_states) <- generation_df$taxa[generation_df$clade == cs]

      .clade_probabilities(clade_states = trait_states)
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
}

microbenchmark(f1_datatable(), f2_baseR(), times = 10)

## Profiz datatable

profvis({
  generation_df = as.data.table(clades)
  lapply(generations, function(g) {
    generation_df = generation_df[generation == g]
    clade_sets <- unique(generation_df$clade)

    clade_result <- sapply(clade_sets, function(cs) {
      trait_states <- generation_df$trait[generation_df$clade == cs]
      names(trait_states) <- generation_df$taxa[generation_df$clade == cs]

      .clade_probabilities(clade_states = trait_states)
    })

    # Combine clade result with generation and clade information
    x = data.frame(
      generation = g,
      clade = clade_sets,
      numerator = unlist(clade_result[1, ]),
      denominator = unlist(clade_result[2, ]),
      stringsAsFactors = FALSE
    )
  })
}, rerun = TRUE)

## Clade probabilities is slow

profvis({.clade_probabilities(trait_states)}, rerun = TRUE)