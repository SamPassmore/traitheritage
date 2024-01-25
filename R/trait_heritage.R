#' Internal function: Calculate probabilities of matching states within a clade
#'
#' @param clade_states A named vector of trait values
#'
#' @return numerator: the number of matching pairs within a clade
#' @return denominator: The number of taxa pairs within a clade
#' @export
#'
.clade_probabilities = function(clade_states) {

  if(length(clade_states) > 1) {
    x = as.numeric(factor(clade_states))
    dd = 1 - c(dist(x, method = "manhattan"))
    numerator = sum(dd)
    denominator = length(dd)
  } else {
    numerator = 0 # if there is only 1 individual in a clade, then the probability is 0
    denominator = 0
  }
  list(numerator = numerator, denominator = denominator)
}

#' Internal Function: Identify clades at all generations
#'
#' @param tree a phylogenetic tree
#' @param n_generations A number showing the number of generations to calculate for the tree. User should calculate this based on their branch lengths
#'
#' @return clades: a data.frame showing the clades for each taxa within each generation.
#' @export
#'
.slice_tree = function(tree, n_generations){
  # Calculate tree depth
  max_tree_depth = max(ape::node.depth.edgelength(tree)[1:ape::Ntip(tree)]) # allows for non-ultrametric trees
  root_depth = 0 # assumes that taxa start at 0. We make the first cut just above zoer

  bin_size = max_tree_depth / n_generations

  cuts = c(cumsum(
    rep(bin_size, n_generations)
  ))

  # Identify taxa within each clade for each generation
  clades = lapply(cuts[-length(cuts)],
                  function(c) {
                    ## get.clades assumes that taxa start at zero, and cuts go backwards from there.
                    .clades = phyloregion::get_clades(tree, cut = c)
                    names(.clades) = seq_along(.clades)
                    ## Stack clades
                    clade_df = stack(.clades)
                    colnames(clade_df) = c("taxa", "clade")

                    clade_df
                  })
  names(clades) = paste0("g_", round(cuts[-length(cuts)], 2))

  clades_df = purrr::map_df(clades, ~as.data.frame(.x), .id="generation")
  clades_df$clade = as.numeric(clades_df$clade) # it is useful later for this to be numeric

  # Return
  clades_df
}


## Main function
#' Calculate the heritage of a trait along a single phylogeny
#'
#' @param tree a single phylogeny
#' @param trait a vector of traits, named for each taxa
#' @param n_generations the number of times to cut a tree and calculate the probability of trait heritage
#'
#' @return a dataframe containing the probability of shared languages within each generation
#' @export
#'
trait_heritage = function(tree, trait, n_generations){
  ## Add argument tests to the function
  if(any(is.na(trait))) stop("No NA trait values are allowed. ")

  # 1. Calculate tree cuts
  clades = .slice_tree(tree, n_generations)

  # Add trait to splits
  clades = dplyr::full_join(clades, data.frame(taxa = names(trait), trait = trait), by = "taxa")

  # 2. For each generation, calculate the probability of a shared trait within each clade
  generations = unique(clades$generation)

  # Use lapply for generations and sapply for clade_sets
  generation_df = as.data.table(clades)
  output = lapply(generations, function(g) {
    generation_df = generation_df[generation_df$generation == g]
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
}

