#' Internal function: Calculate probabilities of matching states within a clade
#'
#' @param clade_states A named vector of trait values
#'
#' @return numerator: the number of matching pairs within a clade
#' @return denominator: The number of taxa pairs within a clade
#' @export
#'
.clade_probabilities <- function(TS) {
  N <- sum(choose(table(TS), 2))
  D <- choose(length(TS), 2)
  return(list(numerator = N, denominator = D))
}

#' Internal Function: Identify clades at all generations
#'
#' @param tree a phylogenetic tree
#' @param generation_time A number showing the number of generations to calculate for the tree. User should calculate this based on their branch lengths
#'
#' @return clades: a data.frame showing the clades for each taxa within each generation.
#' @export
.slice_tree = function(tree, generation_time){
  max_tree_depth = max(ape::node.depth.edgelength(tree)[1:ape::Ntip(tree)]) # allows for non-ultrametric trees
  cuts = seq(from = 0, to = max_tree_depth, by = generation_time)

  # don't calculate any cuts for values the same as the depth of the tree
  cuts = cuts[cuts != max_tree_depth]

  # Identify taxa within each clade for each generation
  clades = lapply(cuts,
                  function(cc) {
                    ## get.clades assumes that taxa start at zero, and cuts go backwards from there.
                    .clades = phyloregion::get_clades(tree, cut = cc)
                    names(.clades) = seq_along(.clades)
                    ## Stack clades
                    clade_df = stack(.clades)
                    colnames(clade_df) = c("taxa", "clade")

                    clade_df
                  })
  names(clades) = paste0("g_", round(cuts, 2))

  clades_df = purrr::map_df(clades, ~as.data.frame(.x), .id="generation")
  clades_df$clade = as.numeric(clades_df$clade) # it is useful later for this to be numeric, but this is just a choice.

  # Return
  clades_df
}


## Main function
#' Calculate the heritage of a trait along a single phylogeny
#'
#' @param tree a single phylogeny
#' @param trait a vector of traits, named for each taxa
#' @param generation_time the number of times to cut a tree and calculate the probability of trait heritage
#'
#' @return a dataframe containing the probability of shared languages within each generation
#' @export
#'
trait_heritage = function(tree, trait, generation_time){
  ## Add argument tests to the function
  if(any(is.na(trait))) stop("No NA trait values are allowed. ")

  # Trait names must match taxa labels
  if(!all(names(trait) %in% tree$tip.label)) stop("Some tips have no matching trait. Make sure all tips have a trait.")

  # 1. Calculate tree cuts
  clades = .slice_tree(tree, generation_time)

  # Add trait to splits
  clades = dplyr::full_join(clades, data.frame(taxa = names(trait), trait = trait), by = "taxa")

  # Use lapply for generations and sapply for clade_sets
  data.table::setDT(clades)

  output <- clades[, .clade_probabilities(trait), by = c("generation", "clade")]
  output = output[, .(numerator_sum = sum(numerator), denominator_sum = sum(denominator)), by = "generation"]
  output[, clade_probability := numerator_sum / denominator_sum]

  # Return
  output
}

