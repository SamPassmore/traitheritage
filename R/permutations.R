#' Calculate the true probabilities of shared traits and the permuted probability
#'
#' @param tree a single phylogeny
#' @param trait a vector of traits, named for each taxa
#' @param generation_time the time between generations to cut the phylogeny at
#' @param n_permutations the number of permutations to run as a comparison (default = 100)
#'
#' @return a list containing the estimated probability, and the permuted probabilities.
#' @export
#'
permute_trait_heritage = function(tree, trait, generation_time, n_permutations = 100){
  # True data results
  true = trait_heritage(tree, trait, generation_time)

  # Permuted results
  permuted_list = lapply(1:n_permutations, function(i){
    p_trait = trait
    names(p_trait) = names(p_trait)[sample(1:length(p_trait))]
    trait_heritage(tree, p_trait, generation_time)
  })
  names(permuted_list) = paste0("p_", 1:n_permutations)

  permuted = do.call(rbind, permuted_list)
  permuted$iteration = lapply(strsplit(row.names(permuted), "\\."), '[[', 1)
  rownames(permuted) = NULL

  list(true = true, permuted = permuted)
}

p.clade_probabilities <- function(TS) {
  N <- sum(choose(table(TS), 2))
  D <- choose(length(TS), 2)
  return(list(numerator = N, denominator = D))
}

p.trait_heritage = function(tree, trait, generation_time){
  ## Add argument tests to the function
  if(any(is.na(trait))) stop("No NA trait values are allowed. ")

  # Trait names must match taxa labels
  if(!all(names(trait) %in% tree$tip.label)) stop("Some tips have no matching trait. Make sure all tips have a trait.")

  # 1. Calculate tree cuts
  clades = .slice_tree(tree, generation_time)

  # Add trait to splits
  clades = dplyr::full_join(clades, data.frame(taxa = names(trait), trait = trait), by = "taxa")

  # Use DT for fast calculations.
  data.table::setDT(clades)

  output = clades[, p.clade_probabilities(trait), by = c("generation", "clade")]

  # Return
  output
}


