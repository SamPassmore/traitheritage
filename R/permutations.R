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
permute_trait_heritage = function(tree, trait, generation_time, state = NULL, n_permutations = 100){
  # True data results
  true = trait_heritage(tree, trait, generation_time, state = state)

  # Permuted results
  permuted_list = lapply(1:n_permutations, function(i){
    p_trait = trait
    names(p_trait) = names(p_trait)[sample(1:length(p_trait))]
    trait_heritage(tree, p_trait, generation_time, state = state)
  })
  names(permuted_list) = paste0("p_", 1:n_permutations)

  permuted = do.call(rbind, permuted_list)
  permuted$iteration = lapply(strsplit(row.names(permuted), "\\."), '[[', 1)
  rownames(permuted) = NULL

  list(true = true, permuted = permuted)
}
