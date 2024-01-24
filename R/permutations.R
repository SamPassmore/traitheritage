#' Title
#'
#' @param tree a single phylogeny
#' @param trait a vector of traits, named for each taxa
#' @param n_generations the number of times to cut a tree and calculate the probability of trait heritage
#' @param n_permutations the number of permutations to run as a comparison (default = 100)
#'
#' @return a list containing the estimated probability, and the permuted probabilities.
#' @export
#'
permute_trait_heritage = function(tree, trait, n_generations, n_permutations = 100){
  # True data results
  true = trait_heritage(tree, trait, n_generations)

  # Permuted results
  permuted_list = lapply(1:n_permutations, function(i){
    p_trait = trait
    names(p_trait) = names(p_trait)[sample(1:length(p_trait))]
    trait_heritage(tree, p_trait, n_generations)
  })
  names(permuted_list) = paste0("p_", 1:n_permutations)

  permuted = do.call(rbind, permuted_list)
  permuted$iteration = lapply(strsplit(row.names(permuted), "\\."), '[[', 1)
  rownames(permuted) = NULL

  list(true = true, permuted = permuted)
}
