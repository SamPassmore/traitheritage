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

  # Permuted results
  permuted_list = lapply(1:n_permutations, function(i){
    p_trait = trait
    names(p_trait) = names(p_trait)[sample(1:length(p_trait))]
    trait_heritage(tree, p_trait, generation_time)
  })

  by_trait = lapply(permuted_list, "[[", 1)
  p_summary = lapply(permuted_list, "[[", 2)

  names(by_trait) = paste0("p_", 1:n_permutations)
  names(p_summary) = paste0("p_", 1:n_permutations)

  permuted_bytrait = data.table::rbindlist(by_trait, idcol = "iteration")
  permuted_summary = data.table::rbindlist(p_summary, idcol = "iteration")

  return(list(by_trait = permuted_bytrait, summary = permuted_summary))
}
