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
permute_trait_heritage = function(tree, trait, generation_time, n_permutations = 100, seed = 937){
  set.seed(seed)
  # Permuted results
  permuted_list = lapply(seq_len(n_permutations), function(i){
    p_trait = trait
    names(p_trait) = names(p_trait)[sample(seq_along(p_trait))]
    trait_heritage(tree, p_trait, generation_time)
  })

  by_trait = lapply(permuted_list, "[[", 1)
  p_summary = lapply(permuted_list, "[[", 2)

  names(by_trait) = paste0("p_", seq_len(n_permutations))
  names(p_summary) = paste0("p_", seq_len(n_permutations))

  permuted_bytrait = data.table::rbindlist(by_trait, idcol = "iteration")
  permuted_summary = data.table::rbindlist(p_summary, idcol = "iteration")

  return(list(by_trait = permuted_bytrait, summary = permuted_summary))
}

permute_trait_heritage_specific = function(tree, trait, generation_time, n_permutations = 100, condition = NULL, seed = 937){
  set.seed(seed)
  # Permuted results
  permuted_list = lapply(seq_len(n_permutations), function(i){
    p_trait = trait
    names(p_trait) = names(p_trait)[sample(seq_along(p_trait))]
    trait_heritage_specific(tree, p_trait, generation_time, condition = condition)
  })

  by_trait = lapply(permuted_list, "[[", 1)
  p_summary = lapply(permuted_list, "[[", 2)

  names(by_trait) = paste0("p_", seq_len(n_permutations))
  names(p_summary) = paste0("p_", seq_len(n_permutations))

  permuted_bytrait = data.table::rbindlist(by_trait, idcol = "iteration")
  permuted_summary = data.table::rbindlist(p_summary, idcol = "iteration")

  return(list(by_trait = permuted_bytrait, summary = permuted_summary))
}
