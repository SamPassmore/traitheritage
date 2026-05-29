#' Calculate the true probabilities of shared traits and the permuted probability
#'
#' @param tree a single phylogeny
#' @param trait a vector of traits, named for each taxa
#' @param generation_time the time between generations to cut the phylogeny at
#' @param n_permutations the number of permutations to run as a comparison (default = 100)
#'
#' @return a list containing the estimated probability, and the permuted probabilities.
#' @export
#' @examples
#' tree <- ape::read.tree(text = "(tA,(tB,(tC,tD)));")
#' tree <- ape::compute.brlen(tree)
#' trait <- c("b", "a", "a", "a")
#' names(trait) <- tree$tip.label
#' permute_trait_heritage(tree, trait, generation_time = 0.2, n_permutations = 5)
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

#' Calculate the true and permuted probabilities for a specific trait condition
#'
#' @param tree a single phylogeny
#' @param trait a vector of traits, named for each taxa
#' @param generation_time the time between generations to cut the phylogeny at
#' @param n_permutations the number of permutations to run as a comparison (default = 100)
#' @param condition the trait state to condition on
#' @param seed random seed for reproducibility (default = 937)
#'
#' @return a list containing the estimated probability and the permuted probabilities.
#' @export
#' @examples
#' tree <- ape::read.tree(text = "((tA:0.2,tX:0.2):0.8,(tB:0.65,(tC:0.3,tD:0.3):0.35):0.35);")
#' trait <- c("b", "b", "b", "a", "a")
#' names(trait) <- tree$tip.label
#' permute_trait_heritage_specific(tree, trait, generation_time = 0.2,
#'   n_permutations = 5, condition = "a")
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
