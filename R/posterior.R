### Calculate the posterior probability of a tree

#' Calculate the probability of trait heritage across a posterior sample of trees
#'
#' @param trees A posterior sample of trees, of class multiPhylo
#' @param trait A named vector of trait values. Names should match taxa
#' @param n_generations A numeric value indicating how many equal cuts to make in the tree
#'
#' @return a Dataframe showing the probability of a shared trait for each generation, witihn each tree
#' @export
#'

posterior_trait_heritage = function(trees, trait, n_generations){

  if(!"multiPhylo" %in% class(trees)) stop("Trees should be of class multiPhylo. If you only have one tree, use the trait_heritage function.")

  output = lapply(trees, function(t){
    trait_heritage(tree = t, trait = trait, n_generations = n_generations)
  })
  names(output) = 1:length(trees)

  dplyr::bind_rows(output, .id = "posterior")

}
