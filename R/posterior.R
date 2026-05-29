#' Calculate the probability of trait heritage across a posterior sample of trees
#'
#' @param trees A posterior sample of trees, of class multiPhylo
#' @param trait A named vector of trait values. Names should match taxa
#' @param generation_time A numeric value indicating how many equal cuts to
#'   make in the tree
#'
#' @return a data frame showing the probability of a shared trait for each
#'   generation, within each tree
#' @export
#' @examples
#' tree <- ape::read.tree(text = "(tA,(tB,(tC,tD)));")
#' tree <- ape::compute.brlen(tree)
#' trees <- c(tree, tree, tree)
#' trait <- c("b", "a", "a", "a")
#' names(trait) <- tree$tip.label
#' posterior_trait_heritage(trees, trait, generation_time = 0.2)
posterior_trait_heritage <- function(trees, trait, generation_time) {

  if (!"multiPhylo" %in% class(trees)) {
    stop(
      "Trees should be of class multiPhylo. ",
      "If you only have one tree, use the trait_heritage function."
    )
  }

  output <- lapply(trees, function(t) {
    trait_heritage(tree = t, trait = trait, generation_time = generation_time)$summary
  })
  names(output) <- seq_along(trees)

  dplyr::bind_rows(output, .id = "posterior")
}
