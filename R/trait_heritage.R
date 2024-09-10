#' Internal function: Calculate probabilities of matching states within a clade
#'
#' @param clade_states A named vector of trait values
#' @param state the state to condition on
#'
#' @return numerator: the number of matching pairs within a clade
#' @return denominator: The number of taxa pairs within a clade
#' @export
#'
.clade_probabilities <- function(state) {

  N <- choose(table(state), 2)
  D <- choose(length(state), 2)

  return(data.table(state = names(N), numerator = c(N), denominator = D))
}

#' Internal function: Cuts a phylogeny at a given point
#'
#' @param tree is a dated phylogenetic tree with branch lengths stored as a phylo object (as in the ape package).
#' @param cut the slice time
#' @param k number of slices
#'
#' @description
#' This function has been adapted from phyloregion, but included separately here to avoid unnecessary dependencies from that package
#' Please cite Daru, B. H., Karunarathne P., & Schliep K. (2020), phyloregion: R package for biogeographic regionalization and macroecology. Methods in Ecology and Evolution, 11: 1483-1491. doi:10.1111/2041-210X.13478
#' if using this function.
#'
#' @return
#' @export
#'
#' @examples
.get_clades = function(tree, cut = NULL, k = NULL){
  nh <- ape::node.depth.edgelength(tree)
  nh <- max(nh) - nh
  # if (!is.null(k)) {
  #   if (k >= Ntip(tree))
  #     return(as.list(tree$tip.label))
  #   if (k == 1)
  #     return(list(tree$tip.label))
  #   kids <- lengths(phangorn::Descendants(tree, type = "children"))
  #   kids[kids > 0] <- kids[kids > 0] - 1L
  #   tmp <- 1
  #   eps <- 1e-08
  #   ordered_nh <- order(nh, decreasing = TRUE)
  #   i <- 1
  #   while (tmp < k) {
  #     j <- ordered_nh[i]
  #     cut <- nh[j] - eps
  #     tmp <- tmp + kids[j]
  #     i <- i + 1
  #   }
  # }
  ind <- which((nh[tree$edge[, 1]] > cut) &
                 (nh[tree$edge[, 2]] <= cut))
  desc <- phangorn::Descendants(tree)
  res <- desc[tree$edge[ind, 2]]
  lapply(res, function(res, tips) tips[res], tree$tip.label)
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
                    .clades = .get_clades(tree, cut = cc)
                    names(.clades) = seq_along(.clades)
                    ## Stack clades
                    clade_df = stack(.clades)
                    colnames(clade_df) = c("taxa", "clade")

                    clade_df
                  })
  names(clades) = paste0("g_", cuts)

  clades_df = data.table::rbindlist(clades, idcol = "generation")
  clades_df$clade = as.numeric(clades_df$clade) # it is useful later for this to be numeric, but this is just a choice.

  # Return
 return(data.table::setDT(clades_df))
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

  if(is.null(names(trait))) stop("trait must have names that match the taxa. Ensure trait is a named vector.")

  # Trait names must match taxa labels
  if(!all(names(trait) %in% tree$tip.label)) stop("Some tips have no matching trait. Make sure all tips have a trait.")

  # 1. Calculate tree cuts
  clades = .slice_tree(tree, generation_time)

  # Add trait to splits
  clades = dplyr::full_join(clades, data.frame(taxa = names(trait), trait = trait), by = "taxa")

  # Use DT for fast calculations.
  data.table::setDT(clades)

  output = clades[, .clade_probabilities(trait), by = c("generation", "clade")]

  # Return
  output
}

