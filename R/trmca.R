#' Calculate Time since Most Recent Common Ancestor for one tree
#'
#'@description
#'This function is for internal use only.
#'
#'
#' @param tree a single phylogeny of class phylo
#' @param states a named vector of states to calculate clades from. Names must match tips.
#'
#' @return a data.frame containing information on the identified clades, the state they contain, and the tree they are from.
#' @export
#'
.tmrca_onetree = function(tree, states){
  # Find all clades
  all.sub.trees = caper::clade.members.list(phy = tree, tip.labels = TRUE)

  # Find all clades where all trait states are common
  common_clades_idx = sapply(all.sub.trees, function(x) length(unique(states[x])) == 1)
  common_clades = all.sub.trees[common_clades_idx]

  ## Remove all clades that are subsets of each other
  # First create comparisons between all clades
  comparisons = expand.grid(names(common_clades), names(common_clades))
  # remove all self-comparisons
  comparisons = comparisons[comparisons[,1] != comparisons[,2],]
  # Then identify which clades are contained within others, and one is longer than the other
  clades_idx = apply(comparisons, 1, function(v)
    all(common_clades[[v[1]]] %in% common_clades[[v[2]]])
  )
  subclades = comparisons[clades_idx,]
  subclades = subclades[!subclades[,2] %in% subclades[,1],]

  clades_ofinterest = names(common_clades)[!names(common_clades) %in% as.character(subclades[,1])]

  # Find the height of all nodes
  node_height = phytools::nodeHeights(tree)
  # There is a possibility of switching this function to an ape function to reduce dependencies
  # node_height = ape::node.depth.edgelength(tree)
  tree_height = max(node_height)
  nh = data.frame(Node = as.character(c(tree$edge[, 1], tree$edge[, 2])),
                  Date = c(
                    round(tree_height - node_height[, 1], 3),
                    round(tree_height - node_height[, 2], 3)
                  ))
  nh = nh[!duplicated(nh$Node),]

  dates = dplyr::left_join(data.frame(Node = clades_ofinterest), nh, by = "Node")
  dates$n_tips = sapply(common_clades[as.character(dates$Node)], length)

  # Tips must be sorted, because they may arise in different orders, but the same clade, in different trees
  dates$tips = sapply(common_clades[as.character(dates$Node)], function(x){
    x = sort(x, na.last = TRUE)
    paste(x, collapse = "; ")
  })

  dates$state = sapply(common_clades[as.character(dates$Node)], function(x) unique(states[x]))

  dates
}

#' Calculate the Time since Most Recent Common Ancestor across a tree posterior
#'
#' @param trees a phylogeny or posterior of phylogeny of class phylo or multiphylo
#' @param states a named vector of states to calculate clades from. Names must match tips.
#'
#' @return a data.frame containing the dates of the MRCA of state clades for each tree in the posterior
#'
#'@export

tmrca = function(trees, states, progress = TRUE){

  if(!class(trees) %in% c("phylo", "multiPhylo")){
    stop("trees must be a phylo or multiPhylo class")
  }

  if(is.null(names(states))){ # or if all names are not in the tree tips
    stop("states must have names matching the tips of the trees")
  }

  if(is.null(names(trees))){
    names(trees) = seq_along(trees)
  }

  if(progress){
    nodes_dates = pbapply::pblapply(trees, function(t) .tmrca_onetree(t, states = states))
    node_dates = dplyr::bind_rows(nodes_dates, .id = "TREE")
  } else {
    nodes_dates = lapply(trees, function(t) .tmrca_onetree(t, states = states))
    node_dates = dplyr::bind_rows(nodes_dates, .id = "TREE")
  }

  node_dates
}

