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
    p.trait_heritage(tree, p_trait, generation_time)
  })
  names(permuted_list) = paste0("p_", 1:n_permutations)

  permuted = do.call(rbind, permuted_list)
  permuted$iteration = lapply(strsplit(row.names(permuted), "\\."), '[[', 1)
  rownames(permuted) = NULL

  list(true = true, permuted = permuted)
}

p.trait_heritage = function(tree, trait, generation_time) {
  ## Add argument tests to the function
  if (any(is.na(trait)))
    stop("No NA trait values are allowed. ")

  # Trait names must match taxa labels
  if (!all(names(trait) %in% tree$tip.label))
    stop("Some tips have no matching trait. Make sure all tips have a trait.")

  # clades sets at each node
  descendants = phangorn::Descendants(tree)
  names(descendants) = seq_along(descendants)
  # clades with more than one taxa
  desc_multi = descendants[sapply(descendants, length) > 1]

  desc_pairs = lapply(desc_multi, function(x)
    data.table(t(combn(x, 2))))
  dp_df = data.table::rbindlist(desc_pairs, idcol = "node")
  #dp_df[,idx := do.call(paste, c(.SD, sep = " ")), .SDcols = c("V1", "V2")]

  # Create alias names for taxa
  ref = data.table(taxa = tree$tip.label, ind = as.numeric(as.factor(tree$tip.label)))

  trait = data.table(taxa = names(trait), trait = trait)
  ref[trait, on = "taxa", trait := trait]

  dp_df = dp_df[ref, on = "V1 == ind", `:=`(trait1 = trait, row_num = ind)][ref, on = "V2 == ind", `:=`(trait2 = trait, col_num = ind)]

  dp_df[, paired := trait1 == trait2, ]

  node_dt = dp_df[, list(numerator = sum(paired),
                         denominator_sum = length(paired)), by = node][, clade_probability := numerator / denominator_sum, ]

  # Identify which clades are under a ceratin time point
  max_tree_depth = max(ape::node.depth.edgelength(tree)[1:ape::Ntip(tree)]) # allows for non-ultrametric trees
  cuts = head(seq(0, max_tree_depth, by = generation_time), -1)
  nh <- ape::node.depth.edgelength(tree)
  nh <- max(nh) - nh

  probs = lapply(cuts, function(cut) {
    ind <- which((nh[tree$edge[, 1]] > cut) &
                   (nh[tree$edge[, 2]] <= cut))
    res = names(descendants[tree$edge[ind, 2]])
    ndt = node_dt[node %in% res]
    oo = data.table(
      numerator_sum = sum(ndt$numerator),
      denominator_sum = sum(ndt$denominator_sum)
    )
    oo[, clade_probability := numerator_sum / denominator_sum, ]
  })
  names(probs) = paste0("g_", cuts)

  pp = data.table::rbindlist(probs, idcol = "generation")

  # Change all NA values to 0.
  # I assume the only time this arises is in singleton clades
  pp[is.na(pp)] = 0

  return(pp)
}


