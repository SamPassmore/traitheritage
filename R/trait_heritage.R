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

  if(any(is.na(trait))) stop("No NA trait values are allowed. ")
  if(is.null(names(trait))) stop("trait must have names that match the taxa. Ensure trait is a named vector.")
  if(!all(names(trait) %in% tree$tip.label)) stop("Some tips have no matching trait. Make sure all tips have a trait.")

  # clades sets at each node
  descendants = phangorn::Descendants(tree)
  names(descendants) = seq_along(descendants)
  # clades with more than one taxa
  desc_multi = descendants[sapply(descendants, length) > 1]
  desc_pairs = lapply(desc_multi, function(x) data.table(t(combn(x, 2))))
  dp_df = data.table::rbindlist(desc_pairs, idcol= "node")

  # make a reference table for taxa id to speed up taxa matching
  ref = data.table(taxa = tree$tip.label, ind = as.numeric(as.factor(tree$tip.label)))
  trait_dt = data.table(taxa = names(trait), trait = trait)
  ref = ref[trait_dt, on = "taxa"]

  # merge in trait data
  dp_df = merge.data.table(dp_df, ref, by.x = "V1", by.y = "ind", all.x = TRUE)
  dp_df = merge.data.table(dp_df, ref, by.x = "V2", by.y = "ind", all.x = TRUE)

  # identify shared traits
  dp_df[, trait := trait.x == trait.y]
  dp_df[, trait_named := ifelse(trait.x == trait.y, trait.x, "DIFFERENT")]

  # Identify which clades are under a certain time point
  max_tree_depth = max(ape::node.depth.edgelength(tree)[1:ape::Ntip(tree)]) # allows for non-ultrametric trees
  cuts = seq(generation_time, max_tree_depth, by = generation_time)

  # Full results table to fill in
  result = data.table(
    generation = rep(cuts, each = length(unique(trait))),
    state = unique(trait)
  )

  # get node ages
  nh <- ape::node.depth.edgelength(tree)
  # make the root 0
  nh <- max(nh) - nh

  # add times to dp_df
  nh_dt = data.table(node = as.character(1:(length(tree$edge.length) + 1)), time = max(nh) - ape::node.depth.edgelength(tree))
  dp_df = merge.data.table(dp_df, nh_dt, by = "node", all.x = TRUE)

  # Calculate shared traits by node times and order by time
  numerator = dp_df[trait == TRUE, .(numerator_sum = .N), by = c("time", "trait_named")][order(time, decreasing = FALSE)]
  denominator = dp_df[, .(denominator_sum = .N), by = c("time")]

  # get start end times for nodes and desired cuts
  node_times = numerator[, .(start = c(0, time[-.N]), end = time)]
  setkey(node_times, start, end)
  cuts_dt = data.table(start = cuts, end = cuts)
  setkey(cuts_dt, start, end)

  cuts_nodes = data.table::foverlaps(y = node_times, x = cuts_dt, type = "within")

  # Create probability table
  node_probs = merge.data.table(cuts_nodes, numerator, by.x = "start", by.y = "time", all = TRUE, allow.cartesian = TRUE)
  node_denom = merge.data.table(cuts_nodes, denominator, by.x = "start", by.y = "time", all = TRUE, allow.cartesian = TRUE)

  probs = merge.data.table(result, node_probs, by.x = c("generation", "state"), by.y = c("i.start", "trait_named"),
                   all.x = TRUE)
  probs = merge.data.table(probs, node_denom, by.x = "generation", by.y = "i.start", all.x = TRUE, allow.cartesian = TRUE)
  probs[, clade_probability := numerator_sum / denominator_sum,]

  # Change all NA values to 0.
  # I assume the only time this arises is in singleton clades
  probs[is.na(probs)] = 0

  ## make summary
  summary = probs[,.(numerator_sum = sum(numerator_sum), denominator_sum = sum(denominator_sum)),
                  by = "generation"][,clade_probability := numerator_sum / denominator_sum]

  return(list(
    ## Results by each level of the trait
    by_trait = probs[,c("generation", "state", "numerator_sum", "denominator_sum", "clade_probability")],
    ## Summary of results by generation
    summary = summary)
  )
}


