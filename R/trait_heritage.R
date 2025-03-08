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

  # Remove duplicated pairs in the phylogenetic hierarchy
  dp_df = dp_df[!duplicated(dp_df, by = c("V1", "V2"), fromLast = TRUE)]

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
    state = as.character(unique(trait))
  )

  # get node ages
  nh <- ape::node.depth.edgelength(tree)
  # make the root 0
  nh <- max(nh) - nh

  # add times to dp_df
  nh_dt = data.table(node = as.character(1:(length(tree$edge.length) + 1)), time = max(nh) - ape::node.depth.edgelength(tree))
  dp_df = merge.data.table(dp_df, nh_dt, by = "node", all.x = TRUE)

  # Calculate shared traits by node times and order by time
  numerator = dp_df[trait == TRUE, .(numerator_node = .N, time = first(time)), by = c("node", "trait_named")][order(time, decreasing = FALSE)]
  numerator[,numerator_sum := cumsum(numerator_node), by = c("trait_named")]

  denominator = dp_df[, .(denominator_node = .N, time = first(time)), by = "node"][order(time, decreasing = FALSE)]
  denominator[,denominator_sum := cumsum(denominator_node)]

  # Full node table
  node_table = expand.grid(node = as.character((length(tree$tip.label)+1):(2*length(tree$tip.label)-1)),
                           trait_named = as.character(unique(trait)),
                           stringsAsFactors = FALSE)
  setDT(node_table, key = c("node", "trait_named"))

  node_table = merge.data.table(node_table, denominator[,.(node, denominator_sum, time)], by = c("node"), all.x = TRUE, allow.cartesian = TRUE)
  node_table = merge.data.table(node_table, numerator[,.(node, trait_named, numerator_sum)], by = c("node", "trait_named"), all.x = TRUE)
  node_table = node_table[order(time),]
  node_table[, numerator_sum := nafill(numerator_sum, "locf"), by = c("trait_named")]

  # get start end times for nodes and desired cuts
  # node_times = numerator[, .(start = c(0, time[-.N]), end = time)]
  # nh_dt_u = unique(nh_dt[,-c("node")])[order(time)]
  # node_times = nh_dt_u[, .(start = c(0, time[-.N]), end = time)]
  # setkey(node_times, start, end)
  # cuts_dt = data.table(start = cuts, end = cuts)
  # setkey(cuts_dt, start, end)

  # node_times = dp_df[order(time),.(end = unique(time), node = unique(node))]
  node_times = unique( dp_df[order(time),list(time, node),] )
  setnames(node_times, c("time", "node"), c("end","node"))
  node_times[,start := c(0, end[-.N])] # add a small amount to start so that intervals are separated
  node_times = node_times[-1,]
  # node_times = round(node_times, 2) # avoids problems with rounding error and comparisons
  setkey(node_times, start, end)
  cuts_dt = data.table(start = cuts, end = cuts)
  setkey(cuts_dt, start, end)

  cuts_nodes = data.table::foverlaps(y = node_times, x = cuts_dt, type = "within")

  # special case for root node
  cuts_nodes[.N, start := cuts_nodes[.N,"end"]]

  # Create probability table
  # node_probs = merge.data.table(cuts_nodes, numerator, by.x = "start", by.y = "time", all = TRUE, allow.cartesian = TRUE)
  # node_probs = merge.data.table(node_probs, denominator, by.x = "start", by.y = "time", all = TRUE, allow.cartesian = TRUE)

  cut_fraction = merge.data.table(cuts_nodes, node_table, by = "node", all = TRUE, allow.cartesian = TRUE)

  probs = merge.data.table(result, cut_fraction, by.x = c("generation", "state"), by.y = c("i.start", "trait_named"),
                   all.x = TRUE)

  probs[, clade_probability := numerator_sum / denominator_sum,]

  # Change all NA values to 0.
  # I assume the only time this arises is in singleton clades
  probs[is.na(probs)] = 0

  # subset to columns of interest
  probs = probs[,c("generation", "state", "numerator_sum", "denominator_sum", "clade_probability")]

  ## make summary
  summary = probs[,.(numerator_sum = sum(numerator_sum), denominator_sum = first(denominator_sum)),
                  by = "generation"][,clade_probability := numerator_sum / denominator_sum]

  return(list(
    ## Results by each level of the trait
    by_trait = probs,
    ## Summary of results by generation
    summary = summary)
  )
}


