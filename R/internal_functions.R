# Internal functions

.get_hierarchy = function(tree, .DescendantsType = "tips"){
  ### First identify the hierachical relationship of the tree
  # clades sets at each node
  descendants = phangorn::Descendants(tree, type = .DescendantsType)
  names(descendants) = seq_along(descendants)

  # Deal with clades that contain more than one taxa
  desc_multi = descendants[sapply(descendants, length) > 1]
  desc_pairs = lapply(desc_multi, function(x) data.table(t(combn(x, 2))))
  dp_df = data.table::rbindlist(desc_pairs, idcol= "node")

  # Remove duplicated pairs in the phylogenetic hierarchy (i.e. reversed)
  dp_df = dp_df[!duplicated(dp_df, by = c("V1", "V2"), fromLast = TRUE)]

  dp_df
}

.extrapolate_results = function(tree, dp_df, trait, generation_time){

  max_tree_depth = max(ape::node.depth.edgelength(tree)[1:ape::Ntip(tree)])
  cuts = seq(generation_time, max_tree_depth, by = generation_time)

  # Full results table to fill in
  result = data.table(
    generation = rep(cuts, each = length(unique(trait))),
    state = as.character(unique(trait))
  )
  setkey(result, generation, state)

  nh <- ape::node.depth.edgelength(tree)
  # make the root 0
  nh <- max(nh) - nh

  # add times to dp_df
  nh_dt = data.table(node = as.character(1:(length(tree$edge.length) + 1)), time = nh)
  dp_df = merge.data.table(dp_df, nh_dt, by = "node", all.x = TRUE)

  # Calculate shared traits by node times and order by time
  numerator = dp_df[trait == TRUE, .(numerator_node = .N, time = first(time)), by = c("node", "trait_named")][order(time, decreasing = FALSE)]
  numerator[,numerator_sum := cumsum(numerator_node), by = c("trait_named")]
  numerator[,trait_named := as.character(trait_named)]

  denominator = dp_df[, .(denominator_node = .N, time = first(time)), by = "node"][order(time, decreasing = FALSE)]
  denominator[,denominator_sum := cumsum(denominator_node)]

  # Full node table
  node_table = expand.grid(node = as.character((length(tree$tip.label)+1):(2*length(tree$tip.label)-1)),
                           trait_named = as.character(unique(trait)),
                           stringsAsFactors = FALSE)
  setDT(node_table, key = c("node", "trait_named"))

  # merge in numerator and denominator
  node_table = merge.data.table(node_table, denominator[,.(node, denominator_sum, time)], by = c("node"), all.x = TRUE, allow.cartesian = TRUE)
  node_table = merge.data.table(node_table, numerator[,.(node, trait_named, numerator_sum)], by = c("node", "trait_named"), all.x = TRUE)
  node_table = node_table[order(time),]
  node_table[, numerator_sum := nafill(numerator_sum, "locf"), by = c("trait_named")]
  # Identify cuts and convert to numeric
  node_table[, time.bin := cut(time, cuts, labels = cuts[-1])]
  node_table[, time.bin := as.numeric(levels(time.bin))[time.bin]]
  setkey(node_table, time)

  # merge the node table to the result table
  result = merge(result, node_table, by.x = c("generation", "state"), by.y = c("time.bin", "trait_named"), all = TRUE)
  # Fill in the missing data, with the previous data (no node changes)
  result = result[, numerator_sum := nafill(numerator_sum, "locf"), by = c("generation", "state")]
  result = result[, `:=` (numerator_sum = nafill(numerator_sum, "locf"),
                          denominator_sum = nafill(denominator_sum, "locf")),
                  by = c("state")]

  result[, clade_probability := numerator_sum / denominator_sum,]

  # Change all NA values to 0.
  # I assume the only time this arises is in singleton clades
  result[is.na(result)] = 0

  # subset to columns of interest
  result[,c("generation", "state", "numerator_sum", "denominator_sum", "clade_probability")]
}
