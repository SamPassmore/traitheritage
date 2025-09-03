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

.extrapolate_results = function(tree, dp_df, trait, generation_time, condition = NULL, tolerance = 2){

  max_tree_depth = max(ape::node.depth.edgelength(tree)[1:ape::Ntip(tree)])
  cuts = seq(generation_time, max_tree_depth, by = generation_time)

  # Add node times to dp_df
  dp_df = get_time(dp_df, tree)

  if(is.null(condition)){
    # Calculate shared traits by node times and order by time
    numerator = dp_df[, .(numerator_node = .N, time = first(time)), by = c("node", "trait_named")][order(time, decreasing = FALSE)]
    numerator[,numerator_sum := cumsum(numerator_node), by = c("trait_named")]
    numerator[,trait_named := as.character(trait_named)]

    # Create results table for later
    result = data.table(
      generation = rep(cuts, each = length(unique(trait))),
      state = as.character(unique(trait))
    )
    setkey(result, "generation", "state")

    # Calculate denominator
    denominator = dp_df[, .(denominator_node = .N, time = first(time)), by = c("node")][order(time, decreasing = FALSE)]
    denominator[,denominator_sum := cumsum(denominator_node)]

    # Calculate Node table
    node_table = expand.grid(node = as.character((length(tree$tip.label)+1):(2*length(tree$tip.label)-1)),
                             trait_named = as.character(unique(trait)),
                             stringsAsFactors = FALSE)
  } else {
    # Create results table for later
    result = data.table(
      generation = rep(cuts, each = length(unique(condition))),
      state = as.character(condition)
    )
    setkey(result, "generation", "state")

    numdenom = custom_counter(dp_df, tree, condition = condition)
    numdenom[, trait_named := as.character(trait_named)]
    numerator = numdenom[,.(node, time, numerator_sum, trait_named)]
    denominator = numdenom[,.(node, time, denominator_sum, trait_named)]

    # make the node table
    node_table = expand.grid(node = as.character((length(tree$tip.label)+1):(2*length(tree$tip.label)-1)),
                             trait_named = as.character(condition),
                             stringsAsFactors = FALSE)
  }

  setDT(node_table, key = c("node", "trait_named"))

  # merge in numerator and denominator
  node_table = merge.data.table(node_table, denominator[,.(node, denominator_sum, time)], by = c("node"), all.x = TRUE, allow.cartesian = TRUE)
  node_table = merge.data.table(node_table, numerator[,.(node, trait_named, numerator_sum)], by = c("node", "trait_named"), all.x = TRUE)
  node_table = node_table[order(time),]
  node_table[, numerator_sum := nafill(numerator_sum, "locf"), by = c("trait_named")]
  # Identify cuts and convert to numeric
  node_table[, time.bin := cut(time, c(-Inf, cuts), labels = cuts)]
  node_table[, time.bin := as.character(levels(time.bin))[time.bin]]
  node_table[, trait_named := as.character(trait_named)]
  setkey(node_table, time)

  # merge the node table to the result table
  if(is.null(condition)){
    result = merge(result[, generation := as.character(generation)], node_table, by.x = c("generation", "state"), by.y = c("time.bin", "trait_named"), all = TRUE)
  } else {
    result = merge(result[,.(generation, state)][, generation := as.character(generation)], node_table[,.SD[.N], by = time.bin][,.(time.bin, numerator_sum, denominator_sum)],
                   by.x = c("generation"), by.y = c("time.bin"), all = TRUE)
  }
  result = result[, generation := as.numeric(generation)][order(generation)]

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

custom_counter = function(dp_df, tree, condition, coevolution = FALSE){
  ## if a new pair includes the trait of interest then calculate the cumulative proability of that node and all descendant nodes
  if(coevolution){
    dd = dp_df[, .(node, time, trait.x, trait.y, taxa.x, taxa.y, lang_trait, dist_trait)][order(time)]
  } else {
    dd = dp_df[, .(node, time, trait.x, trait.y, taxa.x, taxa.y)][order(time)]
  }

  dd[, tp := paste0(taxa.x, taxa.y),]
  nodes = unique(dd$node)
  denominator = list()
  numerator = list()
  for(i in 1:length(nodes)){
    nn = nodes[i]
    ss = dd[node == nn,]
    if(any(ss$trait.x == condition, ss$trait.y == condition)){
      descendants = phangorn::Descendants(tree, nn, type = "all")
      denominator[[i]] = unlist(dd[node %in% c(nn, descendants),.(tp)])
      if(coevolution){
        numerator[[i]] = list(
          lang_dist = unlist(dd[node %in% c(nn, descendants) & lang_trait & dist_trait,.(tp)]),
          lang_nodist = unlist(dd[node %in% c(nn, descendants) & lang_trait & !dist_trait,.(tp)]),
          nolang_dist = unlist(dd[node %in% c(nn, descendants) & !lang_trait & dist_trait,.(tp)]),
          nolang_nodist = unlist(dd[node %in% c(nn, descendants) & !lang_trait & !dist_trait,.(tp)])
        )
      } else {
        numerator[[i]] = unlist(dd[node %in% c(nn, descendants) & trait.x == trait.y & trait.x == condition,.(tp)])
      }
    } else {
      denominator[[i]] = c()
      numerator[[i]] = c()
    }
  }
  denominator_sum = sapply(seq_along(denominator), function(i) {
    length(unique(unlist(denominator[1:i])))
  })

  if(coevolution){
    subnames = names(numerator[[1]])

    numerator_sum = lapply(subnames, function(s) {
      sapply(seq_along(numerator), function(i) {
        length(unique(unlist(lapply(numerator[1:i], `[[`, s))))
      })
    })
    names(numerator_sum) = subnames

    out = data.table(
      node = nodes,
      time = unique(dd$time),
      denominator_sum = denominator_sum,
      trait_named = condition
    ) # cd is clade denominator

    out = cbind(out, data.frame(numerator_sum))

  } else {
    numerator_sum = sapply(seq_along(numerator), function(i) {
      length(unique(unlist(numerator[1:i])))
    })

    out = data.table(
      node = nodes,
      time = unique(dd$time),
      numerator_sum = numerator_sum,
      denominator_sum = denominator_sum,
      trait_named = condition
    ) # cd is clade denominator
  }


  out
}

get_time = function(dp_df, tree){

  nh <- ape::node.depth.edgelength(tree)
  # make the root 0
  nh <- max(nh) - nh

  # add times to dp_df
  nh_dt = data.table(node = as.character(1:(length(tree$edge.length) + 1)), time = nh)
  # Output
  merge.data.table(dp_df, nh_dt, by = "node", all.x = TRUE)
}
