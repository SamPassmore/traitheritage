## trait co-evolution

paste_sort = function(row_num, col_num){
  apply(cbind(row_num, col_num), 1, function(x) paste(sort(x), collapse=" "))
}

#' Calculate the relative probability of a shared trait and distance
#'
#' @param tree a phylogenetic tree of classe phylo
#' @param trait a numeric or character trait for each taxa
#' @param distance_matrix a n by n matrix of distances between all taxa
#' @param generation_time the number of generations to calculate the probaiblity for across the phylogeny
#' @param cut_off The distance cut off value (used to discretize distance)
#'
#' @return a list. time_df shows the probability for each paired category over time. pairs_df shows the data for each viable pair in the dataset.
#' @export
#'
trait_coevolution = function(tree, trait, distance_matrix, generation_time, cut_off){

  if(any(is.na(distance_matrix))) stop("Taxa with missing values should be removed from the analysis and the tree")
  if(!all(tree$tip.label %in% rownames(distance_matrix))) stop("All taxa must match with a row in the distance maxtrix. Ensure row and column names are set.")
  if(!all(tree$tip.label %in% colnames(distance_matrix))) stop("All taxa must match with a row in the distance maxtrix. Ensure row and column names are set.")
  if(max(ape::node.depth.edgelength(tree))/generation_time <= 2) stop("You must make more than one cut in the tree.")

  # # clades sets at each node
  dp_df = .get_hierarchy(tree, .DescendantsType = "all")

  # For this function specifically, we need the object descendants (which we usually calculate in the interal function)
  descendants = phangorn::Descendants(tree, type = "all")
  names(descendants) = seq_along(descendants)

  # dp_df[,idx := do.call(paste_sort, c(.SD, sep = " ")), .SDcols = c("V1", "V2")]
  dp_df[,idx := paste(pmin(V1, V2), pmax(V1, V2), sep = " ")]
  # Create alias names for taxa
  ref = data.table(taxa = tree$tip.label, ind = as.numeric(as.factor(tree$tip.label)))

  trait = data.table(taxa = names(trait), trait = trait)
  ref[trait, on = "taxa", trait := trait]


  # long distance
  n = nrow(distance_matrix)
  s = seq_len(n) - 1L
  nms = dimnames(distance_matrix)
  dm = data.frame(value = distance_matrix[sequence(s, seq.int(1L, length(distance_matrix), n))],
                  taxa.x = gl(n, 1L, labels = nms[[1L]])[sequence(s, 1L)],
                  taxa.y = rep.int(gl(n, 1L, labels = nms[[2L]]), s))
  dm = as.data.table(dm)
  # dm = dm[ref, on = "row == taxa"][ref, on = "col == taxa"]
  dm = dm[ref, on = "taxa.x == taxa", `:=`(trait.x = trait, row_num = ind)][ref, on = "taxa.y == taxa", `:=`(trait.y = trait, col_num = ind)]

  dm[, idx:= do.call(paste_sort, .SD), .SDcols= c("row_num", "col_num")]

  # Calculate distance cut-off
  dm[, dist_trait := value < cut_off,]
  # calculate trait similarity
  dm[, lang_trait := trait.x == trait.y,]

  dp_df = dp_df[dm, on = "idx"]

  ## Possible conditions
  ## D & L, ND & L, D & NL, ND & NL
  node_dt = dp_df[, list(lang_dist = sum(lang_trait & dist_trait),
                         lang_nodist = sum(lang_trait == TRUE & dist_trait == FALSE),
                         nolang_dist = sum(lang_trait == FALSE & dist_trait == TRUE),
                         nolang_nodist = sum(lang_trait == FALSE & dist_trait == FALSE),
                         denominator_sum = length(lang_trait)),
                         by = node]
  node_dt[, `:=`(p_lang_dist = lang_dist / denominator_sum,
                 p_lang_nodist = lang_nodist / denominator_sum,
                 p_nolang_dist = nolang_dist / denominator_sum,
                 p_nolang_nodist = nolang_nodist / denominator_sum)]

  # Identify which clades are under a ceratin time point
  max_tree_depth = max(ape::node.depth.edgelength(tree)[1:ape::Ntip(tree)]) # allows for non-ultrametric trees
  # cuts = head(seq(0, max_tree_depth, by = generation_time), -1)
  cuts = seq(0, max_tree_depth + generation_time, by = generation_time)
  nh <- ape::node.depth.edgelength(tree)
  nh <- max(nh) - nh

  probs = lapply(cuts, function(cut){
    ind <- which((nh[tree$edge[, 1]] > cut) &
                   (nh[tree$edge[, 2]] <= cut))

    # res = names(descendants[tree$edge[ind, 2]])
    # This is a hack to get all nodes and children nodes under a cut.
    # Since we rely on a table on nodes, including taxa shouldn't cause an issue.
    res = c(names(descendants[tree$edge[ind, 2]]), unlist(descendants[tree$edge[ind, 2]]))

    # special case for max
    if(cut >= max(nh)){
      res = node_dt$node
    }

    ndt = node_dt[node %in% res]
    oo = data.table(lang_dist = sum(ndt$lang_dist),
                    lang_nodist = sum(ndt$lang_nodist),
                    nolang_dist = sum(ndt$nolang_dist),
                    nolang_nodist = sum(ndt$nolang_nodist),
                    denominator_sum = sum(ndt$denominator_sum))
    oo[, `:=`(p_lang_dist = lang_dist / denominator_sum,
              p_lang_nodist = lang_nodist / denominator_sum,
              p_nolang_dist = nolang_dist / denominator_sum,
              p_nolang_nodist = nolang_nodist / denominator_sum)]
  })
  names(probs) = paste0("g_", cuts)

  pp = data.table::rbindlist(probs, idcol = "generation")

  # Change all NA values to 0.
  # I assume the only time this arises is in singleton clades
  pp[is.na(pp)] = 0

  return(list(time_df = pp, pairs_df = dp_df))
}


trait_coevolution_permutation =  function(tree, trait, distance_matrix, generation_time, cut_off, n_permutation){
  set.seed(seed)
  # Permuted results
  permuted_list = lapply(1:n_permutations, function(i){
    p_trait = trait
    names(p_trait) = names(p_trait)[sample(1:length(p_trait))]
    trait_coevolution(tree, p_trait, distance_matrix, generation_time, cut_off)
  })

  by_trait = lapply(permuted_list, "[[", 1)
  p_summary = lapply(permuted_list, "[[", 2)

  names(by_trait) = paste0("p_", 1:n_permutations)
  names(p_summary) = paste0("p_", 1:n_permutations)

  permuted_bytrait = data.table::rbindlist(by_trait, idcol = "iteration")
  permuted_summary = data.table::rbindlist(p_summary, idcol = "iteration")

  return(list(by_trait = permuted_bytrait, summary = permuted_summary))
  }


trait_coevolution_specific = function(tree, trait, distance_matrix, generation_time, cut_off, condition){

  if(any(is.na(distance_matrix))) stop("Taxa with missing values should be removed from the analysis and the tree")
  if(!all(tree$tip.label %in% rownames(distance_matrix))) stop("All taxa must match with a row in the distance maxtrix. Ensure row and column names are set.")
  if(!all(tree$tip.label %in% colnames(distance_matrix))) stop("All taxa must match with a row in the distance maxtrix. Ensure row and column names are set.")
  if(max(ape::node.depth.edgelength(tree))/generation_time <= 2) stop("You must make more than one cut in the tree.")


  max_tree_depth = max(ape::node.depth.edgelength(tree)[1:ape::Ntip(tree)]) # allows for non-ultrametric trees
  cuts = seq(0, max_tree_depth + generation_time, by = generation_time)

  result = data.table(
    generation = rep(cuts, each = length(unique(condition))),
    state = as.character(condition)
  )
  setkey(result, "generation", "state")

  # # clades sets at each node
  dp_df = .get_hierarchy(tree, .DescendantsType = "all")

  # For this function specifically, we need the object descendants (which we usually calculate in the interal function)
  descendants = phangorn::Descendants(tree, type = "all")
  names(descendants) = seq_along(descendants)

  # dp_df[,idx := do.call(paste_sort, c(.SD, sep = " ")), .SDcols = c("V1", "V2")]
  dp_df[,idx := paste(pmin(V1, V2), pmax(V1, V2), sep = " ")]
  # Create alias names for taxa
  ref = data.table(taxa = tree$tip.label, ind = as.numeric(as.factor(tree$tip.label)))

  trait = data.table(taxa = names(trait), trait = trait)
  ref[trait, on = "taxa", trait := trait]


  # long distance
  n = nrow(distance_matrix)
  s = seq_len(n) - 1L
  nms = dimnames(distance_matrix)
  dm = data.frame(value = distance_matrix[sequence(s, seq.int(1L, length(distance_matrix), n))],
                  taxa.x = gl(n, 1L, labels = nms[[1L]])[sequence(s, 1L)],
                  taxa.y = rep.int(gl(n, 1L, labels = nms[[2L]]), s))
  dm = as.data.table(dm)
  # dm = dm[ref, on = "row == taxa"][ref, on = "col == taxa"]
  dm = dm[ref, on = "taxa.x == taxa", `:=`(trait.x = trait, row_num = ind)][ref, on = "taxa.y == taxa", `:=`(trait.y = trait, col_num = ind)]

  dm[, idx:= do.call(paste_sort, .SD), .SDcols= c("row_num", "col_num")]

  # Calculate distance cut-off
  dm[, dist_trait := value < cut_off,]
  # calculate trait similarity
  dm[, lang_trait := trait.x == trait.y & trait.x == condition,]

  dp_df = dp_df[dm, on = "idx"]
  dp_df = get_time(dp_df, tree)

  node_dt = custom_counter(dp_df, tree, condition, coevolution = TRUE)

  node_dt[, `:=`(p_lang_dist = lang_dist / denominator_sum,
                 p_lang_nodist = lang_nodist / denominator_sum,
                 p_nolang_dist = nolang_dist / denominator_sum,
                 p_nolang_nodist = nolang_nodist / denominator_sum)]

  node_dt[, time.bin := cut(time, c(-Inf, cuts), labels = cuts)]
  node_dt[, time.bin := as.character(levels(time.bin))[time.bin]]
  node_dt[, trait_named := as.character(trait_named)]
  setkey(node_dt, time)

  result = merge(
    result[, .(generation, state)][, generation := as.character(generation)],
    node_dt[, .SD[.N], by = time.bin][, .(time.bin,
                                          lang_dist,
                                          lang_nodist,
                                          nolang_dist,
                                          nolang_nodist,
                                          denominator_sum)],
    by.x = c("generation"),
    by.y = c("time.bin"),
    all = TRUE
  )

  result = result[, generation := as.numeric(generation)][order(generation)]

  # Fill in the missing data, with the previous data (no node changes)
  result = result[, `:=` (lang_dist = nafill(lang_dist, "locf"),
                          lang_nodist = nafill(lang_nodist, "locf"),
                          nolang_dist = nafill(nolang_dist, "locf"),
                          nolang_nodist = nafill(nolang_nodist, "locf"),
                          denominator_sum = nafill(denominator_sum, "locf"))]

  # Calculate the probbilities
  result = result[, `:=`(
    p_lang_dist = lang_dist / denominator_sum,
    p_lang_nodist = lang_nodist / denominator_sum,
    p_nolang_dist = nolang_dist / denominator_sum,
    p_nolang_nodist = nolang_nodist / denominator_sum
  )
  ]

  # Change all NA values to 0.
  # I assume the only time this arises is in singleton clades
  result[is.na(result)] = 0

  return(list(time_df = result, pairs_df = dp_df))
}


trait_coevolution_permutation_specific = function(tree,
                                                  trait,
                                                  distance_matrix,
                                                  generation_time,
                                                  cut_off,
                                                  condition,
                                                  n_permutations,
                                                  seed = 9872) {

  set.seed(seed)
  # Permuted results
  permuted_list = lapply(1:n_permutations, function(i){
    p_trait = trait
    names(p_trait) = names(p_trait)[sample(1:length(p_trait))]
    trait_coevolution_specific(tree, p_trait, distance_matrix, generation_time, cut_off, condition)
  })

  by_trait = lapply(permuted_list, "[[", 1)
  p_summary = lapply(permuted_list, "[[", 2)

  names(by_trait) = paste0("p_", 1:n_permutations)
  names(p_summary) = paste0("p_", 1:n_permutations)

  permuted_bytrait = data.table::rbindlist(by_trait, idcol = "iteration")
  permuted_summary = data.table::rbindlist(p_summary, idcol = "iteration")

  return(list(by_trait = permuted_bytrait, summary = permuted_summary))
}
