## trait co-evolution

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

  # clades sets at each node
  descendants = phangorn::Descendants(tree)
  names(descendants) = seq_along(descendants)
  # clades with more than one taxa
  desc_multi = descendants[sapply(descendants, length) > 1]

  desc_pairs = lapply(desc_multi, function(x) data.table(t(combn(x, 2))))
  dp_df = data.table::rbindlist(desc_pairs, idcol= "node")
  dp_df[,idx := do.call(paste, c(.SD, sep = " ")), .SDcols = c("V1", "V2")]

  # Create alias names for taxa
  ref = data.table(taxa = tree$tip.label, ind = as.numeric(as.factor(tree$tip.label)))

  trait = data.table(taxa = names(trait), trait = trait)
  ref[trait, on = "taxa", trait := trait]

  # long distance
  n = nrow(distance_matrix)
  s = seq_len(n) - 1L
  nms = dimnames(distance_matrix)
  dm = data.frame(value = distance_matrix[sequence(s, seq.int(1L, length(distance_matrix), n))],
                  row = gl(n, 1L, labels = nms[[1L]])[sequence(s, 1L)],
                  col = rep.int(gl(n, 1L, labels = nms[[2L]]), s))
  dm = as.data.table(dm)
  # dm = dm[ref, on = "row == taxa"][ref, on = "col == taxa"]
  dm = dm[ref, on = "row == taxa", `:=`(trait1 = trait, row_num = ind)][ref, on = "col == taxa", `:=`(trait2 = trait, col_num = ind)]

  paste_sort = function(row_num, col_num){
    apply(cbind(row_num, col_num), 1, function(x) paste(sort(x), collapse=" "))
  }

  dm[, idx:= do.call(paste_sort, .SD), .SDcols= c("row_num", "col_num")]

  # Calculate distance cut-off
  dm[, dist_trait := value < cut_off,]
  # calculate trait similarity
  dm[, lang_trait := trait1 == trait2,]

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
  cuts = head(seq(0, max_tree_depth, by = generation_time), -1)
  nh <- ape::node.depth.edgelength(tree)
  nh <- max(nh) - nh

  probs = lapply(cuts, function(cut){
    ind <- which((nh[tree$edge[, 1]] > cut) &
                   (nh[tree$edge[, 2]] <= cut))
    res = names(descendants[tree$edge[ind, 2]])
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
