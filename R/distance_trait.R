# ## Distance function

distance_trait_heritage = function(tree, distance_matrix, generation_time, cut_off){
  if(any(is.na(distance_matrix))) stop("Taxa with missing values should be removed from the analysis and the tree")
  if(!all(tree$tip.label %in% rownames(distance_matrix))) stop("All taxa must match with a row in the distance maxtrix. Ensure row and column names are set.")
  if(!all(tree$tip.label %in% colnames(distance_matrix))) stop("All taxa must match with a row in the distance maxtrix. Ensure row and column names are set.")
  if(max(ape::node.depth.edgelength(tree))/generation_time <= 2) stop("You must make more than one cut in the tree.")

  # clades sets at each node
  # descendants = phangorn::Descendants(tree)
  # names(descendants) = seq_along(descendants)
  # # clades with more than one taxa
  # desc_multi = descendants[sapply(descendants, length) > 1]
  # desc_pairs = lapply(desc_multi, function(x) data.table(t(combn(x, 2))))
  # dp_df = data.table::rbindlist(desc_pairs, idcol= "node")
  #
  # # Remove duplicated pairs in the phylogenetic hierarchy
  # dp_df = dp_df[!duplicated(dp_df, by = c("V1", "V2"), fromLast = TRUE)]

  dp_df = .get_hierarchy(tree)

  dp_df = dp_df[, idx := paste_sort2(V1, V2), by = seq_len(nrow(dp_df))]

  # make a reference table for taxa id to speed up taxa matching
  ref = data.table(taxa = tree$tip.label, ind = as.numeric(as.factor(tree$tip.label)))

  # merge in trait data
  dp_df = merge.data.table(dp_df, ref, by.x = "V1", by.y = "ind", all.x = TRUE)
  dp_df = merge.data.table(dp_df, ref, by.x = "V2", by.y = "ind", all.x = TRUE)

  # long distance
  n = nrow(distance_matrix)
  s = seq_len(n) - 1L
  nms = dimnames(distance_matrix)
  dm = data.table(value = distance_matrix[sequence(s, seq.int(1L, length(distance_matrix), n))],
                  row = gl(n, 1L, labels = nms[[1L]])[sequence(s, 1L)],
                  col = rep.int(gl(n, 1L, labels = nms[[2L]]), s))
  dm = dm[ref, on = "row == taxa"][ref, on = "col == taxa"]
  # sometimes this join procudes rows with NAs if a taxa is in ref but not dm
  # remove those
  dm = dm[complete.cases(value)]

  # dm[, idx:= do.call(paste_sort, .SD), .SDcols= c("ind", "i.ind")]
  dm[, idx := paste_sort2(ind, i.ind), by = seq_len(nrow(dm))]

  # Calculate cut-off
  dm[, trait := value < cut_off,]
  dm[, trait_named := ifelse(trait, 1, 0)]


  # join distance cut off to nodes table
  dp_df = dp_df[dm, on = "idx"]
  dp_df[,c("value", "row", "col", "ind", "i.ind") := NULL]

  # Identify which clades are under a certain time point
  trait = as.numeric(dp_df$trait)
  result = .extrapolate_results(tree, dp_df, trait, generation_time)

  # max_tree_depth = max(ape::node.depth.edgelength(tree)[1:ape::Ntip(tree)]) # allows for non-ultrametric trees
  # cuts = seq(generation_time, max_tree_depth, by = generation_time)
  #
  # # Full results table to fill in
  # result = data.table(
  #   generation = rep(cuts, each = 2), # there are only two states in the distance function
  #   state = as.character(c(0, 1))
  # )
  #
  # # get node ages
  # nh <- ape::node.depth.edgelength(tree)
  # # make the root 0
  # nh <- max(nh) - nh
  #
  # # add times to dp_df
  # nh_dt = data.table(node = as.character(1:(length(tree$edge.length) + 1)), time = max(nh) - ape::node.depth.edgelength(tree))
  # dp_df = merge.data.table(dp_df, nh_dt, by = "node", all.x = TRUE)
  #
  # # Calculate shared traits by node times and order by time
  # numerator = dp_df[trait == TRUE, .(numerator_node = .N, time = first(time)), by = c("node", "trait_named")][order(time, decreasing = FALSE)]
  # numerator[,numerator_sum := cumsum(numerator_node), by = c("trait_named")]
  # # numerator[,trait_named := as.character(trait_named)]
  #
  # denominator = dp_df[, .(denominator_node = .N, time = first(time)), by = "node"][order(time, decreasing = FALSE)]
  # denominator[,denominator_sum := cumsum(denominator_node)]
  #
  # # Full node table
  # node_table = expand.grid(node = as.character((length(tree$tip.label)+1):(2*length(tree$tip.label)-1)),
  #                          trait_named = as.character(as.numeric(unique(dp_df$trait))),
  #                          stringsAsFactors = FALSE)
  # setDT(node_table, key = c("node", "trait_named"))
  #
  # # merge in numerator and denominator
  # node_table = merge.data.table(node_table, denominator[,.(node, denominator_sum, time)], by = c("node"), all.x = TRUE, allow.cartesian = TRUE)
  # node_table = merge.data.table(node_table, numerator[,.(node, trait_named, numerator_sum)], by = c("node", "trait_named"), all.x = TRUE)
  # node_table = node_table[order(time),]
  # node_table[, numerator_sum := nafill(numerator_sum, "locf"), by = c("trait_named")]
  # # Identify cuts and convert to numeric
  # node_table[, time.bin := cut(time, c(0, cuts), labels = c(0, cuts)[-1])]
  # node_table[, time.bin := as.numeric(levels(time.bin))[time.bin]]
  # setkey(node_table, time)
  #
  # # merge the node table to the result table
  # result = merge(result, node_table, by.x = c("generation", "state"), by.y = c("time.bin", "trait_named"), all = TRUE)
  # # Fill in the missing data, with the previous data (no node changes)
  # result = result[, numerator_sum := nafill(numerator_sum, "locf"), by = c("generation", "state")]
  # result = result[, `:=` (numerator_sum = nafill(numerator_sum, "locf"),
  #                         denominator_sum = nafill(denominator_sum, "locf")),
  #                 by = c("state")]
  # ## fill in empty sets with zero
  # result[][is.na(time), time := 0]
  #
  # # ensure matches to the oldest cut when there are multiple nodes between cuts.
  #
  # result = result[, .SD[which.max(time)], by = c("generation", "state")]
  # # specsial case for root node
  # # cuts_nodes[.N, start := cuts_nodes[.N,"end"]]
  #
  # # Create probability table
  # # node_probs = merge.data.table(cuts_nodes, numerator, by.x = "start", by.y = "time", all = TRUE, allow.cartesian = TRUE)
  # # node_probs = merge.data.table(node_probs, denominator, by.x = "start", by.y = "time", all = TRUE, allow.cartesian = TRUE)
  #
  # # cut_fraction = merge.data.table(cuts_nodes, node_table, by = "node", all = TRUE, allow.cartesian = TRUE)
  #
  # # probs = merge.data.table(result, cut_fraction, by.x = c("generation", "state"), by.y = c("i.start", "trait_named"),
  # #                          all.x = TRUE)
  # probs = result
  # probs[, clade_probability := numerator_sum / denominator_sum,]
  #
  # # Change all NA values to 0.
  # # I assume the only time this arises is in singleton clades
  # probs[is.na(probs)] = 0
  #
  # # subset to columns of interest
  # probs = probs[,c("generation", "state", "numerator_sum", "denominator_sum", "clade_probability")]

  ## make summary
  summary = result[,.(numerator_sum = sum(numerator_sum), denominator_sum = first(denominator_sum)),
                  by = "generation"][,clade_probability := numerator_sum / denominator_sum]

  return(
    ## Summary of results by generation
    summary
  )
}

# https://stackoverflow.com/questions/39005958/r-how-to-get-row-column-subscripts-of-matched-elements-from-a-distance-matri
finv <- function (k, dist_obj) {
  if (!inherits(dist_obj, "dist"))
    stop("please provide a 'dist' object")
  n <- attr(dist_obj, "Size")
  valid <- (k >= 1) & (k <= n * (n - 1) / 2)
  k_valid <- k[valid]
  j <- rep.int(NA_real_, length(k))
  j[valid] <-
    floor(((2 * n + 1) - sqrt((2 * n - 1) ^ 2 - 8 * (k_valid - 1))) / 2)
  i <- j + k - (2 * n - j) * (j - 1) / 2
  data.frame(i = i, j = j)
}

get.prob <- function(cl.i, T1, T2) {
  A <- cl.i[T1]
  B <- cl.i[T2]
  nc <- ncol(cl.i)
  D0 <- apply(cl.i[, 2:nc], 2, function(x) choose(table(x), 2))
  D <- sapply(D0, sum)
  N <- colSums(A[, -1] == B[, -1])
  return(list(numerator = N, denominator = D))
}

paste_sort = function(ind, i.ind){
  apply(cbind(ind, i.ind), 1, function(x) paste(sort(x), collapse=" "))
}

paste_sort2 = function(x, y){
  paste(sort(c(x, y)), collapse=" ")
}
