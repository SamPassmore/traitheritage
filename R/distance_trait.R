# ## Distance function

distance_trait_heritage = function(tree, distance_matrix, generation_time, cut_off){
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

  # long distance
  n = nrow(distance_matrix)
  s = seq_len(n) - 1L
  nms = dimnames(distance_matrix)
  dm = data.frame(value = distance_matrix[sequence(s, seq.int(1L, length(distance_matrix), n))],
                  row = gl(n, 1L, labels = nms[[1L]])[sequence(s, 1L)],
                  col = rep.int(gl(n, 1L, labels = nms[[2L]]), s))
  dm = as.data.table(dm)
  dm = dm[ref, on = "row == taxa"][ref, on = "col == taxa"]

  paste_sort = function(ind, i.ind){
    apply(cbind(ind, i.ind), 1, function(x) paste(sort(x), collapse=" "))
  }

  dm[, idx:= do.call(paste_sort, .SD), .SDcols= c("ind", "i.ind")]

  # Calculate cut-off
  dm[, trait := value < cut_off,]

  dp_df[dm, on = "idx",  trait := i.trait]

  node_dt = dp_df[, list(numerator_sum = sum(trait),
                         denominator_sum = length(trait)), by = node]
  node_dt[, clade_probability := numerator_sum / denominator_sum]

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
    oo = data.table(numerator_sum = sum(ndt$numerator_sum),
                    denominator_sum = sum(ndt$denominator_sum))
    oo[,clade_probability := numerator_sum / denominator_sum]
  })
  names(probs) = paste0("g_", cuts)

  pp = data.table::rbindlist(probs, idcol = "generation")

  # Change all NA values to 0.
  # I assume the only time this arises is in singleton clades
  pp[is.na(pp)] = 0

  return(pp)
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
