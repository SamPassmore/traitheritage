# ## Distance function

distance_trait_heritage = function(tree, distance_matrix, generation_time, cut_off){
  if(any(is.na(distance_matrix))) stop("Taxa with missing values should be removed from the analysis and the tree")
  if(!all(tree$tip.label %in% rownames(distance_matrix))) stop("All taxa must match with a row in the distance maxtrix. Ensure row and column names are set.")
  if(!all(tree$tip.label %in% colnames(distance_matrix))) stop("All taxa must match with a row in the distance maxtrix. Ensure row and column names are set.")
  if(max(ape::node.depth.edgelength(tree))/generation_time <= 2) stop("You must make more than one cut in the tree.")

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
