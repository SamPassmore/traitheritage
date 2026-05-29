#' Calculate the probability of taxa being within a distance threshold
#'
#'@description
#'This function is still under development.
#'
#' @param tree a phylogenetic tree of class phylo
#' @param distance_matrix a n by n matrix of distances between all taxa
#' @param generation_time the number of generations to calculate the probability for across the phylogeny
#' @param cut_off the distance cut-off value; pairs below this value are treated as "close"
#'
#' @return a data.table showing the probability of taxa being within the
#'   distance threshold at each generation
#' @export
#' @examples
#' tree <- ape::read.tree(text = "(tA,(tB,(tC,tD)));")
#' tree <- ape::compute.brlen(tree)
#' coords <- matrix(c(1,2, 1,2, 5,6, 6,5), byrow = TRUE, ncol = 2,
#'   dimnames = list(tree$tip.label, c("X", "Y")))
#' distance_matrix <- as.matrix(dist(coords))
#' distance_trait_heritage(tree, distance_matrix, generation_time = 0.2, cut_off = 5)
distance_trait_heritage = function(tree, distance_matrix, generation_time, cut_off){
  if(any(is.na(distance_matrix))) stop("Taxa with missing values should be removed from the analysis and the tree")
  if(!all(tree$tip.label %in% rownames(distance_matrix))) stop("All taxa must match with a row in the distance matrix. Ensure row and column names are set.")
  if(!all(tree$tip.label %in% colnames(distance_matrix))) stop("All taxa must match with a row in the distance matrix. Ensure row and column names are set.")
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

paste_sort2 = function(x, y){
  paste(sort(c(x, y)), collapse=" ")
}
