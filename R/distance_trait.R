# ## Distance function

distance_trait_heritage = function(tree, distance_matrix, generation_time, cut_off){

  if(any(is.na(distance_matrix))) stop("Taxa with missing values should be removed from the analysis and the tree")
  if(!all(tree$tip.label %in% rownames(distance_matrix))) stop("All taxa must match with a row in the distance maxtrix. Ensure row and column names are set.")
  if(!all(tree$tip.label %in% colnames(distance_matrix))) stop("All taxa must match with a row in the distance maxtrix. Ensure row and column names are set.")

  # Names to numbers
  ref = data.table(taxa = tree$tip.label, ind = as.numeric(as.factor(tree$tip.label)))

  # 1. Calculate tree cuts
  clades = .slice_tree(tree, generation_time)
  clades = clades[ref, on = .(taxa)]
  clades = data.table(clades)

  # 2. Calculate trait distance
  ## This solution is taken from https://stackoverflow.com/a/76731093/1544746
  ## Which notes that it will not work for R < 4.0.0
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

  # Paste indexes together for subsetting
  dm[, idx:= do.call(paste_sort, .SD), .SDcols= c("ind", "i.ind")]

  # Calculate cut-off
  dm[, trait := value < cut_off,]

  # 3. For each generation, calculate the probability of a shared trait within each clade
  generations = unique(clades$generation)

  # Use lapply for generations and sapply for clade_sets
  generation_df = as.data.table(clades)
  # set key on data.table for speeding up subsetting
  # setkey(generation_df, generation)

  ## all pairs by clade and generation
  generation_df = generation_df[, .(idx = list(ind)), by = .(generation, clade)][lengths(idx) > 1]
  generation_df = generation_df[, .(idx = sapply(idx, function(x) combn(x, 2, paste0, collapse = " "))), by = .(generation, clade)]

  # how to join where I have filtered out all clades where there is only one taxa?
  generation_df[dm[,.(idx, trait)], on = "idx", trait := i.trait]

  output_dt = generation_df[, list(numerator_sum = sum(trait),
                               denominator_sum = length(trait)), by = generation]
  output_dt[, clade_probability := numerator_sum / denominator_sum]

  generations_dt = data.table(generation = generations)
  oo = output_dt[generations_dt, on = "generation"]

  # Change all NA values to 0.
  # I assume the only time this arises is in singleton clades
  oo[is.na(oo)] = 0

  return(oo)
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

distance_trait_heritage2 = function(tree, generation_time, distance_matrix, cut_off){
  require(foreach)
  # 1. Calculate tree cuts
  clades = .slice_tree(tree, generation_time)

  # 2. Identify which pairs are under the cut-off
  taxa.numeric = as.numeric(factor(tree$tip.label))
  tm_df = data.table(taxa.numeric, tree$tip.label)

  dist.keep <- which(as.dist(distance_matrix) < cut_off)

  arr.ind <- finv(dist.keep, as.dist(distance_matrix))
  dist.dt <- data.table::data.table(taxa1 = taxa.numeric[arr.ind$i],
                                    taxa2 = taxa.numeric[arr.ind$j])

  data.table::setDT(clades)
  # clades[, taxa.numeric := as.numeric(gsub("t", "", taxa))]

  # create chunks of 100 gens
  cU <- clades[, unique(generation)]
  cN <- length(cU)
  ctz <- seq(100, cN + 100, 100)
  if (max(ctz) != cN)
    ctz <- c(ctz, cN)
  cl.chunk <- foreach::foreach (c.u = ctz) %do% {
    c.l <- c.u - 99
    IND <- cU[c.l:c.u]
    data.table::dcast(clades[generation %in% IND[!is.na(IND)]],
                      formula =  taxa.numeric ~ generation,
                      value.var = "clade")
  }

  xx <- foreach::foreach(i = cl.chunk, j = 1:length(cl.chunk)) %do% {
    print(j)
    get.prob(cl.i = i,
             T1 = dist.dt$taxa1,
             T2 = dist.dt$taxa2)
  }

  output <- clades[1:10, get.prob(trait, dist.dt), by = c("generation", "clade")]

  output[, .(clade_probability = sum(numerator) / sum(denominator)), by = "generation"]

  # Return
  return(output)
}
