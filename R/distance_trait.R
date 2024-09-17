# ## Distance function

distance_trait_heritage = function(tree, distance_matrix, generation_time, cut_off){

  if(any(is.na(distance_matrix))) stop("Taxa with missing values should be removed from the analysis and the tree")
  if(!all(tree$tip.label %in% rownames(distance_matrix))) stop("All taxa must match with a row in the distance maxtrix. Ensure row and column names are set.")
  if(!all(tree$tip.label %in% colnames(distance_matrix))) stop("All taxa must match with a row in the distance maxtrix. Ensure row and column names are set.")

  # 1. Calculate tree cuts
  clades = .slice_tree(tree, generation_time)

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

  dm$trait = dm$value < cut_off

  # 3. For each generation, calculate the probability of a shared trait within each clade
  generations = unique(clades$generation)

  # Use lapply for generations and sapply for clade_sets
  generation_df = as.data.table(clades)
  # set key on data.table for speeding up subsetting
  setkey(generation_df, generation)

  output = lapply(generations, function(g) {
    g_df = generation_df[generation == g]

    clade_sets = unique(generation_df$clade)

    clade_result <- sapply(clade_sets, function(cs) {
      clade_taxa = g_df$taxa[g_df$clade == cs]

      trait = dm$trait[dm$row %in% clade_taxa & dm$col %in% clade_taxa]

      if (length(trait) >= 1) {
        numerator = sum(trait)
        denominator = length(trait)
      } else {
        numerator = 0
        denominator = 0
      }
      list(numerator = numerator, denominator = denominator)
    })

    # Combine clade result with generation and clade information
    data.frame(
      generation = g,
      clade = clade_sets,
      numerator = unlist(clade_result[1, ]),
      denominator = unlist(clade_result[2, ]),
      stringsAsFactors = FALSE
    )
  })
  output = do.call(rbind, output)

  # 3. Calculate the probability of a shared trait for each generation
  output_dt = data.table::data.table(output)
  output_dt = output_dt[, list(numerator_sum = sum(numerator),
                               denominator_sum = sum(denominator)), by = generation]
  output_dt[, clade_probability := numerator_sum / denominator_sum]

  # Return
  data.frame(output_dt)
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
