# ## This function will calculate the within and between group probability for a given pair of traits
# library(ape)
# library(data.table)
# # fake data
# tree = ape::read.tree(text = "((t1, t2, t3, t4, t5),
#                                 (t6, t7, t8, t9, t10),
#                                 (t11, t12, t13, t14, t15));")
# tree = ape::compute.brlen(tree)
#
# # Traits can be a, b, or c
# # trait = factor(c("a", "a", "a", "b", "c", "b", "b", "a", "c", "c", "a", "b", "a", "b", "c"))
# trait = factor(c("a", "a", "a", "a", "a", "b", "b", "b", "b", "b", "c", "c", "c", "c", "c"))
# # trait = factor(c("a", "b", "c", "d", "e", "a", "b", "c", "d", "e", "a", "b", "c", "d", "e"))
# trait_levels = levels(trait)
# names(trait) = tree$tip.label
#
#
# ### Within group variance
# within_variance = function(){}
#
# ### Between group variance
# between_variance = function(){}

## Main function
#' Calculate within and between variance for a trait in a tree
#'
#' @param tree
#' @param generation_time
#' @param trait
#'
#' @return
#' @export
#'
#' @examples
within_between = function(tree, generation_time, trait){

  if(is.null(names(trait)) | !all(tree$tip.label %in% names(trait))){
    stop("the trait object must have names that match the tree tips")
  }

  # make traits a data.table
  trait_dt = data.table(taxa = names(trait), trait = trait)

  # Identify clades at each point in time
  clades = .slice_tree(tree, generation_time = generation_time)

  # Join traits
  clades = clades[trait_dt, on = 'taxa']

  # clade_traits = tapply(clades$taxa, list(clades$generation, clades$clade), function(x) trait[x])
  # clade_traits = t(sapply(clades, function(cc) trait[cc$taxa]))

  #### Within group variance ####
  ### Numerator: Number of observed, within clade pairs of traits

  # The count of each trait within each clade
  counts = clades[, .N, by = .(generation, clade, trait)]
  # counts = apply(clade_traits, 1, function(ct) table(factor(ct, levels = trait_levels)))
  # The number of observed pairs in each clade, summed across clades
  observed_pairs = counts[,op := choose(N, 2)]
  observed_pairs = observed_pairs[, list(op = sum(op)), by = .(generation, trait)]

  ## Denominator:: Number of possible pairs across clades
  ## Number of possible pairs across all clades
  tc = table(trait) # count of traits
  # tp = combn(tc, 2, sum)
  tp = choose(tc, 2)
  tp_dt = data.table(trait = names(tp), tp = c(tp))

  ## Within group Variability
  # (v_within = op / tp)
  observed_pairs = observed_pairs[tp_dt, on = .(trait)]
  v_within = observed_pairs[,v_within := op / tp]

  #### Between group variance ####
  # Numerator: number of observed dissimilar pairs within each clade, summed
  crossed_clades = data.table(table(clades$generation, clades$clade, clades$trait))
  colnames(crossed_clades) = c("generation", "clade", "trait", "count")

  combn_w = function(N, trait){
    if(length(N) > 1){
      combos = combn(N, 2, FUN = prod)
      nmes = combn(trait, 2, FUN=paste, collapse='')
      names(combos) = nmes
    } else {
      NA
    }
    combos
  }

  # Add observed pairs by pair type (ie. 1-2 vs 1-3)
  crossed_clades = crossed_clades[, as.list(combn_w(count, trait)), by = .(generation, clade)]

  # calculate the total number of observed pairs within clades
  numerator = crossed_clades[, .SD, .SDcols = -c("clade")]
  numerator = numerator[, as.list(colSums(.SD)), by = "generation"]

  v_btw_n = melt(numerator, id.vars = "generation", variable.name = "trait_pairs", value.name = "numerator")

  # Denominator
  # count all traits
  cats = table(trait)
  # Number of possible dissimilar pairs across clades
  denominator = combn(cats, 2, FUN = Reduce, f = "*")
  names(denominator) = combn(names(cats), 2, FUN=paste, collapse='')
  v_btw_d = data.table(trait_pairs = names(denominator), denominator = denominator)

  v_btw = v_btw_n[v_btw_d, on = "trait_pairs"]

  v_between = v_btw[, v_between := numerator / denominator, ]

  # Return
  list(v_within = v_within, v_between = v_between)
}



