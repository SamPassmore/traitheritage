#' Calculate the heritage of a trait along a single phylogeny
#'
#' @param tree a single phylogeny
#' @param trait a vector of traits, named for each taxa
#' @param generation_time how frequently the trait probability should be calculated. Must be less than the height of the phylogeny.
#'
#' @return a list containing two dataframes. by_trait contains the probability calculation at each timestep for each level of the trait. summary shows the probability across trait values for each time step.
#' @export
#'
trait_heritage = function(tree, trait, generation_time){

  ## Tests
  if(any(is.na(trait))) stop("No NA trait values are allowed.")
  if(is.null(names(trait))) stop("trait must have names that match the taxa. Ensure trait is a named vector.")
  if(!all(names(trait) %in% tree$tip.label)) stop("Some tips have no matching trait. Make sure all tips have a trait.")

  dp_df = .get_hierarchy(tree)

  # make a reference table of taxa id to speed up taxa matching
  ref = data.table(taxa = tree$tip.label, ind = as.numeric(as.factor(tree$tip.label)))
  trait_dt = data.table(taxa = names(trait), trait = trait)
  ref = ref[trait_dt, on = "taxa"]

  # merge reference table into the trait data
  dp_df = merge.data.table(dp_df, ref, by.x = "V1", by.y = "ind", all.x = TRUE)
  dp_df = merge.data.table(dp_df, ref, by.x = "V2", by.y = "ind", all.x = TRUE)

  # identify shared traits
  dp_df[, trait := trait.x == trait.y]
  dp_df[, trait_named := ifelse(trait.x == trait.y, trait.x, "DIFF")]

  # Identify which clades are under a certain time point
  # allows for non-ultrametric trees

  result = .extrapolate_results(tree, dp_df, trait, generation_time)

  ## make summary
  summary = result[,.(numerator_sum = sum(numerator_sum), denominator_sum = first(denominator_sum)),
                  by = "generation"][,clade_probability := numerator_sum / denominator_sum]

  return(list(
    ## Results by each level of the trait
    by_trait = result,
    ## Summary of results by generation
    summary = summary)
  )
}


