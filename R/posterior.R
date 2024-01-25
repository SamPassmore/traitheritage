### Calculate the posterior probability of a tree

posterior_trait_heritage = function(trees, trait, n_generations){

  if(!"multiPhylo" %in% class(trees)) stop("Trees should be of class multiPhylo. If you only have one tree, use the trait_heritage function.")

  output = lapply(trees, function(t){
    trait_heritage(tree = t, trait = trait, n_generations = n_generations)
  })
  names(output) = 1:length(trees)

  dplyr::bind_rows(output, .id = "posterior")

}
