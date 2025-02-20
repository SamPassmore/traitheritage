# Trait Heritage: An R Package for calculating the heritage of traits on a phylogeny

<!-- badges: start -->
  [![Codecov test coverage](https://codecov.io/gh/SamPassmore/traitheritage/branch/main/graph/badge.svg)](https://app.codecov.io/gh/SamPassmore/traitheritage?branch=main)
  <!-- badges: end -->

## Installation

You can install traitheritage directly from GitHub using `devtools`

```r
devtools::install_github("SamPassmore/traitheritage", dependencies = TRUE)
```

## Usage

Trait heritage currently has three main functions (and other function in development). They are:

1. `traitheritage` : Calculate the probability of taxa sharing a trait over the life of a phylogeny
2. `distance_trait_heritage` : Calculate the probability of taxa being within a distance metric over the life of the phylogeny
3. `trait_coevolution` : Calculate the probability of taxa sharing a trait and being within a distance metric over the life of the phylogeny

## Examples

### `traitheritage`

Starting with a phylogeny, and a trait for each taxa of the phylogeny, we can calculate the probability that taxa share the same trait over the life of the phylogeny. In this case there are four taxa, three of which share a trait (B, C, and, D all have trait a), and one taxa with a different trait value (A has a trait value of b). 

```r
tree = ape::read.tree(text = "(tA,(tB,(tC,tD)));")
tree = ape::compute.brlen(tree)

trait = c("b", "a", "a", "a")
  names(trait) = tree$tip.label
```

![image](https://github.com/user-attachments/assets/881808d8-a3ad-4ff1-ba94-3efeed085263)


In this example, because the tree has a maximum height of 1, we calculate the probabilty every 0.2 steps from the taxa to the root of the tree. This value should reflect the branch lengths of whatever tree is being analysed. 
At each interval, we identify all pairs in the phylogeny that are connected, and calculate the number of pairs that share a trait divided by the total number of pairs. For example, if we look at the cut at 0.4, only tC and tD are connected as a genetic pair, meaning we have total one pair. tC and tD also share a trait (a), and so the probability of a shared trait is 1 / 1. Moving up to the 0.8 cut, we see there are now phylogenetic relationship between tB, tC, and tD, giving three total pairs. Again all traits are shared by all taxa, so the probabilty of a shared trait is again 100%. 

To calculate this we would use the following code:

```r
out = trait_heritage(tree, trait, generation_time = 0.2)
```

`out` contains a list of two dataframes: `summary`, and `by_trait1. Summary contains the average probabilty of a shared trait at each time-point and across all traits. by_trait calculates the probability of sharing each trait across the life of the phylogeny. 

This function also comes with a permutation equivalent, which provides a permuted baseline for hypothesis testing (`permute_trait_heritage`).

### `distance_trait_heritage`

Continuing with the same phylogenetic example, we might additionally know the location of each taxa, and be interested how the relationship between a secondary distance measure and the phylogenetic distance, over the life of the phylogeny. 

```r
distances = matrix(
    c(1, 2,
      2, 1,
      5, 6,
      6, 5),
    byrow = TRUE,
    ncol = 2,
    dimnames = list(t$tip.label, c("X", "Y"))
  )
  distance_matrix = as.matrix(dist(distances))
```
The function calculates distance probabilites within and beyond a certain cut-off. In this fictional example, we choose a euclidean distance of 2, but the appropriate measure will depend on what the distance metric is and the specifics of the example. This is calculated as:

```r
  result = distance_trait_heritage(
    t,
    distance_matrix = distance_matrix,
    generation_time = 0.2,
    cut_off = 2
  )
```
As with the previous function, this outputs a list of two objects, one summarising the average probabiilty over the life of the phylogeny, and one giving the probability by each trait.

### `trait_coevolution`

Finally, if we know both the trait values and the distances between taxa, we might want to know whether these two probabilities are related. This is tested in a quasi-experimental way, by dividing the sample into four categories: those who share a trait and are closer than the user defined distance, those who do not share a trait and are closer than the defined distance, those who share a trait and are further away than the user defined distance, and those who neither share a trait, or are closer than the user defined distance. 

```r
trait_coevolution(t, trait, distance_matrix, 0.2, 2)
```
This function has a slightly different output: a list of two `data.frames`, one which tells us the average probabiliity of shared distance over the life of the phylogeny, and one which shows us the relationships between all pairs, and which of the four categories that pair inhabits. This second object allows detailed evaluation of the data being used to calculate probabilities. 
