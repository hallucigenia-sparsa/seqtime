#' @title Remove lowest abundance species
#'
#' @description Removes a number of lowest total abundace species from a dataset generated with generateDataSet.
#' Note that the function uses the normalized dataset to decide which species to remove, but returns the absolute dataset.
#'
#' @param dataset Abundance dataset
#' @param removeN Number of species to remove
#' @return Dataset without N lowest abundance species
#' @examples
#' klemm = generateA(N=10, type="klemm", c=0.5)
#' dataset = generateDataSet(100, klemm)
#' dataset = removeLowAbundance(dataset, 10)
#' @export

# Takes generated dataset and removes (lowest abundance) species
# Simulates effects of removing rare OTUs
removeLowAbundance = function(dataset, removeN){
  dataset2 = seqtime::normalize(dataset)
  vec = vector(mode = "numeric")
  for (i in 1:removeN){
    vec = rowSums(dataset2)
    min = which.min(vec)
    dataset = dataset[-min,]
  }
  return(dataset)
}

