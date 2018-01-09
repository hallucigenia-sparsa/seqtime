#' @title Generate growth changes induced by environment
#'
#' @description Generate a matrix of growth changes for
#' different environmental conditions. This vector can be supplied to generateDataSet to introduce environmental perturbation in the dataset.
#' Note that the effect of environmental parameters is already summed, so only one growth change is provided per species per condition.
#'
#' @param species Number of species
#' @param env.factors Number of environmental factors
#' @param conditions Number of environmental conditions
#' @param strength Strength of environmental factors
#' @return Matrix of growth changes, 1 column per condition.
#' @examples
#' klemm = generateA(N=10, type="klemm", c=0.5)
#' env = envGrowthChanges(species = 10, env.factors=2, conditions=2, strength=0.5)
#' dataset = generateDataSet(100, klemm, env.matrix = env, perturb.count = c(50, 50))
#' @export

envGrowthChanges = function(species, env.factors=2, conditions=2, strength){
  env.matrix = matrix(nrow = species, ncol = env.factors)
  growth.matrix = matrix(nrow = species, ncol = conditions)
  for (i in 1:env.factors){
    env.matrix[,i] = rnorm(species, sd=strength)
  }
  for (j in 1:conditions){
    env.state = rnorm(env.factors)
    out = env.matrix*env.state
    growth.matrix[,j] = rowSums(out)
  }
  return(growth.matrix)
}
