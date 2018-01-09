#' @title Generate a dataset
#'
#' @description Generate an abundance dataset,
#' with or without environmental perturbance, by using generalized Lotka-Volterra
#'
#' @param samples Number of samples
#' @param matrix Interaction matrix (generated from generateA.R)
#' @param env.matrix Growth rate changes induced by environment; 1 column per environmental condition. Generate this matrix with envGrowthChanges.R
#' @param perturb.count Number of samples per environmental condition. Sum should be equal to total number of samples.
#' @param count Total number of individuals in dataset
#' @param mode Mode for generateAbundances; default value samples counts from Poisson distribution with lambda count/N
#' @return The abundance dataset
#' @examples
#' klemm = generateA(N=10, type="klemm", c=0.5)
#' env = envGrowthChanges(species = 100, env.factors=2, conditions=3, strength=0.5)
#' dataset = generateDataSet(100, klemm, env.matrix = env, perturb.count = c(50, 50))
#' @export
generateDataSet = function(samples, matrix, env.matrix=NULL, perturb.count=NULL, count=1000, mode=4){
  if (is.null(env.matrix)){
    dataset = generateSubset(samples, matrix, count = count, perturb = NULL, mode=mode)
  }
  else {
    if (length(env.matrix[1,]) != length(perturb.count)){
      stop("Please add sample count per environmental condition")
    }
    if (sum(perturb.count) != samples){
      stop("Total number of samples should equal samples per condition.")
    }
    else {
      per = perturbation(growthchanges = env.matrix[,1])
      dataset = generateSubset(perturb.count[1], matrix, perturb = per, count=count, mode=mode)
      for (i in 2:length(env.matrix[1,])){
        per = perturbation(growthchanges = env.matrix[,i])
        subset = generateSubset(perturb.count[i], matrix, perturb = per, count=count, mode=mode)
        dataset = cbind(dataset, subset)
      }
    }
  }
  colnames(dataset) = c(seq(1:samples))
  vec = vector(mode="character")
  for (k in 1:length(matrix[,1])){
    name = paste("sp", k, sep="")
    vec = c(vec, name)
  }
  rownames(dataset) = vec
  return(dataset)
}

generateSubset = function(samples, matrix, perturb=NULL, count, mode){
  N = length(matrix[,1])
  y = generateAbundances(N, count=count, mode=mode)
  dataset = matrix(nrow = N, ncol = samples)
  if (is.null(perturb) == TRUE){
    for (i in 1:samples){
      series = glv(N, y=y, matrix)
      dataset[,i] = series[,1001]
    }
  }
  else {
    for (i in 1:samples){
      series = glv(N, matrix, y=y, perturb = perturb)
      dataset[,i] = series[,1001]
    }
  }
  dataset[dataset < 0] = 0.00000001 # glv produces very small negative values, these are set to a pseudocount value
  return(dataset)
}
