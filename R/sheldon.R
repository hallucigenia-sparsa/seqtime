#' Compute evenness using Sheldon's index
#'
#' @references A.L. Sheldon 1969. Equitability indices: dependence on the species count. Ecology, 50, 466-467.
#'
#' Note that the N2N1 mode results in evenness smaller than 1 for equal taxon
#' probabilities.
#'
#' @params x a vector
#' @params correction whether or not to apply the correction described in Alatalo, Oikos 37, 199-204, 1981
#' @params N2N1 whether to compute Sheldon's evenness as the ratio of e raised to the power of H (H = Shannon diversity) and Simpson's diversity
#' @return Sheldon's evenness
#' @export
##############################################################################
sheldon<-function(x, correction = TRUE, N2N1 = FALSE){
  H = vegan::diversity(x, index="shannon")
  if(N2N1){
    simpson = vegan::diversity(x, index="simpson")
    numerator = 1/simpson
    denominator = exp(1)^H
  }else{
    numerator = exp(1)^H
    denominator = vegan::specnumber(x)
  }
  if(correction){
    numerator = numerator - 1
    denominator = denominator - 1
  }
  S = numerator/denominator
  S
}

##############################################################################
# Compute evenness using Pielou's index
#
# Arguments:
# x = matrix-like object or vector
#
# Value:
# pielou's evenness
#
##############################################################################
pielou<-function(x){
  H = vegan::diversity(x, index="shannon")
  J = H/log(vegan::specnumber(x))
  J
}
