#' @title Compute evenness using Sheldon's index
#'
#' @description Sheldon's index is defined as \eqn{S=\frac{e^H}{N}}, where H is the Shannon diversity and N the species number.
#' It ranges from 0 to 1, where 1 signifies a perfectly even abundance distribution.
#'
#' @references A.L. Sheldon 1969. Equitability indices: dependence on the species count. Ecology, 50, 466-467.
#' @references C Heip 1974. A new index measuring evenness. J. mar. biol. Ass. UK 54, 555-557.
#'
#' @details Note that the N2N1 mode results in evenness smaller than 1 for equal taxon probabilities.
#'
#' @param x a vector of species abundances
#' @param correction whether or not to apply the correction described in Alatalo, Oikos 37, 199-204, 1981
#' @param N2N1 whether to compute Sheldon's evenness as the ratio of e raised to the power of H (H = Shannon diversity) and Simpson's diversity
#' @return Sheldon's evenness
#' @examples
#' N=50
#' # uneven species abundance distribution
#' v = round(generateAbundances(N, mode=5))
#' sheldon(v)
#' # perfectly even species abundance distribution
#' sheldon(rep(N,N))
#' # very uneven species abundance distribution
#' sheldon(c(rep(1,N),1000))
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
