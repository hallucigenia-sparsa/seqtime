#' @importFrom FGN HurstK
#' @importFrom gtools rdirichlet
#' @importFrom tuneR noise
#' @importFrom vegan rrarefy
#' @importFrom vegan radfit
#' @importFrom vegan bstick
#' @importFrom vegan vegdist
#' @importFrom vegan capscale
#' @importFrom deSolve lsoda
#' @importFrom geigen gqz
#' @importFrom MASS ginv
#' @importFrom stinepack stinterp
#' @importFrom untb untb
#' @importFrom untb as.count
#' @importFrom untb species.table
#' @importFrom grDevices colorRampPalette
#' @importFrom grDevices dev.off
#' @importFrom grDevices hsv
#' @importFrom grDevices pdf
#' @importFrom grDevices rgb
#' @importFrom stats cor
#' @importFrom stats cor.test
#' @importFrom stats cov
#' @importFrom stats lm
#' @importFrom stats aov
#' @importFrom stats TukeyHSD
#' @importFrom stats lowess
#' @importFrom stats median
#' @importFrom stats pf
#' @importFrom stats quantile
#' @importFrom stats rbinom
#' @importFrom stats rexp
#' @importFrom stats princomp
#' @importFrom stats rgamma
#' @importFrom stats rlnorm
#' @importFrom stats complete.cases
#' @importFrom stats sd
#' @importFrom stats na.omit
#' @importFrom stats rmultinom
#' @importFrom stats rnbinom
#' @importFrom stats rnorm
#' @importFrom stats rpois
#' @importFrom stats runif
#' @importFrom stats sigma
#' @importFrom stats spline
#' @importFrom stats var
#' @importFrom stats wilcox.test
#' @importFrom stats smooth.spline
#' @importFrom stats predict
#' @importFrom utils data
#' @importFrom utils read.table
#' @import igraph
#' @import graphics

.onAttach <- function(lib, pkg)
{
  packageStartupMessage('\nseqtime: Time Series Analysis of Sequencing Data')
}




