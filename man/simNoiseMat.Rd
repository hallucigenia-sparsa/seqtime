% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simNoiseMat.R
\name{simNoiseMat}
\alias{simNoiseMat}
\title{Simulate Noise}
\usage{
simNoiseMat(samples = 100, noisetypes = list(white = 2, pink = 5, brown
  = 3), brown.type = "bm", mean = 5, sd = 1, r = 2,
  stayDead = FALSE)
}
\arguments{
\item{samples}{the number of samples to be drawn}

\item{noisetypes}{a list specifying the number of rows to be generated for each noise type}

\item{brown.type}{the algorithm to generate brown noise, options: bm (Brown motion) and gbm (geometrical Brown motion; requires sde package)}

\item{mean}{parameter for white and brown noise (type bm)}

\item{sd}{parameter for white and brown noise (type bm)}

\item{r}{parameter for brown noise (type gbm)}

\item{stayDead}{once a time series hits zero or below, it stays at zero}
}
\value{
a matrix where each row represents a distribution following a specified noise type
}
\description{
Simulate a matrix composed of rows with different noise types.
}
\details{
White noise is generated with rnorm and brown noise with cumsum(rnorm).
Pink noise is generated with the noise function from the tuneR package.
Brown noise can optionally be generated with geometric Brown motion via the
GBM function in the sde package. Option stayDead allows simulating extinction
events; i.e. once a value of zero or below has been hit, the row only takes zero values.
}
\examples{
  # plot the power spectrum of pink noise
  ps <- powerspec(simNoiseMat(samples=100,noisetypes=list(white=0,pink=1,brown=0))[1,], plot=TRUE)
  # plot a matrix with mixed noise types
  mat <- simNoiseMat(samples=300,noisetypes=list(white=0,pink=30,brown=20))
  tsplot(mat,header="simulated noise")
  # Taylor law of brown noise
  t <- taylor(simNoiseMat(samples=300,noisetypes=list(white=0,pink=0,brown=50)),type="taylor")
  # Check generation of pink noise
  i <- identifyNoisetypes(simNoiseMat(samples=500,noisetypes=list(white=0,pink=10,brown=0)))
}
