## ----setup, include=FALSE------------------------------------------------
# Global options
library(knitr)
opts_chunk$set(fig.path="figure_noise_simulations/")

## ----noise_simu, echo=TRUE, message=FALSE, warning=FALSE, fig.height = 10, fig.width = 6----
library(tuneR)

# White noise
w <- tuneR::noise(kind = c("white"))

# Brown noise is integrated white noise
# (ie. random walk)
# Use same time series length as in the other series
b <- cumsum(rnorm(length(w@left)))

# Pink noise
p <- tuneR::noise(kind = c("pink"))

# Visualize
par(mfrow=c(3,1))
plot(w,main="white noise")
plot(b,main="brown noise")
plot(p,main="pink noise")

## ----noise_simu2, echo=TRUE, message=FALSE, warning=FALSE----------------
library(pracma)
Hwhite <- hurstexp(w@left, d = 128)
Hbrown <- hurstexp(b, d = 128)
Hpink <- hurstexp(p@left, d = 128)

