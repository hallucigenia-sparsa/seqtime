---
title: "Noise simulation examples"
date: "2017-04-19"
output: 
  rmarkdown::html_vignette:
    toc: true
vignette: >
  %\VignetteIndexEntry{seqtime examples: noise simulation}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteDepends{Cairo}
  %\VignetteEncoding{UTF-8}
  \usepackage[utf8]{inputenc}
---

Simulating noise types (following [this](http://stackoverflow.com/questions/8697567/how-to-simulate-pink-noise-in-r))





```r
require(tuneR)

# White noise
w <- noise(kind = c("white"))

# Brown noise is integrated white noise
# (ie. random walk)
# Use same time series length as in the other series
b <- cumsum(rnorm(length(w@left)))

# Pink noise
p <- noise(kind = c("pink"))

# Visualize
par(mfrow=c(3,1))
plot(w,main="white noise")
plot(b,main="brown noise")
plot(p,main="pink noise")
```

![plot of chunk noise_simu](figure_noise_simulations/noise_simu-1.png)

Estimating Hurst exponent for the noises


```r
library(pracma)
Hwhite <- hurstexp(w@left, d = 128)
```

```
## Simple R/S Hurst estimation:         0.4938473 
## Corrected R over S Hurst exponent:   0.5059308 
## Empirical Hurst exponent:            0.4945896 
## Corrected empirical Hurst exponent:  0.4796675 
## Theoretical Hurst exponent:          0.5151584
```

```r
Hbrown <- hurstexp(b, d = 128)
```

```
## Simple R/S Hurst estimation:         0.9172282 
## Corrected R over S Hurst exponent:   1.007118 
## Empirical Hurst exponent:            0.9977148 
## Corrected empirical Hurst exponent:  0.9937812 
## Theoretical Hurst exponent:          0.5151584
```

```r
Hpink <- hurstexp(p@left, d = 128)
```

```
## Simple R/S Hurst estimation:         0.788127 
## Corrected R over S Hurst exponent:   0.8799208 
## Empirical Hurst exponent:            0.907866 
## Corrected empirical Hurst exponent:  0.9018946 
## Theoretical Hurst exponent:          0.5151584
```



