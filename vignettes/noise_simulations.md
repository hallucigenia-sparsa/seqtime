---
title: "Noise simulation examples"
date: "2017-04-20"
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
```

![plot of chunk noise_simu](figure_noise_simulations/noise_simu-1.png)

Estimating Hurst exponent for the noises


```r
library(pracma)
Hwhite <- hurstexp(w@left, d = 128)
```

```
## Simple R/S Hurst estimation:         0.5055966 
## Corrected R over S Hurst exponent:   0.5251013 
## Empirical Hurst exponent:            0.5240906 
## Corrected empirical Hurst exponent:  0.5085169 
## Theoretical Hurst exponent:          0.5155387
```

```r
Hbrown <- hurstexp(b, d = 128)
```

```
## Simple R/S Hurst estimation:         0.9168978 
## Corrected R over S Hurst exponent:   0.9977773 
## Empirical Hurst exponent:            0.9948558 
## Corrected empirical Hurst exponent:  0.9908147 
## Theoretical Hurst exponent:          0.5155387
```

```r
Hpink <- hurstexp(p@left, d = 128)
```

```
## Simple R/S Hurst estimation:         0.8274734 
## Corrected R over S Hurst exponent:   0.9322969 
## Empirical Hurst exponent:            0.9669187 
## Corrected empirical Hurst exponent:  0.9607199 
## Theoretical Hurst exponent:          0.5155387
```



