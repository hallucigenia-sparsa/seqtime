---
title: "Digital History in Finland, 1488-1828"
author: "Leo Lahti"
date: "December 9, 2015"
output: markdown_document
---

Simulating noise types (following [this](http://stackoverflow.com/questions/8697567/how-to-simulate-pink-noise-in-r))


```r
require(tuneR)

# White noise
w <- noise(kind = c("white"))

# Pink noise
p <- noise(kind = c("pink"))

# Visualize
par(mfrow=c(2,1))
plot(w,main="white noise")
plot(p,main="pink noise")
```

![plot of chunk noise_simu](figure/noise_simu-1.png) 

Estimating Hurst exponent for the noises


```r
library(pracma)
Hwhite <- hurstexp(w@left, d = 128)
```

```
## Simple R/S Hurst estimation:         0.4991999 
## Corrected R over S Hurst exponent:   0.518042 
## Empirical Hurst exponent:            0.5133844 
## Corrected empirical Hurst exponent:  0.4980776 
## Theoretical Hurst exponent:          0.5151584
```

```r
Hpink <- hurstexp(p@left, d = 128)
```

```
## Simple R/S Hurst estimation:         0.862905 
## Corrected R over S Hurst exponent:   0.9364855 
## Empirical Hurst exponent:            0.9399019 
## Corrected empirical Hurst exponent:  0.9339377 
## Theoretical Hurst exponent:          0.5151584
```

