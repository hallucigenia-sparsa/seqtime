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

# Brown noise is integrated white noise
# (ie. random walk)
b <- cumsum(w)
```

```
## Error in eval(expr, envir, enclos): no method for coercing this S4 class to a vector
```

```r
# Pink noise
p <- noise(kind = c("pink"))

# Visualize
par(mfrow=c(3,1))
plot(w,main="white noise")
plot(b,main="brown noise")
```

```
## Error in plot(b, main = "brown noise"): error in evaluating the argument 'x' in selecting a method for function 'plot': Error: object 'b' not found
```

```r
plot(p,main="pink noise")
```

![plot of chunk noise_simu](figure/noise_simu-1.png) 

Estimating Hurst exponent for the noises


```r
library(pracma)
Hwhite <- hurstexp(w@left, d = 128)
```

```
## Simple R/S Hurst estimation:         0.5358669 
## Corrected R over S Hurst exponent:   0.5454923 
## Empirical Hurst exponent:            0.5230614 
## Corrected empirical Hurst exponent:  0.5080329 
## Theoretical Hurst exponent:          0.5151584
```

```r
Hbrown <- hurstexp(b, d = 128)
```

```
## Error in stopifnot(is.numeric(x), is.numeric(d)): object 'b' not found
```

```r
Hpink <- hurstexp(p@left, d = 128)
```

```
## Simple R/S Hurst estimation:         0.8012564 
## Corrected R over S Hurst exponent:   0.888189 
## Empirical Hurst exponent:            0.889189 
## Corrected empirical Hurst exponent:  0.8832739 
## Theoretical Hurst exponent:          0.5151584
```



