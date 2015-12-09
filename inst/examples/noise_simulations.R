# Simulate noise types
# http://stackoverflow.com/questions/8697567/how-to-simulate-pink-noise-in-r
require(tuneR)
w <- noise(kind = c("white"))
p <- noise(kind = c("pink"))
par(mfrow=c(2,1))
plot(w,main="white noise")
plot(p,main="pink noise")

# Hurst exponent of white noise ?
library(pracma)
Hwhite <- hurstexp(w@left, d = 128)
Hpink <- hurstexp(p@left, d = 128)



