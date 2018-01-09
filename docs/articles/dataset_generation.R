## ----setup, include=FALSE------------------------------------------------
library(knitr)

## ----message=FALSE, warning=FALSE----------------------------------------
library(seqtime)
library(ggplot2)
library(reshape2)

## ------------------------------------------------------------------------
N = 50
S = 40
A = generateA(N, "klemm", pep=10, c =0.05)
plotA(A, header="Klemm-Eguiluz interaction matrix")

