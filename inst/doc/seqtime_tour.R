## ----setup, include=FALSE------------------------------------------------
# Global options
library(knitr)
opts_chunk$set(fig.path="figure_seqtime_tour/")

## ---- message=FALSE, warning=FALSE---------------------------------------
library(seqtime)

## ------------------------------------------------------------------------
N=50
A=generateA(N, c=0.1, d=-1)
A=modifyA(A,perc=70,strength="uniform",mode="negpercent")
# Generate a matrix using the algorithm by Klemm and Eguiluz to simulate a species network with a more realistic modular and scale-free structure. This takes a couple of minutes to complete.
#A=generateA(N, type="klemm", c=0.1)

## ---- fig.height = 6, fig.width = 6--------------------------------------
y=round(generateAbundances(N,mode=5))
names(y)=c(1:length(y))
barplot(y,main="Initial species abundances",xlab="Species",ylab="Abundance")

## ---- fig.height = 6, fig.width = 6--------------------------------------
# convert initial abundances in proportions (y/sum(y)) and run without a noise term (sigma=-1)
out.ricker=ricker(N,A=A,y=(y/sum(y)),K=rep(0.1,N), sigma=-1,tend=500)
tsplot(out.ricker,header="Ricker")

## ---- fig.height = 6, fig.width = 6--------------------------------------
# the pseudo-count allows to take the logarithm of species hat went extinct
ricker.taylor=seqtime::taylor(out.ricker, pseudo=0.0001, col="green", type="taylor")

## ---- fig.height = 6, fig.width = 6--------------------------------------
ricker.noise=identifyNoisetypes(out.ricker, abund.threshold = 0)
plotNoisetypes(ricker.noise)

## ---- fig.height = 6, fig.width = 6--------------------------------------
out.soi=soi(N, I=500, A=A, m.vector=y, tend=100)
tsplot(out.soi,header="SOI")

## ---- fig.height = 6, fig.width = 6--------------------------------------
soi.taylor=seqtime::taylor(out.soi, pseudo=0.0001, col="blue", type="taylor")

## ---- fig.height = 6, fig.width = 6--------------------------------------
soi.noise=identifyNoisetypes(out.soi,abund.threshold=0)
plotNoisetypes(soi.noise)

## ---- fig.height = 6, fig.width = 6--------------------------------------
out.hubbell=simHubbell(N=N, M=N,I=1500,d=N, m.vector=(y/sum(y)), m=0.1, tskip=500, tend=1000)
tsplot(out.hubbell,header="Hubbell")

## ---- fig.height = 6, fig.width = 6--------------------------------------
hubbell.taylor=seqtime::taylor(out.hubbell, pseudo=0.0001, col="blue", type="taylor")

## ---- fig.height = 6, fig.width = 6--------------------------------------
hubbell.noise=identifyNoisetypes(out.hubbell,abund.threshold=0)
plotNoisetypes(hubbell.noise)

## ---- fig.height = 6, fig.width = 6--------------------------------------
dm.uneven=simCountMat(N,samples=100,mode=5,k=0.05)
tsplot(dm.uneven,header="Dirichlet-Multinomial")

## ---- fig.height = 6, fig.width = 6--------------------------------------
dm.uneven.taylor=seqtime::taylor(dm.uneven, pseudo=0.0001, col="orange", type="taylor", header="Dirichlet-Multinomial")

## ---- fig.height = 6, fig.width = 6--------------------------------------
dm.uneven.noise=identifyNoisetypes(dm.uneven,abund.threshold=0)
plotNoisetypes(dm.uneven.noise)

## ---- fig.height = 6, fig.width = 6--------------------------------------
dm.even=simCountMat(N,samples=100,mode=1)
dm.even.taylor=seqtime::taylor(dm.even, pseudo=0.0001, col="orange", type="taylor", header="Even Dirichlet-Multinomial")

