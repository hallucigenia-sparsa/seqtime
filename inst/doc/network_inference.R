## ----setup, include=FALSE------------------------------------------------
# Global options
library(knitr)
opts_chunk$set(fig.path="figure_network_inference/")

## ---- message=FALSE, warning=FALSE---------------------------------------
library(seqtime)

## ------------------------------------------------------------------------
N=20
A=generateA(N, c=0.1)
rownames(A)=c(1:N)
colnames(A)=rownames(A)

## ---- fig.height = 6, fig.width = 6--------------------------------------
plotA(A, header="Known interaction matrix")

## ---- fig.height = 6, fig.width = 6--------------------------------------
network=plotA(A,method="network")
# use igraph's function tkplot for manual layout of the network

## ---- fig.height = 6, fig.width = 6--------------------------------------
out.ricker=ricker(N,A=A)
tsplot(out.ricker,header="Ricker")

## ------------------------------------------------------------------------
Aest=limits(out.ricker)$Aest

## ---- fig.height = 6, fig.width = 6--------------------------------------
par(mfrow=c(1,2))
plotA(A,header="known")
plotA(Aest,header="inferred")
par(mfrow=c(1,1))

## ------------------------------------------------------------------------
crossCor=cor(A,Aest)
mean(diag(crossCor), na.rm=TRUE)

## ---- fig.height = 6, fig.width = 6--------------------------------------
limitsqual=limitsQuality(out.ricker,A=Aest,plot=TRUE)

## ---- fig.height = 6, fig.width = 6--------------------------------------
out.hubbell=simHubbell(N=N, M=N,I=1500,d=N, m=0.1, tskip=500, tend=1000)
tsplot(out.hubbell,header="Hubbell")

## ------------------------------------------------------------------------
Aesth=limits(out.hubbell)$Aest

## ---- fig.height = 6, fig.width = 6--------------------------------------
limitsqualh=limitsQuality(out.hubbell,A=Aesth, plot=TRUE)

