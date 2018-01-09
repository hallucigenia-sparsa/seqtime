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

## ------------------------------------------------------------------------
dataset = generateDataSet(S, A)
dataset = normalize(dataset)
dataset = melt(dataset)
colnames(dataset) = c("Species", "Sample", "Abundance")
ggplot(data=dataset, aes(x=dataset$Sample, y=dataset$Abundance, width=1)) + geom_bar(aes(y = dataset$Abundance, x= dataset$Sample, fill=dataset$Species), data=dataset, stat="identity", show.legend=F) + theme(aspect.ratio=.4) + theme_classic()+ ylab("Relative abundance") + xlab("Sample")

## ------------------------------------------------------------------------
env = envGrowthChanges(N, strength=0.8)
dataset = generateDataSet(S, A, env.matrix=env, perturb.count=c(20,20))
dataset = normalize(dataset)
dataset = melt(dataset)
colnames(dataset) = c("Species", "Sample", "Abundance")
ggplot(data=dataset, aes(x=dataset$Sample, y=dataset$Abundance, width=1)) + geom_bar(aes(y = dataset$Abundance, x= dataset$Sample, fill=dataset$Species), data=dataset, stat="identity", show.legend=F) + theme(aspect.ratio=.4) + theme_classic()+ ylab("Relative abundance") + xlab("Sample")

