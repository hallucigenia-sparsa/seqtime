% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/perturbation.R
\name{binaryToPerturb}
\alias{binaryToPerturb}
\title{Given a binary vector, build a perturbation object}
\usage{
binaryToPerturb(bin.metadata.item, false.value = "no",
  true.value = "yes")
}
\arguments{
\item{bin.metadata.item}{a binary metadata item with false and true values}

\item{false.value}{the string/number used to encode false values}

\item{true.value}{the string/number used to encode true values}
}
\description{
True.values are counted as perturbed time points,
false.values as non-perturbed ones.
}
