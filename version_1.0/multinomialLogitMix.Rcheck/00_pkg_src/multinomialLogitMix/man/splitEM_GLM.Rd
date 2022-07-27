\name{splitEM_GLM}
\alias{splitEM_GLM}
\title{
Split-small EM scheme.
}
\description{
Split two randomly selected clusters based on a model with one component smaller than the current one. This procedure is repeated within a small-EM scheme. The best split is chose to initialize the model. 
}
\usage{
splitEM_GLM(y, x, K, smallerModel, tsplit = 10, maxIter = 20, 
	emthreshold = 1e-08, maxNR = 5, nCores, 
	split = TRUE, R0, method)
}

\arguments{
  \item{y}{
y
}
  \item{x}{
x
}
  \item{K}{
k
}
  \item{smallerModel}{
smla
}
  \item{tsplit}{
tsp
}
  \item{maxIter}{
max
}
  \item{emthreshold}{
thr
}
  \item{maxNR}{
maxn
}
  \item{nCores}{
nc
}
  \item{split}{
spi
}
  \item{R0}{
ro
}
  \item{method}{
meth
}
}
\value{
val
}
\references{
Papastamoulis, P., Martin-Magniette, M. L., and Maugis-Rabusseau, C. (2016). On the estimation of mixtures of Poisson regression models with large number of components. Computational Statistics & Data Analysis, 93, 97-106.
}
\author{
Panagiotis Papastamoulis
}
