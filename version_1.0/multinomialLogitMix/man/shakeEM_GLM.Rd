\name{shakeEM_GLM}
\alias{shakeEM_GLM}
\title{
Shake-small EM
}
\description{
Assume that there are
at least two clusters in the fitted model. We randomly select 2 of them and propose to randomly re-allocate the assigned observations within those 2 clusters.
}
\usage{
shakeEM_GLM(y, x, K, equalModel, tsplit = 10, maxIter = 20, 
	emthreshold = 1e-08, maxNR = 5, nCores, 
	split = TRUE, R0, method)
}
\arguments{
  \item{y}{
y
}
  \item{x}{
X
}
  \item{K}{
K
}
  \item{equalModel}{
eq
}
  \item{tsplit}{
tsplit
}
  \item{maxIter}{
maxiter
}
  \item{emthreshold}{
em
}
  \item{maxNR}{
max
}
  \item{nCores}{
nc
}
  \item{split}{
spl
}
  \item{R0}{
ro
}
  \item{method}{
met
}
}
\value{
valu}
\author{
Panagiotis Papastamoulis}

