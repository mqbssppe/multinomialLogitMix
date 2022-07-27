\name{simulate_multinomial_data}
\alias{simulate_multinomial_data}
\title{
Synthetic data generator
}
\description{
This function simulates data from mixture of multinomial logistic regression models.
}
\usage{
simulate_multinomial_data(K, p, D, n, size = 20, prob = 0.025, betaTrue = NULL)
}
\arguments{
  \item{K}{
Number of clusters.
}
  \item{p}{
Number of covariates, including constant.
}
\item{D}{Number of multinomial categories.}
\item{n}{Number of data points to simulate.}
\item{size}{Negative Binomial parameter (number of successes). Default: 20.}
\item{prob}{Negative Binomial parameter (probability of success). Default: 0.025.}
\item{betaTrue}{An array which contains the true values of the logit coefficients per cluster. Default: randomly generated values.}
}
\value{
\item{count_data}{matrix of data counts.}
\item{design_matrix}{design matrix.}
\item{clustering}{Ground-truth partition of the data.}
}
\author{
Panagiotis Papastamoulis
}


