\name{GPtseries}
\alias{GPtseries-package}
\alias{GPtseries}
\docType{package}
\title{
multiple step ahead time series prediction using GP prior.
}
\usage{
MGP(dat,N,Q,num_gen,lag,max,min,Pc,Pm,x_test,M)
}
\arguments{
\item{dat}{[input dat for training]}

\item{N}{[number of training data used in the algorithm]}

\item{Q}{[number of samples per generation in GA]}

\item{num_gen}{[number of generations in GA]}

\item{lag}{[number of time steps used for each training data]}

\item{max}{[vector of maximum value to be searched in GA]}

\item{min}{[vector of minimum value to be searched in GA]}

\item{Pc}{[probability of a crossover in GA]}

\item{Pm}{[probability of a mutation in GA]}

\item{x_test}{[test input data]}

\item{M}{[M-step ahead prediction]}
}
\value{
[M-step ahead prediction for the test input data x_test]
}
\description{
Inference is done using the genetic algorithm.
}
\details{
\tabular{ll}{
Package: \tab GPtseries\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-03-10\cr
}
}
\author{
Xiaoyu Lu

Maintainer: Xiaoyu Lu <xiaoyu.lu@new.ox.ac.uk>
}

\keyword{ package }

\examples{
res = MGP(dat,N,Q,num_gen,lag,max,min,Pc,Pm,x_test,M)
}
