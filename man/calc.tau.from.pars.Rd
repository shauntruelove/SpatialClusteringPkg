% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster_estimates_source.R
\name{calc.tau.from.pars}
\alias{calc.tau.from.pars}
\title{Calculate \eqn{tau(r)} for specific \code{r} distances, given known \eqn{tau(r)} function parameters}
\usage{
calc.tau.from.pars(r, theta, lambda = 0.5, psi = 1)
}
\arguments{
\item{r}{Spatial distance from an indivual.}

\item{theta}{the maximum probability ratio of one susceptible individual being in contact with a susceptible at a distance of 0 (r=0).}

\item{lambda}{rate of decay of the exponential distribution. Describes the dispersion of the \eqn{tau(r)} function.}

\item{psi}{a rescaling value, allowing the distribution to average to 1 across the entire range of \code{r}.}
}
\value{
\eqn{tau(r)} for specific \code{r} distances.
}
\description{
Calculate \eqn{tau(r)} for specific \code{r} distances, given known \eqn{tau(r)} function parameters
}
\examples{

}
\author{
Shaun Truelove
}
