% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster_estimates_source.R
\name{calc.R.analytic.set.tau}
\alias{calc.R.analytic.set.tau}
\title{Calculate R analytically from tau parameters, using integration.}
\usage{
calc.R.analytic.set.tau(v, VE = 1, cbeta = c(10, 15, 20), theta, lambda,
  psi, upper.lim, ...)
}
\arguments{
\item{v}{vaccination coverage.}

\item{VE}{vaccine efficacy. Default is VE=1, assuming we are looking at successful vaccination only (i.e. critical immunization threshold).}

\item{cbeta}{equivalent to \eqn{R[0]}, represents contact number per infectious period (c) x transmission probability per contact (beta).}

\item{theta}{the maximum probability ratio of one susceptible individual being in contact with a susceptible at a distance of 0 (r=0)}

\item{lambda}{rate of decay of the exponential distribution. Describes the dispersion of the \eqn{tau(r)} function}

\item{psi}{a rescaling value, allowing the distribution to average to 1 across the entire range of \code{r}.}

\item{upper.lim}{upper limit of integration.}

\item{...}{additional arguments to be passed to \code{calc.integ.fn}.}
}
\value{
\eqn{R} effective, adjusted for spatial clustering of susceptibility
}
\description{
Calculate R analytically from tau parameters, using integration.
}
\examples{

}
\author{
Shaun Truelove
}
