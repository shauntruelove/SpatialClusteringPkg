
#require(compiler)



#'
#' Calculate Traditional, Homogeneous R Effective from VE, v, and R0
#'
#' @param Bc Bc = R0 for the particular agent.
#' @param VE vaccine efficacy. Default is VE=1, assuming we are looking at successful vaccination only (i.e. critical immunization threshold).
#' @param v vaccination coverage. vaccination coverage.
#'
#' @return traditional, homogeneous \eqn{R} effective estimate.
#'
#' @author Shaun Truelove
#'
#' @export
#'
#' @examples
#'
calc.R <- function(Bc, VE, v){
    res1 <- 1 - (v %*% t(VE))
    colnames(res1) <- VE
    row.names(res1) <- v
    res <- lapply(X=Bc, FUN=function(X, res1){X*res1}, res1=res1)
    return(res)
}



# Tau Distribution --------------------------------------------------------

#'
#' Calculate exponential \eqn{tau(r)} for specific value of \eqn{r}, assuming knowledge of maximum contact probability ratio
#'
#' @param r Spatial distance from an indivual.
#' @param tau.level Maximum relative probability of another susceptible individual, occurring at distance r=0.
#' For example if the overall probability of an individual being susceptible is 0.10 (\code{v=90%}), if \code{tau.level=1.5},
#' the average maximum probability is 0.15, due to clustering.
#' @param lambda rate of decay of the exponential distribution. Describes the dispersion of the \eqn{tau(r)} function.
#' @param max.x maximum distance value for which to calculate the tau vector.
#'
#' @return \eqn{tau(r)} for specific \code{r} distances.
#'
#' @author Shaun Truelove
#'
#' @export
#'
#' @examples
#'
calc.tau <- compiler::cmpfun(function(r, tau.level=1.5, lambda=.5, max.x=50){

    theta <- tau.level-1
    xtmp <- seq(0, max.x, by=.0001)
    mean_y1tmp <- mean(theta*exp(-lambda*xtmp))
    b <- 1 - mean_y1tmp
    theta1 <- (theta + mean_y1tmp)

    return(theta1 * exp(-lambda*r) + b)
})


#'
#' Calculate \eqn{tau(r)} for specific \code{r} distances, given known \eqn{tau(r)} function parameters
#'
#' @param r Spatial distance from an indivual.
#' @param theta the maximum probability ratio of one susceptible individual being in contact with a susceptible at a distance of 0 (r=0).
#' @param lambda rate of decay of the exponential distribution. Describes the dispersion of the \eqn{tau(r)} function.
#' @param psi a rescaling value, allowing the distribution to average to 1 across the entire range of \code{r}.
#'
#' @return \eqn{tau(r)} for specific \code{r} distances.
#'
#' @author Shaun Truelove
#'
#' @export
#'
#' @examples
#'
calc.tau.from.pars <- compiler::cmpfun(function(r, theta, lambda=.5, psi=1){
    return(theta * exp(-lambda*r) + psi)
})



#'
#' Get \eqn{tau(r)} exponential function parameters for use in other functions
#'
#' @param tau.level Maximum relative probability of another susceptible individual, occurring at distance r=0.
#' For example if the overall probability of an individual being susceptible is 0.10 (v=90%),
#' if tau.level=1.5, the average maximum probability is 0.15, due to clustering.
#' @param lambda rate of decay of the exponential distribution. Describes the dispersion of the \eqn{tau(r)} function.
#' @param max.x maximum distance value for which to calculate the tau vector.
#' @param r.delim \code{r} distance window size (in km) for discrete distance ranges.
#'
#' @return list of \eqn{tau(r)} function parameters:
#'   \item{theta}{maximum probability ratio of one susceptible individual being in contact with a susceptible at a distance of 0 (r=0).}
#'   \item{lambda}{rate of decay of the exponential distribution. Describes the dispersion of the \eqn{tau(r)} function.}
#'   \item{psi}{a rescaling value, allowing the distribution to average to 1 across the entire range of \code{r}.}
#'
#' @author Shaun Truelove
#'
#' @export
#'
#' @examples
#'
get.tau.params <- compiler::cmpfun(function(tau.level=1.5, lambda=.5, max.x=50, r.delim=1){

    theta <- tau.level-1
    xtmp <- seq(0, max.x, by=r.delim)
    mean_y1tmp <- mean(theta*exp(-lambda*xtmp))
    psi <- 1 - mean_y1tmp
    theta1 <- (theta + mean_y1tmp)

    return(list(theta=theta1, lambda=lambda, psi=psi))
})


# get.tau.params <- function(upper.lim, rdelin=1, theta, lambda=0.5){
#     x <- seq(0, upper.lim, by=rdelin)
#     ytmp <- theta*exp(-lambda*x)
#     b <- 1 - mean(ytmp)
#     theta1 <- (theta + mean(ytmp))
#     return(list(theta1, b))
# }



#'
#' Calculate tau(r) for r midpoints
#'
#' @param tau.level Maximum relative probability of another susceptible individual, occurring at distance r=0.
#' For example if the overall probability of an individual being susceptible is 0.10 (v=90%),
#' if tau.level=1.5, the average maximum probability is 0.15, due to clustering.
#' @param lambda rate of decay of the exponential distribution. Describes the dispersion of the \eqn{tau(r)} function
#' @param max.x maximum distance value for which to calculate the tau vector.
#' @param r.delim \code{r} distance window size (in km) for discrete distance ranges.
#'
#' @return vector of \eqn{tau(r)} for discrete \code{r} distance window ranges.
#'
#' @author Shaun Truelove
#'
#' @export
#'
#' @examples
#'
calc.tau.distrib <- compiler::cmpfun(function(tau.level=1.5, lambda=0.5, max.x=50, r.delim=0.5){

    r.max=seq(r.delim, max.x, r.delim);
    r.min=r.max-r.max[1];
    r.mid=r.min+(r.max-r.min)/2

    return(calc.tau(r.mid, tau.level=tau.level, lambda=lambda, max.x=max.x))
})


#'
#' Function to produce analytic \eqn{tau(r)} data for plots
#'
#' @param upper.lim upper limit of integration.
#' @param r.delim \code{r} distance window size (in km) for discrete distance ranges.
#' @param theta the maximum probability ratio of one susceptible individual being in contact with a susceptible at a distance of 0 (r=0)
#' @param lambda rate of decay of the exponential distribution. Describes the dispersion of the \eqn{tau(r)} function
#'
#' @return vector of \eqn{tau(r)} for discrete \code{r} distance window ranges.
#'
#' @author Shaun Truelove
#'
#' @export
#'
#' @examples
#'
calc.tau.analytic <- function(theta, lambda=0.5, max.x=50, r.delim=1){

    x <- seq(0, max.x, by=r.delim)
    ytmp <- theta*exp(-lambda*x)
    psi <- 1 - mean(ytmp)
    theta1 <- (theta + mean(ytmp))

    tau <- theta1 * exp(-lambda*x) + psi
    return(tau)
}


#'
#' Integrand function for calculating \eqn{phi} using integration
#'
#' @param r Spatial distance from an indivual
#' @param gx.funct.name a string; function for representing \eqn{g(r)}, the contact distance probability distribution.
#' The gamma distribution function (i.e. \code{dgamma()}, input as \code{gx.funct.name='dgamma'}. Custom functions, such as splines can also be inputted.)
#' @param theta the maximum probability ratio of one susceptible individual being in contact with a susceptible at a distance of 0 (r=0)
#' @param lambda rate of decay of the exponential distribution. Describes the dispersion of the \eqn{tau(r)} function
#' @param psi a rescaling value, allowing the distribution to average to 1 across the entire range of \code{r}.
#' @param ... additional arguments to be passed to function described by \code{gx.funct.name}.
#'
#' @return vector of \eqn{tau(r) * g(r)} at all values of \code{r}.
#'
#' @author Shaun Truelove
#'
#' @examples
#'
calc.integ.fn <- compiler::cmpfun(function(r, gx.funct.name='dgamma', theta=theta, lambda=lambda, psi=psi, ...){
    r.tmp <- r
    if (length(r.tmp)==0)  r.tmp <- 0

    gx.funct <- get(gx.funct.name)
    gx <- gx.funct(r.tmp,...)

    tau <- calc.tau.from.pars(r=r, theta=theta, lambda=lambda, psi=psi)
    if (length(tau)==0)  tau <- 0

    return(gx*tau)
})




# Phi - Clustering Adjustment Factor -------------------------


#' Analytic Calculation of Phi (closed form integration of tau(x)*g(x))
#'
#' tau(r): theta*exp(-lambda*r)+psi
#' g(r): f(r, alpha, beta) ~ gamma(shape=alpha, rate=beta)
#'
#'
#' @param theta the maximum probability ratio of one susceptible individual being in contact with a susceptible at a distance of 0 (r=0).
#' @param lambda rate of decay of the exponential distribution. Describes the dispersion of the \eqn{tau(r)} function.
#' @param psi a rescaling value, allowing the distribution to average to 1 across the entire range of \code{r}.
#' @param alpha shape parameter from gamma distribution.
#' @param beta rate parameter from gamma distribution. Scale parameter = 1/rate.
#'
#' @return \eqn{phi}, the clustering adjustement factor.
#'
#' @author Shaun Truelove
#'
#' @export
#'
#' @examples
#'
calc.phi.cf <- compiler::cmpfun(function(theta, lambda, psi, alpha, beta){
    return(psi + theta*(beta^alpha / (beta+lambda)^alpha))
})


#'
#' Analytic Calculation of Phi (integral of \eqn{tau(x)*g(x)})
#'
#' @param upper.lim upper limit of integration.
#' @param lower.lim lower limit of integration.
#' @param theta the maximum probability ratio of one susceptible individual being in contact with a susceptible at a distance of 0 (r=0)
#' @param lambda rate of decay of the exponential distribution. Describes the dispersion of the \eqn{tau(r)} function
#' @param subdivisions the maximum number of subintervals in the integration.
#' @param by.x r discretization value.
#' @param ... additional arguments to be passed to \code{calc.integ.fn}.
#'
#' @return \eqn{phi}, the clustering adjustement factor.
#'
#' @author Shaun Truelove
#'
#' @export
#'
#' @examples
#'
calc.phi.analytic <- compiler::cmpfun(function(theta, lambda, upper.lim, lower.lim=0.00001, subdivisions=1000L, by.x=0.001, ...){

    xtmp <- seq(0, upper.lim, by=by.x)
    y1tmp <- theta*exp(-lambda*xtmp)
    psi <- 1 - mean(y1tmp)
    theta1 <- (theta + mean(y1tmp))

    phi <- stats::integrate(calc.integ.fn, lower=lower.lim, upper=upper.lim,
                     theta=theta1, lambda=lambda, psi=psi,
                     subdivisions=subdivisions, stop.on.error=FALSE, ...)
    return(phi$value)
})



#'
#' Calculate \eqn{phi} analytically using parametric functions for \eqn{tau(r)} and \eqn{g(r)}
#'
#' @param upper.lim upper limit of integration.
#' @param lower.lim lower limit of integration.
#' @param theta the maximum probability ratio of one susceptible individual being in contact with a susceptible at a distance of 0 (r=0)
#' @param lambda rate of decay of the exponential distribution. Describes the dispersion of the \eqn{tau(r)} function
#' @param psi a rescaling value, allowing the distribution to average to 1 across the entire range of \code{r}.
#' @param subdivisions the maximum number of subintervals in the integration.
#' @param ... additional arguments to be passed to \code{calc.integ.fn}.
#'
#' @return \eqn{phi}, the clustering adjustement factor.
#'
#' @author Shaun Truelove
#'
#' @export
#'
#' @examples
#'
calc.phi.analytic.set.tau <- compiler::cmpfun(function(theta, lambda, psi, upper.lim, lower.lim=0.000001, subdivisions=1000L, ...){
    phi <- stats::integrate(calc.integ.fn, lower=lower.lim, upper=upper.lim,
                     theta=theta, lambda=lambda, psi=psi,
                     subdivisions=subdivisions, stop.on.error=FALSE, ...)
    return(phi$value)
})




#---------- Critical Vaccination Threshold ----------------------------------------

### Function using Closed Form - Analytic
#'
#' Calculate Vc analytically using the closed form solution
#'
#' THis function uses the closed-form solution of the integral of the exponential \code{tau(r)} function and gamma \code{g(r)} function.
#'
#' @param VE vaccine efficacy. Default is VE=1, assuming we are looking at successful vaccination only (i.e. critical immunization threshold).
#' @param R0 basic reproductive number. Default is \eqn{R[0]=15}, for measles.
#' @param theta the maximum probability ratio of one susceptible individual being in contact with a susceptible at a distance of 0 (r=0)
#' @param lambda rate of decay of the exponential distribution. Describes the dispersion of the \eqn{tau(r)} function
#' @param psi a rescaling value, allowing the distribution to average to 1 across the entire range of \code{r}.
#' @param alpha shape parameter from gamma distribution.
#' @param beta rate parameter from gamma distribution. Scale parameter = 1/rate.
#'
#' @return \eqn{V[c]}, adjusted for spatial clustering of susceptibility
#'
#' @author Shaun Truelove
#'
#' @export
#'
#' @examples
#'
calc.Vc.cf <- compiler::cmpfun(function(VE=1, R0=15, theta, lambda, psi, alpha, beta){
    phi <- calc.phi.cf(theta, lambda, psi, alpha, beta)
    res <- sapply(R0, function(x) {(1/VE) * (1 - (1 / (x*phi))) } )
    return(res)
})


#'
#' Calculate Vc analytically, using integration
#'
#' @param VE vaccine efficacy. Default is VE=1, assuming we are looking at successful vaccination only (i.e. critical immunization threshold).
#' @param R0 basic reproductive number. Default is \eqn{R[0]=15}, for measles.
#' @param upper.lim upper limit of integration.
#' @param lower.lim lower limit of integration.
#' @param theta the maximum probability ratio of one susceptible individual being in contact with a susceptible at a distance of 0 (r=0)
#' @param lambda rate of decay of the exponential distribution. Describes the dispersion of the \eqn{tau(r)} function
#' @param subdivisions the maximum number of subintervals in the integration.
#' @param r.delim r discretization value.
#' @param ... additional arguments to be passed to \code{calc.integ.fn}.
#'
#' @return \eqn{V[c]}, adjusted for spatial clustering of susceptibility
#'
#' @author Shaun Truelove
#'
#' @export
#'
#' @examples
#'
calc.Vc.analytic <- compiler::cmpfun(function(VE=1, R0=15, theta, lambda, upper.lim, lower.lim=0.00001, subdivisions=1000L, r.delim=0.001,...){

    xtmp <- seq(0, upper.lim, by=r.delim)
    y1tmp <- theta*exp(-lambda*xtmp)
    b <- 1 - mean(y1tmp)
    theta1 <- (theta + mean(y1tmp))

    phi <- stats::integrate(calc.integ.fn, lower=lower.lim, upper=upper.lim,
                     theta=theta1, lambda=lambda, b=b,
                     subdivisions=subdivisions, stop.on.error=FALSE, ...)
    res <- sapply(R0, function(x) {(1/VE) * (1 - (1 / (x*phi$value))) } )
    return(res)
})


#'
#' Calculate Vc analytically, using integration and previously refined \code{tau(r)} parameters
#'
#' @param VE vaccine efficacy. Default is VE=1, assuming we are looking at successful vaccination only (i.e. critical immunization threshold).
#' @param R0 basic reproductive number. Default is \eqn{R[0]=15}, for measles.
#' @param upper.lim upper limit of integration.
#' @param theta the maximum probability ratio of one susceptible individual being in contact with a susceptible at a distance of 0 (r=0).
#' @param lambda rate of decay of the exponential distribution. Describes the dispersion of the \eqn{tau(r)} function
#' @param psi a rescaling value, allowing the distribution to average to 1 across the entire range of \code{r}.
#' @param ... additional arguments to be passed to \code{calc.integ.fn}.
#'
#' @return \eqn{V[c]}, adjusted for spatial clustering of susceptibility
#'
#' @author Shaun Truelove
#'
#' @export
#'
#' @examples
#'
calc.Vc.analytic.params <- compiler::cmpfun(function(VE=1, R0=15, theta, lambda, psi, upper.lim,...){

    phi <- integrate(calc.integ.fn, lower.lim=0, upper=upper.lim,
                     theta=theta, lambda=lambda, psi=psi,
                     subdivisions=1000L, stop.on.error=FALSE, ...)
    res <- sapply(R0, function(x) {(1/VE) * (1 - (1 / (x*phi$value))) } )
    return(res)
})


#' Calculate Vc analytically, from phi
#'
#' @param VE vaccine efficacy. Default is VE=1, assuming we are looking at successful vaccination only (i.e. critical immunization threshold).
#' @param R0 basic reproductive number. Default is \eqn{R[0]=15}, for measles.
#' @param phi Clustering adjustment factor, \eqn{phi}. This is the relative probability of contact between susceptibles resulting from heterogeneity, compared to a homogenious population. Clustering adjustment factor, \eqn{phi}. This is the relative probability of contact between susceptibles resulting from heterogeneity, compared to a homogenious population.
#'
#' @return \eqn{V[c]}, adjusted for spatial clustering of susceptibility
#'
#' @author Shaun Truelove
#'
#' @export
#'
#' @examples
calc.Vc.analytic.phi <- compiler::cmpfun(function(VE=1, R0=15, phi){
    res <- (1/VE) * (1 - (1 / (R0*phi)))
    return(res)
})




#'
#' Calculate R analytically with the closed-form equation, assuming an exponential form of \code{tau(r)} and a gamma distribution for \code{g(r)}
#'
#' @param v vaccination coverage.
#' @param VE vaccine efficacy. Default is VE=1, assuming we are looking at successful vaccination only (i.e. critical immunization threshold).
#' @param cbeta equivalent to \eqn{R[0]}, represents contact number per infectious period (c) x transmission probability per contact (beta).
#' @param upper.lim upper limit of integration.
#' @param theta the maximum probability ratio of one susceptible individual being in contact with a susceptible at a distance of 0 (r=0)
#' @param lambda rate of decay of the exponential distribution. Describes the dispersion of the \eqn{tau(r)} function
#' @param psi a rescaling value, allowing the distribution to average to 1 across the entire range of \code{r}.
#' @param alpha shape parameter from gamma distribution.
#' @param beta rate parameter from gamma distribution. Scale parameter = 1/rate.
#'
#' @return \eqn{R} effective, adjusted for spatial clustering of susceptibility
#'
#' @author Shaun Truelove
#'
#' @export
#'
#' @examples
#'
calc.R.cf <- compiler::cmpfun(function( v, VE=1, cbeta=c(10, 15, 20), theta, lambda, psi, alpha, beta){
    phi <- calc.phi.cf(theta, lambda, psi, alpha, beta)
    res <- sapply(cbeta, function(x) phi*x*(1-(v*VE)))
    return(res)
})


#'
#' Calculate R analytically from known phi
#'
#' @param phi Clustering adjustment factor, \eqn{phi}. This is the relative probability of contact between susceptibles resulting from heterogeneity, compared to a homogenious population.
#' @param v vaccination coverage.
#' @param VE vaccine efficacy. Default is VE=1, assuming we are looking at successful vaccination only (i.e. critical immunization threshold).
#' @param cbeta equivalent to \eqn{R[0]}, represents contact number per infectious period (c) x transmission probability per contact (beta).
#'
#' @return \eqn{R} effective, adjusted for spatial clustering of susceptibility
#'
#' @author Shaun Truelove
#'
#' @examples
calc.R.analytic.phi <- function(phi, v, VE=1, cbeta=c(10, 15, 20)){
    res <- sapply(cbeta, function(x) phi*x*(1-(v*VE)))
    return(res)
}


#'
#' Calculate R analytically from tau parameters, using integration.
#'
#' @param v vaccination coverage.
#' @param VE vaccine efficacy. Default is VE=1, assuming we are looking at successful vaccination only (i.e. critical immunization threshold).
#' @param cbeta equivalent to \eqn{R[0]}, represents contact number per infectious period (c) x transmission probability per contact (beta).
#' @param upper.lim upper limit of integration.
#' @param theta the maximum probability ratio of one susceptible individual being in contact with a susceptible at a distance of 0 (r=0)
#' @param lambda rate of decay of the exponential distribution. Describes the dispersion of the \eqn{tau(r)} function
#' @param psi a rescaling value, allowing the distribution to average to 1 across the entire range of \code{r}.
#' @param ... additional arguments to be passed to \code{calc.integ.fn}.
#'
#' @return \eqn{R} effective, adjusted for spatial clustering of susceptibility
#'
#' @author Shaun Truelove
#'
#' @export
#'
#' @examples
#'
calc.R.analytic.set.tau <- compiler::cmpfun(function( v, VE=1, cbeta=c(10, 15, 20), theta, lambda, psi, upper.lim, ...){
    phi <- stats::integrate(calc.integ.fn, lower=0, upper=upper.lim,
                     theta=theta, lambda=lambda, psi=psi,
                     subdivisions=500L, stop.on.error=FALSE, ...)
    res <- sapply(cbeta, function(x) phi$value*x*(1-(v*VE)))
    return(res)
})


#'
#' Calculate R analytically, using integration.
#'
#' @param v vaccination coverage.
#' @param VE vaccine efficacy. Default is VE=1, assuming we are looking at successful vaccination only (i.e. critical immunization threshold).
#' @param cbeta equivalent to \eqn{R[0]}, represents contact number per infectious period (c) x transmission probability per contact (beta).
#' @param upper.lim upper limit of integration.
#' @param theta the maximum probability ratio of one susceptible individual being in contact with a susceptible at a distance of 0 (r=0)
#' @param lambda rate of decay of the exponential distribution. Describes the dispersion of the \eqn{tau(r)} function
#' @param subdivisions the maximum number of subintervals in the integration.
#' @param r.delim r discretization value.
#' @param ... additional arguments to be passed to \code{calc.integ.fn}.
#'
#' @return \eqn{R} effective, adjusted for spatial clustering of susceptibility
#'
#' @author Shaun Truelove
#'
#' @export
#'
#' @examples
#'
calc.R.analytic <- compiler::cmpfun(function( v, VE=1, cbeta=c(10, 15, 20), theta, lambda, upper.lim, subdivisions=500L, r.delim=.001, ...){

    xtmp <- seq(0, upper.lim, by=r.delim)
    y1tmp <- theta*exp(-lambda*xtmp)
    psi <- 1 - mean(y1tmp)
    theta1 <- (theta + mean(y1tmp))

    phi <- stats::integrate(calc.integ.fn, lower=0, upper=upper.lim,
                     theta=theta1, lambda=lambda, psi=psi,
                     subdivisions=subdivisions, stop.on.error=FALSE, ...)
    res <- sapply(cbeta, function(x) phi$value*x*(1-(v*VE)))
    return(res)
})




# Functions Using Splines -------------------------------------------------



#'
#' Integrand function to feed to the R effective and Phi functions
#'
#' @param r Spatial distance from an indivual
#' @param tau.fun.spline spline function defined for \eqn{tau(r)} function.
#' @param gx.funct.name a string; function for representing \eqn{g(r)}, the contact distance probability distribution.
#' The gamma distribution function (i.e. \code{dgamma()}, input as \code{gx.funct.name='dgamma'}. Custom functions, such as splines can also be inputted.).
#' @param ... additional arguments to be passed to \code{gx.funct}.
#'
#' @return vector of \eqn{tau(r) * g(r)} at all values of \code{r}.
#'
#' @author Shaun Truelove
#'
#' @export
#'
#' @examples
#'
calc.integ.fn.spline <- function(r, tau.fun.spline, gx.funct.name='dgamma', ...){
    r.tmp <- r
    if (length(r.tmp)==0)  r.tmp <- 0
    # Get g(x) function and calculate the value for r
    gx.funct <- get(gx.funct.name)
    gx <- gx.funct(r.tmp, ...)
    # calculate the value for tau, given the spline
    tau <- predict(tau.fun.spline, r.tmp)$y
    if (length(tau)==0)  tau <- 0
    return(gx*tau)
}


#'
#' Calculate \eqn{phi} using spline function for \eqn{tau(r)}
#'
#' @param tau.fun.spline spline function defined for \eqn{tau(r)} function.
#' @param gx.funct.name a string; function for representing \eqn{g(r)}, the contact distance probability distribution.
#' The gamma distribution function (i.e. \code{dgamma()}, input as \code{gx.funct.name='dgamma'}. Custom functions, such as splines can also be inputted.).
#' @param upper.lim upper limit of integration.
#' @param lower.lim lower limit of integration.
#' @param subdivisions the maximum number of subintervals in the integration.
#' @param ... additional arguments to be passed to \code{calc.intTeg.fn.spline}.
#'
#' @return \eqn{phi}, the clustering adjustement factor.
#'
#' @author Shaun Truelove
#'
#' @export
#'
#' @examples
#'
calc.phi.spline <- compiler::cmpfun(function(tau.fun.spline, gx.funct.name='dgamma', upper.lim, lower.lim=0.000001, subdivisions=1000L, ...){
    phi <- stats::integrate(calc.intTeg.fn.spline, lower=lower.lim, upper=upper.lim, tau.fun.spline, gx.funct.name,
                     subdivisions=subdivisions, stop.on.error=FALSE, ...)
    return(phi$value)
})


#'
#' Calculate \eqn{phi} using spline function for \eqn{tau(r)}
#'
#' @param v vaccination coverage.
#' @param VE vaccine efficacy. Default is VE=1, assuming we are looking at successful vaccination only (i.e. critical immunization threshold).
#' @param cbeta equivalent to \eqn{R[0]}, represents contact number per infectious period (c) x transmission probability per contact (beta).
#' @param tau.fun.spline spline function defined for \eqn{tau(r)} function.
#' @param gx.funct.name a string; function for representing \eqn{g(r)}, the contact distance probability distribution.
#' The gamma distribution function (i.e. \code{dgamma()}, input as \code{gx.funct.name='dgamma'}. Custom functions, such as splines can also be inputted.).
#' @param upper.lim upper limit of integration.
#' @param lower.lim lower limit of integration.
#' @param subdivisions the maximum number of subintervals in the integration.
#' @param ... additional arguments to be passed to \code{calc.integ.fn.spline}.
#'
#' @return \eqn{R} effective, adjusted for spatial clustering of susceptibility
#'
#' @author Shaun Truelove
#'
#' @export
#'
#' @examples
#'
calc.R.spline <- function(v, VE=1, cbeta=c(12, 15, 18), tau.fun.spline, gx.funct.name='dgamma', upper.lim, lower.lim=0.000001, subdivisions=1000L, ...){
    phi <- stats::integrate(calc.integ.fn.spline, lower=lower.lim, upper=upper.lim, tau.fun.spline, gx.funct.name,
                     subdivisions=subdivisions, stop.on.error=FALSE, ...)
    res <- sapply(cbeta, function(x) phi$value*x*(1-(v*VE)))
    return(res)
}




# BASIC REPRODUCTIVE NUMBER -----------------------------------------------

# These are functions for back-calculation of \eqn{R[0]}. These are useful for diagnostics, sensitivity analyses, and graphing.


#'
#' Calculate \eqn{R[0]} from \eqn{V[c]}
#'
#' Back-calculate the \eqn{R[0]} from known \eqn{V[c]} and \eqn{phi}
#'
#' @param VE vaccine efficacy. Default is VE=1, assuming we are looking at successful vaccination only.
#' @param Vc Critical vaccination threshold (\eqn{V[c]}), a.k.a critical immunity threshold (\eqn{I[c]}) when \eqn{VE=1}.
#' @param phi Clustering adjustment factor, \eqn{phi}. This is the relative probability of contact between susceptibles resulting from heterogeneity, compared to a homogenious population. Clustering adjustment factor, \eqn{phi}. This is the relative probability of contact between susceptibles resulting from heterogeneity, compared to a homogenious population.
#'
#' @return \eqn{R[0]}
#'
#' @author Shaun Truelove
#'
#' @export
#'
#' @examples
#'
calc.R0.from.Vc <- function(VE=1, Vc=.95, phi){
    res <- (-phi * ((Vc*VE)-1))^-1
    return(res)
}


#'
#' Calculate \eqn{R[0]} from \eqn{V[c]}
#'
#' Analytic version using integration to back-calculate \eqn{R[0]} from \eqn{V[c]} for diagnostic and graphical purposes.
#' Here the \eqn{g(r)} function must be specified along with its parameters (i.e. a gamma distribution is not assumed).
#' An exponential distribution for \eqn{tau(r)} is assumed.
#'
#' @param VE vaccine efficacy. Default is VE=1, assuming we are looking at successful vaccination only (i.e. critical immunization threshold).
#' @param v vaccination coverage.
#' @param upper.lim upper limit of integration.
#' @param theta the maximum probability ratio of one susceptible individual being in contact with a susceptible at a distance of 0 (r=0).
#' @param lambda rate of decay of the exponential distribution. Describes the dispersion of the \eqn{tau(r)} function.
#' @param ... additional arguments to be passed to \code{calc.integ.fn}.
#'
#' @return \eqn{R[0]}
#'
#' @author Shaun Truelove
#'
#' @export
#'
#' @examples
#'
calc.R0.from.Vc.analytic <- compiler::cmpfun(function(VE=1, v=.95, theta, lambda, upper.lim, ...){

    xtmp <- seq(0, upper.lim, by=.001)
    y1tmp <- theta*exp(-lambda*xtmp)
    psi <- 1 - mean(y1tmp)
    theta1 <- (theta + mean(y1tmp))

    phi <- stats::integrate(calc.integ.fn, lower=0, upper=upper.lim,
                     theta=theta1, lambda=lambda, psi=psi,
                     subdivisions=500L, stop.on.error=FALSE, ...)

    res <- sapply(v, function(x) { -((x-1) * VE * phi$value)^-1 } ) ## x=v
    return(res)
})




