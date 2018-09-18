##' Age perturb proxy data
##'
##' This function age perturbs input proxy data given on a reference chronology
##' by mapping the proxy values onto a perturbed chronology.
##'
##' The mapping is performed by index matching after rounding of the perturbed
##' age values to integer values (while perturbed ages from \code{\link{BAM}}
##' are integer, non-integer ages can result from tranferring to a Brownian
##' Bridge process (see \code{\link{BrownianBridge}})).
##' @param t.pert numeric vector with the age perturbed chronology to map the
##' proxy data onto.
##' @param t.orig numeric vector with the original reference chronology of the
##' data in \code{x}.
##' @param x numeric vector of the proxy data to perturb.
##' @return Numeric vector with the age perturbed proxy data of the same length
##' as \code{x}.
##' @author Thomas Münch
##' @export
PerturbData <- function(t.pert, t.orig, x) {

    i      <- match(round(t.pert), t.orig)
    x.pert <- x[i]

    return(x.pert)

}

