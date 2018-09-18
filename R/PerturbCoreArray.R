##' Age perturb proxy core array
##'
##' \code{PerturbCoreArray} age perturbs the proxy data of a similarly dated
##' core array.
##'
##' The function allows one to model the age uncertainty of proxy records with
##' several age control points (ACP) where the age uncertainty evolves according
##' to the model of Comboul et al. (2014) between the age control points but is
##' forced back to zero time uncertainty at the ACPs following a Brownian Bridge
##' process.
##'
##' Per default, the start age (i.e. the youngest age at the core top)
##' is assumed to be the first ACP, thus, if not set explicitly in the vector
##' \code{acp}, the youngest age is added as its first element. You can specify
##' an arbitrary number of additional ACPs. Between each pair of ACPs (starting
##' at the core top), \code{ns} number of constrained age perturbation
##' realisations are performed following the Comboul et al. (2014) model and the
##' Brownian Bridge concept. If the last ACP equals the oldest age in \code{t}
##' (i.e. the last proxy data point), the last core segment also follows a
##' Brownian Bridge process. Alternatively, if the last ACP is younger than the
##' final age, \code{NA} is added as the last element to the vector \code{acp}
##' which results in an unconstrained age perturbation process for the last core
##' segment.
##'
##' ACPs lying outside the time interval defined by \code{t} will be removed
##' from the vector \code{acp} with a warning.
##'
##' The actual age perturbations are performed by calling
##' \code{\link{PerturbCoreSegment}} for each consecutive time interval between
##' two ACPs.
##' @param X an nt * nc array of the proxy data from the core array, where the
##' number of rows (\code{nt}) is the length of the cores (i.e. the length of
##' the proxy data) and the number of columns (\code{nc}) the number of cores in
##' the array.
##' @param t numeric vector providing the original layer-counted
##' chronology. Must equal the length of the proxy data, thus \code{dim(X)[1]}.
##' @param acp numeric vector of age control points where the age uncertainty is
##' assumed to be zero. Per default, a two-element vector where the first
##' element is the start age (\code{t[1]}) and the second element \code{NA}
##' resulting in an unconstrained process; see details for handling more age
##' control points.
##' @inheritParams BAM
##' @return A list of two components:
##' \describe{
##' \item{Tp:}{an nt * nc * ns array of age perturbation realisations for the
##' core array.}
##' \item{Xp:}{an nt * nc * ns array of the corresponding age pertubed proxy
##' data.}
##' }
##'
##' However, the dimension of the arrays depend on the values of \code{nc}
##' (number of cores) and \code{ns} (number of realisations). If one or both
##' values equal 1, the output arrays are converted to the lowest possile
##' dimensionality.
##' @author Thomas Münch
##' @references Comboul, M., Emile-Geay, J., Evans, M. N., Mirnateghi, N., Cobb,
##' K. M. and Thompson, D. M.: A probabilistic model of chronological errors in
##' layer-counted climate proxies: applications to annually banded coral
##' archives, Clim. Past, 10(2), 825-841, doi: 10.5194/cp-10-825-2014, 2014.
##' @seealso \code{\link{BAM}}, \code{\link{BrownianBridge}},
##' \code{\link{PerturbCoreSegment}}
##' @export
PerturbCoreArray <- function(X, t, acp = c(t[1], NA), ns = 1000,
                             model = "poisson", rate = 0.05, resize = 1) {

    nt <- dim(X)[1]
    nc <- dim(X)[2]

    if (nt != length(t)) {
        stop("Incompatible lengths of time vector and data.")
    }

    i <- match(acp, t)
    if (!all(!is.na(i))) {
        warning("Some ACPs lie outside of given time vector; will be clipped.")
        acp <- acp[-which(is.na(i))]
    }

    # always treat youngest age as first age control point
    if (acp[1] != t[1]) {
        acp <- c(t[1], acp)
    }

    # treat oldest age as 'NA' age control point
    # if not equal to actual last acp
    if (acp[length(acp)] != t[nt] & !is.na(acp[length(acp)])) {
        acp <- c(acp, NA)
    }

    Xp <- array(data = NA, dim = c(nt, nc, ns))
    Tp <- array(data = NA, dim = c(nt, nc, ns))

    for (i in 1 : (length(acp) - 1)) {

        acp1 <- acp[i]
        acp2 <- acp[i + 1]
        pert <- PerturbCoreSegment(X, t, Xp, Tp, acp1, acp2, ns = ns,
                                   model = model, rate = rate, resize = resize)
        
        Xp <- pert$Xp
        Tp <- pert$Tp

    }

    # simplify arrays to simplest format
    Tp <- drop(Tp)
    Xp <- drop(Xp)

    res <- list()
    res$Tp <- Tp
    res$Xp <- Xp

    return(res)

}

        
