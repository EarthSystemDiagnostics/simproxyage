##' Age perturb proxy core segments
##'
##' \code{PerturbCoreSegment} age perturbs the proxy data of an individual
##' segment of a similarly dated core array, i.e. the section of the cores that
##' lies between two specified age control points.
##'
##' \code{PerturbCoreSegment}'s intended use is in combination with
##' \code{PerturbCoreArray} only rather than as a standalone function (though
##' possible).
##' @param Xp array of dimension \code{c(dim(X), ns)}, where \code{ns} is the
##' number of age perturbation realisations, to store the perturbed proxy
##' data. To be provided with the calling function.
##' @param Tp array of the same dimension as \code{Xp} to store the age
##' perturbation realisations. To be provided with the calling function.
##' @param acp1 the first age control point, i.e. the start age at the top of
##' the core segment (per default \code{t[1]}).
##' @param acp2 the second age control point at the end of the core segment;
##' \code{NA} (the default) results in an unconstrained age perturbation
##' process.
##' @param nc integer number of cores in the core array; per default set to the
##' number of columns in \code{X}.
##' @inheritParams PerturbCoreArray
##' @return A list of two components:
##' \describe{
##' \item{Tp:}{The array of age perturbation realisations filled with the
##' realisations corresponding to the time interval between \code{acp1} and
##' \code{acp2}.}
##' \item{Xp:}{The array of age perturbed proxy data filled with the
##' perturbation realisations corresponding to the time interval between
##' \code{acp1} and \code{acp2}.}
##' }
##' @author Thomas Münch
##' @references Comboul, M., Emile-Geay, J., Evans, M. N., Mirnateghi, N., Cobb,
##' K. M. and Thompson, D. M.: A probabilistic model of chronological errors in
##' layer-counted climate proxies: applications to annually banded coral
##' archives, Clim. Past, 10(2), 825-841, doi: 10.5194/cp-10-825-2014, 2014.
##' @seealso \code{\link{PerturbCoreArray}}
##' @export
PerturbCoreSegment <- function(X, t, Xp = array(dim = c(dim(X), ns)),
                               Tp = Xp, acp1 = t[1], acp2 = NA, nc = dim(X)[2],
                               ns = 1000, model = "poisson", rate = 0.05,
                               resize = 1) {

    t.final <- acp2
    if (is.na(t.final)) {
        t.final <- t[length(t)]
    }
    
    i <- which(t <= acp1 & t >= t.final)
    n.seg <- length(i)

    t.in <- t[i]
    x.in <- X[i, , drop = FALSE]

    bam <- BAM(t.in, ns = ns, nc = nc, model = model,
               rate = rate, resize = resize)

    bam.Tp <- bam$Tp    
    if (!is.na(acp2)) {

        for (i.c in 1 : nc) {
            
            bam.Tp[, i.c, ] <- apply(bam.Tp[, i.c, , drop = FALSE], c(2, 3),
                                     BrownianBridge, bam$tp, t.final)
        }
    }

    bam.Xp <- array(dim = dim(bam.Tp))
    
    for (i.c in 1 : nc) {

        bam.Xp[, i.c, ] <- apply(bam.Tp[, i.c, , drop = FALSE], c(2, 3),
                                 PerturbData, bam$tp, x.in[, i.c])
    }

    Xp[i, , ] <- bam.Xp[1 : n.seg, , ]
    Tp[i, , ] <- bam.Tp[1 : n.seg, , ]

    res <- list()
    res$Tp <- Tp
    res$Xp <- Xp

    return(res)

}


    
