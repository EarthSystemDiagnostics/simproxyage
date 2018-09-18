##' Constrained Brownian Bridge process
##'
##' Transfer an unconstrained age perturbation process into a constrained
##' perturbation process following the Brownian Bridge concept, assuming an age
##' control point where the age uncertainty is taken to be zero.
##' @param x numeric vector of the perturbed ages from the unconstrained
##' process; must be of the same length as \code{t}.
##' @param t numeric vector with the original (reference) ages.
##' @param acp the age control point where the age uncertainty is assumed zero.
##' @return Numeric vector of the same length as \code{t} with the perturbed
##' ages following a Brownian Bridge process.
##' @author Thomas Münch
##' @examples
##' t <- 100 : 0
##' ns <- 100
##' bam <- BAM(t, ns = ns)
##' bb <- bam$Tp
##' bb[, 1, ] <- apply(bam$Tp[, 1, , drop = FALSE], c(2, 3), BrownianBridge,
##'                    bam$tp, acp = 1)
##' plot(t, t, type = "n", xlim = c(100, 0), ylim = c(100, 0),
##'      xlab = "Original time", ylab = "Perturbed time")
##' for (i in 1 : ns) {
##' lines(bam$tp, bb[, 1, i], col = "grey", lwd = 1)
##' }
##' abline(0, 1, lwd = 2)

##' @export
BrownianBridge <- function(x, t, acp) {

    t.i   <- seq(t) - 1
    t.acp <- which(t == acp)

    bb <- x - t
    bb <- bb - t.i / (t.acp - 1) * bb[t.acp]
    bb <- bb + t

    return(bb)

}

