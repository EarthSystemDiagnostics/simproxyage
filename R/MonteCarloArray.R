##' Monte Carlo simulation of age-perturbed core array
##'
##' For a core array of \code{nc} proxy records, this wrapper function creates
##' \code{nc} identical surrogate time series, independently age-perturbs
##' each series following the model in Comboul et al. (2014) and averages the
##' age-perturbed time series. The process is repeated \code{ns} times. For a
##' common signal recorded by the core array, the approach allows one to study
##' the effect of independent age errors on this signal upon averaging across
##' the core array.
##'
##' Note the principle difference between using \code{PerturbCoreArray} alone
##' and using \code{MonteCarloArray}: while \code{PerturbCoreArray} can also
##' realise \code{ns} age perturbation realisations, these realisations are
##' applied to the same data set in the cores. Here, \code{ns} simulations are
##' run over different surrogate data sets of the cores (however, for the same
##' data in each core) and each simulation run corresponds to one set of age
##' perturbation realisations across the cores. In other words,
##' \code{MonteCarloArray} runs \code{ns} simulations where for each simulation
##' \code{PerturbCoreArray} is called setting the \code{ns} parameter in the
##' function call there to 1.
##' @param t numeric vector of integer values providing a reference chronology
##' for the age perturbations (starting with the youngest age)
##' @param acp numeric vector of age control points where the age uncertainty is
##' assumed to be zero. Per default, a two-element vector where the first
##' element is the start age (\code{t[1]}) and the second element \code{NA}
##' resulting in an unconstrained age perturbation process; see
##' \code{?PerturbCoreArray} for details.
##' @param nt the length of the records (i.e. the number of data points) to
##' simulate; per default set to \code{length(t)}.
##' @param nc the number of cores in the modelled core array
##' @param ns the number of Monte Carlo simulations
##' @inheritParams BAM
##' @param surrogate.fun the random number generator to use for creating the
##' surrogate time series; per default the base \code{R} white noise generator.
##' @param ... additional parameters passed to \code{surrogate.fun}
##' @return a list with components \code{input} and \code{stack}, where both are
##' \code{nt * ns} arrays:
##' \itemize{
##' \item \code{input} contains the \code{ns} realisations of
##' the original surrogate time series
##' \item \code{stack} contains
##' the \code{ns} averages across the \code{nc} individually age-perturbed
##' surrogate time series
##' }
##' @author Thomas Münch
##' @seealso \code{\link{PerturbCoreArray}}
##' @export
MonteCarloArray <- function(t = 100 : 1, acp = c(t[1], NA), nt = length(t),
                            nc = 1, ns = 100, model = "poisson",
                            rate = 0.05, resize = 1,
                            surrogate.fun = rnorm, ...) {

    stacks <- array(data = NA, dim = c(nt, ns))
    input  <- array(data = NA, dim = c(nt, ns))
    
    for (i in 1 : ns) {
        
        # surrogate data to perturb (= nc times the same noise time series)
        X <- array(surrogate.fun(nt, ...), dim = c(nt, nc))
        input[, i] <- X[, 1]

        bam <- PerturbCoreArray(X, t, acp = acp, ns = 1,
                                model = model, rate = rate, resize = resize)

        stacks[, i] <- apply(bam$Xp, 1, mean)

    }

    res <- list()
    res$input <-  input
    res$stacks <- stacks
    
    return(res)

}

