##' Generate age perturbations of layer-counted proxy data
##'
##' This function creates an ensemble of perturbed ages of a given reference
##' chronology that is based on layer counting ("banded age model").
##'
##' Ensemble of age perturbations can be realised assuming a single proxy record,
##' or a set of similarly layer-counted proxy data (i.e. a certain number of
##' "cores"). The random age perturbations can either follow a Poisson process
##' or a Bernoulli process; they are realised by randomly removing or doubling
##' time points of the original age model which can result in realisations that
##' are longer (shorter) than the original chronology; in such a case, the
##' output data can be cropped (extended) to match the original length depending
##' on the option specified via \code{resize}. Note that the perturbed ages are
##' automatically flipped to range from most recent to oldest measurements when
##' the input ages are given in increasing order. The code is a modified version
##' of the MATLAB function \code{BAM_simul_perturb.m} provided with the
##' publication  Comboul et al. (2014); please see there for a detailed
##' description of the method.
##' @param t numeric vector providing the original layer-counted chronology.
##' @param ns integer number of age perturbation realisations for each core.
##' @param nc integer number of cores for which an ensemble of age perturbations
##' shall be realised.
##' @param model name string of the random process to use for perturbing the age
##' model; must be either "poisson" (the default) or "bernoulli"; see Comboul et
##' al. (2014) for details on the two models.
##' @param rate numeric vector of probability rate(s) that an age band is
##' perturbed; you can specify a vector of two rates where the first entry is
##' the probability for a missing band and the second entry the probability for
##' a double-counting of a band. If only a single value is specified (per
##' default 0.05), symmetric perturbations are assumed.
##' @param resize the resizing option in case of shorter/longer than original
##' time axes: 0 = do not resize, -1 = resize to shortest realisation, 1 =
##' resize to longest realisation (default).
##' @return A list of two elements:
##' \describe{
##' \item{tp:}{the old chronology cropped/extended depending on the resizing
##' option (thus tp = t for resize = 0).}
##' \item{Tp:}{an array of dimension length(\code{t}) * \code{nc} * \code{ns} with
##' the \code{ns} realisations of perturbed ages for each core \code{nc}.}
##' }
##' @author Maud Comboul; translated and adopted from the original MATLAB code
##' by Thomas Münch
##' @references Comboul, M., Emile-Geay, J., Evans, M. N., Mirnateghi, N., Cobb,
##' K. M. and Thompson, D. M.: A probabilistic model of chronological errors in
##' layer-counted climate proxies: applications to annually banded coral
##' archives, Clim. Past, 10(2), 825-841, doi: 10.5194/cp-10-825-2014, 2014.
##' @examples
##' t <- 100 : 0
##' ns <- 100
##' bam <- BAM(t, ns = ns)
##' plot(t, t, type = "n", xlim = c(100, 0), ylim = c(100, 0),
##'      xlab = "Original time", ylab = "Perturbed time")
##' for (i in 1 : ns) {
##' lines(bam$tp, bam$Tp[, 1, i], col = "grey", lwd = 1)
##' }
##' abline(0, 1, lwd = 2)
##' @export
BAM <- function(t, ns = 1000, nc = 1, model = "poisson",
                rate = 0.05, resize = 1) {

    if (length(t) == 0) {
        stop("Invalid input time vector.")
    }
    
    if (mean(diff(t)) > 0) {
        isflipped <- 1
        t <- rev(t)
    } else {
        isflipped <- 0
    }
    
    n <- length(t)
    p <- nc

    if (length(rate) == 1) {
        rate <- rep(rate, 2)
    }


    # Generate an ensemble of time perturbation models

    time.model <- array(data = 1, dim = c(n, p, ns))
    if (model == "poisson") {
        
        # apply poisson model
        for (nn in 1 : ns) {

            # poisson model for missing bands
            num_event_mis <- rpois(p, rate[1] * n)
            # poisson model for layers counted multiple times
            num_event_dbl <- rpois(p, rate[2] * n)

            for (ii in 1 : p) {
                # place events uniformly on {2,...,n}
                jumps <- sample(n-1, num_event_mis[ii], replace = FALSE) + 1
                # remove 1 at jump locations
                time.model[jumps, ii, nn] <- time.model[jumps, ii, nn] - 1
                # add 1 at jump locations
                jumps <- sample(n-1, num_event_dbl[ii], replace = FALSE) + 1
                time.model[jumps, ii, nn] <- time.model[jumps, ii, nn] + 1
            }
        }

    } else if (model == "bernoulli") {

        # apply bernoulli model
        for (nn in 1 : ns){

            # binomial model for missing bands
            time.model[, , nn] <- time.model[, , nn] -
                array(data = rbinom(n = n * p, size = 1, prob = rate[1]),
                      dim = c(n, p))
            # binomial model for layers counted multiple times
            time.model[, , nn] <- time.model[, , nn] +
                array(data = rbinom(n = n * p, size = 1, prob = rate[2]),
                      dim = c(n, p))
        }

    } else {
        
        stop(paste("Unknown age model;",
                   "acceptable inputs are ''poisson'' and ''bernoulli''."))
    }


    
    # Expand length of time vector if resizing is required
    
    if (resize) {
        
        t_ext <- ceiling(2 * rate[2] * n)
        tn <- n + t_ext

        dt <- t[2] - t[1]
        time_ext <- seq(t[length(t)] + dt, t[length(t)] + t_ext * dt, dt)
        tp <- c(t, time_ext)
        
    } else {
        
        tn <- n
        tp <- t
    }

    
    # Create the pertubed chronologies
    
    Tmax <- 0
    Tmin <- n
    Tp <- array(dim = c(tn, p, ns))
    
    for (nn in 1 : ns) {

        for (ii in 1 : p) {
            xcount <- 1
            Xcount <- 1
            tt <- 1
            while (tt < (n + 1)) {

                if (time.model[tt, ii, nn] == 0) {
                    # remove band
                    Xcount <- min(Xcount + 1, tn)
                } else if (time.model[tt, ii, nn] == 2) {
                    # insert double band
                    Tp[xcount, ii, nn] <- tp[Xcount]
                    xcount <- min(tn, xcount + 1)
                }

                Tp[xcount, ii, nn] <- tp[Xcount]
                xcount <- min(tn, xcount + 1)
                Xcount <- min(tn, Xcount + 1)
                tt <- tt + 1
            }
            
            k <- which(!is.na(Tp[, ii, nn]))
            k <- length(k)
            if (k > Tmax) {
                Tmax <- k
            }
            if (k < Tmin) {
                Tmin <- k
            }
        }
    }


    # Crop output size to shortest sequence
    if (resize == -1) {
        n <- Tmin
    }
    
    # Expand output size to longest (non-NA) sequence
    if (resize == 1) {
        n <- Tmax
    }
    
    Tp <- Tp[1 : n, , , drop = FALSE]
    tp <- tp[1 : n]

    if (isflipped == 1) {
        
        for (nn in 1 : ns) {
            Tp[, , nn] <- apply(Tp[, , nn, drop = FALSE], 2, rev)
        }
        
        tp <- rev(tp)
    }

    res <- list()
    res$tp <- tp
    res$Tp <- Tp
    
    return(res)

}


