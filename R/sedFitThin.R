sedFitThin <-
function(s, e=s*0.2, z=2.5, nsamp=100, alpha=2, beta=1.5,
                       wl=c(250, 350, 500), sc.start=1.e-6, T.start=50) {

    # Function takes Herschel-SPIRE photometry and fits optically-thin
    # greybody function for a single-component temperature and galaxy
    # luminosity, as described in Blain, Barnard & Chapman 2003,
    # MNRAS 338, 733.  Conversion from observed to rest frame is from
    # equation (24) in Hogg 2000, astro-ph 9905116v4.
    #
    # Function generates nsamp realizations of observed flux densities s
    # with standard devaitions e for error analysis.  Results are returned
    # in a matrix with nsamp rows and 3 colums: temperature, luminosity from
    # greybody alone, and luminosity from a greybody with smoothly-joined
    # power law to short wavelengths (see Blain et al. 2003).

    # A. Harris, 2013.1.5, 2013.1.14

    # fit routine needs numeric values, is touchy about lists from data tables
    s <- as.numeric(s)
    e <- as.numeric(e)
    z <- as.numeric(z)
    alpha <- as.numeric(alpha)
    beta <- as.numeric(beta)
    wl <- as.numeric(wl)

    # generate data for monte-carlo analysis; first row gets
    # center-of-error values
    l.nue <- matrix(nrow=nsamp, ncol=length(s))
    for (i in 1:length(s)) l.nue[,i] <- rnorm(nsamp, s[i], e[i])
    l.nue[1,] <- s

    # rest frame band centers in Hz
    nue <- 3e5/wl*(1+z)

    # observed frame flux density to rest frame luminosity density conversion
    scaleFactor <- 4*pi*(D.L(z)*3.08567758e22)^2/(1+z)  # in SI units
    scaleFactor <- scaleFactor*1.e9  # GHz to Hz
    scaleFactor <- scaleFactor*1e-26  # Jy to SI units
    scaleFactor <- scaleFactor/3.839e26  # W to L_sun

    # define functions
    # optically thin greybody
    otGreybody <- function(nu, T, beta, sc=1) {
                           # nu in GHz, T in K, beta and sc unitless
                           sc*nu^(3+beta)/(exp(0.04801449*nu/T) - 1)
                       }
    # high frequency tail
    hfTail <- function(nu, alpha) nu^-alpha

    # matrix for results
    results <- matrix(nrow=nsamp, ncol=5)
    # fit control
    control <- nls.control(maxiter=20, warnOnly=TRUE)
    # Frequency range for luminosities: 1000um = 300 GHz; 8um = 37500 GHz
    nu.low <- 300
    nu.high <- 37500

    # fit model to rest frame luminosities
    for (i in 1:nsamp) {
        fd <- l.nue[i,]
        fit <- nls(fd ~ otGreybody(nue, T, beta, sc),
                   control=control,
                   start=list(sc=sc.start, T=T.start))
        cfit <- coef(fit)
        if (fit$convInfo$isConv) {
            # temperature
            results[i, 1] <- cfit[2]
            results[i, 4] <- cfit[1]
            # high-frequency tail
            nu.t <- (3 + beta + alpha)*20.82705*cfit[2]
            val.t <- otGreybody(nu=nu.t, T=cfit[2], beta=beta, sc=cfit[1])
            results[i,5] <- nu.t

            # luminosity of grey body and hybrid greybody+power law
            results[i, 2] <- integrate(otGreybody, nu.low, nu.high,
                                       T=cfit[2], beta=beta, sc=cfit[1])$value
            results[i, 3] <- integrate(otGreybody, nu.low, nu.t,
                                       T=cfit[2], beta=beta, sc=cfit[1])$value
            results[i, 3] <-  results[i, 3] +
                    (val.t*nu.t^alpha/(1-alpha))*
                        (nu.high^(1-alpha) - nu.t^(1-alpha))
        }
    }
    # scale luminosities to L_sun
    results[,2] <- results[,2]*scaleFactor
    results[,3] <- results[,3]*scaleFactor
    results[,4] <- results[,4]*scaleFactor

    # some statistics
    td <- mean(results[,1], na.rm=TRUE)
    e.td <- sd(results[,1], na.rm=TRUE)
    lum.gb <- mean(results[,2], na.rm=TRUE)
    e.lum.gb <- sd(results[,2], na.rm=TRUE)
    lum.gbpl <- mean(results[,3], na.rm=TRUE)
    e.lum.gbpl <- sd(results[,3], na.rm=TRUE)

    # return
    success <- 1 - length(which(is.na(results[,1])))/nsamp
    out <- list(td=td, e.td=e.td, lum.gb=lum.gb, e.lum.gb=e.lum.gb,
                lum.gbpl=lum.gbpl, e.lum.gbpl=e.lum.gbpl,
                scaleFactor=scaleFactor, success=success, results=results)
    class(out) <- 'sedfit'
    out
}
