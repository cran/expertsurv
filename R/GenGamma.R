
## Log-gamma or generalized gamma distribution (parameterisation as in Farewell and Prentice, Technometrics 1977)
### FIXME value for x = 0 
dgengamma <- function(x, mu=0, sigma=1, Q, log=FALSE) {
    check_numeric(x=x, mu=mu, sigma=sigma, Q=Q)
    dgengamma_work(x, mu, sigma, Q, log)
}


pgengamma <- function(q, mu=0, sigma=1, Q, lower.tail = TRUE, log.p = FALSE) {
    check_numeric(q=q, mu=mu, sigma=sigma, Q=Q)
    pgengamma_work(q, mu, sigma, Q, lower.tail, log.p)
}



Hgengamma <- function(x, mu=0, sigma=1, Q)
{
    -pgengamma(q=x, mu=mu, sigma=sigma, Q=Q, lower.tail=FALSE, log.p=TRUE)
}

##' 

hgengamma <- function(x, mu=0, sigma=1, Q)
{
  logdens <- dgengamma(x = x, mu = mu, sigma = sigma, Q = Q, log=TRUE)
  logsurv <- pgengamma(q = x, mu = mu, sigma = sigma, Q = Q, lower.tail = FALSE, log.p=TRUE)
  loghaz <- logdens - logsurv
  exp(loghaz)
}


qgengamma <- function(p, mu=0, sigma=1, Q, lower.tail = TRUE, log.p = FALSE)
{
    d <- dbase("gengamma", lower.tail=lower.tail, log=log.p, p=p, mu=mu, sigma=sigma, Q=Q)
    for (i in seq_along(d)) assign(names(d)[i], d[[i]])
    p[Q<0] <- 1 - p[Q<0]
    ret[ind] <- numeric(sum(ind))
    ret[ind][Q==0] <- qlnorm(p[Q==0], mu[Q==0], sigma[Q==0])
    qn0 <- Q!=0
    p <- p[qn0]; mu <- mu[qn0]; sigma <- sigma[qn0]; Q <- Q[qn0]
    ret[ind][qn0] <- exp(mu + sigma*(log(Q^2*qgamma(p, 1/Q^2, 1)) / Q))
    ret
}



rgengamma <- function(n, mu=0, sigma=1, Q) {
    r <- rbase("gengamma", n=n, mu=mu, sigma=sigma, Q=Q)
    for (i in seq_along(r)) assign(names(r)[i], r[[i]])
    ret[ind][Q==0] <- rlnorm(n, mu, sigma)
    qn0 <- Q!=0
    if (any(qn0)) {
        mu <- mu[qn0]; sigma <- sigma[qn0]; Q <- Q[qn0]
        w <- log(Q^2*rgamma(n, 1/Q^2, 1)) / Q
        ret[ind][qn0] <- exp(mu + sigma*w)
    }
    ret
}


rmst_gengamma = function(t, mu=0, sigma=1, Q, start=0){
  rmst_generic(pgengamma, t, start=start, mu=mu, sigma=sigma, Q=Q)
}


mean_gengamma = function(mu=0, sigma=1, Q){
  rmst_generic(pgengamma, Inf, start=0, mu=mu, sigma=sigma, Q=Q)
}

### FIXME limiting value for x=0:  0 if bk >1, 1 if b=k=1, ... ? 


dgengamma.orig <- function(x, shape, scale=1, k, log=FALSE){
    d <- dbase("gengamma.orig", log=log, x=x, shape=shape, scale=scale, k=k)
    for (i in seq_along(d)) assign(names(d)[i], d[[i]])
    logdens <- log(shape) - lgamma(k) + (shape*k - 1)*log(x) - shape*k*log(scale) - (x/scale)^shape
    ret[ind] <- if (log) logdens else exp(logdens)
    ret
}



pgengamma.orig <- function(q, shape, scale=1, k, lower.tail = TRUE, log.p = FALSE) {
    d <- dbase("gengamma.orig", lower.tail=lower.tail, log=log.p, q=q, shape=shape, scale=scale, k=k)
    for (i in seq_along(d)) assign(names(d)[i], d[[i]])
    y <- log(q)
    w <- (y - log(scale))*shape
    prob <- pgamma(exp(w), shape=k)
    if (!lower.tail) prob <- 1 - prob
    if (log.p) prob <- log(prob)
    ret[ind] <- prob
    ret
}



Hgengamma.orig <- function(x, shape, scale=1, k)
{
    -log(pgengamma.orig(q=x, shape=shape, scale=scale, k=k, lower.tail=FALSE))
}



hgengamma.orig <- function(x, shape, scale=1, k)
{
  logdens <- dgengamma.orig(x = x, shape = shape, scale = scale, k = k, log=TRUE)
  logsurv <- pgengamma.orig(q = x, shape = shape, scale = scale, k = k, lower.tail = FALSE, log.p=TRUE)
  loghaz <- logdens - logsurv
  exp(loghaz)
}



qgengamma.orig <- function(p, shape, scale=1, k, lower.tail = TRUE, log.p = FALSE)
{
    d <- dbase("gengamma.orig", lower.tail=lower.tail, log=log.p, p=p, shape=shape, scale=scale, k=k)
    for (i in seq_along(d)) assign(names(d)[i], d[[i]])
    w <- log(qgamma(p, shape=k))
    y <- w / shape  + log(scale)
    ret[ind] <- exp(y)
    ret
}


rgengamma.orig <- function(n, shape, scale=1, k) {
    r <- rbase("gengamma.orig", n=n, shape=shape, scale=scale, k=k)
    for (i in seq_along(r)) assign(names(r)[i], r[[i]])
    w <- log(rgamma(n, shape=k))
    y <- w / shape  + log(scale)
    ret[ind] <- exp(y)
    ret
}



rmst_gengamma.orig = function(t, shape, scale=1, k, start=0){
  rmst_generic(pgengamma.orig, t, start=start, shape=shape, scale=scale, k=k)
}


mean_gengamma.orig = function(shape, scale=1, k){
  rmst_generic(pgengamma.orig, Inf, start=0, shape=shape, scale=scale, k=k)
}

check.gengamma.orig <- function(shape, scale, k){
    ret <- rep(TRUE, length(shape))
    if (missing(shape)) stop("shape parameter \"shape\" not given")
    if (missing(k)) stop("shape parameter \"k\" not given")
    if (any(!is.na(shape) & shape < 0)) {
        warning("Negative shape parameter \"shape\"")
        ret[!is.na(shape) & shape < 0] <- FALSE
    }
    if (any(!is.na(scale) & scale < 0)) {
        warning("Negative scale parameter");
        ret[!is.na(scale) & scale < 0] <- FALSE
    }
    if (any(!is.na(k) & k < 0)) {
        warning("Negative shape parameter \"k\"");
        ret[!is.na(k) & k < 0] <- FALSE
    }
    ret
}
