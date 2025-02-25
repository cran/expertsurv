##' Summaries of fitted flexible survival models
##'
##' Return fitted survival, cumulative hazard or hazard at a series of times
##' from a fitted \code{\link{flexsurvreg}} or \code{\link{flexsurvspline}}
##' model.
##'
##' Time-dependent covariates are not currently supported.  The covariate
##' values are assumed to be constant through time for each fitted curve.
##'
##' @param object Output from \code{\link{flexsurvreg}} or
##' \code{\link{flexsurvspline}}, representing a fitted survival model object.
##'
##' @param newdata Data frame containing covariate values to produce fitted
##' values for.  Or a list that can be coerced to such a data frame.  There
##' must be a column for every covariate in the model formula, and one row for
##' every combination of covariates the fitted values are wanted for.  These
##' are in the same format as the original data, with factors as a single
##' variable, not 0/1 contrasts.
##'
##' If this is omitted, if there are any continuous covariates, then a single
##' summary is provided with all covariates set to their mean values in the
##' data - for categorical covariates, the means of the 0/1 indicator variables
##' are taken.  If there are only factor covariates in the model, then all
##' distinct groups are used by default.
##'
##' @param X Alternative way of defining covariate values to produce fitted
##' values for.  Since version 0.4, \code{newdata} is an easier way that
##' doesn't require the user to create factor contrasts, but \code{X} has been
##' kept for backwards compatibility.
##'
##' Columns of \code{X} represent different covariates, and rows represent
##' multiple combinations of covariate values.  For example
##' \code{matrix(c(1,2),nrow=2)} if there is only one covariate in the model,
##' and we want survival for covariate values of 1 and 2.  A vector can also be
##' supplied if just one combination of covariates is needed.
##'
##' For ``factor'' (categorical) covariates, the values of the contrasts
##' representing factor levels (as returned by the \code{\link{contrasts}}
##' function) should be used.  For example, for a covariate \code{agegroup}
##' specified as an unordered factor with levels \code{20-29, 30-39, 40-49,
##' 50-59}, and baseline level \code{20-29}, there are three contrasts.  To
##' return summaries for groups \code{20-29} and \code{40-49}, supply \code{X =
##' rbind(c(0,0,0), c(0,1,0))}, since all contrasts are zero for the baseline
##' level, and the second contrast is ``turned on'' for the third level
##' \code{40-49}.
##'
##' @param type \code{"survival"} for survival probabilities.
##'
##' \code{"cumhaz"} for cumulative hazards.
##'
##' \code{"hazard"} for hazards.
##'
##' \code{"rmst"} for restricted mean survival.
##'
##' \code{"mean"} for mean survival.
##'
##' \code{"median"} for median survival (alternative to \code{type="quantile"} with \code{quantiles=0.5}).
##'
##' \code{"quantile"} for quantiles of the survival time distribution.
##'
##' \code{"link"} for the fitted value of the location parameter (i.e. the "linear predictor" but on the natural scale of the parameter, not on the log scale)
##'
##' Ignored if \code{"fn"} is specified.
##'
##' @param fn Custom function of the parameters to summarise against time.
##' This has optional first two arguments \code{t} representing time, and
##' \code{start} representing left-truncation points, and any remaining
##' arguments must be parameters of the distribution.  It should be vectorised, and
##' return a vector corresponding to the vectors given by \code{t}, \code{start} and
##' the parameter vectors.
##'
##' @param t Times to calculate fitted values for. By default, these are the
##' sorted unique observation (including censoring) times in the data - for
##' left-truncated datasets these are the "stop" times.
##'
##' @param quantiles If \code{type="quantile"}, this specifies the quantiles of the survival time distribution to return estimates for.
##'
##'
##' @param cross If \code{TRUE} (the default) then summaries are calculated for all combinations of times
##' specified in \code{t} and covariate vectors specifed in \code{newdata}.
##'
##' If \code{FALSE},
##' then the times \code{t} should be of length equal to the number of rows in \code{newdata},
##' and one summary is produced for each row of \code{newdata} paired with the corresponding
##' element of \code{t}. This is used, e.g. when determining Cox-Snell residuals.
##'
##' @param ci Set to \code{FALSE} to omit confidence intervals.
##'
##' @param se Set to \code{TRUE} to include standard errors.
##'
##' @param B Number of simulations from the normal asymptotic distribution of
##' the estimates used to calculate confidence intervals or standard errors.
##' Decrease for greater
##' speed at the expense of accuracy, or set \code{B=0} to turn off calculation
##' of CIs and SEs.
##'
##' @param cl Width of symmetric confidence intervals, relative to 1.
##'
##' @param tidy If \code{TRUE}, then the results are returned as a tidy data
##' frame instead of a list.  This can help with using the \pkg{ggplot2}
##' package to compare summaries for different covariate values.
##'
##' @param na.action Function determining what should be done with missing values in \code{newdata}.  If \code{na.pass} (the default) then summaries of \code{NA} are produced for missing covariate values.  If \code{na.omit}, then missing values are dropped, the behaviour of \code{summary.flexsurvreg} before \code{flexsurv} version 1.2.
##'
##' @param start NULL
##'
##' @param ... Further arguments passed to or from other methods.  Currently
##' unused.
##'
##' @return If \code{tidy=FALSE}, a list with one component for each unique
##' covariate value (if there are only categorical covariates) or one component
##' (if there are no covariates or any continuous covariates).  Each of these
##' components is a matrix with one row for each time in \code{t}, giving the
##' estimated survival (or cumulative hazard, or hazard) and 95\% confidence
##' limits.  These list components are named with the covariate names and
##' values which define them.
##'
##' If \code{tidy=TRUE}, a data frame is returned instead.  This is formed by
##' stacking the above list components, with additional columns to identify the
##' covariate values that each block corresponds to.
##'
##' If there are multiple summaries, an additional list component named
##' \code{X} contains a matrix with the exact values of contrasts (dummy
##' covariates) defining each summary.
##'
##' The \code{\link{plot.flexsurvreg}} function can be used to quickly plot
##' these model-based summaries against empirical summaries such as
##' Kaplan-Meier curves, to diagnose model fit.
##'
##' Confidence intervals are obtained by sampling randomly from the asymptotic
##' normal distribution of the maximum likelihood estimates and then taking quantiles
##' (see, e.g. Mandel (2013)).
##'
##' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
##'
##' @seealso \code{\link{flexsurvreg}}, \code{\link{flexsurvspline}}.
##'
##' @references Mandel, M. (2013). "Simulation based confidence intervals for
##' functions with complicated derivatives." The American Statistician (in
##' press).
##'
##' @keywords models
##' @noRd
summary.flexsurvreg <- function(object, newdata=NULL, X=NULL, type="survival",
                                    fn=NULL, t=NULL, quantiles=0.5, start=0, cross=TRUE,
                                    ci=TRUE, se=FALSE, B=1000, cl=0.95,
                                    tidy=FALSE, na.action=na.pass, ...){
 x <- object
 X <- newdata_to_X(x, newdata, X, na.action)
 type <- match.arg(type, c("survival","cumhaz","hazard","rmst","mean","median", "quantile","link"))
 if (is.null(fn)) {
  fn <- summary_fns(x, type)
 }
 fn <- expand.summfn.args(fn)
 if (type=="link")
   x$aux$"(location)" <- x$dlist$location
 args <- xt_to_fnargs(x, X, t, quantiles, start, type, cross)

 est <- do.call(fn, args)
 if (type %in% c("median","mean"))
     res <- data.frame(est=est)
 else if (type=="quantile")
     res <- data.frame(quantile = args$t, est=est)
 else res <- data.frame(time = args$t, est=est)
 if (ci || se){
     res.ci <- cisumm.flexsurvreg(x, args$t, args$start, attr(args, "X"), fn=fn, B=B, cl=cl)
     if (ci) {
         res$lcl <- res.ci[1,]
         res$ucl <-  res.ci[2,]
     }
     if (se) res$se <-  res.ci[3,]
 }
 nodata <- is.null(attr(args, "newdata"))
 if (x$ncovs > 0 && !nodata) {
     res <- cbind(res, attr(args, "newdata"))
 }
 rownames(res) <- NULL
 if (!tidy){ ## For backwards compatibility
     if (x$ncovs == 0 || nodata) res <- list(res)
     else {
         nd <- attr(args, "newdata.orig")
         covnames_untidy <- apply(as.data.frame(nd), 1,
                                  function(x)paste0(names(nd), "=", x, collapse=","))
         nx <- attr(args, "nx")
         nt <- attr(args, "nt")
         resdf <- res[,setdiff(names(res), colnames(nd)),drop=FALSE]
         res <- split(resdf, rep(1:nx, each=nt))
         names(res) <- covnames_untidy
         res <- lapply(res, function(x)setNames(as.data.frame(x), names(resdf)))
         for (i in seq_along(res)) row.names(res[[i]]) <- NULL
     }
 }
 if (x$ncovs > 0) attr(res,"X") <- X
 class(res) <- c("summary.flexsurvreg",class(res))
 res
}


newdata_to_X <- function(x, newdata=NULL, X=NULL, na.action=na.pass){
    if (!is.null(X)) {
        X <- as.matrix(X)
        if (!is.matrix(X) || (is.matrix(X) && ncol(X) != x$ncoveffs)) {
            plural <- if (x$ncoveffs > 1) "s" else ""
            stop("expected X to be a matrix with ", x$ncoveffs, " column", plural, " or a vector with ", x$ncoveffs, " element", plural)
        }
        attr(X, "newdata") <- as.data.frame(X)
    }
    else if (is.null(newdata)){
        if (is.null(x[["data"]]))
            stop("`newdata` should be supplied if the data have been removed from the model object")
        Xraw <- model.frame(x)[,unique(attr(model.frame(x),"covnames.orig")),drop=FALSE]
        isfac <- sapply(Xraw, function(x){is.factor(x) || is.character(x)})
        if (is.vector(X)) X <- matrix(X, nrow=1)
        if (x$ncovs > 0 && is.null(X)) {
            ## if any continuous covariates, calculate fitted survival for "average" covariate value
            if (!all(isfac)){
                nd <- colMeans(model.matrix(x))
                X <- matrix(nd ,nrow=1, dimnames=list(NULL,names(nd)))
                attr(X, "newdata") <- as.data.frame(X)
            }
            ## else calculate for all different factor groupings
            else {
                X <- unique(model.matrix(x))
                ## build names like "COVA=value1,COVB=value2"
                nam <- as.matrix(unique(Xraw))
                for (i in 1:ncol(nam)) nam[,i] <- paste(colnames(nam)[i], nam[,i], sep="=")
                rownames(X) <- apply(nam, 1, paste, collapse=",")
                attr(X, "newdata") <- unique(Xraw)
            }
        }
        else if (is.null(X)) X <- as.matrix(0, nrow=1, ncol=max(x$ncoveffs,1))
        else {
            attr(X, "newdata") <- X
            colnames(attr(X, "newdata")) <- colnames(model.matrix(x))
        }
    }
    else
        X <- form.model.matrix(x, as.data.frame(newdata), na.action=na.action)
    X
}



xt_to_fnargs <- function(x, X, t, quantiles=0.5, start=0, type="survival", cross=TRUE){
 tstart <- summfn_to_tstart(x, type, t, quantiles, start)
 t <- tstart$t
 nd <- ndorig <- attr(X, "newdata")
 nt <- length(t)
 nx <- nrow(X)
 if (!cross){
  if (nt != nrow(X)){
   stop(sprintf("length(t)=%s, should equal nrow(X)=%s", nt, nrow(X)))
  }
 } else {
  tstart$t <- rep(t, nx)
  tstart$start <- rep(tstart$start, nx)
  X <- X[rep(1:nx, each=nt),,drop=FALSE]
  nd <- nd[rep(1:nx, each=nt),,drop=FALSE]
 }
 pbase <- x$res.t[x$dlist$pars,"est"]
 beta <- if (x$ncovs==0) 0 else x$res[x$covpars,"est"]
 basepars.mat <- add.covs(x, pbase, beta, X, transform=FALSE)
 basepars <- as.list(as.data.frame(basepars.mat))
 fnargs <- c(tstart, basepars)
 for (j in seq_along(x$aux)){
   fnargs[[names(x$aux)[j]]] <- x$aux[[j]]
 }
 attr(fnargs, "newdata") <- nd
 attr(fnargs, "nx") <- nx
 attr(fnargs, "nt") <- nt
 attr(fnargs, "newdata.orig") <- ndorig
 attr(fnargs, "X") <- X
 fnargs
}


summfn_to_tstart <- function(x, type="survival", t=NULL, quantiles=0.5, start=0){
  nodata_msg <- "prediction times `t` should be defined explicitly if the data are not included in the model object"
  if(type == "mean"){
    if(!is.null(t))
      warning("Mean selected, but time specified.  For restricted mean, set type to 'rmst'.")
    # Type = mean same as RMST w/ time = Inf
    t <- rep_len(Inf,length(start))
  }
  else if(type == "median"){
    if(!is.null(t)) warning("Median selected, but time specified.")
    t <- rep_len(0.5,length(start))
  }
  else if(type == "link"){
    if(!is.null(t)) warning("`link` selected, but time specified.")
    t <- rep_len(0,length(start))
  }
  else if(type == "quantile"){
    t <- quantiles
    if((any(t<0) | any(t>1))){
      stop("Quantiles should not be less than 0 or greater than 1")
    }
    maxlen <- max(length(t), length(start))
    t <- rep_len(t,maxlen)
  }
  else if(type == "rmst"){
      if (is.null(x[["data"]]))
          stop(nodata_msg)
      if (is.null(t))
          t <- max(x$data$Y[,"time1"])
  }
  else if (is.null(t)){
      if (is.null(x[["data"]]))
          stop(nodata_msg)
      t <- sort(unique(x$data$Y[,"stop"]))
  }
  if (length(start)==1)
    start <- rep_len(start, length(t))
  else if (length(start) != length(t))
    stop("length of \"start\" is ",length(start),". Should be 1, or length of \"t\" which is ",length(t))
  list(t=t, start=start)
}


cisumm.flexsurvreg <- function(x, t, start, X, fn, B=1000, cl=0.95) {
    if (all(is.na(x$cov)) || (B==0))
        ret <- array(NA, dim=c(2, length(t)))
    else {
      sim <- normboot.flexsurvreg(x, B, X=X, tidy=TRUE)
      args <- list(t = rep(t, each=B),
                   start = rep(start, each=B))
      args <- c(args, as.list(as.data.frame(sim))[x$dlist$pars])
      for (j in seq_along(x$aux))
          args[[names(x$aux)[j]]] <- x$aux[[j]]
      ret <- do.call(fn, args)
      ret <- matrix(ret, nrow=B)
      retci <- apply(ret, 2, function(x)quantile(x, c((1-cl)/2, 1 - (1-cl)/2), na.rm=TRUE))
      retse <- apply(ret, 2, sd)
      ret <- rbind(retci, retse)
    }
    ret
}

#' @noRd
summary_fns <- function(x, type){
    switch(type,   # TODO warn for clashing arguments in dfns
           "survival" = function(t,start,...) {
               ret <- (1 - x$dfns$p(t,...))/(1 - x$dfns$p(start,...))
               ret[t<start] <- 1 # prob[t<start] was previously 0
               ret
           },
           "median" = function(start,...) {
             start_p = x$dfns$p(start,...)
             med_from_start = start_p + (1 - start_p)/2
             ret = x$dfns$q(med_from_start,...)
           },
           "quantile" = function(t=0.5, start,...) {
             start_p = x$dfns$p(start,...)
             qu_from_start = start_p + (1 - start_p)* t
             ret = x$dfns$q(qu_from_start,...)
           },
           "hazard" = function(t,start,...) {
               ret <- x$dfns$h(t,...) * (1 - x$dfns$p(start,...))
               ret[t<start] <- 0
               ret
           },
           "cumhaz" = function(t,start,...) {
               ret <- x$dfns$H(t,...) - x$dfns$H(start,...)
               ret[t<start] <- 0
               ret
           },
           "rmst" = function(t,start,...) x$dfns$rmst(t,start=start, ...),
           "mean" = function(t,start,...) x$dfns$mean(...),
           "link" = function(...){
               args <- list(...)
               args[[args$"(location)"]]
           }
    )
}

##' @exportS3Method NULL
print.summary.flexsurvreg <- function(x, ...){
    if (!inherits(x, "data.frame")){
        for (i in seq_along(x)){
            cat(names(x)[i], "\n")
            print(x[[i]])
            if (i<length(x)) cat("\n")
        }
    } else print.data.frame(x)
}


add.covs <- function(x, pars, beta, X, transform=FALSE){  ## TODO option to transform on input
    nres <- nrow(X)
    if (!is.matrix(pars)) pars <- matrix(pars, nrow=nres, ncol=length(pars), byrow=TRUE)
    if (!is.matrix(beta)) beta <- matrix(beta, nrow=1)
    for (j in seq_along(x$dlist$pars)){
        covinds <- x$mx[[x$dlist$pars[j]]]
        if (length(covinds) > 0){
            pars[,j] <- pars[,j] + beta[,covinds] %*% t(X[,covinds,drop=FALSE])
        }
        if (!transform)
            pars[,j] <- x$dlist$inv.transforms[[j]](pars[,j])
    }
    colnames(pars) <- x$dlist$pars
    pars
}

## Draw B samples from multivariate normal distribution of baseline
## parameter estimators, for given covariate values



##' Simulate from the asymptotic normal distribution of parameter estimates.
##'
##' Produce a matrix of alternative parameter estimates under sampling
##' uncertainty, at covariate values supplied by the user.  Used by
##' \code{\link{summary.flexsurvreg}} for obtaining confidence intervals around
##' functions of parameters.
##'
##'
##' @param x A fitted model from \code{\link{flexsurvreg}} (or \code{\link{flexsurvspline}}).
##'
##' @param B Number of samples.
##'
##' @param newdata Data frame or list containing the covariate values to
##' evaluate the parameters at.  If there are covariates in the model, at least
##' one of \code{newdata} or \code{X} must be supplied, unless \code{raw=TRUE}.
##'
##' @param X Alternative (less convenient) format for covariate values: a
##' matrix with one row, with one column for each covariate or factor contrast.
##' Formed from all the "model matrices", one for each named parameter of the
##' distribution, with intercepts excluded, \code{cbind}ed together.
##'
##' @param transform \code{TRUE} if the results should be transformed to the
##' real-line scale, typically by log if the parameter is defined as positive.
##' The default \code{FALSE} returns parameters on the natural scale.
##'
##' @param raw Return samples of the baseline parameters and the covariate
##' effects, rather than the default of adjusting the baseline parameters for
##' covariates.
##'
##' @param rawsim allows input of raw samples from a previous run of
##' \code{normboot.flexsurvreg}. This is useful if running
##' \code{normboot.flexsurvreg} multiple time on the same dataset but with
##' counterfactual contrasts, e.g. treat =0 vs. treat  =1.
##' Used in \code{standsurv.flexsurvreg}.
##'
##' @param tidy If \code{FALSE} (the default) then
##' a list is returned.  If \code{TRUE} a data frame is returned, consisting
##' of the list elements \code{rbind}ed together, with integer variables
##' labelling the covariate number and simulation replicate number.
##'
##' @return If \code{newdata} includes only one covariate combination, a matrix
##' will be returned with \code{B} rows, and one column for each named
##' parameter of the survival distribution.
##'
##' If more than one covariate combination is requested (e.g. \code{newdata} is
##' a data frame with more than one row), then a list of matrices will be
##' returned, one for each covariate combination.
##' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
##' @seealso \code{\link{summary.flexsurvreg}}
##' @references Mandel, M. (2013). "Simulation based confidence intervals for
##' functions with complicated derivatives." The American Statistician (in
##' press).
##' @keywords models
##' @examples
##'
##'     fite <- expertsurv:::flexsurvreg(Surv(futime, fustat) ~ age, data = ovarian, dist="exp")
##'     expertsurv:::normboot.flexsurvreg(fite, B=10, newdata=list(age=50))
##'     expertsurv:::normboot.flexsurvreg(fite, B=10, X=matrix(50,nrow=1))
##'     expertsurv:::normboot.flexsurvreg(fite, B=10, newdata=list(age=0))  ## closer to...
##'     fite$res
##' @noRd
normboot.flexsurvreg <- function (x, B, newdata = NULL, X = NULL, transform = FALSE, 
          raw = FALSE, tidy = FALSE, rawsim = NULL,MLE = FALSE){
  if (x$ncovs > 0 && !raw) {
    if (is.null(X)) {
      if (is.null(newdata)) 
        stop("neither \"newdata\" nor \"X\" supplied")
      X <- form.model.matrix(x, as.data.frame(newdata))
    }
  }
  else X <- as.matrix(0, nrow = 1, ncol = 1)
  sim <- matrix(nrow = B, ncol = nrow(x$res))
  colnames(sim) <- rownames(x$res)
  if (is.na(x$cov[1])) 
    stop("Covariance matrix not available from non-converged model")
  if (is.null(rawsim)) {
    
    if(MLE){
      sim[, x$optpars] <- matrix(rep(x$opt$par, B),nrow = B, byrow = TRUE)  
    }else{
      sim[, x$optpars] <- rmvnorm(B, x$opt$par, x$cov) 
    }
    sim[, x$fixedpars] <- rep(x$res.t[x$fixedpars, "est"], 
                              each = B)
    rawsim <- sim
  }
  else {
    sim <- rawsim
  }
  if (x$ncovs > 0 && !raw) {
    beta <- sim[, x$covpars, drop = FALSE]
    if (nrow(X) == 1) {
      res <- sim[, x$dlist$pars, drop = FALSE]
      res <- add.covs(x = x, pars = res, beta = beta, X = X, 
                      transform = transform)
    }
    else {
      res <- vector(nrow(X), mode = "list")
      for (i in 1:nrow(X)) {
        res[[i]] <- sim[, x$dlist$pars, drop = FALSE]
        res[[i]] <- add.covs(x = x, pars = res[[i]], 
                             beta = beta, X = X[i, , drop = FALSE], transform = transform)
      }
    }
  }
  else {
    res <- sim
    if (!transform) {
      for (j in seq_along(x$dlist$pars)) {
        res[, j] <- x$dlist$inv.transforms[[j]](res[, 
                                                    j])
      }
    }
  }
  if (tidy && is.list(res)) {
    res <- cbind(covno = rep(1:nrow(X), each = B), repno = rep(1:B, 
                                                               nrow(X)), do.call("rbind", res))
    res <- as.data.frame(res)
  }
  attr(res, "X") <- X
  attr(res, "rawsim") <- rawsim
  res
}

### Compute CIs for survival, cumulative hazard, hazard, or user
### defined function, at supplied times t and covariates X, using
### random sample of size B from the assumed MVN distribution of MLEs.

normbootfn.flexsurvreg <- function(x, t, start, newdata=NULL, X=NULL, fn, B, rawsim=NULL){
    sim <- normboot.flexsurvreg(x, B, newdata=newdata, X=X, rawsim=rawsim)
    X <- attr(sim, "X")
    if (!is.list(sim)) sim <- list(sim)
    ret <- array(NA_real_, dim=c(nrow(X), B, length(t)))
    fncall0 <- list(t,start)
    for (j in seq_along(x$aux))
      fncall0[[names(x$aux)[j]]] <- x$aux[[j]]
    for (k in 1:nrow(X)){
        for (i in seq_len(B)) {
          fncall <- c(fncall0, lapply(sim[[k]][i,seq_along(x$dlist$pars)], function(z) z))
          ret[k,i,] <- do.call(fn, fncall)
        }
    }
    if (nrow(X)==1) ret[1,,,drop=FALSE] else ret
}
