#' Engine for Probabilistic Sensitivity Analysis on the survival curves
#' 
#' Creates the survival curves for the fitted model(s) - Original code from ``survHE``
#' 
#' 
#' @param fit the result of the call to the \code{fit.models} function,
#' containing the model fitting (and other relevant information)
#' @param mod the index of the model. Default value is 1, but the user can
#' choose which model fit to visualise, if the call to fit.models has a vector
#' argument for distr (so many models are fitted & stored in the same object)
#' @param t the time vector to be used for the estimation of the survival curve
#' @param newdata a list (of lists), specifiying the values of the covariates
#' at which the computation is performed. For example
#' \code{list(list(arm=0),list(arm=1))} will create two survival curves, one
#' obtained by setting the covariate \code{arm} to the value 0 and the other by
#' setting it to the value 1. In line with \code{flexsurv} notation, the user
#' needs to either specify the value for *all* the covariates or for none (in
#' which case, \code{newdata=NULL}, which is the default). If some value is
#' specified and at least one of the covariates is continuous, then a single
#' survival curve will be computed in correspondence of the average values of
#' all the covariates (including the factors, which in this case are expanded
#' into indicators).
#' @param nsim The number of simulations from the distribution of the survival
#' curves. Default at \code{nsim=1}, in which case uses the point estimate for
#' the relevant distributional parameters and computes the resulting survival
#' curve
#' @param ...  Additional options
#' @return A list with survival times for the fitted models
#' @author Gianluca Baio
#' @seealso  \code{psa.plot} (for example)
#' @keywords Survival models Bootstrap Probabilistic sensitivity analysis
#' @references 
#' \insertRef{Baio.2020}{expertsurv}
#' @export make.surv
make.surv <- function(fit,mod=1,t=NULL,newdata=NULL,nsim=1,...) {
  ## Creates the survival curves for the fitted model(s)
  # fit = the result of the call to the fit.models function, containing the model fitting (and other relevant information)
  # mod = the index of the model. Default value is 1, but the user can choose which model fit to visualise, 
  #     if the call to fit.models has a vector argument for distr (so many models are fitted & stored in the same object)
  # t = the time framework to be used for the estimation of the survival curve
  # newdata = a list (of lists), specifiying the values of the covariates at which the computation is performed. For example
  #           'list(list(arm=0),list(arm=1))' will create two survival curves, one obtained by setting the covariate 'arm'
  #           to the value 0 and the other by setting it to the value 1. In line with 'flexsurv' notation, the user needs
  #           to either specify the value for *all* the covariates or for none (in which case, 'newdata=NULL', which is the
  #           default). If some value is specified and at least one of the covariates is continuous, then a single survival
  #           curve will be computed in correspondence of the average values of all the covariates (including the factors, 
  #           which in this case are expanded into indicators). The order of the variables in the list *must* be the same
  #           as in the formula used for the model
  # nsim = the number of simulations from the distribution of the survival curves. Default at nsim=1, in which case
  #          uses the point estimate for the relevant distributional parameters and computes the resulting survival curve
  # ... = additional options
  
  # Defines list with optional parameters
  exArgs <- list(...)
  
  # The Poly-Weibull/HMC model is special and needs a custom function
  if(fit$misc$model_name[mod]=="pow") {
    pwstuff=make_surv_pw(fit=fit,mod=mod,t=t,newdata=newdata,nsim=nsim,exArgs)
    sim=pwstuff$sim
    mat=pwstuff$mat
    X=pwstuff$X
    t=pwstuff$t
  } else {
    # Any other model/method has a streamlined process using helper functions
    # Extracts the model object and the data from the survHE output
    m <- fit$models[[mod]]
    data <- fit$misc$data
    # Create a vector of times, if the user hasn't provided one, based on the observed data
    if(is.null(t)) {
      t <- sort(unique(fit$misc$km$time))
    }
    # By default uses the mean to compute summary statistics (eg for the case when nsim=1)
    # but the user can specify the median, which works better for very skewed parameters
    # which happens for example when trying to fit a model with too many parameters, which
    # results in a huge uncertainty, blowing up the estimates and making even the mean survival
    # curve impossible to compute
    if(exists("summary_stat",exArgs)){summary_stat=exArgs$summary_stat} else {summary_stat="mean"}

    # Makes sure the distribution name(s) vector is in a useable format
    dist <- fit$misc$model_name[mod]
    
    # Now creates the profile of covariates for which to compute the survival curves
    X <- make_profile_surv(fit$misc$formula,data,newdata)
    
    # Draws a sample of nsim simulations from the distribution of the model parameters
    sim <- do.call(paste0("make_sim_",fit$method),
                   args=list(m=m,t=t,X=X,nsim=nsim,newdata=newdata,dist=dist,data=data,
                             formula=fit$misc$formula,summary_stat=summary_stat)
    )
    # Computes the survival curves - first in matrix form with all the simulations
    # Needs to add more inputs for the case of hmc/rps
    if(fit$method=="hmc" & dist=="rps") {
      exArgs$data.stan <- fit$misc$data.stan[[mod]]
      t[t==0] <- min(0.00001,min(t[t>0]))
    }
    if(fit$method=="mle" & dist=="rps") {
      exArgs$knots=fit$models[[mod]]$knots
    }
    mat <- do.call(compute_surv_curve,
                   args=list(sim=sim,exArgs=exArgs,nsim=nsim,dist=dist,t=t,method=fit$method,X=X) 
    )
  }
  
  # Finally computes the actual survival curves, in summary forms
  if (nsim==1) {
    # If nsim=1 then only save the point estimates of the survival curve
    S <-lapply(mat,function(x) x %>% mutate(S=rowMeans(select(.,contains("S")))) %>% select(t,S))
  } else {
    # If nsim>1 then also give the lower and upper quartile of the underlying distribution
    # avoids no visibile binding for global variable
    . = low = upp = NULL
    S <- lapply(mat,function(x) {
      x %>% mutate(S=rowMeans(select(.,contains("S"))),
                   low=(apply(x %>% select(contains("S")),1,quantile,.025)),
                   upp=(apply(x %>% select(contains("S")),1,quantile,.975))
                   ) %>% select(t,S,low,upp)
    })
  }
  
  # Formats the output
  list(S=S,sim=sim,nsim=nsim,mat=mat,des.mat=X,times=t)
}
