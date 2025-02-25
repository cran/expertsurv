
#' Elicit survival judgements interactively and estimate survival models
#'
#' Opens up a web browser (using the shiny package), from which you can specify
#' judgements and fit distributions for multiple timepoints and experts.
#' Plots of the fitted density functions are provided overlayed on the survival data (where appropriate).
#'
#' Once the elicitation is complete the analysis can be run.
#' Click "Download R objects" to download the ``expertsurv`` object generated from the analysis.
#' Click "Download report" to generate a report including plots and parameter values for the parametric survival models.
#' For detailed instructions use \code{browseVignettes("expertsurv")}
#' @return
#' If "Download R objects" selected an \code{expertsurv} object containing the results of the elicitation and analysis. The object includes:
#' \itemize{
#'   \item \code{model}: The fitted survival models.
#'   \item \code{parameters}: The estimated parameters for each model.
#'	}
#' If "Download report" selected either .html, .pdf or .docx document is downloaded with:
#' \itemize{
#'   \item \code{plot}: Plots of the fitted survival overlayed on the Kaplan Meier data and plot of expert opinion as probability distributions.
#'   \item \code{summary}: A summary report of the analysis including goodness of fit and parameter values.
#'}
#' @param compile_mods list of compiled stan models generated by \link[=compile_stan]{compile_stan}. Supplying the compiled stan models will greatly speed up the computation of the Bayesian analysis (otherwise each time a stan model is run it will be compiled (and not reused between runs)).
#' @export
#' @author Philip Cooney <phcooney@@tcd.ie>
#' @examples
#' if (interactive()) {
#' elicit_surv()
#' }
elicit_surv <- function (compile_mods = NULL){
  
  
  # Set compile_mods in a global environment variable
  .GlobalEnv$compile_mods <- compile_mods
  
  
  library("expertsurv", character.only = TRUE)
  
  required_packages <- c("shiny", "shinyWidgets", "shinycssloaders", "shinyjs", "shinyMatrix", "shinybusy")
  
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  
  if (length(missing_packages) > 0) {
    stop("You need to install the following R packages to run the application: ", paste(missing_packages, collapse = ", "))
  }
  
  app_dir <- system.file("app/app.R", package = "expertsurv")
  if (app_dir == "") {
    stop("Could not find Shiny app directory. Try re-installing `expertsurv`.", call. = FALSE)
  }
  
  
  shiny::runApp(app_dir, display.mode = "normal",launch.browser = TRUE)
  
}
# 
# tmpfun <- get("elicit_surv", envir = asNamespace("expertsurv"))
# environment(elicit_surv) <- environment(tmpfun)
# attributes(elicit_surv) <- attributes(tmpfun)  
# assignInNamespace("elicit_surv", elicit_surv, ns="expertsurv")
