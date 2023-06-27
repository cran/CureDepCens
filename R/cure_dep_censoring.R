#---------------------------------------------

#' Cure Dependent Censoring model
#' @aliases cure_dep_censoring
#' @export
#' @description cure_dep_censoring can be used to fit survival data with cure fraction and dependent censoring. It can also be utilized to take into account informative censoring.
#' @param formula an object of class "formula": should be used as 'time ~ cure covariates | informative covariates'.
#' @param data a data frame, list or environment containing the variables.
#' @param delta_t Indicator function of the event of interest.
#' @param delta_c Indicator function of the dependent censoring.
#' @param ident Cluster variable.
#' @param Num_intervals Number of intervals of the time grid (mep only).
#' @param dist distribution to be used in the model adjustment, specifies the marginal distribution of times (must be either weibull or mep).
#' @details This function estimates the parameters of the Piecewise exponential model (dist = "mep") or Weibull model (dist = "weibull") with cure rate and dependent censoring, considering the frailty model to estimate the clusters variability and a parameter that captures the dependence between failure and dependent censoring times.
#' @return cure_dep_censoring returns an object of class "dcensoring" containing the results of the fitted models.
#' An object of class "dcensoring" is a list containing at least the following components:
#' \itemize{
#'   \item \code{param_est} a vector containing estimated parameters (dependency parameter, regression coefficients associated with the cure rate, regression coefficients associated with dependent censoring times, and time distribution parameters (Weibull or piecewise exponential)).
#'   \item \code{stde} a vector containing the estimated standard errors of the estimated parameters vector.
#'   \item \code{crit} a vector containing the information criteria, Akaike's information criterion (AIC), Bayesian information criterion (BIC), Hannan-Quinn information criterion (HQ), calculated according to Louis, T. A. (1982).
#'   \item \code{pvalue} p-value of the estimated parameters vector.
#'   \item \code{n} number of observations in the dataset.
#'   \item \code{p} number of covariates associated with the cure fraction.
#'   \item \code{q} number of covariates associated with the dependent censoring times (informative censoring times or competitive risk times).
#'   \item \code{formula} formula used in the function call.
#'   \item \code{terms} the terms object used, containing the covariates associated with the cure fraction and with the dependent censoring times.
#'   \item \code{labels1} labels of the covariates associated with the cure fraction.
#'   \item \code{labels2} labels of the covariates associated with the dependent censoring times.
#'   \item \code{risco_a_T} a vector containing the cumulative baseline hazar of failure times.
#'   \item \code{risco_a_C} a vector containing the cumulative baseline hazar of dependent censoring times.
#'   \item \code{bi} a matrix containing the generated frailties, one of the outputs of the function cure_dep_censoring, in which the individuals are in the rows and the Monte Carlo replicas in the columns.
#'   \item \code{X_Cure} a matrix of variables associated with the cure fraction.
#'   \item \code{X_C} a matrix of variables associated with the dependent censoring times.
#'   \item \code{time} a vector of the observable times.
#' }
#' @examples
#' \donttest{
#' library(CureDepCens)
#'
#' delta_t = ifelse(Dogs_MimicData$cens==1,1,0)
#' delta_c = ifelse(Dogs_MimicData$cens==2,1,0)
#'
#' fit <- cure_dep_censoring(formula = time ~ x1_cure + x2_cure | x_c1 + x_c2,
#'                           data = Dogs_MimicData,
#'                           delta_t = delta_t,
#'                           delta_c = delta_c,
#'                           ident = Dogs_MimicData$ident,
#'                           dist = "mep")
#'}
cure_dep_censoring <- function(formula, data, delta_t, delta_c, ident, dist = c("weibull", "mep"), Num_intervals = 3){

  dist <- match.arg(dist)

  switch(dist,
         "weibull" = model_Weibull_dep(formula=formula, data=data, delta_t=delta_t, delta_c=delta_c, ident=ident),
         "mep" = model_MEP_dep(formula=formula, data=data, delta_t=delta_t, delta_c=delta_c, ident=ident, Num_intervals = 3))
}
