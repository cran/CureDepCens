#' Plot the survival function
#'
#' @aliases plot_cure
#' @export
#' @description This graph helps to visualize the survival function.
#' @param object an object of the class "dcensoring".
#' @param scenario which defines the scenario in the graph (t: failure times, c: dependent censoring times, or both).
#' @details In order to smooth the line presented in the graph, we used the 'lowess' function. So, it can result in a non-monotonous survival function.
#' @return a survival function graph of the fitted model.
#'
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
#'
#' plot_cure(fit, scenario = "t")
#'}
#'
plot_cure <- function(object, scenario  = c("t", "c", "both")){

  scenario  <- match.arg(scenario)

  bmax <- object$bmax
  #Caso MEP
  if (is.null(bmax) == FALSE){
    switch(scenario ,
           "t" = invisible(plot_mep_t(object)),
           "c" = invisible(plot_mep_c(object)),
           "both"= invisible(plot_mep_t(object) + plot_mep_c(object)))
  }

  #Caso Weibull
  else{
    switch(scenario ,
           "t" = invisible(plot_weibull_t(object)),
           "c" = invisible(plot_weibull_c(object)),
           "both"= invisible(plot_weibull_t(object) + plot_weibull_c(object)))
  }
}
