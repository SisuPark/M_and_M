## "In order to avoid overfitting we will select a model based on the AIC/BIC
## criterion. The selection is simply finding the lowest value where in general AIC
## allows slightly more complex models compared to BIC."
##
## from
## https://www.r-bloggers.com/an-exercise-in-non-linear-modeling/

##' A simple function for updating the formula and extracting
##' the information criteria
##'
##' @param no A number that is used together with the add_var_str
##' @param fit A regression fit that is used for the update
##' @param rm_var The variable that is to be substituted
##' @param add_var_str A sprintf() string that accepts the no
##'  parameter for each update
##' @param ic_fn The information criteria function (AIC/BIC)

getInfCrit <- function(no, fit, rm_var, add_var_str, ic_fn)
{
  new_var <- sprintf(add_var_str, no)
  updt_frml <- as.formula(sprintf(".~.-%s+%s", rm_var, new_var))
  ret <- ic_fn(update(fit, updt_frml))
  names(ret) <- new_var
  return(ret)
}
