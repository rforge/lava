##' @export
vcov.lvmfit <- function(object,...) {
  res <- object$vcov
  if ("lvm.missing"%in%class(object)) {
    resnames <- names(pars(object))
  } else {
    resnames <- coef(Model(object),fix=FALSE, mean=object$control$meanstructure)
  }
  colnames(res) <- rownames(res) <- resnames
  return(res)
}

##' @export
vcov.multigroupfit <- function(object,...) {
  res <- object$vcov
  colnames(res) <- rownames(res) <- object$model$name
  return(res)
}

