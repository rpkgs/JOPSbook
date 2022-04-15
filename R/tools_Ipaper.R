#' @export 
listk <- function(...) {
  cols <- as.list(substitute(list(...)))[-1]
  vars <- names(cols)
  Id_noname <- if (is.null(vars)) {
    seq_along(cols)
  } else {
    which(vars == "")
  }
  if (length(Id_noname) > 0) {
    vars[Id_noname] <- sapply(cols[Id_noname], deparse)
  }
  x <- setNames(list(...), vars)
  return(x)
}

#' @export
"%||%" <- function(a, b) {
  if (!is.null(a)) a else b
}

#' @export 
require2 <- function(pkg, ...) {
  pkgname = deparse(substitute(pkg))
  if (!require(pkgname, character.only = TRUE)){
    pak::pkg_install(pkgname)
    require(pkgname, character.only = TRUE)
  }
}
