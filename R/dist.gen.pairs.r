#' @export
dist.gen.pairs <- function(x,y,...) dist.gen(t(cbind(x,y)),...)
