#' @export
dist.gen <- function(x,method="euclidean", ...){
	if ( method %in% c("spearman","pearson","kendall") ) as.dist( 1 - cor(t(x),method=method,...) ) else dist(x,method=method,...)
}
