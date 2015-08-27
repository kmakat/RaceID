#' @export
KmeansCBI <- function (data, krange, k = NULL, scaling = FALSE, runs = 1, 
    criterion = "ch", method="euclidean",...) 
{
    if (!is.null(k)) 
        krange <- k
    if (!identical(scaling, FALSE)) 
        sdata <- scale(data, center = TRUE, scale = scaling)
    else sdata <- data
    c1 <- Kmeansruns(sdata, krange, runs = runs, criterion = criterion, method = method,
        ...)
    partition <- c1$cluster
    cl <- list()
    nc <- krange
    for (i in 1:nc) cl[[i]] <- partition == i
    out <- list(result = c1, nc = nc, clusterlist = cl, partition = partition, 
        clustermethod = "kmeans")
    out
}